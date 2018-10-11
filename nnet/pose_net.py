import re

import tensorflow as tf
import tensorflow.contrib.slim as slim
from tensorflow.contrib.slim.nets import resnet_v1

from dataset.pose_dataset import Batch
from nnet import losses


net_funcs = {'resnet_50': resnet_v1.resnet_v1_50,
             'resnet_101': resnet_v1.resnet_v1_101}


def prediction_layer(cfg, input, name, num_outputs):
    with slim.arg_scope([slim.conv2d, slim.conv2d_transpose], padding='SAME',
                        activation_fn=None, normalizer_fn=None,
                        weights_regularizer=slim.l2_regularizer(cfg.weight_decay)):
        with tf.variable_scope(name):
            pred = slim.conv2d_transpose(input, num_outputs,
                                         kernel_size=[3, 3], stride=2,
                                         scope='block4')
            return pred


def get_batch_spec(cfg):
    num_joints = cfg.num_joints
    batch_size = cfg.batch_size
    batch_spec = {
        Batch.inputs: [batch_size, None, None, 3],
        Batch.part_score_targets: [batch_size, None, None, num_joints],
        Batch.part_score_weights: [batch_size, None, None, num_joints]
    }
    if cfg.location_refinement:
        batch_spec[Batch.locref_targets] = [batch_size, None, None, num_joints * 2]
        batch_spec[Batch.locref_mask] = [batch_size, None, None, num_joints * 2]
    if cfg.pairwise_predict:
        batch_spec[Batch.pairwise_targets] = [batch_size, None, None, num_joints * (num_joints - 1) * 2]
        batch_spec[Batch.pairwise_mask] = [batch_size, None, None, num_joints * (num_joints - 1) * 2]
    return batch_spec


class PoseNet:
    def __init__(self, cfg):
        self.cfg = cfg

    def extract_features(self, inputs):
        net_fun = net_funcs[self.cfg.net_type]

        mean = tf.constant(self.cfg.mean_pixel,
                           dtype=tf.float32, shape=[1, 1, 1, 3], name='img_mean')
        im_centered = inputs - mean

        with slim.arg_scope(resnet_v1.resnet_arg_scope()):
            net, end_points = net_fun(im_centered, global_pool=False, output_stride=16, is_training=False)

        return net, end_points

    def prediction_layers(self, features, end_points, reuse=None, no_interm=False, scope='pose'):
        cfg = self.cfg

        num_layers = re.findall("resnet_([0-9]*)", cfg.net_type)[0]
        layer_name = 'resnet_v1_{}'.format(num_layers) + '/block{}/unit_{}/bottleneck_v1'

        out = {}
        with tf.variable_scope(scope, reuse=reuse):
            out['part_pred'] = prediction_layer(cfg, features, 'part_pred',
                                                cfg.num_joints)
            if cfg.location_refinement:
                out['locref'] = prediction_layer(cfg, features, 'locref_pred',
                                                 cfg.num_joints * 2)
            if cfg.pairwise_predict:
                out['pairwise_pred'] = prediction_layer(cfg, features, 'pairwise_pred',
                                                       cfg.num_joints * (cfg.num_joints - 1) * 2)
            if cfg.intermediate_supervision and not no_interm:
                interm_name = layer_name.format(3, cfg.intermediate_supervision_layer)
                block_interm_out = end_points[interm_name]
                out['part_pred_interm'] = prediction_layer(cfg, block_interm_out,
                                                           'intermediate_supervision',
                                                           cfg.num_joints)

        return out

    def get_net(self, inputs):
        net, end_points = self.extract_features(inputs)
        return self.prediction_layers(net, end_points)

    def test(self, inputs):
        heads = self.get_net(inputs)
        return self.add_test_layers(heads)

    def add_test_layers(self, heads):
        prob = tf.sigmoid(heads['part_pred'])
        outputs = {'part_prob': prob}
        if self.cfg.location_refinement:
            outputs['locref'] = heads['locref']
        if self.cfg.pairwise_predict:
            outputs['pairwise_pred'] = heads['pairwise_pred']
        return outputs

    def part_detection_loss(self, heads, batch, locref, pairwise, intermediate):
        cfg = self.cfg

        weigh_part_predictions = cfg.weigh_part_predictions
        part_score_weights = batch[Batch.part_score_weights] if weigh_part_predictions else 1.0

        def add_part_loss(pred_layer):
            return tf.losses.sigmoid_cross_entropy(batch[Batch.part_score_targets],
                                                   heads[pred_layer],
                                                   part_score_weights)

        loss = {}
        loss['part_loss'] = add_part_loss('part_pred')
        total_loss = loss['part_loss']
        if intermediate:
            loss['part_loss_interm'] = add_part_loss('part_pred_interm')
            total_loss = total_loss + loss['part_loss_interm']

        if locref:
            locref_pred = heads['locref']
            locref_targets = batch[Batch.locref_targets]
            locref_weights = batch[Batch.locref_mask]

            loss_func = losses.huber_loss if cfg.locref_huber_loss else tf.losses.mean_squared_error
            loss['locref_loss'] = cfg.locref_loss_weight * loss_func(locref_targets, locref_pred, locref_weights)
            total_loss = total_loss + loss['locref_loss']

        if pairwise:
            pairwise_pred = heads['pairwise_pred']
            pairwise_targets = batch[Batch.pairwise_targets]
            pairwise_weights = batch[Batch.pairwise_mask]

            loss_func = losses.huber_loss if cfg.pairwise_huber_loss else tf.losses.mean_squared_error
            loss['pairwise_loss'] = cfg.pairwise_loss_weight * loss_func(pairwise_targets, pairwise_pred,
                                                                         pairwise_weights)
            total_loss = total_loss + loss['pairwise_loss']

        # loss['total_loss'] = slim.losses.get_total_loss(add_regularization_losses=params.regularize)
        loss['total_loss'] = total_loss
        return loss

    def train(self, batch):
        cfg = self.cfg

        intermediate = cfg.intermediate_supervision
        locref = cfg.location_refinement
        pairwise = cfg.pairwise_predict

        heads = self.get_net(batch[Batch.inputs])
        return self.part_detection_loss(heads, batch, locref, pairwise, intermediate)
