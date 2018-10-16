import logging
import threading

import tensorflow as tf
import tensorflow.contrib.slim as slim

from ..config import load_config, cfg_from_file
from ..dataset.factory import create as create_dataset
from ..nnet.net_factory import pose_net
from ..nnet.pose_net import get_batch_spec
from ..util.logging import setup_logging


class LearningRate(object):
    def __init__(self, cfg):
        self.steps = cfg.multi_step
        self.current_step = 0

    def get_lr(self, iteration):
        lr = self.steps[self.current_step][0]
        if iteration == self.steps[self.current_step][1]:
            self.current_step += 1

        return lr


def setup_preloading(batch_spec):
    placeholders = {name: tf.placeholder(tf.float32, shape=spec) for (name, spec) in batch_spec.items()}
    names = placeholders.keys()
    placeholders_list = list(placeholders.values())

    QUEUE_SIZE = 20

    q = tf.FIFOQueue(QUEUE_SIZE, [tf.float32]*len(batch_spec))
    enqueue_op = q.enqueue(placeholders_list)
    batch_list = q.dequeue()

    batch = {}
    for idx, name in enumerate(names):
        batch[name] = batch_list[idx]
        batch[name].set_shape(batch_spec[name])
    return batch, enqueue_op, placeholders


def load_and_enqueue(sess, enqueue_op, coord, dataset, placeholders):
    while not coord.should_stop():
        batch_np = dataset.next_batch()
        food = {pl: batch_np[name] for (name, pl) in placeholders.items()}
        sess.run(enqueue_op, feed_dict=food)


def start_preloading(sess, enqueue_op, dataset, placeholders):
    coord = tf.train.Coordinator()

    t = threading.Thread(target=load_and_enqueue,
                         args=(sess, enqueue_op, coord, dataset, placeholders))
    t.start()

    return coord, t


def get_optimizer(loss_op, cfg):
    learning_rate = tf.placeholder(tf.float32, shape=[])

    if cfg.optimizer == "sgd":
        optimizer = tf.train.MomentumOptimizer(learning_rate=learning_rate, momentum=0.9)
    elif cfg.optimizer == "adam":
        optimizer = tf.train.AdamOptimizer(cfg.adam_lr)
    else:
        raise ValueError('unknown optimizer {}'.format(cfg.optimizer))
    train_op = slim.learning.create_train_op(loss_op, optimizer)

    return learning_rate, train_op


def train(cfg_filename=None, dataset_filename=None,
    max_to_keep=5, memfrac=None, disable_autotune=False,
    choose_gpu=None):
    """Train on a dataset
    
    cfg_filename : path to a configuration file
        If None, it will be automatically determined by load_config, which
        will look for pose_cfg.yaml either in the current directory or in
        os.environ("POSE_PARAM_PATH")
    
    dataset_filename : path to a dataset, or None
        If not None, this will override the parameter `dataset` 
        in the configuration file.
    
    max_to_keep : int or None
        Number of snapshots to keep on disk
        See documentation at tensorflow.train.Saver
    
    memfrac : float or None
        Fraction of memory to allocate on the GPU
        Specify a value between 0.0 and 1.0. If None, all available memory
        will be used. See additional documentation at tensorflow.ConfigProto
    
    disable_autotune : bool
        If True, set the environment variable CUDNN_USE_AUTOTUNE to '0'
    
    choose_gpu : string, int, or None
        If not None, set the environment variable CUDA_VISIBLE_DEVICES to this
        It will be converted to a string first
    """
    # Set environment variables
    if disable_autotune:
        os.environ['CUDNN_USE_AUTOTUNE'] = '0'
    
    if choose_gpu is not None:
        os.environ['CUDA_VISIBLE_DEVICES'] = str(choose_gpu)
    
    # Set up logging
    setup_logging()

    # Load configuration options
    if cfg_filename is None:
        # Use the old load_config method
        cfg = load_config()
    else:
        # Load the configuration options from the specified file
        cfg = cfg_from_file(cfg_filename)
    
    # Optionally set the dataset in cfg
    if dataset_filename is not None:
        cfg.dataset = dataset_filename
    
    # Create the dataset
    dataset = create_dataset(cfg)

    batch_spec = get_batch_spec(cfg)
    batch, enqueue_op, placeholders = setup_preloading(batch_spec)

    losses = pose_net(cfg).train(batch)
    total_loss = losses['total_loss']

    for k, t in losses.items():
        tf.summary.scalar(k, t)
    merged_summaries = tf.summary.merge_all()

    variables_to_restore = slim.get_variables_to_restore(include=["resnet_v1"])
    restorer = tf.train.Saver(variables_to_restore)
    
    # Initialize a Saver
    saver = tf.train.Saver(max_to_keep=max_to_keep)

    # ConfigProto for tensorflow session
    tfcfg = tf.ConfigProto()
    if memfrac is not None:
        tfcfg.gpu_options.per_process_gpu_memory_fraction = memfrac
    
    # Initialize the tensorflow session
    sess = tf.Session(config=tfcfg)
    
    coord, thread = start_preloading(sess, enqueue_op, dataset, placeholders)

    train_writer = tf.summary.FileWriter(cfg.log_dir, sess.graph)

    learning_rate, train_op = get_optimizer(total_loss, cfg)

    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())

    # Restore variables from disk.
    restorer.restore(sess, cfg.init_weights)

    max_iter = int(cfg.multi_step[-1][1])

    display_iters = cfg.display_iters
    cum_loss = 0.0
    lr_gen = LearningRate(cfg)

    for it in range(max_iter+1):
        current_lr = lr_gen.get_lr(it)
        [_, loss_val, summary] = sess.run([train_op, total_loss, merged_summaries],
                                          feed_dict={learning_rate: current_lr})
        cum_loss += loss_val
        train_writer.add_summary(summary, it)

        if it % display_iters == 0:
            average_loss = cum_loss / display_iters
            cum_loss = 0.0
            logging.info("iteration: {} loss: {} lr: {}"
                         .format(it, "{0:.4f}".format(average_loss), current_lr))

        # Save snapshot
        if (it % cfg.save_iters == 0 and it != 0) or it == max_iter:
            model_name = cfg.snapshot_prefix
            saver.save(sess, model_name, global_step=it)

    sess.close()
    coord.request_stop()
    coord.join([thread])

