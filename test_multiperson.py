import argparse
import logging
import os

import numpy as np
import scipy.io
import scipy.ndimage
import json
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

from config import load_config
from dataset.factory import create as create_dataset
from dataset.pose_dataset import Batch
from util.mscoco_util import pose_predict_with_gt_segm
from nnet.predict import *
from util import visualize
from multiperson.detections import extract_detections
from multiperson.predict import SpatialModel, eval_graph, get_person_conf_multicut
from multiperson.visualize import PersonDraw, visualize_detections

import matplotlib.pyplot as plt


def test_net(visualise, cache_scoremaps, development):
    logging.basicConfig(level=logging.INFO)

    cfg = load_config()
    dataset = create_dataset(cfg)
    dataset.set_shuffle(False)

    sm = SpatialModel(cfg)
    sm.load()

    draw_multi = PersonDraw()

    from_cache = "cached_scoremaps" in cfg
    if not from_cache:
        sess, inputs, outputs = setup_pose_prediction(cfg)

    if cache_scoremaps:
        out_dir = cfg.scoremap_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    pairwise_stats = dataset.pairwise_stats
    num_images = dataset.num_images if not development else min(10, dataset.num_images)
    coco_results = []

    for k in range(num_images):
        print('processing image {}/{}'.format(k, num_images-1))

        batch = dataset.next_batch()

        cache_name = "{}.mat".format(batch[Batch.data_item].coco_id)

        if not from_cache:
            outputs_np = sess.run(outputs, feed_dict={inputs: batch[Batch.inputs]})
            scmap, locref, pairwise_diff = extract_cnn_output(outputs_np, cfg, pairwise_stats)

            if cache_scoremaps:
                if visualise:
                    img = np.squeeze(batch[Batch.inputs]).astype('uint8')
                    pose = argmax_pose_predict(scmap, locref, cfg.stride)
                    arrows = argmax_arrows_predict(scmap, locref, pairwise_diff, cfg.stride)
                    visualize.show_arrows(cfg, img, pose, arrows)
                    visualize.waitforbuttonpress()
                    continue

                out_fn = os.path.join(out_dir, cache_name)
                dict = {'scoremaps': scmap.astype('float32'),
                        'locreg_pred': locref.astype('float32'),
                        'pairwise_diff': pairwise_diff.astype('float32')}
                scipy.io.savemat(out_fn, mdict=dict)
                continue
        else:
            #cache_name = '1.mat'
            full_fn = os.path.join(cfg.cached_scoremaps, cache_name)
            mlab = scipy.io.loadmat(full_fn)
            scmap = mlab["scoremaps"]
            locref = mlab["locreg_pred"]
            pairwise_diff = mlab["pairwise_diff"]

        detections = extract_detections(cfg, scmap, locref, pairwise_diff)
        unLab, pos_array, unary_array, pwidx_array, pw_array = eval_graph(sm, detections)
        person_conf_multi = get_person_conf_multicut(sm, unLab, unary_array, pos_array)

        if visualise:
            img = np.squeeze(batch[Batch.inputs]).astype('uint8')
            #visualize.show_heatmaps(cfg, img, scmap, pose)

            """
            # visualize part detections after NMS
            visim_dets = visualize_detections(cfg, img, detections)
            plt.imshow(visim_dets)
            plt.show()
            visualize.waitforbuttonpress()
            """

#            """
            visim_multi = img.copy()
            draw_multi.draw(visim_multi, dataset, person_conf_multi)

            plt.imshow(visim_multi)
            plt.show()
            visualize.waitforbuttonpress()
#            """


        if cfg.use_gt_segm:
            coco_img_results = pose_predict_with_gt_segm(scmap, locref, cfg.stride, batch[Batch.data_item].gt_segm,
                                                      batch[Batch.data_item].coco_id)
            coco_results += coco_img_results
            if len(coco_img_results):
                dataset.visualize_coco(coco_img_results, batch[Batch.data_item].visibilities)

    if cfg.use_gt_segm:
        with open('predictions_with_segm.json', 'w') as outfile:
            json.dump(coco_results, outfile)

    sess.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--novis', default=False, action='store_true')
    parser.add_argument('--cache', default=False, action='store_true')
    parser.add_argument('--dev', default=False, action='store_true')
    args, unparsed = parser.parse_known_args()

    test_net(not args.novis, args.cache, args.dev)
