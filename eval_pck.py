import argparse

import numpy as np
from numpy import array as arr
import scipy.io as sio

from config import load_config
from dataset.factory import create as dataset_create


def enclosing_rect(points):
    xs = points[:, 0]
    ys = points[:, 1]
    return np.array([np.amin(xs), np.amin(ys), np.amax(xs), np.amax(ys)])


def rect_size(rect):
    return np.array([rect[2]-rect[0], rect[3]-rect[1]])


def print_results(pck, cfg):
    str = ""
    for heading in (cfg.all_joints_names + ["total"]):
        str += " & " + heading
    print(str)

    str = ""
    all_joint_ids = cfg.all_joints + [np.arange(cfg.num_joints)]
    for j_ids in all_joint_ids:
        j_ids_np = arr(j_ids)
        pck_av = np.mean(pck[j_ids_np])
        str += " & {0:.1f}".format(pck_av)
    print(str)


def eval_pck(cfg):
    dataset = dataset_create(cfg)
    filename = 'predictions.mat'
    pred = sio.loadmat(filename)

    joints = pred['joints']
    pck_ratio_thresh = cfg.pck_threshold

    num_joints = cfg.num_joints
    num_images = joints.shape[1]

    pred_joints = np.zeros((num_images, num_joints, 2))
    gt_joints = np.zeros((num_images, num_joints, 2))
    pck_thresh = np.zeros((num_images, 1))
    gt_present_joints = np.zeros((num_images, num_joints))

    for k in range(num_images):
        pred = joints[0,k]
        gt = dataset.data[k].joints[0]
        if gt.shape[0] == 0:
            continue
        gt_joint_ids = gt[:, 0].astype('int32')
        rect = enclosing_rect(gt[:,1:3])
        pck_thresh[k] = pck_ratio_thresh*np.amax(rect_size(rect))

        gt_present_joints[k, gt_joint_ids] = 1
        gt_joints[k, gt_joint_ids, :] = gt[:,1:3]
        pred_joints[k, :, :] = pred[:,0:2]

    dists = np.sqrt(np.sum((pred_joints - gt_joints)**2, axis=2))
    correct = dists <= pck_thresh

    num_all = np.sum(gt_present_joints, axis=0)

    num_correct = np.zeros((num_joints, ))
    for j_id in range(num_joints):
        num_correct[j_id] = np.sum(correct[gt_present_joints[:,j_id] == 1, j_id], axis=0)

    pck = num_correct/num_all*100.0

    print_results(pck, cfg)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args, unparsed = parser.parse_known_args()

    cfg = load_config()

    eval_pck(cfg)