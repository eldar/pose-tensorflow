import math
import os
import sys
from collections import namedtuple

import numpy as np

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + '/../lib/nms_cython/')
from nms_grid import nms_grid

Detections = namedtuple('Detections', ['coord', 'coord_grid', 'conf', 'pairwise'])


def pos_from_grid_raw(cfg, gridpos):
    return gridpos * cfg.stride + 0.5 * cfg.stride


def pos_from_gridpos_offset(cfg, gridpos, pred_offset):
    return pos_from_grid_raw(cfg, gridpos) + pred_offset


def make_nms_grid(nms_radius):
    nms_radius = math.ceil(nms_radius)
    dist_grid = np.zeros([2 * nms_radius + 1, 2 * nms_radius + 1], dtype=np.uint8)
    for yidx in range(dist_grid.shape[0]):
        for xidx in range(dist_grid.shape[1]):
            if (yidx - nms_radius) ** 2 + (xidx - nms_radius) ** 2 <= nms_radius ** 2:
                dist_grid[yidx][xidx] = 1
    return dist_grid


def extract_detections(cfg, scmap, locref, pairwise_diff):
    num_joints = cfg.num_joints
    num_pairwise_relations = pairwise_diff.shape[2]

    # get dist_grid
    dist_grid = make_nms_grid(cfg.nms_radius)

    unProb = [None] * num_joints
    unPos = [None] * num_joints
    unPos_grid = [None] * num_joints
    pairwiseDiff = [None] * num_joints

    # apply nms
    for p_idx in range(num_joints):
        # IMPORTANT, as C++ function expects row-major
        prob_map = np.ascontiguousarray(scmap[:, :, p_idx])
        # print(prob_map.flags) has to be C_CONTIGUOUS

        dets = nms_grid(prob_map, dist_grid, cfg.det_min_score)

        cur_prob = np.zeros([len(dets), 1], dtype=np.float64)
        cur_pos = np.zeros([len(dets), 2], dtype=np.float64)
        cur_pos_grid = np.zeros([len(dets), 2], dtype=np.float64)
        cur_pairwise = np.zeros([len(dets), num_pairwise_relations, 2], dtype=np.float64)

        for idx, didx in enumerate(dets):
            ix = didx % scmap.shape[1]
            iy = didx // scmap.shape[1]

            cur_prob[idx, 0] = scmap[iy, ix, p_idx]
            cur_pos_grid[idx, :] = pos_from_grid_raw(cfg, np.array([ix, iy]))
            cur_pos[idx, :] = cur_pos_grid[idx, :] + locref[iy, ix, p_idx, :]
            cur_pairwise[idx, :, :] = pairwise_diff[iy, ix, :, :]

        unProb[p_idx] = cur_prob
        unPos[p_idx] = cur_pos
        unPos_grid[p_idx] = cur_pos_grid
        pairwiseDiff[p_idx] = cur_pairwise

    return Detections(coord=unPos, coord_grid=unPos_grid, conf=unProb, pairwise=pairwiseDiff)
