import math
import numpy as np

import scipy.spatial

import matplotlib.pyplot as plt

import munkres

from util.visualize import check_point, _npcircle
from util import visualize


min_match_dist = 200
marker_size = 5

draw_conf_min_count = 3


def get_ref_points(person_conf):
    avg_conf = np.sum(person_conf, axis=1) / person_conf.shape[1]

    # last points is tip of the head -> use it as reference
    ref_points = person_conf[:, -1, :]

    # use average of other points if head tip is missing
    emptyidx = (np.sum(ref_points, axis=1) == 0)
    ref_points[emptyidx, :] = avg_conf[emptyidx, :]

    return ref_points


class PersonDraw:
    def __init__(self):
        self.mk = munkres.Munkres()

        self.prev_person_conf = np.zeros([0, 1])
        self.prev_color_assignment = None

        # generated colors from http://tools.medialab.sciences-po.fr/iwanthue/
        track_colors_str = ["#F5591E",
                            "#3870FB",
                            "#FE5DB0",
                            "#B4A691",
                            "#43053F",
                            "#3475B1",
                            "#642612",
                            "#B3B43D",
                            "#DD9BFE",
                            "#28948D",
                            "#E99D53",
                            "#012B46",
                            "#9D2DA3",
                            "#04220A",
                            "#62CB22",
                            "#EE8F91",
                            "#D71638",
                            "#00613A",
                            "#318918",
                            "#B770FF",
                            "#82C091",
                            "#6C1333",
                            "#973405",
                            "#B19CB2",
                            "#F6267B",
                            "#284489",
                            "#97BF17",
                            "#3B899C",
                            "#931813",
                            "#FA76B6"]

        self.track_colors = [(int(s[1:3], 16), int(s[3:5], 16), int(s[5:7], 16)) for s in track_colors_str]

    def draw(self, visim, dataset, person_conf):
        minx = 2 * marker_size
        miny = 2 * marker_size
        maxx = visim.shape[1] - 2 * marker_size
        maxy = visim.shape[0] - 2 * marker_size

        num_people = person_conf.shape[0]
        color_assignment = dict()

        # MA: assign same color to matching body configurations
        if self.prev_person_conf.shape[0] > 0 and person_conf.shape[0] > 0:
            ref_points = get_ref_points(person_conf)
            prev_ref_points = get_ref_points(self.prev_person_conf)

            # MA: this munkres implementation assumes that num(rows) >= num(columns)
            if person_conf.shape[0] <= self.prev_person_conf.shape[0]:
                cost_matrix = scipy.spatial.distance.cdist(ref_points, prev_ref_points)
            else:
                cost_matrix = scipy.spatial.distance.cdist(prev_ref_points, ref_points)

            assert (cost_matrix.shape[0] <= cost_matrix.shape[1])

            conf_assign = self.mk.compute(cost_matrix)

            if person_conf.shape[0] > self.prev_person_conf.shape[0]:
                conf_assign = [(idx2, idx1) for idx1, idx2 in conf_assign]
                cost_matrix = cost_matrix.T

            for pidx1, pidx2 in conf_assign:
                if cost_matrix[pidx1][pidx2] < min_match_dist:
                    color_assignment[pidx1] = self.prev_color_assignment[pidx2]

        print("#tracked objects:", len(color_assignment))

        free_coloridx = sorted(list(set(range(len(self.track_colors))).difference(set(color_assignment.values()))),
                               reverse=True)

        for pidx in range(num_people):
            # color_idx = pidx % len(self.track_colors)
            if pidx in color_assignment:
                color_idx = color_assignment[pidx]
            else:
                if len(free_coloridx) > 0:
                    color_idx = free_coloridx[-1]
                    free_coloridx = free_coloridx[:-1]
                else:
                    color_idx = np.random.randint(len(self.track_colors))

                color_assignment[pidx] = color_idx

            assert (color_idx < len(self.track_colors))

            if np.sum(person_conf[pidx, :, 0] > 0) < draw_conf_min_count:
                continue

            for kidx1, kidx2 in dataset.get_pose_segments():
                p1 = (int(math.floor(person_conf[pidx, kidx1, 0])), int(math.floor(person_conf[pidx, kidx1, 1])))
                p2 = (int(math.floor(person_conf[pidx, kidx2, 0])), int(math.floor(person_conf[pidx, kidx2, 1])))

                if check_point(p1[0], p1[1], minx, miny, maxx, maxy) and check_point(p2[0], p2[1], minx, miny, maxx,
                                                                                     maxy):
                    color = np.array(self.track_colors[color_idx][::-1], dtype=np.float64) / 255.0
                    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], marker='o', linestyle='solid', linewidth=2.0, color=color)


        self.prev_person_conf = person_conf
        self.prev_color_assignment = color_assignment



keypoint_colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255), (0, 255, 255), (255, 0, 255), (255, 255, 0),
                   (255, 0, 0), (0, 255, 0), (0, 0, 255), (0, 255, 255), (255, 0, 255), (255, 255, 0),
                   (255, 0, 0), (0, 255, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255)]

def visualize_detections(cfg, img, detections):
    vis_scale = 1.0
    marker_size = 4

    minx = 2 * marker_size
    miny = 2 * marker_size
    maxx = img.shape[1] - 2 * marker_size
    maxy = img.shape[0] - 2 * marker_size

    unPos = detections.coord
    joints_to_visualise = range(cfg.num_joints)
    visim_dets = img.copy()
    for pidx in joints_to_visualise:
        for didx in range(unPos[pidx].shape[0]):
            cur_x = unPos[pidx][didx, 0] * vis_scale
            cur_y = unPos[pidx][didx, 1] * vis_scale

            # / cfg.global_scale

            if check_point(cur_x, cur_y, minx, miny, maxx, maxy):
                _npcircle(visim_dets,
                          cur_x, cur_y,
                          marker_size,
                          keypoint_colors[pidx])
    return visim_dets
