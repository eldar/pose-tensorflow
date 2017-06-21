import math
import time
import sys
import os
from collections import namedtuple

import numpy as np
import scipy.io as sio

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + '/../lib/multicut_cython/')
from multicut import solve_nl_lmp

from dataset.pose_dataset import get_pairwise_index


def logit_transform(p):
    p = np.minimum(np.maximum(p, 1e-7), 1.0 - 1e-7)
    return np.log((1-p)/p)


def eval_graph(sm, detections):

    time_start = time.time()

    unary_prob = detections.conf
    coordinates = detections.coord
    pairwise_regr = detections.pairwise

    cidx_list = range(0, sm.num_keypoints)
    #cidx_list = [0, 1, 2]

    unary_counts = []

    for idx, cidx in enumerate(cidx_list):
        unary_counts.append(unary_prob[cidx].shape[0])
    
    num_unary = sum(unary_counts)
    num_pairwise = 0

    for idx1, cidx1 in enumerate(cidx_list):
        for cidx2 in cidx_list[idx1:]:
            if cidx1 == cidx2:
                num_pairwise += unary_prob[cidx1].shape[0]*(unary_prob[cidx1].shape[0] - 1) // 2
            else:
                num_pairwise += unary_prob[cidx1].shape[0]*unary_prob[cidx2].shape[0]

    pos_array = np.zeros([num_unary, 2], dtype=np.float64)
    unary_array = np.zeros([num_unary, 1], dtype=np.float64)
    pw_array = np.zeros([num_pairwise, 1], dtype=np.float64)
    pwidx_array = np.zeros([num_pairwise, 2], dtype=np.uint16)

    firstidx = 0
    firstidx_list = []
    for idx, cidx in enumerate(cidx_list):
        lastidx = firstidx + unary_prob[cidx].shape[0]
        unary_array[firstidx:lastidx] = unary_prob[cidx]
        pos_array[firstidx:lastidx] = coordinates[cidx]

        firstidx_list.append(firstidx)
        firstidx = lastidx

    firstidx = 0
    for idx1, cidx1 in enumerate(cidx_list):
        for idx2 in range(idx1, len(cidx_list)):

            if coordinates[cidx1].shape[0] > 0:
                cidx2 = cidx_list[idx2]

                if coordinates[cidx2].shape[0] > 0:

                    if not sm.need_this_pairwise(cidx1, cidx2):
                        continue

                    cur_prob, ptidx = sm.eval(cidx1, cidx2, detections)
                    lastidx = firstidx + cur_prob.shape[0]

                    ptidx[:, 0] += firstidx_list[idx1]
                    ptidx[:, 1] += firstidx_list[idx2]

                    pw_array[firstidx:lastidx, 0] = cur_prob
                    pwidx_array[firstidx:lastidx, :] = ptidx

                    firstidx = lastidx

    is_sparse_graph = True
    solver_type = False
    do_suppression = True
    logit_in_solver = False

    if unary_array.shape[0] > 0:

        unary_array_solver = unary_array if logit_in_solver else logit_transform(unary_array)
        pw_array_solver = pw_array if logit_in_solver else logit_transform(pw_array)

        time_start = time.time()

        res = solve_nl_lmp(unary_array_solver, pwidx_array, pw_array_solver,
                          is_sparse_graph, solver_type, do_suppression, logit_in_solver)

        unLab = np.array(res, dtype=np.uint64)

        firstidx = 0
        for cidx in cidx_list:
            lastidx = firstidx + unary_prob[cidx].shape[0]
            unLab[firstidx:lastidx, 0] = cidx
            firstidx = lastidx

    else:
        unLab = np.array([])

    return unLab, pos_array, unary_array, pwidx_array, pw_array

def get_person_conf_single(sm, unProb, pos_array, pwidx_array, pw_array):

    assert(len(unProb) == sm.num_keypoints)
    assert(pwidx_array.shape[0] == pw_array.shape[0])
    assert(pw_array.shape[1] == 1)

    det_type_idx = []
    
    firstidx = 0
    for pidx in range(len(unProb)):
        lastidx = firstidx + unProb[pidx].shape[0]

        curidx = np.array([False]*pos_array.shape[0])
        curidx[firstidx:lastidx] = True
        
        firstidx = lastidx
        det_type_idx.append(curidx)

    # num people == number of heads 
    head_idx = 13
    num_people = unProb[head_idx].shape[0]

    connect_graph = dict()
    connect_graph[13] = (12,)
    connect_graph[12] = (8, 9, 2, 3)
    connect_graph[8] = (7,)
    connect_graph[7] = (6,)
    connect_graph[9] = (10,)
    connect_graph[10] = (11,)
    connect_graph[2] = (1,)
    connect_graph[1] = (0,)
    connect_graph[3] = (4,)
    connect_graph[4] = (5,)

    person_conf = np.zeros([num_people, sm.num_keypoints, 2])
    SearchNode = namedtuple('SearchNode', ['pidx', 'kidx', 'hidx'])
    
    search_queue = []

    for pidx, hidx in enumerate(np.flatnonzero(det_type_idx[head_idx])):
        search_queue.append(SearchNode(pidx=pidx, kidx=head_idx, hidx=hidx))
        person_conf[pidx, head_idx, :] = pos_array[hidx, :]

    assert(len(search_queue) == num_people)

    while len(search_queue) > 0:
        node = search_queue.pop()

        pidx = node.pidx
        kidx = node.kidx
        hidx = node.hidx

        # find the closes match for current part 
        if kidx in connect_graph:

            for kidx2 in connect_graph[kidx]:
                best_hidx = None
                best_pw = None
                
                # search all pairwise with compatible type 
                for idx in range(pwidx_array.shape[0]):
                    if pwidx_array[idx, 0] == hidx or pwidx_array[idx, 1] == hidx:
                        idx2 = np.flatnonzero(pwidx_array[idx, :] != hidx)[0]
                        other_hidx = pwidx_array[idx, idx2]

                        if det_type_idx[kidx2][other_hidx]:
                            if pw_array[idx] > best_pw:
                                best_hidx = other_hidx
                                best_pw = pw_array[idx]
                        

                if best_pw > 0.5:
                    person_conf[pidx, kidx2, :] = pos_array[best_hidx, :]
                    search_queue.append(SearchNode(pidx = pidx, kidx = kidx2, hidx = best_hidx))

    return person_conf

    

def get_person_conf_multicut(sm, unLab, unary_array, pos_array):
    if unLab.shape[0] > 0:
        num_people = int(np.max(unLab[:, 1])) + 1
    else:
        num_people = 0

    person_conf = np.zeros([num_people, sm.num_keypoints, 2])
    sum_prob = np.zeros([num_people, sm.num_keypoints, 1])

    # compute weighted average of keypoints of the same time 
    for didx in range(unLab.shape[0]):
        kidx = unLab[didx,0]
        pidx = unLab[didx,1]

        person_conf[pidx, kidx, :] += pos_array[didx, :]*unary_array[didx]
        sum_prob[pidx, kidx] += unary_array[didx]

    print("num_people: ", num_people)

    for pidx in range(num_people):
        for kidx in range(sm.num_keypoints):
            if sum_prob[pidx, kidx] > 0:
                person_conf[pidx, kidx, :] /= sum_prob[pidx, kidx]

    return person_conf


def compute_angle(deltaX, deltaY):
    angle = np.arctan2(deltaY,deltaX)

    # atan2 is between [-pi, pi]
    # wrapMinusPiPifast(angle)

    assert((angle > math.pi).sum() == 0)
    assert((angle < -math.pi).sum() == 0)

    return angle
    #assert(all(angle <= math.pi) && all(angle >= -math.pi))


def wrap_angle(a):
    larger = a > math.pi
    smaller = a < -math.pi
    a[larger]  = a[larger] - 2*math.pi
    a[smaller] = a[smaller] + 2*math.pi

    return a


def compute_features(delta_real, delta_predicted):
    a = compute_angle(delta_real[:, 0], delta_real[:, 1])
    a_forward = compute_angle(delta_predicted[:, 0], delta_predicted[:, 1])
    delta1 = delta_real - delta_predicted

    dist = np.linalg.norm(delta1, axis=1, ord=2)
    abs_a = np.abs(wrap_angle(a - a_forward))
    return dist, abs_a


class SpatialModel:

    def __init__(self, cfg):
        self.num_keypoints = cfg.num_joints
        num_keypoints = cfg.num_joints
        self.cfg = cfg

        self.graph_dict = dict()

        self.same_part_pw_coef = 0.2

        self.X_min = [[None]*num_keypoints for idx in range(num_keypoints)]
        self.X_max = [[None]*num_keypoints for idx in range(num_keypoints)]
        self.w = [[None]*num_keypoints for idx in range(num_keypoints)]

    def load(self):
        for cidx1 in range(self.num_keypoints):
            for cidx2 in range(cidx1+1, self.num_keypoints):
                model_name  = "{}/spatial_model_cidx_{}_{}.mat".format(self.cfg.pairwise_model_dir, cidx1 + 1, cidx2 + 1)
                #print "loading:", model_name
                if not os.path.isfile(model_name):
                    continue
    
                spatial_model = sio.loadmat(model_name)

                self.X_max[cidx1][cidx2] = spatial_model['spatial_model']['training_opts'][0][0][0]['X_max'][0][0]
                self.X_min[cidx1][cidx2] = spatial_model['spatial_model']['training_opts'][0][0][0]['X_min'][0][0]
                self.w[cidx1][cidx2] = spatial_model['spatial_model']['log_reg'][0][0][0]['w'][0][0][:]

                #self.pair_class[cidx1][cidx2] = PairClassDistAngle(w=w, X_min=X_min, X_max=X_max)
                
        if not self.cfg.tensorflow_pairwise_order:
            tmpval = sio.loadmat(self.cfg.pairwise_stats_fn)
            self.graph = tmpval['graph']
            self.pairwise_means = tmpval['means']
            self.pairwise_std_devs = tmpval['std_devs']

            for gidx in range(self.graph.shape[0]):
                cidx1 = self.graph[gidx, 0]
                cidx2 = self.graph[gidx, 1]
                self.graph_dict[(cidx1, cidx2)] = gidx

    def get_fwd_bwd_index(self, cidx1, cidx2):
        if self.cfg.tensorflow_pairwise_order:
            fwd_idx = get_pairwise_index(cidx1, cidx2, self.cfg.num_joints)
            bwd_idx = get_pairwise_index(cidx2, cidx1, self.cfg.num_joints)
        else:
            fwd_idx = self.graph_dict[(cidx1+1, cidx2+1)]
            bwd_idx = self.graph_dict[(cidx2+1, cidx1+1)]
        return fwd_idx, bwd_idx

    def need_this_pairwise(self, cidx1, cidx2):
        if cidx1 == cidx2:
            return True
        sparse_graph = self.cfg.sparse_graph
        return not sparse_graph or [cidx1, cidx2] in sparse_graph

    def eval(self, cidx1, cidx2, detections):
        unPos = detections.coord

        idx_type1 = np.array(range(unPos[cidx1].shape[0]))
        idx_type2 = np.array(range(unPos[cidx2].shape[0]))

        assert(idx_type1.shape[0] > 0)
        assert(idx_type2.shape[0] > 0)

        num_edges = len(idx_type1) * len(idx_type2)

        tmpidx1, tmpidx2 = np.meshgrid(idx_type1, idx_type2)
        ptidx = np.hstack((tmpidx1.T.reshape((num_edges, 1)), tmpidx2.T.reshape((num_edges, 1))))

        if cidx1 != cidx2:
            cur_prob = self.compute_different_part_pairwise(cidx1, cidx2, detections, ptidx, num_edges)
        else:
            cur_prob = None
            ptidx = ptidx[ptidx[:, 0] < ptidx[:, 1]]

            delta = unPos[cidx2][ptidx[:, 1], :] - unPos[cidx1][ptidx[:, 0], :]
            dists = np.linalg.norm(delta, axis=1, ord=2)

            cur_prob = 1./(1 + np.exp(self.same_part_pw_coef*dists-7.5))

        return cur_prob, ptidx

    def compute_different_part_pairwise(self, cidx1, cidx2, detections, ptidx, num_edges):
        unPos = detections.coord
        unPos_grid = detections.coord_grid
        nextReg = detections.pairwise

        fwd_idx, bwd_idx = self.get_fwd_bwd_index(cidx1, cidx2)
        #print("cidxs ", cidx1, cidx2)
        #print("pairwise index", fwd_idx, bwd_idx)

        assert (ptidx.shape[0] > 0)

        delta_real_forward = unPos[cidx2][ptidx[:, 1], :] - unPos_grid[cidx1][ptidx[:, 0], :]
        delta_real_backward = unPos[cidx1][ptidx[:, 0], :] - unPos_grid[cidx2][ptidx[:, 1], :]

        delta_forward = nextReg[cidx1][ptidx[:, 0], fwd_idx, :]
        delta_backward = nextReg[cidx2][ptidx[:, 1], bwd_idx, :]

        dist1, abs_a1 = compute_features(delta_real_forward, delta_forward)
        dist2, abs_a2 = compute_features(delta_real_backward, delta_backward)

        featAugm = np.hstack((dist1.reshape(num_edges, 1), abs_a1.reshape(num_edges, 1), dist2.reshape(num_edges, 1),
                              abs_a2.reshape(num_edges, 1)))

        # print "features: ", cidx1, cidx2, featAugm

        featAugm = np.hstack((featAugm, np.exp(-featAugm), np.ones((num_edges, 1))))

        featAugm[:, :-1] = (featAugm[:, :-1] - self.X_min[cidx1][cidx2]) / (
            self.X_max[cidx1][cidx2] - self.X_min[cidx1][cidx2])
        cur_prob = 1 / (1 + np.exp(-featAugm.dot(self.w[cidx1][cidx2])))

        return cur_prob
