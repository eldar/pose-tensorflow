import numpy as np

import tensorflow as tf

from nnet.net_factory import pose_net


def setup_pose_prediction(cfg):
    inputs = tf.placeholder(tf.float32, shape=[cfg.batch_size, None, None, 3])

    outputs = pose_net(cfg).test(inputs)

    restorer = tf.train.Saver()

    sess = tf.Session()

    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())

    # Restore variables from disk.
    restorer.restore(sess, cfg.init_weights)

    return sess, inputs, outputs


def extract_cnn_output(outputs_np, cfg, pairwise_stats = None):
    scmap = outputs_np['part_prob']
    scmap = np.squeeze(scmap)
    locref = None
    pairwise_diff = None
    if cfg.location_refinement:
        locref = np.squeeze(outputs_np['locref'])
        shape = locref.shape
        locref = np.reshape(locref, (shape[0], shape[1], -1, 2))
        locref *= cfg.locref_stdev
    if cfg.pairwise_predict:
        pairwise_diff = np.squeeze(outputs_np['pairwise_pred'])
        shape = pairwise_diff.shape
        pairwise_diff = np.reshape(pairwise_diff, (shape[0], shape[1], -1, 2))
        num_joints = cfg.num_joints
        for pair in pairwise_stats:
            pair_id = (num_joints - 1) * pair[0] + pair[1] - int(pair[0] < pair[1])
            pairwise_diff[:, :, pair_id, 0] *= pairwise_stats[pair]["std"][0]
            pairwise_diff[:, :, pair_id, 0] += pairwise_stats[pair]["mean"][0]
            pairwise_diff[:, :, pair_id, 1] *= pairwise_stats[pair]["std"][1]
            pairwise_diff[:, :, pair_id, 1] += pairwise_stats[pair]["mean"][1]
    return scmap, locref, pairwise_diff


def argmax_pose_predict(scmap, offmat, stride):
    """Combine scoremat and offsets to the final pose."""
    num_joints = scmap.shape[2]
    pose = []
    for joint_idx in range(num_joints):
        maxloc = np.unravel_index(np.argmax(scmap[:, :, joint_idx]),
                                  scmap[:, :, joint_idx].shape)
        offset = np.array(offmat[maxloc][joint_idx])[::-1] if offmat is not None else 0
        pos_f8 = (np.array(maxloc).astype('float') * stride + 0.5 * stride +
                  offset)
        pose.append(np.hstack((pos_f8[::-1],
                               [scmap[maxloc][joint_idx]])))
    return np.array(pose)


def argmax_arrows_predict(scmap, offmat, pairwise_diff, stride):
    num_joints = scmap.shape[2]
    arrows = {}
    for joint_idx in range(num_joints):
        maxloc = np.unravel_index(np.argmax(scmap[:, :, joint_idx]),
                                  scmap[:, :, joint_idx].shape)
        offset = np.array(offmat[maxloc][joint_idx])[::-1] if offmat is not None else 0
        pos_f8 = (np.array(maxloc).astype('float') * stride + 0.5 * stride +
                  offset)[::-1]
        for joint_idx_end in range(num_joints):
            if joint_idx_end != joint_idx:
                pair_id = (num_joints - 1) * joint_idx + joint_idx_end - int(joint_idx < joint_idx_end)
                difference = np.array(pairwise_diff[maxloc][pair_id])[::-1] if pairwise_diff is not None else 0
                pos_f8_end = (np.array(maxloc).astype('float') * stride + 0.5 * stride + difference)[::-1]
                arrows[(joint_idx, joint_idx_end)] = (pos_f8, pos_f8_end)

    return arrows