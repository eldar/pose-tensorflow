import numpy as np
import scipy.io


def multi_dim_argmax(arr):
    """ Getting argmax over the first 2 axes """
    m, n, i, j = arr.shape
    arr = np.reshape(arr, (m * n, i, j))
    return np.unravel_index(np.argmax(arr, 0), (m, n))


def interweave_matrices(x, y, z):
    """Combine matrices by concatenating their cols: x.col(1),y.col(1),z.col(1) ... x.col(n),y.col(n),z.col(n) """
    num_joints = x.shape[1]
    id_x = (np.arange(0, num_joints, 0.5) + 1).astype('int')
    id_y = (np.arange(0, num_joints, 0.5) + 0.5).astype('int')
    id_z = (np.arange(0, num_joints, 0.5) + 0).astype('int')
    x = np.insert(x, id_x, 0, axis=1)
    y = np.insert(y, id_y, 0, axis=1)
    z = np.insert(z, id_z, 0, axis=1)
    return x + y + z


def pose_predict_with_gt_segm(scmap, offmat, stride, gt_segm, coco_id):
    """Generate all poses in an image using segmentations"""
    if gt_segm.size == 0:
        img_keypoints = []
    else:
        num_persons = gt_segm.shape[2]
        num_joints = scmap.shape[2]
        init_w = gt_segm.shape[1]
        init_h = gt_segm.shape[0]

        upsized_w = scmap.shape[1] * stride
        upsized_h = scmap.shape[0] * stride

        diff_w_l = int((upsized_w - init_w) // 2)
        diff_w_r = int(upsized_w - init_w - diff_w_l)

        diff_h_u = int((upsized_h - init_h) // 2)
        diff_h_d = int(upsized_h - init_h - diff_h_u)

        padded_gt_segm = np.pad(gt_segm, ((diff_h_u, diff_h_d), (diff_w_l, diff_w_r), (0, 0)), 'constant')
        upsized_scmap = scipy.ndimage.zoom(scmap, (stride, stride, 1))

        person_joint_scmap = padded_gt_segm[:, :, :, np.newaxis] * upsized_scmap[:, :, np.newaxis, :]
        upsized_maxloc = multi_dim_argmax(person_joint_scmap)
        maxloc_0 = (upsized_maxloc[0]//stride).astype(int)
        maxloc_1 = (upsized_maxloc[1]//stride).astype(int)
        indices = np.array([np.arange(num_joints)] * num_persons)
        offset = np.moveaxis(offmat[(maxloc_0, maxloc_1, indices)][:,:,::-1], -1, 0) if offmat is not None else 0
        pos_f8 = (np.array((maxloc_0, maxloc_1)).astype('float') * stride + 0.5 * stride + offset)
        v = scmap[(maxloc_0, maxloc_1, indices)]
        img_keypoints = (interweave_matrices(pos_f8[1].astype('int'), pos_f8[0].astype('int'), v)).tolist()

    coco_img_results = []
    for person_keypoints in img_keypoints:
        person_result = {}
        person_result["image_id"] = coco_id
        person_result["category_id"] = 1
        person_result["keypoints"] = person_keypoints
        person_result["score"] = 1
        coco_img_results.append(person_result)

    return coco_img_results

