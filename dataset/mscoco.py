import os
import sys
import skimage.io as io
import matplotlib.pyplot as plt

import numpy as np
from scipy.misc import imresize
import json

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + '/../lib/coco/PythonAPI')

from pycocotools.coco import COCO as COCO
from pycocotools import mask as maskUtils

from dataset.pose_dataset import PoseDataset, DataItem


def get_gt_visibilities(inFile, visibilities):
    with open(inFile) as data_file:
        data = json.load(data_file)

    for person_id in range(len(data)):
        keypoints = data[person_id]["keypoints"]
        keypoints = [visibilities[person_id][i//3] if i % 3 == 2 else int(keypoints[i]) for i in
                     range(len(keypoints))]
        data[person_id]["keypoints"] = keypoints

    with open(inFile, 'w') as data_file:
        json.dump(data, data_file)


class MSCOCO(PoseDataset):
    def __init__(self, cfg):
        cfg.all_joints = [[0], [2, 1], [4, 3], [6, 5], [8, 7],[10, 9], [12, 11], [14, 13], [16, 15]]
        cfg.all_joints_names = ["nose", 'eye', 'ear', 'shoulder', 'elbow', 'hand', 'hip', 'knee', 'foot']
        cfg.num_joints = 17
        super().__init__(cfg)

    def load_dataset(self):
        dataset  = self.cfg.dataset
        dataset_phase = self.cfg.dataset_phase
        dataset_ann = self.cfg.dataset_ann

        # initialize COCO api
        annFile = '%s/annotations/%s_%s.json'%(dataset,dataset_ann,dataset_phase)
        self.coco = COCO(annFile)

        imgIds = self.coco.getImgIds()

        data = []

        # loop through each image
        for imgId in imgIds:
            item = DataItem()

            img = self.coco.loadImgs(imgId)[0]
            item.im_path = "%s/images/%s/%s"%(dataset, dataset_phase, img["file_name"])
            item.im_size = [3, img["height"], img["width"]]
            item.coco_id = imgId
            annIds = self.coco.getAnnIds(imgIds=img['id'], iscrowd=False)
            anns = self.coco.loadAnns(annIds)

            all_person_keypoints = []
            masked_persons_RLE = []
            visible_persons_RLE = []
            all_visibilities = []

            # Consider only images with people
            has_people = len(anns) > 0
            if not has_people and self.cfg.coco_only_images_with_people:
                continue

            for ann in anns: # loop through each person
                person_keypoints = []
                visibilities = []
                if ann["num_keypoints"] != 0:
                    for i in range(self.cfg.num_joints):
                        x_coord = ann["keypoints"][3 * i]
                        y_coord = ann["keypoints"][3 * i + 1]
                        visibility = ann["keypoints"][3 * i + 2]
                        visibilities.append(visibility)
                        if visibility != 0: # i.e. if labeled
                            person_keypoints.append([i, x_coord, y_coord])
                    all_person_keypoints.append(np.array(person_keypoints))
                    visible_persons_RLE.append(maskUtils.decode(self.coco.annToRLE(ann)))
                    all_visibilities.append(visibilities)
                if ann["num_keypoints"] == 0:
                    masked_persons_RLE.append(self.coco.annToRLE(ann))

            item.joints = all_person_keypoints
            item.im_neg_mask = maskUtils.merge(masked_persons_RLE)
            if self.cfg.use_gt_segm:
                item.gt_segm = np.moveaxis(np.array(visible_persons_RLE), 0, -1)
                item.visibilities = all_visibilities
            data.append(item)

        self.has_gt = self.cfg.dataset is not "image_info"
        return data


    def compute_scmap_weights(self, scmap_shape, joint_id, data_item):
        size = scmap_shape[0:2]
        scmask = np.ones(size)
        m = maskUtils.decode(data_item.im_neg_mask)
        if m.size:
            scmask = 1.0 - imresize(m, size)
        scmask = np.stack([scmask] * self.cfg.num_joints, axis=-1)
        return scmask

    def get_pose_segments(self):
       return [[0, 1], [0, 2], [1, 3], [2, 4], [5, 7], [6, 8], [7, 9], [8, 10], [11, 13], [12, 14], [13, 15], [14, 16]]

    def visualize_coco(self, coco_img_results, visibilities):
        inFile = "tmp.json"
        with open(inFile, 'w') as outfile:
            json.dump(coco_img_results, outfile)
        get_gt_visibilities(inFile, visibilities)

        # initialize cocoPred api
        cocoPred = self.coco.loadRes(inFile)
        os.remove(inFile)

        imgIds = [coco_img_results[0]["image_id"]]

        for imgId in imgIds:
            img = cocoPred.loadImgs(imgId)[0]
            im_path = "%s/images/%s/%s" % (self.cfg.dataset, self.cfg.dataset_phase, img["file_name"])
            I = io.imread(im_path)

            fig = plt.figure()
            a = fig.add_subplot(2, 2, 1)
            plt.imshow(I)
            a.set_title('Initial Image')

            a = fig.add_subplot(2, 2, 2)
            plt.imshow(I)
            a.set_title('Predicted Keypoints')
            annIds = cocoPred.getAnnIds(imgIds=img['id'])
            anns = cocoPred.loadAnns(annIds)
            cocoPred.showAnns(anns)

            a = fig.add_subplot(2, 2, 3)
            plt.imshow(I)
            a.set_title('GT Keypoints')
            annIds = self.coco.getAnnIds(imgIds=img['id'])
            anns = self.coco.loadAnns(annIds)
            self.coco.showAnns(anns)

            plt.show()
