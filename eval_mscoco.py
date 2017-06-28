import os
import sys

import argparse
import json

from config import load_config

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path + '/lib/coco/PythonAPI')
from pycocotools.coco import COCO as COCO
from pycocotools.cocoeval import COCOeval

def apply_threhsold(inFile, threshold):
    outFile = inFile[:-5] + '-' + str(threshold) + '.json'

    with open(inFile) as data_file:
        data = json.load(data_file)

    for person_id in range(len(data)):
        keypoints = data[person_id]["keypoints"]
        keypoints = [int(keypoints[i] > threshold) if i % 3 == 2 else int(keypoints[i]) for i in range(len(keypoints))]
        data[person_id]["keypoints"] = keypoints

    with open(outFile, 'w') as outfile:
        json.dump(data, outfile)

    return outFile


def eval_init(cfg):
    dataset = cfg.dataset
    dataset_phase = cfg.dataset_phase
    dataset_ann = cfg.dataset_ann
    threshold = 0

    # initialize cocoGT api
    annFile = '%s/annotations/%s_%s.json' % (dataset, dataset_ann, dataset_phase)
    cocoGT = COCO(annFile)

    # initialize cocoPred api
    inFile = "predictions_with_segm.json"
    predFile = apply_threhsold(inFile, threshold)
    cocoPred = cocoGT.loadRes(predFile)

    return cocoGT, cocoPred


def eval_mscoco_with_segm(cocoGT, cocoPred):
    # running evaluation
    cocoEval = COCOeval(cocoGT, cocoPred, "keypoints")
    cocoEval.evaluate()
    cocoEval.accumulate()
    cocoEval.summarize()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args, unparsed = parser.parse_known_args()

    cfg = load_config()

    cocoGT, cocoPred = eval_init(cfg)
    eval_mscoco_with_segm(cocoGT, cocoPred)