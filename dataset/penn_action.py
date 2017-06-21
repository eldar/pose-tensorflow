from functools import reduce
import numpy as np

from dataset.pose_dataset import PoseDataset, Batch


def merge_batch(batches):
    """
    Merges n=len(batches) batches of size 1 into
    one batch of size n
    """
    res = {}
    for key, tensor in batches[0].items():
        elements = [batch[key] for batch in batches]
        if type(tensor) is np.ndarray:
            elements = reduce(lambda x, y: np.concatenate((x, y), axis=0), elements)
        res[key] = elements
    return res


class PennAction(PoseDataset):
    def __init__(self, cfg):
        cfg.all_joints = [[0], [1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12]]
        cfg.all_joints_names = ["head", "shoulder", "elbow", "wrist", "hip", "knee", "ankle"]
        cfg.num_joints = 13
        super().__init__(cfg)
        self.add_extra_fields()

    def add_extra_fields(self):
        dataset = self.raw_data['dataset']
        for i in range(self.num_images):
            raw_item = dataset[0, i]
            item = self.data[i]
            item.seq_id = raw_item[4][0][0]
            item.frame_id = raw_item[5][0][0]

    def mirror_joint_coords(self, joints, image_width):
        joints[:, 1] = image_width - joints[:, 1] + 1  # 1-indexed
        return joints

    def next_batch(self):
        while True:
            imidx, mirror = self.next_training_sample()
            data_item = self.get_training_sample(imidx)

            scale = self.get_scale()
            if not self.is_valid_size(data_item.im_size, scale):
                continue

            if self.cfg.video_batch:
                sequences = self.raw_data['sequences']
                seq_ids = sequences[0, data_item.seq_id][0]
                num_frames = len(seq_ids)
                start_frame = data_item.frame_id
                num_frames_model = self.cfg.batch_size

                if start_frame + num_frames_model - 1 >= num_frames:
                    start_frame = num_frames - num_frames_model

                seq_subset = seq_ids[start_frame:start_frame+num_frames_model]
                data_items = [self.get_training_sample(imidx) for imidx in seq_subset]
                batches = [self.make_batch(item, scale, mirror) for item in data_items]

                batch = merge_batch(batches)
            else:
                batch = self.make_batch(data_item, scale, mirror)

            return batch
