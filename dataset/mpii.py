from dataset.pose_dataset import PoseDataset


class MPII(PoseDataset):
    def __init__(self, cfg):
        cfg.all_joints = [[0, 5], [1, 4], [2, 3], [6, 11], [7, 10], [8, 9], [12], [13]]
        cfg.all_joints_names = ['ankle', 'knee', 'hip', 'wrist', 'elbow', 'shoulder', 'chin', 'forehead']
        cfg.num_joints = 14
        super().__init__(cfg)

    def mirror_joint_coords(self, joints, image_width):
        joints[:, 1] = image_width - joints[:, 1]
        return joints

    def get_pose_segments(self):
       return [[0, 1], [1, 2], [3, 4], [4, 5], [6, 7], [7, 8], [9, 10], [10, 11], [12, 13]]
