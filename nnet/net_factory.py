from . import pose_net

def pose_net(cfg):
    cls = pose_net.PoseNet
    return cls(cfg)
