from nnet.pose_net import PoseNet


def pose_net(cfg):
    if cfg.video:
        from nnet.pose_seq_net import PoseSeqNet
        cls = PoseSeqNet
    else:
        cls = PoseNet
    return cls(cfg)
