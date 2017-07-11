Before training the models you need to download the ImageNet pre-trained ResNet weights

```
$ cd models/pretrained
$ ./download.sh
```

Training parameters are specified in the `pose_cfg.yaml` file.

Here are the dataset specific instructions. 

## Training a model with MPII Pose Dataset (Single Person)


1. Download the dataset from [this page](http://human-pose.mpi-inf.mpg.de/),
both images and annotations. Unpack it to the path `<path_to_dataset>`
to have the following directory structure:

```
<path_to_dataset>/images/*.jpg
<path_to_dataset>/mpii_human_pose_v1_u12_1.mat
```

2. Preprocess dataset (crop and rescale)

```
$ cd matlab/mpii
$ matlab -nodisplay -nosplash

# in matlab execute the following function, be sure to specify the *absolute* path
preprocess_single('<path_to_dataset>')
```

3. Edit the training definition file
`models/mpii/train/pose_cfg.yaml` such that:

```
dataset: `<path_to_dataset>/cropped/dataset.mat`
```

4. Train the model

```
$ cd models/mpii/train/
$ TF_CUDNN_USE_AUTOTUNE=0 CUDA_VISIBLE_DEVICES=0 python3 ../../../train.py
```

## Training a model with MS COCO dataset (Multi-Person)

1. Download [MS COCO](http://mscoco.org/dataset/#download)
train2014 set with keypoint and object instances annotations.

2. Download pairwise statistics:
```
$ cd models/coco
$ ./download_models.sh
```

3. Edit the training definition file
`models/coco/train/pose_cfg.yaml` such that:

```
dataset: `<path_to_mscoco>`
```

4. Train the model:

```
$ cd models/coco/train/
$ TF_CUDNN_USE_AUTOTUNE=0 CUDA_VISIBLE_DEVICES=0 python3 ../../../train.py
```

## Training on your own dataset

If you wish to train keypoint detectors on your own data, first of all
you have to prepare dataset definition file `dataset.mat` (Matlab data format).
We included  an example of such file [here](dataset_example.mat). It must
contain a variable `dataset` that is a struct array, with each entry
corresponding to an image and with the following fields:

1. `image` - path to the image
2. `size` - 1x3 array containing [num_channels, image_height, image_width]. Set num_channels=3 for an RGB image.
3. `joints` - a cell array of nx3 joint annotations, for example:

```matlab
joints = {[ ...
  0,  175,  261; ... % 0-indexed joint ID, X coordinate, Y coordinate
  1,  173,  178; ...
  2,  144,  122; ...
  3,  193,  124; ...
]};
```

You will potentially need to adjust the training definition file `pose_cfg.yaml`.
Here we will describe some of its options:

```yaml
# path to the dataset description file
dataset: /path/to/dataset.mat

# all locations within this distance threshold are considered
# positive training samples for detector
pos_dist_thresh: 17 

# all images in the dataset will be rescaled by the following
# scaling factor to be processed by the CNN. You can select the
# optimal scale by cross-validation
global_scale: 0.80

# During training an image will be randomly scaled within the
# range [scale_jitter_lo; scale_jitter_up] to augment training data,
# We found +/- 15% scale jitter works quite well.
scale_jitter_lo: 0.85
scale_jitter_up: 1.15

# Randomly flips an image horizontally to augment training data
mirror: true

# list of pairs of symmetric joint IDs, for example in this case
# 0 and 5 are IDs for the symmetric parts, and 12 or 13 do not have
# symmetric parts. This is used to do flip training data correctly. 
all_joints = [[0, 5], [1, 4], [2, 3], [6, 11], [7, 10], [8, 9], [12], [13]]

# Type of the CNN to use, currently resnet_101 and resnet_50
# are supported
net_type: resnet_101
init_weights: ../../pretrained/resnet_v1_101.ckpt

# Location refinement parameters (check https://arxiv.org/abs/1511.06645)
location_refinement: true
locref_huber_loss: true
locref_loss_weight: 0.05
locref_stdev: 7.2801

# Enabling this adds additional loss layer in the middle of the ConvNet,
# which helps accuracy.
intermediate_supervision: true
intermediate_supervision_layer: 12

# all images larger with size
# width * height > max_input_size*max_input_size are not used in training.
# Prevents training from crashing with out of memory exception for very
# large images.
max_input_size: 850

# Learning rate schedule for the SGD optimiser. 
multi_step:
- [0.005, 10000]
- [0.02, 430000]
- [0.002, 730000]
- [0.001, 1030000]

# How often display loss
display_iters: 20

# How often to save training snapshot
save_iters: 60000
```

You don't have to crop images such that they all have the same size,
as training is done with `batch_size=1` (batch size larger than 1 is
currently not supported anyway).