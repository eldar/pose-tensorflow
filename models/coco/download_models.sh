#!/bin/sh

curl -O http://datasets.d2.mpi-inf.mpg.de/deepercut-models-tensorflow/coco-resnet-101.data-00000-of-00001
curl -O http://datasets.d2.mpi-inf.mpg.de/deepercut-models-tensorflow/coco-resnet-101.meta
curl -O http://datasets.d2.mpi-inf.mpg.de/deepercut-models-tensorflow/coco-resnet-101.index

curl -O http://datasets.d2.mpi-inf.mpg.de/deepercut-models-tensorflow/pairwise_coco.tar.gz
tar xvzf pairwise_coco.tar.gz
