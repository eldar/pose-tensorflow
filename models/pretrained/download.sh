#!/bin/sh

if command -v curl 1>/dev/null; then
    DOWNLOADER="curl"
else
    DOWNLOADER="wget -qO-"
fi

$DOWNLOADER http://download.tensorflow.org/models/resnet_v1_50_2016_08_28.tar.gz | tar xvz
$DOWNLOADER http://download.tensorflow.org/models/resnet_v1_101_2016_08_28.tar.gz | tar xvz
