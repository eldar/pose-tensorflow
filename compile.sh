#!/bin/sh

cd lib

cd coco/PythonAPI/
python setup.py build_ext --inplace
cd -

cd multicut_cython
python setup.py build_ext --inplace
cd -

cd nms_cython
python setup.py build_ext --inplace
cd -

cd ..
