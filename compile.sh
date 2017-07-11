#!/bin/sh

cd lib

cd coco/PythonAPI/
python3 setup.py build_ext --inplace
cd -

cd multicut_cython
python3 setup.py build_ext --inplace
cd -

cd nms_cython
python3 setup.py build_ext --inplace
cd -

cd ..
