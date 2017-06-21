#!/bin/sh

cd lib

cd coco/PythonAPI/
python_tf setup.py build_ext --inplace
cd -

cd multicut_cython
python_tf setup.py build_ext --inplace
cd -

cd nms_cython
python_tf setup.py build_ext --inplace
cd -

cd ..
