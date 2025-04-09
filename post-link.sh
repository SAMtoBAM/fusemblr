#!/bin/bash
sed -i 's/v, 1000, 10000)/v, 1000, 200000)/g' $CONDA_PREFIX/lib/python3.*/site-packages/flye/main.py