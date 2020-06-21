#!/bin/sh

pip install .
cd tests
python test.py
cd ..
