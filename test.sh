#!/bin/sh

pip install .
cd ziggie/tests
python test.py
cd ..
