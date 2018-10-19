#!/bin/bash

rm -rf api-gen git/docs
make clean
sphinx-apidoc -fMeET -o api-gen ../sp2graph ../sp2graph/**/setup.py
make html
rm -r build/html/_sources
mv build/html git/docs
