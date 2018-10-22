#!/bin/bash

rm -rf api-gen latest
make clean
sphinx-apidoc -fMeET -o api-gen ../sp2graph ../sp2graph/**/setup.py
make html
rm -r build/html/_sources
mv build/html latest
