#!/bin/bash/

python MD25box-write-positionFiles.py
python MD04box-write-positionFiles.py  
python MD10box-write-positionFiles.py

python MD04box-compute2PCF-vbins.py
python MD10box-compute2PCF-vbins.py
python MD25box-compute2PCF-vbins.py

