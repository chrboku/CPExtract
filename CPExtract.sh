#!/bin/sh
export LD_LIBRARY_PATH="/opt/R-4.1-patched/lib64/R/lib:$LD_LIBRARY_PATH"
python2 ./CPExtract.py
