#!/bin/bash

python process_rnaseq.py

R -f normalize_data.r
