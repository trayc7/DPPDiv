#!/bin/bash

## This executes a run under a strict clock using the FBD model for fossil calibration
./dppdiv-seq-sse -in test_seq.dat -tre test_tre.phy -cal test_fos.cal -tga -clok -n 1000000 -sf 100 -pf 1000 -out test_fbd 

## This runs the same analysis but under the prior by returning a constant value for the likelihood
./dppdiv-seq-sse -in test_seq.dat -tre test_tre.phy -cal test_fos.cal -tga -clok -n 1000000 -sf 100 -pf 1000 -out test_fbd.pr -rnp
