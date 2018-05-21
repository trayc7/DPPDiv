#!/bin/bash

## This executes a run under the FBD range model 1 for the estimation of diversification rates from stratigraphic ranges
./dppdiv-seq-sse -fbdrsky -cal ranges.1.cal -int ranges.1.int -rho 1 -out test_fbd_range.m1 -n 1000000 -sf 100

## This runs the same analysis but under the prior by returning a constant value for the likelihood
./dppdiv-seq-sse -fbdrsky -cal ranges.1.cal -int ranges.1.int -rho 1 -out test_fbd_range.pr -n 1000000 -sf 100 -rnp

## This executes a run under the FBD range model 2 for the estimation of diversification rates from stratigraphic ranges
./dppdiv-seq-sse -fbdrsky -cal ranges.1.cal -int ranges.1.int -rho 1 -out test_fbd_range.m2 -n 1000000 -sf 100 -fbdrlk 2

## This executes a run under the FBD range model 3 for the estimation of diversification rates from stratigraphic ranges
./dppdiv-seq-sse -fbdrsky -cal ranges.1.cal -int ranges.1.int -pa ranges.1.pa -rho 1 -out test_fbd_range.m3 -n 1000000 -sf 100 -fbdrlk 3
