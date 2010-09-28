#!/bin/bash
./batchrun.py bp/mix/uniaxi-planar-l6 -p DF
mv extractTaylorDF__00001.txt up-sp-DF.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr-l6 -p DF
mv extractTaylorDF__00001.txt up-pr-DF.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr-l6 -p DF
mv extractTaylorDF__00001.txt up-ob-DF.txt

./batchrun.py bp/mix/biaxi-planar-l6 -p DF
mv extractTaylorDF__00001.txt bp-sp-DF.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr-l6 -p DF
mv extractTaylorDF__00001.txt bp-pr-DF.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr-l6 -p DF
mv extractTaylorDF__00001.txt bp-ob-DF.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-l6 -f
grep '^tau* *=.* .*/.*' fit.log > up-sp-tau.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr-l6 -f
grep '^tau* *=.* .*/.*' fit.log > up-pr-tau.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr-l6 -f
grep '^tau* *=.* .*/.*' fit.log > up-ob-tau.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-l6 -f
grep '^tau* *=.* .*/.*' fit.log > bp-sp-tau.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-prolate-varyr-l6 -f
grep '^tau* *=.* .*/.*' fit.log > bp-pr-tau.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-oblate-varyr-l6 -f
grep '^tau* *=.* .*/.*' fit.log > bp-ob-tau.txt

./batchrun.py bp/mix/uniaxi-planar-l6 -p volumearea00001.txt 
mv extractvolumearea00001.txt up-sp-area.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr-l6 -p volumearea00001.txt 
mv extractvolumearea00001.txt up-pr-area.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr-l6 -p volumearea00001.txt 
mv extractvolumearea00001.txt up-ob-area.txt

./batchrun.py bp/mix/biaxi-planar-l6 -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-sp-area.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr-l6 -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-pr-area.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr-l6 -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-ob-area.txt
