#!/bin/bash
./batchrun.py bp/mix/uniaxi-planar -p DF
mv extractTaylorDF__00001.txt up-sp-DF.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p DF
mv extractTaylorDF__00001.txt up-pr-DF.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p DF
mv extractTaylorDF__00001.txt up-ob-DF.txt

./batchrun.py bp/mix/biaxi-planar -p DF
mv extractTaylorDF__00001.txt bp-sp-DF.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p DF
mv extractTaylorDF__00001.txt bp-pr-DF.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p DF
mv extractTaylorDF__00001.txt bp-ob-DF.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar -f
grep '^tau* *=.* .*/.*' fit.log > up-sp-tau.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > up-pr-tau.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > up-ob-tau.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar -f
grep '^tau* *=.* .*/.*' fit.log > bp-sp-tau.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > bp-pr-tau.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > bp-ob-tau.txt

./batchrun.py bp/mix/uniaxi-planar -E '../../Alex/strainenergy.py' > up-sp-W.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' > up-pr-W.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' > up-ob-W.txt

./batchrun.py bp/mix/biaxi-planar -E '../../Alex/strainenergy.py' > bp-sp-W.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' > bp-pr-W.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' > bp-ob-W.txt

./batchrun.py bp/mix/uniaxi-planar -p volumearea00001.txt 
mv extractvolumearea00001.txt up-sp-area.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt up-pr-area.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt up-ob-area.txt

./batchrun.py bp/mix/biaxi-planar -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-sp-area.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-pr-area.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-ob-area.txt
