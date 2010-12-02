#!/bin/bash
# uniaxial-planar DF
./batchrun.py bp/mix/uniaxi-planar -p DF -o "\$capillary_no\$;1.d-2"
mv extractTaylorDF__00001.txt up-sp-DF-g1.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p DF -o "\$capillary_no\$;1.d-2"
mv extractTaylorDF__00001.txt up-pr-DF-g1.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p DF -o "\$capillary_no\$;1.d-2"
mv extractTaylorDF__00001.txt up-ob-DF-g1.d-2.txt

./batchrun.py bp/mix/uniaxi-planar -p DF -o "\$capillary_no\$;2.d-2"
mv extractTaylorDF__00001.txt up-sp-DF-g2.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p DF -o "\$capillary_no\$;2.d-2"
mv extractTaylorDF__00001.txt up-pr-DF-g2.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p DF -o "\$capillary_no\$;2.d-2"
mv extractTaylorDF__00001.txt up-ob-DF-g2.d-2.txt

./batchrun.py bp/mix/uniaxi-planar -p DF
mv extractTaylorDF__00001.txt up-sp-DF-g4.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p DF
mv extractTaylorDF__00001.txt up-pr-DF-g4.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p DF
mv extractTaylorDF__00001.txt up-ob-DF-g4.d-2.txt

# Biaxial-Planar DF
./batchrun.py bp/mix/biaxi-planar -p DF -o "\$capillary_no\$;1.d-2"
mv extractTaylorDF__00001.txt bp-sp-DF-g1.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p DF -o "\$capillary_no\$;1.d-2"
mv extractTaylorDF__00001.txt bp-pr-DF-g1.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p DF -o "\$capillary_no\$;1.d-2"
mv extractTaylorDF__00001.txt bp-ob-DF-g1.d-2.txt

./batchrun.py bp/mix/biaxi-planar -p DF -o "\$capillary_no\$;2.d-2"
mv extractTaylorDF__00001.txt bp-sp-DF-g2.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p DF -o "\$capillary_no\$;2.d-2"
mv extractTaylorDF__00001.txt bp-pr-DF-g2.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p DF -o "\$capillary_no\$;2.d-2"
mv extractTaylorDF__00001.txt bp-ob-DF-g2.d-2.txt

./batchrun.py bp/mix/biaxi-planar -p DF
mv extractTaylorDF__00001.txt bp-sp-DF-g4.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p DF
mv extractTaylorDF__00001.txt bp-pr-DF-g4.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p DF
mv extractTaylorDF__00001.txt bp-ob-DF-g4.d-2.txt

# uniaxial-planar tau
rm -f fit.log 
./batchrun.py bp/mix/uniaxi-planar -f -o "\$capillary_no\$;1.d-2"
grep '^tau* *=.* .*/.*' fit.log > up-sp-tau-g1.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -f -o "\$capillary_no\$;1.d-2"
grep '^tau* *=.* .*/.*' fit.log > up-pr-tau-g1.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -f -o "\$capillary_no\$;1.d-2"
grep '^tau* *=.* .*/.*' fit.log > up-ob-tau-g1.d-2.txt

rm -f fit.log 
./batchrun.py bp/mix/uniaxi-planar -f -o "\$capillary_no\$;2.d-2"
grep '^tau* *=.* .*/.*' fit.log > up-sp-tau-g2.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -f -o "\$capillary_no\$;2.d-2"
grep '^tau* *=.* .*/.*' fit.log > up-pr-tau-g2.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -f -o "\$capillary_no\$;2.d-2"
grep '^tau* *=.* .*/.*' fit.log > up-ob-tau-g2.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar -f
grep '^tau* *=.* .*/.*' fit.log > up-sp-tau-g4.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > up-pr-tau-g4.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > up-ob-tau-g4.d-2.txt

# biaxial-planar tau
rm -f fit.log 
./batchrun.py bp/mix/biaxi-planar -f -o "\$capillary_no\$;1.d-2"
grep '^tau* *=.* .*/.*' fit.log > bp-sp-tau-g1.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -f -o "\$capillary_no\$;1.d-2"
grep '^tau* *=.* .*/.*' fit.log > bp-pr-tau-g1.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -f -o "\$capillary_no\$;1.d-2"
grep '^tau* *=.* .*/.*' fit.log > bp-ob-tau-g1.d-2.txt

rm -f fit.log 
./batchrun.py bp/mix/biaxi-planar -f -o "\$capillary_no\$;2.d-2"
grep '^tau* *=.* .*/.*' fit.log > bp-sp-tau-g2.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -f -o "\$capillary_no\$;2.d-2"
grep '^tau* *=.* .*/.*' fit.log > bp-pr-tau-g2.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -f -o "\$capillary_no\$;2.d-2"
grep '^tau* *=.* .*/.*' fit.log > bp-ob-tau-g2.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar -f
grep '^tau* *=.* .*/.*' fit.log > bp-sp-tau-g4.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > bp-pr-tau-g4.d-2.txt

rm -f fit.log
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -f
grep '^tau* *=.* .*/.*' fit.log > bp-ob-tau-g4.d-2.txt

# uniaxial-planar strain energy
./batchrun.py bp/mix/uniaxi-planar -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;1.d-2" > up-sp-W-g1.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;1.d-2" > up-pr-W-g1.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;1.d-2" > up-ob-W-g1.d-2.txt

./batchrun.py bp/mix/uniaxi-planar -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;2.d-2" > up-sp-W-g2.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;2.d-2" > up-pr-W-g2.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;2.d-2" > up-ob-W-g2.d-2.txt

./batchrun.py bp/mix/uniaxi-planar -E '../../Alex/strainenergy.py' > up-sp-W-g4.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' > up-pr-W-g4.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' > up-ob-W-g4.d-2.txt

# biaxial-planar strain energy
./batchrun.py bp/mix/biaxi-planar -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;1.d-2" > bp-sp-W-g1.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;1.d-2" > bp-pr-W-g1.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;1.d-2" > bp-ob-W-g1.d-2.txt

./batchrun.py bp/mix/biaxi-planar -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;2.d-2" > bp-sp-W-g2.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;2.d-2" > bp-pr-W-g2.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' -o "\$capillary_no\$;2.d-2" > bp-ob-W-g2.d-2.txt

./batchrun.py bp/mix/biaxi-planar -E '../../Alex/strainenergy.py' > bp-sp-W-g4.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -E '../../Alex/strainenergy.py' > bp-pr-W-g4.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -E '../../Alex/strainenergy.py' > bp-ob-W-g4.d-2.txt

# uniaxial-planar surface area
./batchrun.py bp/mix/uniaxi-planar -p volumearea00001.txt -o "\$capillary_no\$;1.d-2"
mv extractvolumearea00001.txt up-sp-area-g1.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p volumearea00001.txt -o "\$capillary_no\$;1.d-2"
mv extractvolumearea00001.txt up-pr-area-g1.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p volumearea00001.txt -o "\$capillary_no\$;1.d-2"
mv extractvolumearea00001.txt up-ob-area-g1.d-2.txt

./batchrun.py bp/mix/uniaxi-planar -p volumearea00001.txt -o "\$capillary_no\$;2.d-2"
mv extractvolumearea00001.txt up-sp-area-g2.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p volumearea00001.txt -o "\$capillary_no\$;2.d-2"
mv extractvolumearea00001.txt up-pr-area-g2.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p volumearea00001.txt -o "\$capillary_no\$;2.d-2"
mv extractvolumearea00001.txt up-ob-area-g2.d-2.txt

./batchrun.py bp/mix/uniaxi-planar -p volumearea00001.txt 
mv extractvolumearea00001.txt up-sp-area-g4.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-prolate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt up-pr-area-g4.d-2.txt
./batchrun.py bp/mix/uniaxi-planar-oblate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt up-ob-area-g4.d-2.txt

# biaxial-planar surface area
./batchrun.py bp/mix/biaxi-planar -p volumearea00001.txt -o "\$capillary_no\$;1.d-2"
mv extractvolumearea00001.txt bp-sp-area-g1.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p volumearea00001.txt -o "\$capillary_no\$;1.d-2"
mv extractvolumearea00001.txt bp-pr-area-g1.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p volumearea00001.txt -o "\$capillary_no\$;1.d-2"
mv extractvolumearea00001.txt bp-ob-area-g1.d-2.txt

./batchrun.py bp/mix/biaxi-planar -p volumearea00001.txt -o "\$capillary_no\$;2.d-2"
mv extractvolumearea00001.txt bp-sp-area-g2.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p volumearea00001.txt -o "\$capillary_no\$;2.d-2"
mv extractvolumearea00001.txt bp-pr-area-g2.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p volumearea00001.txt -o "\$capillary_no\$;2.d-2"
mv extractvolumearea00001.txt bp-ob-area-g2.d-2.txt

./batchrun.py bp/mix/biaxi-planar -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-sp-area-g4.d-2.txt
./batchrun.py bp/mix/biaxi-planar-prolate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-pr-area-g4.d-2.txt
./batchrun.py bp/mix/biaxi-planar-oblate-varyr -p volumearea00001.txt 
mv extractvolumearea00001.txt bp-ob-area-g4.d-2.txt
