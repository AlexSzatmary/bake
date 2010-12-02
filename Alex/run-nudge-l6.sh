./batchrun.py bp/mix/uniaxi-planar-prolate-nudge-l6 -l
./batchrun.py bp/mix/uniaxi-planar-oblate-nudge-l6 -l
./batchrun.py bp/mix/biaxi-planar-prolate-nudge-l6 -l
./batchrun.py bp/mix/biaxi-planar-oblate-nudge-l6 -l
./batchrun.py bp/mix/uniaxi-planar-prolate-nudge-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.50d0" -o "\$phif\$;0.01d0" -l
./batchrun.py bp/mix/uniaxi-planar-oblate-nudge-l6 -o "\$r\$;1.0d0" -o -o "\$phif\$;0.49d0" -l
./batchrun.py bp/mix/uniaxi-planar-oblate-nudge-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.49d0" -o "\$phif\$;0.00d0" -l
