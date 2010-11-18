./batchrun.py bp/mix/uniaxi-planar-prolate-l6 -R tara
./batchrun.py bp/mix/uniaxi-planar-oblate-l6 -R tara
./batchrun.py bp/mix/biaxi-planar-prolate-l6 -R tara
./batchrun.py bp/mix/biaxi-planar-oblate-l6 -R tara
./batchrun.py bp/mix/uniaxi-planar-prolate-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.50d0" -o "\$phif\$;0.00d0" -R tara
./batchrun.py bp/mix/uniaxi-planar-prolate-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0" -o "\$phif\$;0.25d0;0.50d0" -R tara
./batchrun.py bp/mix/uniaxi-planar-oblate-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.50d0" -o "\$phif\$;0.00d0" -R tara
./batchrun.py bp/mix/uniaxi-planar-oblate-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0" -o "\$phif\$;0.25d0;0.50d0" -R tara

#./batchrun.py bp/mix/uniaxi-planar-prolate-l6 -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.50d0" -R tara
#./batchrun.py bp/mix/uniaxi-planar-oblate-l6 -o "\$r\$;1.0d0" -R tara
