./batchrun.py bp/mix/uniaxi-planar-prolate -e nvec1_0001.txt
./batchrun.py bp/mix/uniaxi-planar-oblate -e nvec1_0001.txt
./batchrun.py bp/mix/biaxi-planar-prolate -e nvec1_0001.txt
./batchrun.py bp/mix/biaxi-planar-oblate -e nvec1_0001.txt
./batchrun.py bp/mix/uniaxi-planar-prolate -o "\$r\$;1.d0" -o "\$thetaf\$;0.d0;0.25d0;0.50d0" -o "\$phif\$;0.d0" -e nvec1_0001.txt
./batchrun.py bp/mix/uniaxi-planar-prolate -o "\$r\$;1.d0" -o "\$thetaf\$;0.d0" -o "\$phif\$;0.25d0;0.50d0" -e nvec1_0001.txt
./batchrun.py bp/mix/uniaxi-planar-oblate -o "\$r\$;1.d0" -o "\$thetaf\$;0.d0;0.25d0;0.5d0" -o "\$phif\$;0.d0" -e nvec1_0001.txt
./batchrun.py bp/mix/uniaxi-planar-oblate -o "\$r\$;1.d0" -o "\$thetaf\$;0.d0" -o "\$phif\$;0.25d0;0.50d0" -e nvec1_0001.txt
