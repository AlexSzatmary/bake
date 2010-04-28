#!/bin/bash
./batchrun.py bp/mix/uniaxi-planar-prolate -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.50d0" -o "\$phif\$;0.00d0" -c "nvec* mv* Tay* caps* minmax* volume* lambda* mean*"
./batchrun.py bp/mix/uniaxi-planar-prolate -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0" -o "\$phif\$;0.25d0;0.50d0" -c "nvec* mv* Tay* caps* minmax* volume* lambda* mean*"
./batchrun.py bp/mix/uniaxi-planar-oblate -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0;0.25d0;0.50d0" -o "\$phif\$;0.00d0" -c "nvec* mv* Tay* caps* minmax* volume* lambda* mean*"
./batchrun.py bp/mix/uniaxi-planar-oblate -o "\$r\$;1.0d0" -o "\$thetaf\$;0.00d0" -o "\$phif\$;0.25d0;0.50d0" -c "nvec* mv* Tay* caps* minmax* volume* lambda* mean*"
