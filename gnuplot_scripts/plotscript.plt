set terminal png
set yr[0:]
set output planar-axiall6b6.d0g1.d-3DF.png
plot "planar-axiall6b6.d0DF.txt" using 2:5 every ::33::43, "analytical_data/BBR_axial_mix.txt" using 1:3 every ::33::43

set output planar-axiall6b6.d0g1.d-2DF.png
plot "planar-axiall6b6.d0.txt" using 2:5 every ::66::76, "analytical_data/BBR_axial_mix.txt" using 1:3 every ::66::76

set output planar-axiall6b6.d0g1.d-1DF.png
plot "planar-axiall6b6.d0.txt" using 2:5 every ::99::108, "analytical_data/BBR_axial_mix.txt" using 1:3 every ::99::108

set output planar-biaxiall6b6.d0g1.d-3DF.png
plot "planar-biaxiall6b6.d0DF.txt" using 2:5 every ::33::43, "analytical_data/BBR_biaxial_mix.txt" using 1:3 every ::33::43

set output planar-biaxiall6b6.d0g1.d-2DF.png
plot "planar-biaxiall6b6.d0.txt" using 2:5 every ::66::76, "analytical_data/BBR_biaxial_mix.txt" using 1:3 every ::66::76

set output planar-biaxiall6b6.d0g1.d-1DF.png
plot "planar-biaxiall6b6.d0.txt" using 2:5 every ::99::108, "analytical_data/BBR_biaxial_mix.txt" using 1:3 every ::99::108
