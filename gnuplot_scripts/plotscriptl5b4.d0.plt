set terminal png
set yr[0:]
set output "../planar-axiall5b4.d0g1.d-2DF.png"
plot "../planar-axiall5b4.d0DF.txt" using 2:5, "../analytical_data/BBR_axial_mix.txt" using 1:3 every ::66::76

set output "../planar-biaxiall5b4.d0g1.d-2DF.png"
plot "../planar-biaxiall5b4.d0DF.txt" using 2:5, "../analytical_data/BBR_biaxial_mix.txt" using 1:3 every ::66::76
