### @plot sp-DF-g0.04
unset key
set xlabel "q"
set ylabel "DF" offset 2,0
set yr[0:]
set terminal postscript eps enhanced color font "Helvetica" 30 size 12, 4.5
set output 'sp-DF-g0.04-big-presentation-big.eps'
set multiplot
set label "uniaxial" at graph 0.05, graph 0.05
set label "planar" at graph 0.85, graph 0.05
set size 0.5, 1
set origin 0., 0.
set lmargin 8
set rmargin 0
plot 'up-sp-DF-g4.d-2.txt' u 2:5 pt 7 ps 2 t '0.04'
unset label
unset ylabel
#set ylabel "DF/Ca"
set yr[0:]
set origin 0.45, 0.
#set title "Figure 7b"
set lmargin 8
set rmargin 2
set format y ""
set label "biaxial" at graph 0.05, graph 0.05
set label "planar" at graph 0.7, graph 0.05
plot 'bp-sp-DF-g4.d-2.txt' u 2:5 pt 7 ps 2 t '0.04'
unset label
unset format
set nomultiplot
unset output
unset size
unset origin

### @plot sp-all-g0.04
#unset key
#set yr[0:]
unset xtics
set output 'sp-all-g0.04-presentation.eps'
set terminal postscript eps enhanced color font "Helvetica" 30 size 12, 4.5
set multiplot
set size 0.25, 0.5
set format x ""
unset xlabel
set ylabel "DF" offset 2,0
set origin 0., 0.5
set lmargin 8
set rmargin 0
plot 'up-sp-DF-g4.d-2.txt' u 2:5 pt 7 ps 2 t '0.04'
# unset ylabel
# set yr[0:]
set origin 0.225, 0.5
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
plot 'bp-sp-DF-g4.d-2.txt' u 2:5 pt 7 ps 2 t '0.04'
unset format
set format x ""

set ylabel "{/Symbol D}A/A (%)" offset 2,0
set origin 0.5, 0.5
set lmargin 8
set rmargin 0
sp_area_0 = 3.57416564718625807E+02
plot 'up-sp-area-g4.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100) pt 7 ps 2 t '0.04'
# unset ylabel
# set yr[0:]
set origin 0.725, 0.5
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
plot 'bp-sp-area-g4.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100) pt 7 ps 2 t '0.04'
unset format

set ylabel "{/Symbol t}" offset 2,0
set origin 0., 0.
set lmargin 8
set rmargin 0
set yr[0:0.25]
set label "uniaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.70, graph -0.06
plot 'up-sp-tau-g4.d-2.txt' u ($0/10.):3 pt 7 ps 2 t '0.04'
unset label
# unset ylabel
set origin 0.225, 0.
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
set label "biaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.70, graph -0.06
plot 'bp-sp-tau-g4.d-2.txt' u ($0/10.):3 pt 7 ps 2 t '0.04'
unset label
unset format
set yr[0:*]

set ylabel "(Strain energy){/Symbol \264}1000" offset 2,0
set origin 0.5, 0.
set lmargin 8
set rmargin 0
set label "uniaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.7, graph -0.06
plot 'up-sp-W-g4.d-2.txt' u ($0/10.):($1*1000.) pt 7 ps 2 t '0.04'
unset label
# unset ylabel
# set yr[0:]
set origin 0.725, 0.
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
set label "biaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.7, graph -0.06
plot 'bp-sp-W-g4.d-2.txt' u ($0/10.):($1*1000.) pt 7 ps 2 t '0.04'
unset label
unset format
set nomultiplot
unset output
unset size
unset origin

### @plot sp-all-varyca
#unset key
#set yr[0:]
set output 'sp-all-varyca-presentation.eps'
set multiplot
set size 0.25, 0.5
set format x ""
unset xlabel
set ylabel "DF/Ca" offset 2,0
set origin 0., 0.5
set lmargin 8
set rmargin 0
plot 'up-sp-DF-g1.d-2.txt' u 2:($5/0.01) pt 1 ps 2 lw 6 t '0.01', 'up-sp-DF-g2.d-2.txt' u 2:($5/0.02) pt 5 ps 2 t '0.02', 'up-sp-DF-g4.d-2.txt' u 2:($5/0.04) pt 7 ps 2 t '0.04'
set origin 0.225, 0.5
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
plot 'bp-sp-DF-g1.d-2.txt' u 2:($5/0.01) pt 1 ps 2 lw 6 t '0.01', 'bp-sp-DF-g2.d-2.txt' u 2:($5/0.02) pt 5 ps 2 t '0.02', 'bp-sp-DF-g4.d-2.txt' u 2:($5/0.04) pt 7 ps 2 t '0.04'
unset format
set format x ""

set ylabel "{/Symbol D}A/(1000 A Ca^2) (a%)" offset 2,0
set origin 0.5, 0.5
set lmargin 8
set rmargin 0
sp_area_0 = 3.57416564718625807E+02
plot 'up-sp-area-g1.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100/0.01**2/1000) pt 1 ps 2 lw 6 t '0.01', 'up-sp-area-g2.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100/0.02**2/1000) pt 5 ps 2 t '0.02', 'up-sp-area-g4.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100/0.04**2/1000) pt 7 ps 2 t '0.04'
# unset ylabel
# set yr[0:]
set origin 0.725, 0.5
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
plot 'bp-sp-area-g1.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100/0.01**2/1000) pt 1 ps 2 lw \
6 t '0.01', 'bp-sp-area-g2.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100/0.02**2/1000) pt 5 ps 2 t '0.02', 'bp-sp-area-g4.d-2.txt' u 2:(($6-sp_area_0)/sp_area_0*100/0.04**2/1000) pt 7 ps 2 t '0.04'
unset format

set ylabel "{/Symbol t}/Ca" offset 2,0
set origin 0., 0.
set lmargin 8
set rmargin 0
set yr[0:.35]
set label "uniaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.7, graph -0.06
plot 'up-sp-tau-g1.d-2.txt' u ($0/10.):3 pt 1 ps 2 lw 6 t '0.01', 'up-sp-tau-g2.d-2.txt' u ($0/10.):3 pt 5 ps 2 t '0.02', 'up-sp-tau-g4.d-2.txt' u ($0/10.):3 pt 7 ps 2 t '0.04'
unset label
set origin 0.225, 0.
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
set label "biaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.7, graph -0.06
plot 'bp-sp-tau-g1.d-2.txt' u ($0/10.):3 pt 1 ps 2 lw 6 t '0.01', 'bp-sp-tau-g2.d-2.txt' u ($0/10.):3 pt 5 ps 2 t '0.02', 'bp-sp-tau-g4.d-2.txt' u ($0/10.):3 pt 7 ps 2 t '0.04'
unset label
unset format
set yr[0:*]

set ylabel "(Strain energy)/Ca^2" offset 1,0
set origin 0.5, 0.
set lmargin 8
set rmargin 0
set label "uniaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.7, graph -0.06
plot 'up-sp-W-g1.d-2.txt' u ($0/10.):($1/0.01**2) pt 1 ps 2 lw 6 t '0.01', \
    'up-sp-W-g2.d-2.txt' u ($0/10.):($1/0.02**2) pt 5 ps 2 t '0.02', \
    'up-sp-W-g4.d-2.txt' u ($0/10.):($1/0.04**2) pt 7 ps 2 t '0.04'
unset label
# set yr[0:]
set origin 0.725, 0.
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
set label "biaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.7, graph -0.06
plot 'bp-sp-W-g1.d-2.txt' u ($0/10.):($1/0.01**2) pt 1 ps 2 lw 6 t '0.01', 'bp-sp-W-g2.d-2.txt' u ($0/10.):($1/0.02**2) pt 5 ps 2 t '0.02', 'bp-sp-W-g4.d-2.txt' u \
($0/10.):($1/0.04**2) pt 7 ps 2 t '0.04'
unset label
unset format
set nomultiplot
unset output
unset size
unset origin

### @plot varyshape-all-g0.04-presentation
#unset key
#set yr[0:]
set output 'varyshape-all-g0.04-presentation.eps'
set multiplot
set size 0.25, 0.5
set format x ""
unset xlabel
set ylabel "DF" offset 2,0
set origin 0., 0.5
set lmargin 8
set rmargin 0
plot 'up-sp-DF.txt' u 2:5 pt 1 ps 2 lw 6 t 'Sphere', 'up-pr-DF.txt' u 2:(($7-$6*2)/($7+$6*2)) pt 5 ps 2 t 'Prolate', 'up-ob-DF.txt' u 2:(($7-$6*2)/($7+$6*2)) pt 7 ps 2 t 'Oblate'
unset label
# set yr[0:]
set origin 0.225, 0.5
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
plot 'bp-sp-DF.txt' u 2:5 pt 1 ps 2 lw 6 t 'Sphere', 'bp-pr-DF.txt' u 2:(($7-$6*2)/($7+$6*2)) pt 5 ps 2 t 'Prolate', 'bp-ob-DF.txt' u 2:(($7-$6*2)/($7+$6*2)) pt 7 ps 2 t 'Oblate'
unset format
set format x ""

set ylabel "{/Symbol D}A/A (%)" offset 2,0
set origin 0.5, 0.5
set lmargin 8
set rmargin 0
sp_area_0 = 3.57416564718625807E+02
pr_area_0 = sp_area_0
ob_area_0 = sp_area_0
#pr_area_0 = 3.05626567372839816E+02
#ob_area_0 = 3.28995251894758781E+02
plot 'up-sp-area.txt' u 2:(($6-sp_area_0)/sp_area_0*100) pt 1 ps 2 lw 6 t 'Sphere','up-pr-area.txt' u 2:(($6-pr_area_0)/pr_area_0*100) pt 5 ps 2 t 'Prolate', 'up-ob-area.txt' u 2:(($6-ob_area_0)/ob_area_0*100) pt 7 ps 2 t 'Oblate'
# unset ylabel
# set yr[0:]
set origin 0.725, 0.5
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
plot 'bp-sp-area.txt' u 2:(($6-sp_area_0)/sp_area_0*100) pt 1 ps 2 lw 6 t 'Sphere', 'bp-pr-area.txt' u 2:(($6-pr_area_0)/pr_area_0*100) pt 5 ps 2 t 'Prolate', 'bp-ob-area.txt' u 2:(($6-ob_area_0)/ob_area_0*100) pt 7 ps 2 t 'Oblate'
unset format

set ylabel "{/Symbol t}" offset 2,0
set origin 0., 0.
set lmargin 8
set rmargin 0
set yr[0:0.25]
set label "uniaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.70, graph -0.06
plot 'up-sp-tau.txt' u ($0/10.):3 pt 1 ps 2 lw 6 t 'Sphere', 'up-pr-tau.txt' u ($0/10.):3 pt 5 ps 2 t 'Prolate', 'up-ob-tau.txt' u \
($0/10.):3 pt 7 ps 2 t 'Oblate'
unset label
set origin 0.225, 0.
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
set label "biaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.70, graph -0.06
plot 'bp-sp-tau.txt' u ($0/10.):3 pt 1 ps 2 lw 6 t 'Sphere', 'bp-pr-tau.txt' u ($0/10.):3 pt 5 ps 2 t 'Prolate', 'bp-ob-tau.txt' u ($0/10.):3 pt 7 ps 2 t 'Oblate'
unset label
unset format
set yr[0:*]

set ylabel "(Strain energy){/Symbol \264}1000" offset 2,0
set origin 0.5, 0.
set lmargin 8
set rmargin 0
set label "uniaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.70, graph -0.06
plot 'up-sp-W.txt' u ($0/10.):($1*1000.) pt 1 ps 2 lw 6 t 'Sphere', 'up-pr-W.txt' u ($0/10.):($1*1000.) pt 5 ps 2 t 'Prolate', 'up-ob-W.txt' u \
($0/10.):($1*1000.) pt 7 ps 2 t 'Oblate'
unset label
# set yr[0:]
set origin 0.725, 0.
set lmargin 8
set rmargin 2
unset ylabel
set format y ""
set label "biaxial" at graph 0.05, graph -0.06
set label "planar" at graph 0.70, graph -0.06
plot 'bp-sp-W.txt' u ($0/10.):($1*1000.) pt 1 ps 2 lw 6 t 'Sphere', 'bp-pr-W.txt' u ($0/10.):($1*1000.) pt 5 ps 2 t 'Prolate', 'bp-ob-W.txt' u ($0/10.):($1*1000.) pt 7 ps 2 t 'Oblate'
unset label
unset format
set nomultiplot
unset output
unset size
unset origin
