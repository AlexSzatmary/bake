# A heads-up display of lots of variables for one run.

set multiplot
set size 0.5, 1./3.
set origin 0., 2./3.
plot './capsulex__00001.txt' u 1:2 w l t 'x', './capsulex__00001.txt' u 1:3 w l t 'y', './capsulex__00001.txt' u 1:4 w l t 'z'
set origin 0., 1./3.
plot 'TaylorDF__00001.txt' w l t 'DF'
set origin 0., 0.
plot 'TaylorDF__00001.txt' u 1:3 w l t 'rmin', 'TaylorDF__00001.txt' u 1:4 w l t 'rmax', 'minmaxx___00001.txt' w l t 'min x', 'minmaxx___00001.txt' u 1:3 w l t 'max x', 'minmaxy___00001.txt' w l t 'min y', 'minmaxy___00001.txt' u 1:3 w l t 'max y', 'minmaxz___00001.txt' w l t 'min z', 'minmaxz___00001.txt' u 1:3 w l t 'max z'
set origin 0.5, 2./3.
plot 'lambda1___00001.txt' w l t 'min(l_1)', 'lambda1___00001.txt' u 1:3 w l t 'max(l_1)', 'lambda2___00001.txt' w l t 'min(l_2)', 'lambda2___00001.txt' u 1:3 w l t 'max(l_2)'
set origin 0.5, 1./3.
plot 'meanfluidv.txt' w l t 'u', 'meanfluidv.txt' u 1:3 w l t 'v', 'meanfluidv.txt' u 1:4 w l t 'w'
set origin 0.5, 0.
plot 'meanforce.txt' w l t 'fx', 'meanforce.txt' u 1:3 w l t 'fy', 'meanforce.txt' u 1:4 w l t 'fz'

set nomultiplot
pause -1 "Type return to continue"
