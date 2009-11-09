f(x) = DFinfty*(1-exp(-x/tau))
DFinfty = 37.5*$capillary_no$
tau = 12.0793*$capillary_no$
fit f(x) "./$wd$/TaylorDF__00001.txt" using ($1/1000):2 via DFinfty, tau
