f(x) = DFinfty*(1-exp(-t/tau))
DFinfty = 37.5*$capillary_no$
tau = 12.0793*$capillary_no$
fit f(x) "$wd$" using ($1/1000):2 via DFinfty, tau
