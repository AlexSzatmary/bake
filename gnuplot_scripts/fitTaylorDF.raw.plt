set t unknown
f(x) = DF_Infty*(1-exp(-x/tau))
DF_Infty = 37.5*$capillary_no$
tau = 12.0793*$capillary_no$**2
R = 1.+$nellip$.
set yr [*:*]
fit f(x) "./$wd$/TaylorDF__00001.txt" using ($1/1000):(((1+R)*$2+1-R)/((1-R)*$2+1+R)) via DF_Infty, tau
set yrange [*:*] writeback
plot "./$wd$/TaylorDF__00001.txt" u ($1/1000):(((((1+R)*$2+1-R)/((1-R)*$2+1+R))-f($1/1000))/DF_Infty)
set yrange restore
show yrange
