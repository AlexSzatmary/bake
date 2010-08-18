subroutine inmv(lcube,radx,h,xfn,elmnew,n1,n2,elmv,imic,irec, rlbs,ltb)
!     Initialize microvilli
  implicit none

  integer, parameter :: nfsize=10242, nfsize2=20480
! mic - number of microvilli; mrec number of microvilli per receptor
  integer, parameter :: mic=252,mrec=50
  integer n1,n2,i,j,k
  real lcube,h,radx
  real xc1, yc1, zc1, xc2, yc2, zc2, dist, distm
  real x1,y1,z1,x2,y2,z2,x3,y3,z3
  integer j1,j2,j3,m1,m2,m3,m4
  integer imic(1:mic),irec(1:mrec,1:mic)
!     elmv - element nearest a given microvillus
  integer elmnew(1:3,1:nfsize2),elmv(1:12,1:mic)
  integer el(1:3),ev(1:12),ltb(1:2,1:mrec,1:mic)

!     cmv - coordinates of microvilli
  real xfn(1:3,1:nfsize),cmv(1:3,1:mic)
!     edr does nothing. rlbs - 
  real edr(1:nfsize2),rlbs(1:3,1:mrec,1:mic)

!     Load positions of microvilli. These are on a sphere of radius 1. They 
!     don't seem to match the exact positions of membrane nodes.
  open(25,file='icos252.sph',status='unknown')
  xc1 = 0.0
  yc1 = 0.0
  zc1 = 0.0
!     Scale the microvilli coordinates
  do i = 1,mic
     read(25,21) cmv(1,i),cmv(2,i),cmv(3,i)
     cmv(1,i) = radx*cmv(1,i)/h
     cmv(2,i) = radx*cmv(2,i)/h
     cmv(3,i) = radx*cmv(3,i)/h
     cmv(1,i) = (cmv(1,i) + 0.5*lcube/h)
     cmv(2,i) = (cmv(2,i) + 0.38*lcube/h)
     cmv(3,i) = (cmv(3,i )+ 0.375*lcube/h)
  end do
  print *, cmv(1, i), cmv(2, i), cmv(3, i)
21 format(4x, e20.13, 3x, e20.13, 3x, e20.13)

!     Zero arrays. In Fortran 90, these could be single lines.
  do j = 1, mic
     imic(j)=0
     do i = 1, mrec
        irec(i, j) = 0
        rlbs(1, i, j) = 0.
        rlbs(2, i, j) = 0.
        rlbs(3, i, j) = 0.
     end do
  end do

!     Zero ltb
  ltb=0

!     Sets up elmv: for a given microvillus, pick the element with the closest
!     centroid.
  do k = 1,mic
     distm = 1.e05
     xc1 = cmv(1,k)
     yc1 = cmv(2,k)
     zc1 = cmv(3,k)
     do i = 1,2*n2
        j1 = elmnew(1,i)
        j2 = elmnew(2,i)
        j3 = elmnew(3,i)
        x1 = xfn(1,j1)
        y1 = xfn(2,j1)
        z1 = xfn(3,j1)
        x2 = xfn(1,j2)
        y2 = xfn(2,j2)
        z2 = xfn(3,j2)
        x3 = xfn(1,j3)
        y3 = xfn(2,j3)
        z3 = xfn(3,j3)
        xc2= (x1+x2+x3)/3.
        yc2= (y1+y2+y3)/3.
        zc2= (z1+z2+z3)/3.
        dist=sqrt((xc1-xc2)**2+(yc1-yc2)**2+(zc1-zc2)**2)
        if(dist <= distm) then
           distm=dist
           elmv(1,k)=j1
           elmv(2,k)=j2
           elmv(3,k)=j3
        endif
     end do
  end do

  print *, dist, distm
  print *, xfn(:, 1)
  print *, cmv(:, 1)
!     Looping over all elements near microvilli and over all elements, period.
!     This seems to make elmv hold the element nearest a given microvilli,
!     as well as the three elements touching that element.
  do k=1,mic
     do m1 = 1, 3
        ev(m1)=elmv(m1,k)
     end do
     do i = 1, 2*n2
        do m2 = 1, 3
           el(m2) = elmnew(m2, i)
        end do
        if((ev(1) == el(1)).and.(ev(2) == el(2)) &
             .and.(ev(3) == el(3))) then

        else
           if(((ev(1) == el(1)).or.(ev(1) == el(2))&
                .or.(ev(1) == el(3)))&
                .and.((ev(2) == el(1)).or.(ev(2) == el(2))&
                .or.(ev(2) == el(3)))) then
              if((ev(1) == el(1)).and.(ev(2) == el(2))) ev(4)=el(3)
              if((ev(1) == el(1)).and.(ev(2) == el(3))) ev(4)=el(2)
              if((ev(1) == el(2)).and.(ev(2) == el(1))) ev(4)=el(3)
              if((ev(1) == el(2)).and.(ev(2) == el(3))) ev(4)=el(1)
              if((ev(1) == el(3)).and.(ev(2) == el(2))) ev(4)=el(1)
              if((ev(1) == el(3)).and.(ev(2) == el(1))) ev(4)=el(2)
           endif
           if(((ev(1) == el(1)).or.(ev(1) == el(2)) &
                .or.(ev(1) == el(3))) &
                .and.((ev(3) == el(1)).or.(ev(3) == el(2))&
                .or.(ev(3) == el(3)))) then
              if((ev(1) == el(1)).and.(ev(3) == el(2))) ev(5)=el(3)
              if((ev(1) == el(1)).and.(ev(3) == el(3))) ev(5)=el(2)
              if((ev(1) == el(2)).and.(ev(3) == el(1))) ev(5)=el(3)
              if((ev(1) == el(2)).and.(ev(3) == el(3))) ev(5)=el(1)
              if((ev(1) == el(3)).and.(ev(3) == el(2))) ev(5)=el(1)
              if((ev(1) == el(3)).and.(ev(3) == el(1))) ev(5)=el(2)
           endif
           if(((ev(2) == el(1)).or.(ev(2) == el(2)) &
                .or.(ev(2) == el(3))) &
                .and.((ev(3) == el(1)).or.(ev(3) == el(2)) &
                .or.(ev(3) == el(3)))) then
              if((ev(2) == el(1)).and.(ev(3) == el(2))) ev(6)=el(3)
              if((ev(2) == el(1)).and.(ev(3) == el(3))) ev(6)=el(2)
              if((ev(2) == el(2)).and.(ev(3) == el(1))) ev(6)=el(3)
              if((ev(2) == el(2)).and.(ev(3) == el(3))) ev(6)=el(1)
              if((ev(2) == el(3)).and.(ev(3) == el(2))) ev(6)=el(1)
              if((ev(2) == el(3)).and.(ev(3) == el(1))) ev(6)=el(2)
           endif
        endif
     end do
     do i=1,2*n2
        do m3=1,3
           el(m3)=elmnew(m3,i)
        end do
        if(((ev(2) /= el(1)).and.(ev(2) /= el(2))) &
             .and.(ev(2) /= el(3))) then
           if((ev(1) == el(1)).and.(ev(4) == el(2))) ev(7)=el(3)
           if((ev(1) == el(1)).and.(ev(4) == el(3))) ev(7)=el(2)
           if((ev(1) == el(2)).and.(ev(4) == el(1))) ev(7)=el(3)
           if((ev(1) == el(2)).and.(ev(4) == el(3))) ev(7)=el(1)
           if((ev(1) == el(3)).and.(ev(4) == el(2))) ev(7)=el(1)
           if((ev(1) == el(3)).and.(ev(4) == el(1))) ev(7)=el(2)
        endif
        if((ev(3) /= el(1)).and.(ev(3) /= el(2)) &
             .and.(ev(3) /= el(3))) then
           if((ev(1) == el(1)).and.(ev(5) == el(2))) ev(8)=el(3)
           if((ev(1) == el(1)).and.(ev(5) == el(3))) ev(8)=el(2)
           if((ev(1) == el(2)).and.(ev(5) == el(1))) ev(8)=el(3)
           if((ev(1) == el(2)).and.(ev(5) == el(3))) ev(8)=el(1)
           if((ev(1) == el(3)).and.(ev(5) == el(2))) ev(8)=el(1)
           if((ev(1) == el(3)).and.(ev(5) == el(1))) ev(8)=el(2)
        endif
        if((ev(1) /= el(1)).and.(ev(1) /= el(2)) &
             .and.(ev(1) /= el(3))) then
           if((ev(2) == el(1)).and.(ev(4) == el(2))) ev(9)=el(3)
           if((ev(2) == el(1)).and.(ev(4) == el(3))) ev(9)=el(2)
           if((ev(2) == el(2)).and.(ev(4) == el(1))) ev(9)=el(3)
           if((ev(2) == el(2)).and.(ev(4) == el(3))) ev(9)=el(1)
           if((ev(2) == el(3)).and.(ev(4) == el(2))) ev(9)=el(1)
           if((ev(2) == el(3)).and.(ev(4) == el(1))) ev(9)=el(2)
        endif
        if((ev(3) /= el(1)).and.(ev(3) /= el(2)) &
             .and.(ev(3) /= el(3))) then
           if((ev(2) == el(1)).and.(ev(6) == el(2))) ev(10)=el(3)
           if((ev(2) == el(1)).and.(ev(6) == el(3))) ev(10)=el(2)
           if((ev(2) == el(2)).and.(ev(6) == el(1))) ev(10)=el(3)
           if((ev(2) == el(2)).and.(ev(6) == el(3))) ev(10)=el(1)
           if((ev(2) == el(3)).and.(ev(6) == el(2))) ev(10)=el(1)
           if((ev(2) == el(3)).and.(ev(6) == el(1))) ev(10)=el(2)
        endif
        if((ev(1) /= el(1)).and.(ev(1) /= el(2)) &
             .and.(ev(1) /= el(3))) then
           if((ev(3) == el(1)).and.(ev(5) == el(2))) ev(11)=el(3)
           if((ev(3) == el(1)).and.(ev(5) == el(3))) ev(11)=el(2)
           if((ev(3) == el(2)).and.(ev(5) == el(1))) ev(11)=el(3)
           if((ev(3) == el(2)).and.(ev(5) == el(3))) ev(11)=el(1)
           if((ev(3) == el(3)).and.(ev(5) == el(2))) ev(11)=el(1)
           if((ev(3) == el(3)).and.(ev(5) == el(1))) ev(11)=el(2)
        endif
        if((ev(2) /= el(1)).and.(ev(2) /= el(2)) &
             .and.(ev(2) /= el(3))) then
           if((ev(3) == el(1)).and.(ev(6) == el(2))) ev(12)=el(3)
           if((ev(3) == el(1)).and.(ev(6) == el(3))) ev(12)=el(2)
           if((ev(3) == el(2)).and.(ev(6) == el(1))) ev(12)=el(3)
           if((ev(3) == el(2)).and.(ev(6) == el(3))) ev(12)=el(1)
           if((ev(3) == el(3)).and.(ev(6) == el(2))) ev(12)=el(1)
           if((ev(3) == el(3)).and.(ev(6) == el(1))) ev(12)=el(2)
        endif
     end do
     do m4 = 1, 12
        elmv(m4, k) = ev(m4)
     end do
  end do

!       do i = 1,12
!          print *, ev(i)
!       end do


  return
end subroutine inmv
!************************************************************
subroutine fmv(klok,mass,length,time,xfn,frc,elmv, &
     imic,irec,rlbs,rlbf,ihist,ltb)
  implicit none

  integer, parameter :: nfsize=10242
  integer, parameter :: mic=252, mrec=50
  real mass,length,time
  real xfn(1:3,1:nfsize),frc(1:3,1:nfsize),rlbf(1:mrec,1:mic)
  real rlbs(1:3,1:mrec,1:mic),cmv(1:3,1:mic),frm(1:3,1:mic)
  integer imic(1:mic),irec(1:mrec,1:mic),elmv(1:12,1:mic),ev(1:12)
  integer ihist(0:50),ltb(1:2,1:mrec,1:mic)
  real a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,nlig,blig
  real rmv,tmx,tmy,tmz,cmx,cmy,cmz,ft1,fx1,fy1,fz1,fhist
  real dms,dism,pb,pran,rm,rl,rk,rk0,fk,fk0,bt,f0
  real sigb,sigt,sige,dt,sub,fx,fy,fz
  integer i,j,j1,j2,j3,imictemp,irectemp,klok,icount,jcount,m

  nlig=10.0
! height of the substrate?
  sub=length*08.
  dt=time
  rmv=0.35
  rl=0.1
  bt=310.*.013807
  fk0=1.0
  rk0=1.0
  sigt=0.99e06
  sigb=1.0e06
  sige=sigb
  dism=0.5
  ft1=0.
  fx1=0.
  fy1=0.
  fz1=0.
  ihist=0
  rlbf=0.

  do j=1,mic
     j1=elmv(1,j)
     j2=elmv(2,j)
     j3=elmv(3,j)
     do m=1,12
        ev(m)=elmv(m,j)
     end do
! Centroid of element
     cmx=length*(xfn(1,j1)+xfn(1,j2)+xfn(1,j3))/3.
     cmy=length*(xfn(2,j1)+xfn(2,j2)+xfn(2,j3))/3.
     cmz=length*(xfn(3,j1)+xfn(3,j2)+xfn(3,j3))/3.
     a1=length*(xfn(1,j2)-xfn(1,j1))
     a2=length*(xfn(2,j2)-xfn(2,j1))
     a3=length*(xfn(3,j2)-xfn(3,j1))
     b1=length*(xfn(1,j3)-xfn(1,j1))
     b2=length*(xfn(2,j3)-xfn(2,j1))
     b3=length*(xfn(3,j3)-xfn(3,j1))
     c1=a2*b3-b2*a3
     c2=b1*a3-a1*b3
     c3=a1*b2-b1*a2
! d_i is the normal to the element closest to microvillus j
     d1=c1/sqrt(c1**2+c2**2+c3**2)
     d2=c2/sqrt(c1**2+c2**2+c3**2)
     d3=c3/sqrt(c1**2+c2**2+c3**2)
! tm_i is the tip of the microvillus.
     tmx=d1*rmv+cmx
     tmy=d2*rmv+cmy
     tmz=d3*rmv+cmz
     fx=0.
     fy=0.
     fz=0.
     frm(1,j)=0.
     frm(2,j)=0.
     frm(3,j)=0.

!count number of bound receptors on microvillus j
     blig=0.
     do i=1,mrec
        if(irec(i,j) == 1) then
           blig=blig+1.
        endif
     end do

     dms = (tmy-sub)

     if(imic(j) == 0) then
        imictemp=0
        irectemp=0
        if(dms <= dism) then
           imictemp=2
           do i=1,mrec
! Probabilistic binding
              call random_number(pran)
              rm=sqrt((tmy-sub)**2)
! If the actual distance between the receptor and ligand is less than the
! unstressed length, the binding constant is k_0 for a single bond, so fk is
! flig*fk0, where flig is the number of free ligands, nlig-blig
              if((rm-rl) <= 0.) then
                 fk=(nlig-blig)*fk0
! If the distance between the receptor and ligand is greater than the
! unstressed length, use the dembo model:
! Dembo, M. 1994. On peeling an adherent cell from a surface. In vol. 24 of 
! series: Lectures on Mathematics in the Life Sciences, Some mathematical 
! problems in biology. American Mathematical Society, Providence, RI. 51–77.
              else
                 fk=(nlig-blig)*fk0*exp(-sigt*(rm-rl)**2/(2.*bt))
              endif
! The probability of breaking in a single timestep is given by this exponential
              pb=1.-exp(-fk*dt)
              if(pb.gt.pran) then
                 imictemp=3
                 irectemp=1
! The point of attachment is the closest point from the microvillus on the
! surface.
                 rlbs(1,i,j)=tmx
                 rlbs(2,i,j)=sub
                 rlbs(3,i,j)=tmz
! I think this is related to the bond lifetime records.
                 ltb(1,i,j)=1
                 ltb(2,i,j)=klok
                 write(216,*)klok,i,j,"binding"
              else
                 irectemp=0
              endif
              irec(i,j)=irectemp
              irectemp=0
           end do
           imic(j)=imictemp
           imictemp=0
        endif
     elseif(imic(j) == 1) then 
        imictemp=0     
        irectemp=0
        if(dms <= dism) then
           imictemp=2
           do i=1,mrec
              if(irec(i,j) == 1) then
                 irectemp=1
!     check probabilistic unbinding
                 call random_number(pran)
                 rm=sqrt((tmx-rlbs(1,i,j))**2+(tmy-rlbs(2,i,j))**2+ &
                      (tmz-rlbs(3,i,j))**2)
                 if(rm <= rl) rk=rk0
                 if(rm > rl) rk=rk0*exp((sige-sigt)*(rm-rl)**2/(2.*bt))
                 pb=1.-exp(-rk*dt)
                 if(pb < pran) then
                    irectemp=1
                    imictemp=3
!     calculate force of bond, add to mic force
                    if(rm > rl) then
                       fx=fx+sige*(rm-rl)*(rlbs(1,i,j)-tmx)/rm
                       fy=fy+sige*(rm-rl)*(rlbs(2,i,j)-tmy)/rm
                       fz=fz+sige*(rm-rl)*(rlbs(3,i,j)-tmz)/rm
                       fhist=sqrt(fx**2+fy**2+fz**2)/10000.
                       rlbf(i,j)=fhist*10.
                       if(fhist < 50.) ihist(int(fhist))=ihist(int(fhist))+1
                       if(fhist.ge.50.) ihist(50)=ihist(50)+1
                    else
                       ihist(0)=ihist(0)+1
                    endif
                 else
                    if(rm <= rl) f0=0.
                    if(rm > rl) f0=sige*(rm-rl)/1000.
                    write(218,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    write(219,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    ltb(1,i,j)=0
                    ltb(2,i,j)=0
                    irectemp=0
                    rlbs(1,i,j)=0.
                    rlbs(2,i,j)=0.
                    rlbs(3,i,j)=0.
                    write(216,*)klok,i,j,"unbinding"
                 endif
              else
!     check probabilistic binding
                 call random_number(pran)
                 rm=sqrt((tmy-sub)**2)
                 if((rm-rl) <= 0.) then
                    fk=(nlig-blig)*fk0
                 else
                    fk=(nlig-blig)*fk0*exp(-sigt*(rm-rl)**2/(2.*bt))
                 endif
                 pb=1.-exp(-fk*dt)
                 if(pb > pran) then
                    irectemp=1
! I think the following line might be wrong. This statement is in the chunk for
! if the microvillus root is in the binding zone, but an imic of 1 means that
! the microvillus root is outside of the binding zone. On the other hand,
! having the wrong imic shouldn't matter too much, especially if it's a 1.
                    imictemp=1
                    rlbs(1,i,j)=tmx
                    rlbs(2,i,j)=sub
                    rlbs(3,i,j)=tmz
                    ltb(1,i,j)=1
                    ltb(2,i,j)=klok
                    write(216,*)klok,i,j,"binding"
                 else
                    irectemp=0
                 endif
              endif
              irec(i,j)=irectemp
              irectemp=0
           end do
        else
           imictemp=0
           irectemp=0
           do i=1,mrec
              if(irec(i,j) == 1) then
                 irectemp=1
!     check probabilistic unbinding
                 call random_number(pran)
                 rm=sqrt((tmx-rlbs(1,i,j))**2+(tmy-rlbs(2,i,j))**2+ &
                      (tmz-rlbs(3,i,j))**2)
                 if(rm <= rl) rk=rk0
                 if(rm > rl) rk=rk0*exp((sige-sigt)*(rm-rl)**2/(2.*bt))
                 pb=1.-exp(-rk*dt)
                 if(pb < pran) then
                    irectemp=1
                    imictemp=1
!     calculate force of bond, add to mic force
                    if(rm > rl) then
                       fx=fx+sige*(rm-rl)*(rlbs(1,i,j)-tmx)/rm
                       fy=fy+sige*(rm-rl)*(rlbs(2,i,j)-tmy)/rm
                       fz=fz+sige*(rm-rl)*(rlbs(3,i,j)-tmz)/rm
                       fhist=sqrt(fx**2+fy**2+fz**2)/10000.
                       rlbf(i,j)=fhist*10.
                       if(fhist < 50.) ihist(int(fhist))=ihist(int(fhist))+1
                       if(fhist.ge.50.) ihist(50)=ihist(50)+1
                    else
                       ihist(0)=ihist(0)+1
                    endif
                 else
                    if(rm <= rl) f0=0.
                    if(rm > rl) f0=sige*(rm-rl)/1000.
                    write(218,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    write(219,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    ltb(1,i,j)=0
                    ltb(2,i,j)=0
                    irectemp=0
                    rlbs(1,i,j)=0.
                    rlbs(2,i,j)=0.
                    rlbs(3,i,j)=0.
                    write(216,*)klok,i,j,"unbinding"
                 endif
                 irec(i,j)=irectemp
                 irectemp=0
              endif
           end do
        endif
        imic(j)=imictemp
        imictemp=0
     elseif(imic(j) == 2) then    
        imictemp=0     
        if(dms <= dism) then
           imictemp=2
           do i=1,mrec
!     check probabilistic binding
              call random_number(pran)
              rm=sqrt((tmy-sub)**2)
              if((rm-rl) <= 0.) then
                 fk=(nlig-blig)*fk0
              else
                 fk=(nlig-blig)*fk0*exp(-sigt*(rm-rl)**2/(2.*bt))
              endif
              pb=1.-exp(-fk*dt)
              if(pb > pran) then
                 imictemp=3
                 irectemp=1
                 rlbs(1,i,j)=tmx
                 rlbs(2,i,j)=sub
                 rlbs(3,i,j)=tmz
                 ltb(1,i,j)=1
                 ltb(2,i,j)=klok
                 write(216,*)klok,i,j,"binding"
              else
                 irectemp=0
              endif
              irec(i,j)=irectemp
              irectemp=0
           end do
        endif
        imic(j)=imictemp
        imictemp=0
     elseif(imic(j) == 3) then 
        imictemp=0     
        irectemp=0
        if(dms <= dism) then
           imictemp=2
           do i=1,mrec
              if(irec(i,j) == 1) then
                 irectemp=1
!     check probabilistic unbinding
                 call random_number(pran)
                 rm=sqrt((tmx-rlbs(1,i,j))**2+(tmy-rlbs(2,i,j))**2+ &
                      (tmz-rlbs(3,i,j))**2)
                 if(rm <= rl) rk=rk0
                 if(rm > rl) rk=rk0*exp((sige-sigt)*(rm-rl)**2/(2.*bt))
                 pb=1.-exp(-rk*dt)
                 if(pb < pran) then
                    irectemp=1
                    imictemp=3
!     calculate force of bond, add to mic force
                    if(rm > rl) then
                       fx=fx+sige*(rm-rl)*(rlbs(1,i,j)-tmx)/rm
                       fy=fy+sige*(rm-rl)*(rlbs(2,i,j)-tmy)/rm
                       fz=fz+sige*(rm-rl)*(rlbs(3,i,j)-tmz)/rm
                       fhist=sqrt(fx**2+fy**2+fz**2)/10000.
                       rlbf(i,j)=fhist*10.
                       if(fhist < 50.) ihist(int(fhist))= &
                            ihist(int(fhist))+1
                       if(fhist >= 50.) ihist(50)=ihist(50)+1
                    else
                       ihist(0)=ihist(0)+1
                    endif
                 else
                    if(rm <= rl) f0=0.
                    if(rm > rl) f0=sige*(rm-rl)/1000.
                    write(218,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    write(219,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    ltb(1,i,j)=0
                    ltb(2,i,j)=0
                    irectemp=0
                    rlbs(1,i,j)=0.
                    rlbs(2,i,j)=0.
                    rlbs(3,i,j)=0.
                    write(216,*)klok,i,j,"unbinding"
                 endif
              else
!     check probabilistic binding
                 call random_number(pran)
                 rm=sqrt((tmy-sub)**2)
                 if((rm-rl) <= 0.) then
                    fk=(nlig-blig)*fk0
                 else
                    fk=(nlig-blig)*fk0*exp(-sigt*(rm-rl)**2/(2.*bt))
                 endif
                 pb=1.-exp(-fk*dt)
                 if(pb > pran) then
                    irectemp=1
                    imictemp=1
                    rlbs(1,i,j)=tmx
                    rlbs(2,i,j)=sub
                    rlbs(3,i,j)=tmz
                    ltb(1,i,j)=1
                    ltb(2,i,j)=klok
                    write(216,*)klok,i,j,"binding"
                 else
                    irectemp=0
                 endif
              endif
              irec(i,j)=irectemp
              irectemp=0
           end do
        else
           imictemp=0
           irectemp=0
           do i=1,mrec
              if(irec(i,j) == 1) then
                 irectemp=1
! check probabilistic unbinding
                 call random_number(pran)
                 rm=sqrt((tmx-rlbs(1,i,j))**2+(tmy-rlbs(2,i,j))**2+ &
                      (tmz-rlbs(3,i,j))**2)
                 if(rm <= rl) rk=rk0
                 if(rm > rl) rk=rk0*exp((sige-sigt)*(rm-rl)**2/(2.*bt))
                 pb=1.-exp(-rk*dt)
                 if(pb < pran) then
                    irectemp=1
                    imictemp=1
! calculate force of bond, add to mic force
                    if(rm > rl) then
                       fx=fx+sige*(rm-rl)*(rlbs(1,i,j)-tmx)/rm
                       fy=fy+sige*(rm-rl)*(rlbs(2,i,j)-tmy)/rm
                       fz=fz+sige*(rm-rl)*(rlbs(3,i,j)-tmz)/rm
                       fhist=sqrt(fx**2+fy**2+fz**2)/10000.
                       rlbf(i,j)=fhist*10.
                       if(fhist < 50.) ihist(int(fhist))= &
                            ihist(int(fhist))+1
                       if(fhist >= 50.) ihist(50)=ihist(50)+1
                    else
                       ihist(0)=ihist(0)+1
                    endif
                 else
                    if(rm <= rl) f0=0.
                    if(rm > rl) f0=sige*(rm-rl)/1000.
                    write(218,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    write(219,252) ltb(2,i,j),klok-ltb(2,i,j),f0
                    ltb(1,i,j)=0
                    ltb(2,i,j)=0
                    irectemp=0
                    rlbs(1,i,j)=0.
                    rlbs(2,i,j)=0.
                    rlbs(3,i,j)=0.
                    write(216,*)klok,i,j,"unbinding"
                 endif
                 irec(i,j)=irectemp
                 irectemp=0
              endif
           end do
        endif
        imic(j)=imictemp
        imictemp=0
     endif
     fx1=fx1+fx/1000.
     fy1=fy1+fy/1000.
     fz1=fz1+fz/1000.
     ft1=ft1+sqrt(fx**2+fy**2+fz**2)/1000.
     ! Scales force back to program units; fostar could be used instead of the
     ! (mass*length/time**2) group.
     frm(1,j)=fx/(mass*length/time**2)
     frm(2,j)=fy/(mass*length/time**2)
     frm(3,j)=fz/(mass*length/time**2)
     do m=1,12
        frc(1,ev(m))=frc(1,ev(m))+frm(1,j)/12.
        frc(2,ev(m))=frc(2,ev(m))+frm(2,j)/12.
        frc(3,ev(m))=frc(3,ev(m))+frm(3,j)/12.
     end do
  end do
  icount=0
  jcount=0
  do j=1,mic
     if(imic(j) == 1) jcount=jcount+1
     if(imic(j) == 3) jcount=jcount+1
     do i=1,mrec
        if(irec(i,j) == 1) icount=icount+1
     end do
  end do
!     Force in pN
  write(210,251)klok,icount,jcount,ft1,fx1,fy1,fz1
251 format(i9,x,i5,x,i5,f10.3,x,f10.3,x,f10.3,x,f10.3)
252 format(i9,x,i9,x,f15.5)
  return
end subroutine fmv
!************************************************************
