subroutine sph(fineness, xfnew, ilmnew)
  implicit none
  integer fineness
  integer, parameter :: nfsize1=12, nfsize2=20
  integer, parameter :: nfsize3=100000, nfsize4=200000
  integer i, j, k, ncount1, ncount2, nmax, level, n1, n2
  double precision, allocatable :: xfold(:,:)
  double precision :: xfnew(:,:)
  integer, allocatable :: ilmold(:,:), ilm(:)
  integer :: ilmnew(:,:)
  integer, allocatable :: ivert(:,:)
  double precision :: par, par1, par2, rd2, rd4, rd6
  allocate(xfold(1:3,1:size(xfnew,2)))
  allocate(ilmold(1:3,1:nfsize4), ilm(1:3))
  allocate(ivert(1:6,1:nfsize3))

  print *, 'meshgenl17', size(xfnew,1), size(xfnew,2), size(ilmnew,1),&
       size(ilmnew, 2)
  xfold = 0.d0
  xfnew = 0.d0
  ilmold = 0
  ilmnew = 0
  ivert = 0
  nmax = nfsize2
  par=(1.d0+(5.d0**0.5d0))/2.d0
  par1=1.d0/((1.d0+(par**2))**0.5d0)
  par2=par/((1.d0+(par**2))**0.5d0)

  xfold(1:3, 1)=(/ par2, par1,   0.d0/)
  xfold(1:3, 2)=(/-par2, par1,   0.d0/)
  xfold(1:3, 3)=(/-par2,-par1,   0.d0/)
  xfold(1:3, 4)=(/ par2,-par1,   0.d0/)
  xfold(1:3, 5)=(/ par1,   0.d0, par2/)
  xfold(1:3, 6)=(/ par1,   0.d0,-par2/)
  xfold(1:3, 7)=(/-par1,   0.d0,-par2/)
  xfold(1:3, 8)=(/-par1,   0.d0, par2/)
  xfold(1:3, 9)=(/   0.d0, par2, par1/)
  xfold(1:3,10)=(/   0.d0,-par2, par1/)
  xfold(1:3,11)=(/   0.d0,-par2,-par1/)
  xfold(1:3,12)=(/   0.d0, par2,-par1/)

  ilmold(1:3, 1)=(/ 5, 9, 8/)
  ilmold(1:3, 2)=(/ 5, 8,10/)
  ilmold(1:3, 3)=(/ 6, 7,12/)
  ilmold(1:3, 4)=(/ 6,11, 7/)
  ilmold(1:3, 5)=(/ 1, 5, 4/)
  ilmold(1:3, 6)=(/ 1, 4, 6/)
  ilmold(1:3, 7)=(/ 3, 8, 2/)
  ilmold(1:3, 8)=(/ 3, 2, 7/)
  ilmold(1:3, 9)=(/ 9, 1,12/)
  ilmold(1:3,10)=(/ 9,12, 2/)
  ilmold(1:3,11)=(/10,11, 4/)
  ilmold(1:3,12)=(/10, 3,11/)
  ilmold(1:3,13)=(/ 9, 5, 1/)
  ilmold(1:3,14)=(/12, 1, 6/)
  ilmold(1:3,15)=(/ 5,10, 4/)
  ilmold(1:3,16)=(/ 6, 4,11/)
  ilmold(1:3,17)=(/ 8, 9, 2/)
  ilmold(1:3,18)=(/ 7, 2,12/)
  ilmold(1:3,19)=(/ 8, 3,10/)
  ilmold(1:3,20)=(/ 7,11, 3/)

  do level=1,fineness
     ncount1=0
     ncount2=0
     xfnew=0.d0
     ilmnew=0
     ivert=0
     do i=1,6
        ivert(i,1)= i
     end do

     do i=1,3
        xfnew(i,1)= xfold(i,ilmold(1,1))
        xfnew(i,3)= xfold(i,ilmold(2,1))
        xfnew(i,5)= xfold(i,ilmold(3,1))
        xfnew(i,2)= (xfnew(i,1)+xfnew(i,3))/2.d0
        xfnew(i,4)= (xfnew(i,3)+xfnew(i,5))/2.d0
        xfnew(i,6)= (xfnew(i,5)+xfnew(i,1))/2.d0
     end do

     rd2= dsqrt(xfnew(1,2)**2.d0+xfnew(2,2)**2.d0+xfnew(3,2)**2.d0)
     rd4= dsqrt(xfnew(1,4)**2.d0+xfnew(2,4)**2.d0+xfnew(3,4)**2.d0)
     rd6= dsqrt(xfnew(1,6)**2.d0+xfnew(2,6)**2.d0+xfnew(3,6)**2.d0)
     do i=1,3
        xfnew(i,2)= xfnew(i,2)/rd2
        xfnew(i,4)= xfnew(i,4)/rd4
        xfnew(i,6)= xfnew(i,6)/rd6
     end do


     ilmnew(1,1)= ivert(1,1) 
     ilmnew(2,1)= ivert(2,1) 
     ilmnew(3,1)= ivert(6,1) 
     ilmnew(1,2)= ivert(3,1) 
     ilmnew(2,2)= ivert(4,1) 
     ilmnew(3,2)= ivert(2,1) 
     ilmnew(1,3)= ivert(5,1) 
     ilmnew(2,3)= ivert(6,1) 
     ilmnew(3,3)= ivert(4,1) 
     ilmnew(1,4)= ivert(2,1) 
     ilmnew(2,4)= ivert(4,1) 
     ilmnew(3,4)= ivert(6,1) 

     ncount1=6
     ncount2=4

     do i=2,nmax
        do j=1,i-1

           ilm=0
           do n1=1,3
              if(ilmold(1,i)==ilmold(n1,j))then
                 ivert(1,i)=ivert(2*n1-1,j)
                 ilm(1)=n1
              endif
              if(ilmold(2,i)==ilmold(n1,j))then
                 ivert(3,i)=ivert(2*n1-1,j)
                 ilm(2)=n1
              endif
              if(ilmold(3,i)==ilmold(n1,j))then
                 ivert(5,i)=ivert(2*n1-1,j)
                 ilm(3)=n1
              endif
           end do

           if((ilm(1)/=0).and.(ilm(2)/=0)) then
              if((ilm(1)+ilm(2))==3) then
                 ivert(2,i)=ivert(2,j)
              elseif((ilm(1)+ilm(2))==5) then
                 ivert(2,i)=ivert(4,j)
              elseif((ilm(1)+ilm(2))==4) then
                 ivert(2,i)=ivert(6,j)
              endif
           endif

           if((ilm(2)/=0).and.(ilm(3)/=0)) then
              if((ilm(2)+ilm(3))==3) then
                 ivert(4,i)=ivert(2,j)
              elseif((ilm(2)+ilm(3))==5) then
                 ivert(4,i)=ivert(4,j)
              elseif((ilm(2)+ilm(3))==4) then
                 ivert(4,i)=ivert(6,j)
              endif
           endif

           if((ilm(1)/=0).and.(ilm(3)/=0)) then
              if((ilm(1)+ilm(3))==3) then
                 ivert(6,i)=ivert(2,j)
              elseif((ilm(1)+ilm(3))==5) then
                 ivert(6,i)=ivert(4,j)
              elseif((ilm(1)+ilm(3))==4) then
                 ivert(6,i)=ivert(6,j)
              endif
           endif

        end do


        if(ivert(2,i)==0) then
           if((ivert(1,i)==0).and.(ivert(3,i)==0)) then
              ivert(1,i)=ncount1+1
              ivert(2,i)=ncount1+2
              ivert(3,i)=ncount1+3

              do k=1,3
                 xfnew(k,ivert(1,i))= xfold(k,ilmold(1,i))
                 xfnew(k,ivert(3,i))= xfold(k,ilmold(2,i))
                 xfnew(k,ivert(2,i))=(xfnew(k,ivert(1,i))+xfnew(k,ivert(3,i)))/2.d0
              end do

              ncount1=ncount1+3
           elseif((ivert(1,i)==0).and.(ivert(3,i)/=0)) then
              ivert(1,i)=ncount1+1
              ivert(2,i)=ncount1+2
              do k=1,3
                 xfnew(k,ivert(1,i))= xfold(k,ilmold(1,i))
                 xfnew(k,ivert(2,i))=(xfnew(k,ivert(1,i))+xfnew(k,ivert(3,i)))/2.d0
              end do
              ncount1=ncount1+2
           elseif((ivert(1,i)/=0).and.(ivert(3,i)==0)) then
              ivert(2,i)=ncount1+1
              ivert(3,i)=ncount1+2
              do k=1,3
                 xfnew(k,ivert(3,i))= xfold(k,ilmold(2,i))
                 xfnew(k,ivert(2,i))=(xfnew(k,ivert(1,i))+xfnew(k,ivert(3,i)))/2.d0
              end do

              ncount1=ncount1+2
           elseif((ivert(1,i)/=0).and.(ivert(3,i)/=0)) then
              ivert(2,i)=ncount1+1
              do k=1,3
                 xfnew(k,ivert(2,i))=(xfnew(k,ivert(1,i))+xfnew(k,ivert(3,i)))/2.d0
              end do
              ncount1=ncount1+1
           endif
        endif

        if(ivert(4,i)==0) then
           if((ivert(3,i)==0).and.(ivert(5,i)==0)) then
              ivert(3,i)=ncount1+1
              ivert(4,i)=ncount1+2
              ivert(5,i)=ncount1+3
              do k=1,3
                 xfnew(k,ivert(3,i))= xfold(k,ilmold(2,i))
                 xfnew(k,ivert(5,i))= xfold(k,ilmold(3,i))
                 xfnew(k,ivert(4,i))=(xfnew(k,ivert(3,i))+xfnew(k,ivert(5,i)))/2.d0
              end do
              ncount1=ncount1+3
           elseif((ivert(3,i)==0).and.(ivert(5,i)/=0)) then
              ivert(3,i)=ncount1+1
              ivert(4,i)=ncount1+2
              do k=1,3
                 xfnew(k,ivert(3,i))= xfold(k,ilmold(2,i))
                 xfnew(k,ivert(4,i))=(xfnew(k,ivert(3,i))+xfnew(k,ivert(5,i)))/2.d0
              end do
              ncount1=ncount1+2
           elseif((ivert(3,i)/=0).and.(ivert(5,i)==0)) then
              ivert(4,i)=ncount1+1
              ivert(5,i)=ncount1+2
              do k=1,3
                 xfnew(k,ivert(5,i))= xfold(k,ilmold(3,i))
                 xfnew(k,ivert(4,i))=(xfnew(k,ivert(3,i))+xfnew(k,ivert(5,i)))/2.d0
              end do

              ncount1=ncount1+2
           elseif((ivert(3,i)/=0).and.(ivert(5,i)/=0)) then
              ivert(4,i)=ncount1+1
              do k=1,3
                 xfnew(k,ivert(4,i))=(xfnew(k,ivert(3,i))+xfnew(k,ivert(5,i)))/2.d0
              end do
              ncount1=ncount1+1
           endif
        endif

        if(ivert(6,i)==0) then
           if((ivert(5,i)==0).and.(ivert(1,i)==0)) then
              ivert(5,i)=ncount1+2
              ivert(6,i)=ncount1+3
              ivert(1,i)=ncount1+1
              do k=1,3
                 xfnew(k,ivert(5,i))= xfold(k,ilmold(3,i))
                 xfnew(k,ivert(1,i))= xfold(k,ilmold(1,i))
                 xfnew(k,ivert(6,i))=(xfnew(k,ivert(5,i))+xfnew(k,ivert(1,i)))/2.d0
              end do

              ncount1=ncount1+3
           elseif((ivert(5,i)==0).and.(ivert(1,i)/=0)) then
              ivert(5,i)=ncount1+1
              ivert(6,i)=ncount1+2
              do k=1,3
                 xfnew(k,ivert(5,i))= xfold(k,ilmold(3,i))
                 xfnew(k,ivert(6,i))=(xfnew(k,ivert(5,i))+xfnew(k,ivert(1,i)))/2.d0
              end do

              ncount1=ncount1+2
           elseif((ivert(5,i)/=0).and.(ivert(1,i)==0)) then
              ivert(6,i)=ncount1+2
              ivert(1,i)=ncount1+1
              do k=1,3
                 xfnew(k,ivert(1,i))= xfold(k,ilmold(1,i))
                 xfnew(k,ivert(6,i))=(xfnew(k,ivert(5,i))+xfnew(k,ivert(1,i)))/2.d0
              end do

              ncount1=ncount1+2
           elseif((ivert(5,i)/=0).and.(ivert(1,i)/=0)) then
              ivert(6,i)=ncount1+1
              do k=1,3
                 xfnew(k,ivert(6,i))=(xfnew(k,ivert(5,i))+xfnew(k,ivert(1,i)))/2.d0
              end do
              ncount1=ncount1+1
           endif
        endif

        rd2= dsqrt(xfnew(1,ivert(2,i))**2.d0 +xfnew(2,ivert(2,i))**2.d0 +xfnew(3,ivert(2,i))**2.d0)
        rd4= dsqrt(xfnew(1,ivert(4,i))**2.d0 +xfnew(2,ivert(4,i))**2.d0 +xfnew(3,ivert(4,i))**2.d0)
        rd6= dsqrt(xfnew(1,ivert(6,i))**2.d0 +xfnew(2,ivert(6,i))**2.d0 +xfnew(3,ivert(6,i))**2.d0)
        do k=1,3
           xfnew(k,ivert(2,i))= xfnew(k,ivert(2,i))/rd2
           xfnew(k,ivert(4,i))= xfnew(k,ivert(4,i))/rd4
           xfnew(k,ivert(6,i))= xfnew(k,ivert(6,i))/rd6
        end do


        ilmnew(1,ncount2+1)= ivert(1,i) 
        ilmnew(2,ncount2+1)= ivert(2,i) 
        ilmnew(3,ncount2+1)= ivert(6,i) 
        ilmnew(1,ncount2+2)= ivert(3,i) 
        ilmnew(2,ncount2+2)= ivert(4,i) 
        ilmnew(3,ncount2+2)= ivert(2,i) 
        ilmnew(1,ncount2+3)= ivert(5,i) 
        ilmnew(2,ncount2+3)= ivert(6,i) 
        ilmnew(3,ncount2+3)= ivert(4,i) 
        ilmnew(1,ncount2+4)= ivert(2,i) 
        ilmnew(2,ncount2+4)= ivert(4,i) 
        ilmnew(3,ncount2+4)= ivert(6,i) 
        ncount2=ncount2+4
     end do

     nmax=ncount2
     do n2=1,ncount1
        xfold(1,n2)=xfnew(1,n2)
        xfold(2,n2)=xfnew(2,n2)
        xfold(3,n2)=xfnew(3,n2)
     end do

     do n2=1,ncount2
        ilmold(1,n2)=ilmnew(1,n2)
        ilmold(2,n2)=ilmnew(2,n2)
        ilmold(3,n2)=ilmnew(3,n2)
     end do
  end do


  open(20,file='sphere.sph',status='unknown')
  print *, 'meshgenl316', ncount1, ncount2
  write(20,300)ncount1,ncount2/2
300 format(7x,i5,2x,i5)

  do i=1,ncount1
     write(20,320)xfnew(1,i),xfnew(2,i),xfnew(3,i)
  end do
320 format(13x,e20.13,e20.13,e20.13)

  do i=1,ncount2
     write(20,340)ilmnew(1,i),ilmnew(2,i),ilmnew(3,i)
  end do
340 format(3(i8))
  close(20)
  !     d=0.d0

  deallocate (xfold,ilmold,ilm,ivert)
end subroutine sph
!*********************************************
subroutine sf(xfn, elmnew, shpint, shpfs)
  implicit none
  integer, parameter :: NFSIZE=10242,NFSIZE2=20480
  integer elmnew(:,:)
  double precision :: xfn(:,:)
  double precision :: xf(3,size(xfn,2))
  double precision :: shpint(:,:), shpfs(:,:)
  double precision :: e11, e12, e13, e21, e22, e23, e31, e32, e33
  double precision :: e1m, e2m
  double precision :: es1, es2, es3
  double precision :: et1, et2, et3, etm
  double precision :: R11, R12, R13, R21, R22, R23, R31, R32, R33
  double precision :: pi, theta
  double precision :: xi1, xi2, xi3, xj1, xj2, xj3, xk1, xk2, xk3
  double precision :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3
  double precision :: x1, x2, x3, y1, y2, y3, z1, z2, z3
  integer n1, n2, i, j1, j2, j3

  open(20,file='sphere.sph',status='unknown')
  !     read in number of nodes (n1) and number of elements (n2)
  !  read(20,100)n1,n2
  n1 = size(xfn, 2)
  n2 = size(elmnew, 2)/2
  !100 format(7x,i5,2x,i5)
  !     read in nodal positions
  !     angle of rotation theta
  pi = 4.0d0*datan(1.0d0)
  theta =  pi/4.0d0
  do i = 1,n1
     !     read(20,101)XFN(1,i),XFN(2,i), XFN(3,i)
     XF(1,i) = XFN(1,i)
     XF(2,i) = XFN(2,i)
     XF(3,i) = XFN(3,i)
     !101  format(13x,e20.13,e20.13,e20.13)
  enddo
  !  do i = 1,(2*n2)
  !     read(20,103)elmnew(1,i),elmnew(2,i),elmnew(3,i)
  !  enddo
  !103 format(3(i8))

  open(21,file='shpfcta.sph',status='unknown')
  open(22,file='shpfctb.sph',status='unknown')
  open(23,file='shpint.sph',status='unknown')
  do i = 1,2*n2
     j1 = elmnew(1,i)
     j2 = elmnew(2,i)
     j3 = elmnew(3,i)
     !     node 1
     X1 = XFN(1,j1)
     Y1 = XFN(2,j1)
     Z1 = XFN(3,j1)
     !     node 2
     X2 = XFN(1,j2)
     Y2 = XFN(2,j2)
     Z2 = XFN(3,j2)
     !     node 3
     X3 = XFN(1,j3)
     Y3 = XFN(2,j3)
     Z3 = XFN(3,j3)

     e1m = dsqrt((-X1 + X2)**2 + (-Y1 + Y2)**2 + (-Z1 + Z2)**2)
     e11 =(-X1 + X2)/e1m
     e12 =(-Y1 + Y2)/e1m
     e13 =(-Z1 + Z2)/e1m

     es1 = (X3 - X1)
     es2 = (Y3 - Y1)
     es3 = (Z3 - Z1)

     et1 =  e13*es2 - e12*es3
     et2 = -e13*es1 + e11*es3
     et3 = -e11*es2 + e12*es1

     etm = dsqrt(et1**2 + et2**2 + et3**2)
     e31 = et1/etm
     e32 = et2/etm
     e33 = et3/etm
     !     
     e21 = -e13*e32 + e12*e33
     e22 =  e13*e31 - e11*e33
     e23 =  e11*e32 - e12*e31
     !     
     e2m = dsqrt(e21**2 + e22**2 + e23**2)
     e21 = e21/e2m
     e22 = e22/e2m
     e23 = e23/e2m
     !     
     R11=e11
     R12=e12
     R13=e13
     R21=e21
     R22=e22
     R23=e23
     R31=e31
     R32=e32
     R33=e33
     !     
     xi1 = 0.d0
     xi2 = 0.d0
     xi3 = 0.d0
     !     
     a1 = X2 - X1
     a2 = Y2 - Y1
     a3 = Z2 - Z1
     !     
     xj1 = R11*a1 + R12*a2 + R13*a3
     xj2 = R21*a1 + R22*a2 + R23*a3
     xj3 = R31*a1 + R32*a2 + R33*a3
     !     
     a1 = X3 - X1
     a2 = Y3 - Y1
     a3 = Z3 - Z1
     !     
     xk1 = R11*a1 + R12*a2 + R13*a3
     xk2 = R21*a1 + R22*a2 + R23*a3
     xk3 = R31*a1 + R32*a2 + R33*a3
     !     
     a1=xj2-xk2
     a2=xk2-xi2
     a3=xi2-xj2

     b1=xk1-xj1
     b2=xi1-xk1
     b3=xj1-xi1

     c1=xj1*xk2-xk1*xj2
     c2=xk1*xi2-xi1*xk2
     c3=xi1*xj2-xj1*xi2

     d1=a1*xi1+b1*xi2+c1
     d2=a2*xj1+b2*xj2+c2
     d3=a3*xk1+b3*xk2+c3

     a1=a1/d1
     a2=a2/d2
     a3=a3/d3

     b1=b1/d1
     b2=b2/d2
     b3=b3/d3

     shpfs(1,i) = a1
     shpfs(2,i) = a2
     shpfs(3,i) = a3
     shpfs(4,i) = b1
     shpfs(5,i) = b2
     shpfs(6,i) = b3
     shpfs(7,i) = d1
     shpint(1,i) = xj1
     shpint(2,i) = xk1
     shpint(3,i) = xk2

     write(21,104) a1,a2,a3,d1
     write(22,105) b1,b2,b3
     write(23,105) xj1,xk1,xk2
104  format(1pe15.8,3x,1pe15.8,3x,1pe15.8,3x,1pe15.8)
105  format(1pe15.8,3x,1pe15.8,3x,1pe15.8)
  enddo
  close(21)
  close(22)
  close(23)
end subroutine sf

