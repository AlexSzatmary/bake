      subroutine inspher(lcube,radx,h,xfi,xfn,elmnew,shpint,shpfs) 

      implicit none

      integer, parameter :: nfsize=10242,nfsize2=20480
      integer n1,n2,i
      real lcube,h,theta,x,y,radx
      real xfn(1:3,1:nfsize)
      real xfi(1:3,1:nfsize)
      integer elmnew(1:3,1:nfsize2)
      real shpint(1:3,1:nfsize2),shpfs(1:7,1:nfsize2)
      real, parameter :: pi = 3.141592654
      theta =  pi/4.0

      open(20,file='sphererad1.out',status='unknown')
      open(21,file='shpfcta.sph',status='unknown')
      open(22,file='shpfctb.sph',status='unknown')
      open(23,file='shpint.sph',status='unknown')
!     !m This n1 and n2 ought to correspond to nfsize and nfsize2
!     and should be reflected throughout the code.
      read(20,100) n1, n2
 100  format(7x,i5,2x,i5)

!     Loads the sphere
      do i = 1, n1
         read(20,101) xfn(1,i),xfn(2,i), xfn(3,i)

!     Rotate & scale sphere
         x = xfn(1,i)*cos(theta) - xfn(2,i)*sin(theta)
         y = xfn(1,i)*sin(theta) + xfn(2,i)*cos(theta)
         xfn(1,i) = radx*x/h
         xfn(2,i) = radx*y/h
         xfn(3,i) = radx*xfn(3,i)/h

!     Position sphere
         xfn(1,i) = (xfn(1,i) + 0.5*lcube/h)
         xfn(2,i) = (xfn(2,i) + 0.380*lcube/h)
         xfn(3,i) = (xfn(3,i) + 0.375*lcube/h)

!     Make xfi identical to xfn
         xfi(1,i) = xfn(1,i)
         xfi(2,i) = xfn(2,i)
         xfi(3,i) = xfn(3,i)
 101     format(13x,e20.13,e20.13,e20.13)
      enddo

!     !m The factor of 2 in this do loop is a relic of an inconsistent
!     standard, so either the code should be made smarter, the mesh gen
!     should be standardized, or this should be clearly documented.

      do i = 1, nfsize2
         read(20,103) elmnew(1,i), elmnew(2,i), elmnew(3,i)
 103     format(3(i8))
!     !m The mesh generator isn't making much sense. All of shpfs
!     ought to be in one file.
         read(21,104) shpfs(1,i), shpfs(2,i), shpfs(3,i), shpfs(7,i)
         read(22,105) shpfs(4,i), shpfs(5,i), shpfs(6,i)
         read(23,105)shpint(1,i), shpint(2,i), shpint(3,i)
      enddo

 104  format(1pe15.8,3x,1pe15.8,3x,1pe15.8,3x,1pe15.8)
 105  format(1pe15.8,3x,1pe15.8,3x,1pe15.8)
 106  format(3(i8))

!     !m This scaling could be done as the file is loaded.
      do i= 1,nfsize2
         shpfs(1,i)=shpfs(1,i)/radx
         shpfs(2,i)=shpfs(2,i)/radx
         shpfs(3,i)=shpfs(3,i)/radx
         shpfs(4,i)=shpfs(4,i)/radx
         shpfs(5,i)=shpfs(5,i)/radx
         shpfs(6,i)=shpfs(6,i)/radx
         shpfs(7,i)=shpfs(7,i)*radx*radx
      enddo

      return
      end subroutine inspher
!************************************************************
      subroutine membnx1(xfi,xfn,frc)
!     Copies xfn from xfi, zeros frc

      implicit none

      integer, parameter :: nfsize=10242
      real xfi(1:3,1:nfsize)
      real xfn(1:3,1:nfsize)
      real frc(1:3,1:nfsize)
      integer i
      do i=1,nfsize
         xfn(1,i)=xfi(1,i)
         xfn(2,i)=xfi(2,i)
         xfn(3,i)=xfi(3,i)
         frc(1,i)=0.
         frc(2,i)=0.
         frc(3,i)=0.
      end do
      return
      end subroutine membnx1
!************************************************************
      subroutine membnx2(n1,n2,clock,xfn,elmnew,shpint,shpfs,frc,&
           h64,fostar,radx)

      implicit none

      interface elm
       subroutine elmfrc(shpfs,ielm,u2,u3,v3,fx1,fy1,fx2,fy2,&
            fx3,fy3,fz1,fz2,fz3,clock)
      implicit none
      integer,  intent (in) :: ielm,clock
      real,dimension(:,:), intent (in) :: shpfs
      real, intent (in) ::u2,u3,v3
      real, intent (out) :: fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3
      end subroutine elmfrc
      end interface elm
      integer nfsize,nfsize2
      parameter(nfsize=10242,nfsize2=20480)
      real lcube,h64,fostar,smarea,smvol,areamn,areamx,area0,vol0
      real xc1,yc1,zc1,dot,vol
      real r11,r12,r13,r21,r22,r23,r31,r32,r33
      real e1m,e2m,e3m,a1,a2,a3,xj1,xj2,xj3
      real e11,e12,e13,e21,e22,e23,e31,e32,e33,et1,et2,et3
      real area,parea,radx
      real x1,x2,x3,y1,y2,y3,z1,z2,z3
      real xcl1,xcl2,xcl3
      real xjl1,xjl2,xjl3,xkl1,xkl2,xkl3
      real fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3 
      real fx11,fx12,fx13,fx21,fx22,fx23,fx31,fx32,fx33
      real fy11,fy12,fy13,fy21,fy22,fy23,fy31,fy32,fy33
      real fz11,fz12,fz13,fz21,fz22,fz23,fz31,fz32,fz33
      real forcex1,forcey1,forcez1
      integer j1,j2,j3,n1,n2,clock
      integer i,icount,iunit1,iunit2
      real xfn(1:3,1:nfsize),frc(1:3,1:nfsize)
      integer elmnew(1:3,1:nfsize2)
      real shpint(1:3,1:nfsize2),shpfs(1:7,1:nfsize2)
      area0 = 12.56131
      vol0 = 4.185322
      smarea = 0.0
      smvol = 0.0
      areamx = -100.
      areamn =  100.
      xc1 = 0.0
      yc1 = 0.0
      zc1 = 0.0
      forcex1=0.
      forcey1=0.
      forcez1=0.

      do i = 1,n1
         frc(1,i) = 0.0
         frc(2,i) = 0.0
         frc(3,i) = 0.0
      end do

      do i = 1,n1
         xc1 = xc1 + xfn(1,i)
         yc1 = yc1 + xfn(2,i)
         zc1 = zc1 + xfn(3,i)
      end do

      xc1 = xc1/float(n1)*h64
      yc1 = yc1/float(n1)*h64
      zc1 = zc1/float(n1)*h64
      if((clock/100)*100 == clock) then
         write(401,222) clock,xc1/h64,yc1/h64,zc1/h64
         write(218,223) -1,-1,-1.
      endif
 222  format(i8,x,e13.6,x,e13.6,x,e13.6)
 223  format(i9,x,i9)

      do i = 1,nfsize2
         j1 = elmnew(1,i)
         j2 = elmnew(2,i)
         j3 = elmnew(3,i)
         x1 = xfn(1,j1)*h64
         y1 = xfn(2,j1)*h64
         z1 = xfn(3,j1)*h64
         x2 = xfn(1,j2)*h64
         y2 = xfn(2,j2)*h64
         z2 = xfn(3,j2)*h64
         x3 = xfn(1,j3)*h64
         y3 = xfn(2,j3)*h64
         z3 = xfn(3,j3)*h64
         e1m = sqrt((-x1 + x2)**2 + (-y1 + y2)**2 + (-z1 + z2)**2)
         r11 =(-x1 + x2)/e1m
         r12 =(-y1 + y2)/e1m
         r13 =(-z1 + z2)/e1m  
         et1 = (x3 - x1)
         et2 = (y3 - y1)
         et3 = (z3 - z1)
         e31 = -r13*et2 + r12*et3
         e32 =  r13*et1 - r11*et3
         e33 =  r11*et2 - r12*et1
         e3m = sqrt(e31**2 + e32**2 + e33**2)
         r31 = e31/e3m
         r32 = e32/e3m
         r33 = e33/e3m
         e21 = -r33*r12 + r32*r13
         e22 =  r33*r11 - r31*r13
         e23 =  r31*r12 - r32*r11
         e2m = sqrt(e21**2 + e22**2 + e23**2)
         r21 = e21/e2m
         r22 = e22/e2m
         r23 = e23/e2m
         a1 = x2 - x1
         a2 = y2 - y1
         a3 = z2 - z1
         xjl1 = r11*a1 + r12*a2 + r13*a3
         xjl2 = r21*a1 + r22*a2 + r23*a3
         xjl3 = r31*a1 + r32*a2 + r33*a3
         a1 = x3 - x1 
         a2 = y3 - y1 
         a3 = z3 - z1 
         xkl1 = r11*a1 + r12*a2 + r13*a3 
         xkl2 = r21*a1 + r22*a2 + r23*a3 
         xkl3 = r31*a1 + r32*a2 + r33*a3 
         a1 = xc1 -x1
         a2 = yc1 -y1
         a3 = zc1 -z1
         xcl1 = r11*a1 + r12*a2 + r13*a3
         xcl2 = r21*a1 + r22*a2 + r23*a3
         xcl3 = r31*a1 + r32*a2 + r33*a3
         dot = 1.0
         if(r33  >  0.0) dot = -1.0
         area = 0.5*xjl1*xkl2
         parea = 2.0*area/shpfs(7,i) 
         vol = area*xcl3/3.0
         smarea = smarea + abs(area)
         smvol = smvol + abs(vol)
         if(parea  >  areamx) areamx = parea
         if(parea .lt. areamn) areamn = parea
         xjl1 = (xjl1 - shpint(1,i))*radx
         xkl1 = (xkl1 - shpint(2,i))*radx
         xkl2 = (xkl2 - shpint(3,i))*radx

         call elmfrc(shpfs,i,xjl1,xkl1,xkl2,fx1,fy1,fx2,fy2,&
              fx3,fy3,fz1,fz2,fz3,clock)
         fx12 =(r11*fx1 + r21*fy1 + r31*fz1)/fostar
         fy12 =(r12*fx1 + r22*fy1 + r32*fz1)/fostar
         fz12 =(r13*fx1 + r23*fy1 + r33*fz1)/fostar
         fx22 =(r11*fx2 + r21*fy2 + r31*fz2)/fostar
         fy22 =(r12*fx2 + r22*fy2 + r32*fz2)/fostar
         fz22 =(r13*fx2 + r23*fy2 + r33*fz2)/fostar
         fx32 =(r11*fx3 + r21*fy3 + r31*fz3)/fostar
         fy32 =(r12*fx3 + r22*fy3 + r32*fz3)/fostar
         fz32 =(r13*fx3 + r23*fy3 + r33*fz3)/fostar
         frc(1,j1) = (frc(1,j1) - fx12)
         frc(2,j1) = (frc(2,j1) - fy12)
         frc(3,j1) = (frc(3,j1) - fz12)
         frc(1,j2) = (frc(1,j2) - fx22)
         frc(2,j2) = (frc(2,j2) - fy22)
         frc(3,j2) = (frc(3,j2) - fz22)
         frc(1,j3) = (frc(1,j3) - fx32)
         frc(2,j3) = (frc(2,j3) - fy32)
         frc(3,j3) = (frc(3,j3) - fz32)
      end do
 100  format(6(e12.5))
      smarea= smarea/area0
      smvol = smvol/vol0

      return
      end subroutine membnx2
!**********************************************************************
      subroutine elmfrc(shpfs, ielm, u2, u3, v3, fx1, fy1, fx2, fy2, &
        fx3, fy3, fz1, fz2, fz3, clock)

      implicit none

      integer, intent (in) :: ielm,clock
      real,dimension(:,:), intent (in) :: shpfs
      real, intent(in) ::u2,u3,v3
      real, intent(out) :: fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3
      integer nfsize,nfsize2
      real u1,v1,v2
      real a0,a1,a2,a3,b1,b2,b3,b,bmu,c,ck
      real sumua,sumub,sumva,sumvb
      real g11,g22,g12,t0,t1,t2,t3,eh,xclock
      real dg11u1,dg11u2,dg11u3,dg11v1,dg11v2,dg11v3
      real dg22u1,dg22u2,dg22u3,dg22v1,dg22v2,dg22v3
      real dg12u1,dg12u2,dg12u3,dg12v1,dg12v2,dg12v3
      real dl1du1,dl1du2,dl1du3,dl2du1,dl2du2,dl2du3
      real dl1dv1,dl1dv2,dl1dv3,dl2dv1,dl2dv2,dl2dv3
      real dwdl1,dwdl2
      real ji1,ji2,jl1,jl2
      u1 = 0.0
      v1 = 0.0
      v2 = 0.0

      a0 = shpfs(7,ielm)
      a1 = shpfs(1,ielm)
      a2 = shpfs(2,ielm)
      a3 = shpfs(3,ielm)
      b1 = shpfs(4,ielm)
      b2 = shpfs(5,ielm)
      b3 = shpfs(6,ielm)

      sumua = u2*a2 + u3*a3
      sumub = u2*b2 + u3*b3
      sumva = v3*a3
      sumvb = v3*b3

      g11 =  1. + 2.*sumua + sumua*sumua + sumva*sumva
      g22 =  1. + 2.*sumvb + sumvb*sumvb + sumub*sumub
      g12 =  sumub + sumub*sumua + sumva + sumvb*sumva  

!     dg11/du
      dg11u1 = 2.*a1 + 2.*sumua*a1
      dg11u2 = 2.*a2 + 2.*sumua*a2
      dg11u3 = 2.*a3 + 2.*sumua*a3
      dg11v1 = 2.*sumva*a1
      dg11v2 = 2.*sumva*a2 
      dg11v3 = 2.*sumva*a3 

!     dg22/du 
      dg22u1 = 2.*sumub*b1 
      dg22u2 = 2.*sumub*b2 
      dg22u3 = 2.*sumub*b3 
      dg22v1 = 2.*b1 + 2.*sumvb*b1 
      dg22v2 = 2.*b2 + 2.*sumvb*b2  
      dg22v3 = 2.*b3 + 2.*sumvb*b3       

!     dg12/du
      dg12u1 = b1 + sumua*b1 + sumub*a1 
      dg12u2 = b2 + sumua*b2 + sumub*a2 
      dg12u3 = b3 + sumua*b3 + sumub*a3 
      dg12v1 = a1 + sumva*b1 + sumvb*a1  
      dg12v2 = a2 + sumva*b2 + sumvb*a2  
      dg12v3 = a3 + sumva*b3 + sumvb*a3

!     dlambda1/du
      t0 = sqrt((g11-g22)**2 + 4.*g12*g12)
      t1 = ((g11 + g22 + t0))
      jl1 = sqrt(0.5*t1)
      jl2 = sqrt(0.5*(g11+g22-t0))
      t2 = dg11u1 - dg22u1 

      if(abs(t0)  >  1.0e-03) then 
         t3  = 0.5/t0*(2.*(g11-g22)*(dg11u1 - dg22u1) &
              + 8.*g12*dg12u1 )
      else       
         t3 = 0.0
      endif

      dl1du1 = 0.5*sqrt(0.5)/sqrt(t1)*(t2+t3)
      dl2du1 = 0.5*sqrt(0.5)/sqrt(t1-2.*t0)*(t2-t3)
      t2 = dg11u2 - dg22u2

      if(abs(t0)  >  1.0e-03) then
         t3  = + 0.5/t0*(2.*(g11-g22)*(dg11u2 - dg22u2) &
              + 8.*g12*dg12u2 )
      else
         t3 =0.0
      endif

      dl1du2 = 0.5*sqrt(0.5)/sqrt(t1)*(t2+t3)
      dl2du2 = 0.5*sqrt(0.5)/sqrt(t1-2.*t0)*(t2-t3)
      t2 = dg11u3 + dg22u3

      if(abs(t0)  >  1.0e-03) then
         t3  = + 0.5/t0*(2.*(g11-g22)*(dg11u3 - dg22u3) &
              + 8.*g12*dg12u3 )
      else
         t3 = 0.0
      endif

      dl1du3 = 0.5*sqrt(0.5)/sqrt(t1)*(t2+t3)
      dl2du3 = 0.5*sqrt(0.5)/sqrt(t1-2.*t0)*(t2-t3)
      t2 = dg11v1 + dg22v1

      if(abs(t0)  >  1.0e-03) then
         t3  = + 0.5/t0*(2.*(g11-g22)*(dg11v1 - dg22v1) &
              + 8.*g12*dg12v1 )
      else
         t3 = 0.0
      endif

      dl1dv1 = 0.5*sqrt(0.5)/sqrt(t1)*(t2+t3)
      dl2dv1 = 0.5*sqrt(0.5)/sqrt(t1-2.*t0)*(t2-t3)
      t2 = dg11v2 + dg22v2

      if(abs(t0)  >  1.0e-03) then
         t3  = + 0.5/t0*(2.*(g11-g22)*(dg11v2 - dg22v2) &
              + 8.*g12*dg12v2 )
      else
         t3 = 0.0
      endif

      dl1dv2 = 0.5*sqrt(0.5)/sqrt(t1)*(t2+t3)
      dl2dv2 = 0.5*sqrt(0.5)/sqrt(t1-2.*t0)*(t2-t3)
      t2 = dg11v3 + dg22v3

      if(abs(t0)  >  1.0e-03) then
         t3  = + 0.5/t0*(2.*(g11-g22)*(dg11v3 - dg22v3) &
              + 8.*g12*dg12v3 )
      else
         t3 = 0.0
      endif

      dl1dv3 = 0.5*sqrt(0.5)/sqrt(t1)*(t2+t3)
      dl2dv3 = 0.5*sqrt(0.5)/sqrt(t1-2.*t0)*(t2-t3)
      b =  9.765699e-03/0.5
      bmu = b
      c = 1.0e+01*b
      ck = c

      if(clock <= 1000) then
         eh = 3.0e+06
      elseif((clock > 1000).and.(clock <= 50000)) then
         xclock=(real(clock-1000))/49000.
         eh=3.0e05*(xclock+10.*(1.-xclock))
      elseif(clock > 50000) then
         eh = 3.0e+05
      endif
      ji1 = jl1*jl1 + jl2*jl2 - 2.
      ji2 = jl1*jl1*jl2*jl2 - 1.
!     strain energy for skalak red blood cell membrane
!     skalak, et al. j. o biophys. 1973
!     
!     dwdl1 = (b/2.)*(ji1*jl1 + jl1 - jl1*jl2*jl2)
!     &     + (c/2.)*ji2*jl1*jl2**2
!     
!     dwdl2 = (b/2.)*(ji1*jl2 + jl2 - jl1*jl1*jl2)
!     &     + (c/2.)*ji2*jl2*jl1**2
!     
!     strain energy for barthes-besiel elastic membrane under
!     small deformations jfm 1981
!     
!     Why is strain energy time dependent?
!
      dwdl1 = -(eh)/(3.*jl1) + (eh*jl1)/3.  &
           + (eh)/(3.*jl1)*log(jl1**2*jl2**2)
!     
      dwdl2 = -(eh)/(3.*jl2) + (eh*jl2)/3.  &
          + (eh)/(3.*jl2)*log(jl1**2*jl2**2)
!     
!     
!     strain energy function from evans and skalak
!     mech and therm of biomembranes 1980 pg 98 (4.8.1)
!     
!     dwdl1 = ck*(jl1*jl2 - 1.0)*jl2
!     & + 0.5*bmu*(jl1*jl1 - jl2*jl2)/(jl1*jl1*jl2)
!     
!     dwdl2 = ck*(jl1*jl2 - 1.0)*jl1 
!     & - 0.5*bmu*(jl1*jl1 - jl2*jl2)/(jl1*jl2*jl2)
!     
!     
!     linearzid evans and skalak strain energy for barthes-besiel invariants under
!     small deformations 
!     
!     dwdl1 = -(bmu)/(jl1) + (bmu*jl1)
!     &          + (ck-bmu)/(jl1)*log(jl1*jl2)
!     
!     dwdl2 = -(bmu)/(jl2) + (bmu*jl2)
!     &           + (ck-bmu)/(jl2)*log(jl1*jl2)
!     
!     calculate nodal forces

      fx1 =(dwdl1*dl1du1 + dwdl2*dl2du1)*0.5*a0
      fx2 =(dwdl1*dl1du2 + dwdl2*dl2du2)*0.5*a0
      fx3 =(dwdl1*dl1du3 + dwdl2*dl2du3)*0.5*a0
      fy1 =(dwdl1*dl1dv1 + dwdl2*dl2dv1)*0.5*a0
      fy2 =(dwdl1*dl1dv2 + dwdl2*dl2dv2)*0.5*a0
      fy3 =(dwdl1*dl1dv3 + dwdl2*dl2dv3)*0.5*a0
      fz1 = 0.0
      fz2 = 0.0
      fz3 = 0.0
      return
      end subroutine elmfrc

