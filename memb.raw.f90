!**********************************************************************
SUBROUTINE INHIST(XFN,FIRSTN,NUMBER,NEXTN)
  IMPLICIT NONE
  INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ,NFSIZE
  INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
  double precision :: FLNGX,FLNGY,FLNGZ
  PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
  PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
  PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
  PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
  PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)

  INTEGER N,II,JJ
  double precision :: xfn(:,:)
  integer firstn(:,:),number(:,:),nextn(:)       

  nfsize = size(xfn, 2)

  DO N=1,NFSIZE
     NEXTN(N)=0
  end do
  DO JJ=1,NGY
     DO II=1,NGX
        FIRSTN(II,JJ)=0
        NUMBER(II,JJ)=0
     end do
  end do
  DO N=1,NFSIZE
     JJ=MOD(INT(XFN(2,N)+FLNGY)-1+NGY,NGY) + 1
     II=MOD(INT(XFN(1,N)+FLNGX)-1+NGX,NGX) + 1
     NEXTN(N)=FIRSTN(II,JJ)
     FIRSTN(II,JJ)=N
     NUMBER(II,JJ)=NUMBER(II,JJ)+1
  end do
  RETURN
END SUBROUTINE INHIST
!************************************************************
subroutine generatecapsule(RAD,H,XFN, elmnew,shpint,shpfs, &
     my_cap_center, my_fineness, cap_i, a_prestress, my_nvec_i)
  IMPLICIT NONE
  interface 
     subroutine sph(fineness, xfnew, ilmnew)
       integer fineness
       double precision :: xfnew(:,:)
       integer :: ilmnew(:,:)
     end subroutine sph

     subroutine sf(xfn, elmnew, shpint, shpfs)
       double precision :: xfn(:,:), rad, h
       integer elmnew(:,:)
       double precision :: shpint(:,:), shpfs(:,:)
     end subroutine sf
  end interface

  integer, parameter :: lxng=$lngx$,lyng=lxng,lzng=lxng
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  double precision, parameter :: flngx=ngx,flngy=ngy,flngz=ngz
  INTEGER i 
  double precision H,pi,x,y,rad, a_prestress
  double precision :: XFN(:,:)
  INTEGER elmnew(:,:)
  double precision :: shpint(:,:), shpfs(:,:) 
  double precision :: my_cap_center(:)
  integer my_fineness, cap_i
  integer my_nvec_i(:)
  character(len=12) strfname

  pi = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits

  !     The sphere comes in as the unit sphere.
  call sph(my_fineness, xfn, elmnew)
  call sf(xfn, elmnew, shpint, shpfs)

  do i = 1, size(xfn, 2)
     XFN(1,i) =RAD*xfn(1,i)
     XFN(2,i) =RAD*xfn(2,i)
     XFN(3,i) =RAD*XFN(3,i)
  end do

  call sf(xfn, elmnew, shpint, shpfs)

  do i = 1, 3
     ! Scale the capsule to the grid size, and position it in the computational
     ! domain
     XFN(i,:) =xfn(i,:)/H + my_cap_center(i)
  enddo

  my_nvec_i(1:3) = maxloc(xfn, 2)
  my_nvec_i(4:6) = minloc(xfn, 2)

  !     This scales the sphere's (non-dimensional) finite element parameters
  !     to real units.
  do i= 1, size(elmnew,2)
     shpfs(1,i)=shpfs(1,i)/(1+a_prestress)
     shpfs(2,i)=shpfs(2,i)/(1+a_prestress)
     shpfs(3,i)=shpfs(3,i)/(1+a_prestress)
     shpfs(4,i)=shpfs(4,i)/(1+a_prestress)
     shpfs(5,i)=shpfs(5,i)/(1+a_prestress)
     shpfs(6,i)=shpfs(6,i)/(1+a_prestress)
! shpfs(7,:) is probably the element in shpfs/shpint that would be most useful
! to someone writing analysis code. This number is twice the area of the
! element. It's twice the area of the element, rather than just the area of the
! element, for historical reasons.
     shpfs(7,i)=shpfs(7,i)*rad*rad/(1+a_prestress)/(1+a_prestress)
     !     The following 3 lines added 5-2-07 due to changes in scaling in
     !     membnx.
     shpint(1,i) = shpint(1,i)/(1+a_prestress)
     shpint(2,i) = shpint(2,i)/(1+a_prestress)
     shpint(3,i) = shpint(3,i)/(1+a_prestress)
  enddo


  write(strfname,'(a4,i4,a4)') 'capx', cap_i, '.txt'
  call padzeros(strfname)
  open(400,file=strfname,status='unknown')
  do i=1,size(xfn,2)
     write(400,'(es24.17,2(x,es24.17))')xfn(1,i),xfn(2,i),xfn(3,i)
  end do
  close(400)

  write(strfname,'(a4,i4,a4)') 'elm_', cap_i, '.txt'
  call padzeros(strfname)
  open(400, file=strfname)
  do i = 1, size(elmnew,2)
     write(400, '(3(i7))') elmnew(1, i), elmnew(2, i), elmnew(3, i)
  end do
  close(400)

  write(strfname,'(a4,i4,a4)') 'shpa', cap_i, '.txt'
  call padzeros(strfname)
  open(400,file=strfname,status='unknown')
  write(strfname,'(a4,i4,a4)') 'shpb', cap_i, '.txt'
  call padzeros(strfname)
  open(401,file=strfname,status='unknown')
  write(strfname,'(a4,i4,a4)') 'shpi', cap_i, '.txt'
  call padzeros(strfname)
  open(402,file=strfname,status='unknown')
  do i = 1, size(shpfs, 2)
     write(400,'(es24.17,3(x,es24.17))') shpfs(1,i), &
          shpfs(2,i), shpfs(3,i), shpfs(7,i)
     write(401,'(es24.17,2(x,es24.17))') shpfs(4,i), shpfs(5,i), &
          shpfs(6,i)
     write(402,'(es24.17,2(x,es24.17))') shpint(1,i), shpint(2,i), &
          shpint(3,i)
  end do
  close(400)
  close(401)
  close(402)

  return
end subroutine generatecapsule
!************************************************************
subroutine MEMBNX(XFN,elmnew,shpint,shpfs,FRC, H,FOSTAR, lambda1, lambda2)
  IMPLICIT NONE
  INTERFACE ELM
     subroutine elmfrc(shpfs,ielm,u2,u3,v3,fx1,fy1,fx2,fy2, &
          fx3,fy3,fz1,fz2,fz3, jl1, jl2)
       IMPLICIT NONE
       INTEGER,  INTENT (IN) :: ielm
       double precision,DIMENSION(:,:), INTENT (IN) :: shpfs
       double precision, INTENT (IN) :: u2,u3,v3
       double precision, INTENT (OUT) :: fx1,fx2,fx3,fy1,fy2,fy3,fz1, &
            fz2,fz3, jl1, jl2
     END subroutine elmfrc
  END INTERFACE

  double precision :: H,FOSTAR
  double precision :: R11,R12,R13,R21,R22,R23,R31,R32,R33
  double precision :: e1m,e2m,e3m,a1,a2,a3
  double precision :: e21,e22,e23,e31,e32,e33,et1,et2, et3
  double precision :: x1,x2,x3,y1,y2,y3,z1,z2,z3
  double precision :: xjl1,xkl1,xkl2
  double precision :: fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3 
  double precision :: fx12,fx22,fx32
  double precision :: fy12,fy22,fy32
  double precision :: fz12,fz22,fz32
  double precision :: jl1, jl2
  double precision :: FORCEX1,FORCEY1,FORCEZ1
  double precision :: FORCEX2,FORCEY2,FORCEZ2
  integer j1,j2,j3
  integer i


  double precision :: XFN(:,:),FRC(:,:)
  double precision :: lambda1(:), lambda2(:)
  INTEGER elmnew(:,:)
  double precision :: shpint(:,:),shpfs(:,:)

  integer nnode, nelm

  nnode = size(xfn,2)
  nelm = size(elmnew,2)

  FORCEX1=0.d0
  FORCEY1=0.d0
  FORCEZ1=0.d0
  FORCEX2=0.d0
  FORCEY2=0.d0
  FORCEZ2=0.d0

  do i = 1, nnode
     FRC(1,i) = 0.0d0
     FRC(2,i) = 0.0d0
     FRC(3,i) = 0.0d0
  end do

  do i = 1, nelm
     j1 = elmnew(1,i)
     j2 = elmnew(2,i)
     j3 = elmnew(3,i)
     !     Changed these lines 5-2-07 so that coordinates would be in real 
     !     units (as opposed to non-dimensional, scaling by cell radius) 
     !     consistently in this subroutine. The original scaling was by 
     !     (h/radx); now it is simply by h. This means that, although xfn is in
     !     program units, x1, y1, ... are in real units.
     !     Node 1
     x1 = xfn(1,j1)*h
     y1 = xfn(2,j1)*h
     z1 = xfn(3,j1)*h
     !     Node 2
     x2 = xfn(1,j2)*h
     y2 = xfn(2,j2)*h
     z2 = xfn(3,j2)*h
     !     Node 3
     x3 = xfn(1,j3)*h
     y3 = xfn(2,j3)*h
     z3 = xfn(3,j3)*h

     e1m = dsqrt((-X1 + X2)**2 + (-Y1 + Y2)**2 + (-Z1 + Z2)**2)

     !     These R's are the rotation matrix entries
     R11 =(-X1 + X2)/e1m
     R12 =(-Y1 + Y2)/e1m
     R13 =(-Z1 + Z2)/e1m  
     et1 = (X3 - X1)
     et2 = (Y3 - Y1)
     et3 = (Z3 - Z1)
     e31 = R12*et3 - R13*et2
     e32 = R13*et1 - R11*et3
     e33 = R11*et2 - R12*et1
     e3m = dsqrt(e31**2 + e32**2 + e33**2)
     R31 = e31/e3m
     R32 = e32/e3m
     R33 = e33/e3m
     e21 = R32*R13 - R33*R12
     e22 = R33*R11 - R31*R13
     e23 = R31*R12 - R32*R11
     e2m = dsqrt(e21**2 + e22**2 + e23**2)
     R21 = e21/e2m
     R22 = e22/e2m
     R23 = e23/e2m
     a1 = X2 - X1
     a2 = Y2 - Y1
     a3 = Z2 - Z1
     xjl1 = R11*a1 + R12*a2 + R13*a3
     !     Changed this line 5-2-07 so that coordinates would be in real units
     !     consistently in this subroutine.
     !     shpint values are subtracted here so elmfrc doesn't get coordinates,
     !     it gets displacements
     xjl1 = xjl1 - shpint(1,i)
     a1 = X3 - X1 
     a2 = Y3 - Y1 
     a3 = Z3 - Z1 
     xkl1 = R11*a1 + R12*a2 + R13*a3 
     xkl2 = R21*a1 + R22*a2 + R23*a3 
     !     Changed the following two lines 5-2-07 so that coordinates would be in
     !     real units consistently in this subroutine.
     xkl1 = xkl1 - shpint(2,i)
     xkl2 = xkl2 - shpint(3,i)

     call elmfrc(shpfs,i,xjl1,xkl1,xkl2,fx1,fy1,fx2,fy2, &
          fx3,fy3,fz1,fz2,fz3, jl1, jl2)
     !     Rotate forces from the 2-d plane to the actual orientation of the
     !     element; scale force from cgs to program units.
     fx12 =(R11*fx1 + R21*fy1 + R31*fz1)/FOSTAR
     fy12 =(R12*fx1 + R22*fy1 + R32*fz1)/FOSTAR
     fz12 =(R13*fx1 + R23*fy1 + R33*fz1)/FOSTAR
     fx22 =(R11*fx2 + R21*fy2 + R31*fz2)/FOSTAR
     fy22 =(R12*fx2 + R22*fy2 + R32*fz2)/FOSTAR
     fz22 =(R13*fx2 + R23*fy2 + R33*fz2)/FOSTAR
     fx32 =(R11*fx3 + R21*fy3 + R31*fz3)/FOSTAR
     fy32 =(R12*fx3 + R22*fy3 + R32*fz3)/FOSTAR
     fz32 =(R13*fx3 + R23*fy3 + R33*fz3)/FOSTAR

     FRC(1,j1) = (FRC(1,j1) - fx12)
     FRC(2,j1) = (FRC(2,j1) - fy12)
     FRC(3,j1) = (FRC(3,j1) - fz12)
     FRC(1,j2) = (FRC(1,j2) - fx22)
     FRC(2,j2) = (FRC(2,j2) - fy22)
     FRC(3,j2) = (FRC(3,j2) - fz22)
     FRC(1,j3) = (FRC(1,j3) - fx32)
     FRC(2,j3) = (FRC(2,j3) - fy32)
     FRC(3,j3) = (FRC(3,j3) - fz32)
     lambda1(i) = jl1
     lambda2(i) = jl2
  end do
  return
end subroutine MEMBNX
!**********************************************************************
subroutine elmfrc(shpfs,ielm,u2,u3,v3,fx1,fy1,fx2,fy2, &
     fx3,fy3,fz1,fz2,fz3, jl1, jl2)
  IMPLICIT NONE
  INTEGER,  INTENT (IN) :: ielm
  double precision, DIMENSION(:,:), INTENT (IN) :: shpfs
  double precision, INTENT (IN) :: u2,u3,v3
  double precision, INTENT (OUT) :: fx1,fx2,fx3,fy1,fy2,fy3,fz1, &
       fz2,fz3
  double precision :: u1,v1,v2
  double precision :: a0,a1,a2,a3,b1,b2,b3
  double precision :: sumua,sumub,sumva,sumvb
  double precision :: g11,g22,g12,t0,t1,t2,t3,eh
  double precision :: dg11u1,dg11u2,dg11u3,dg11v1,dg11v2,dg11v3
  double precision :: dg22u1,dg22u2,dg22u3,dg22v1,dg22v2,dg22v3
  double precision :: dg12u1,dg12u2,dg12u3,dg12v1,dg12v2,dg12v3
  double precision :: dl1du1,dl1du2,dl1du3,dl2du1,dl2du2,dl2du3
  double precision :: dl1dv1,dl1dv2,dl1dv3,dl2dv1,dl2dv2,dl2dv3
  double precision :: dwdl1,dwdl2
  double precision :: ji1,ji2,jl1,jl2

  !     u is displacement in x direction, v in y, in the 2-D coordinate
  !     system with point 1 of the element at the origin, point 2 on the
  !     +x axis, and point 3 in the first quadrant. Thus, the next three
  !     lines are given. None of these variables are used, though, and so
  !     they should be removed.
  u1 = 0.0d0
  v1 = 0.0d0
  v2 = 0.0d0

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

  !     Entries for right Cauchy strain tensor
  g11 = 1.d0 + 2.d0*sumua + sumua*sumua + sumva*sumva
  g22 = 1.d0 + 2.d0*sumvb + sumvb*sumvb + sumub*sumub
  g12 = sumub + sumub*sumua + sumva + sumvb*sumva  

  !     These derivatives let us use strain energy gradients to calculate
  !     force
  !     dg11/du
  dg11u1 = 2.d0*a1 + 2.d0*sumua*a1
  dg11u2 = 2.d0*a2 + 2.d0*sumua*a2
  dg11u3 = 2.d0*a3 + 2.d0*sumua*a3
  dg11v1 = 2.d0*sumva*a1
  dg11v2 = 2.d0*sumva*a2 
  dg11v3 = 2.d0*sumva*a3 

  !     dg22/du 
  dg22u1 = 2.d0*sumub*b1 
  dg22u2 = 2.d0*sumub*b2 
  dg22u3 = 2.d0*sumub*b3 
  dg22v1 = 2.d0*b1 + 2.d0*sumvb*b1 
  dg22v2 = 2.d0*b2 + 2.d0*sumvb*b2  
  dg22v3 = 2.d0*b3 + 2.d0*sumvb*b3       

  !     dg12/du
  dg12u1 = b1 + sumua*b1 + sumub*a1 
  dg12u2 = b2 + sumua*b2 + sumub*a2 
  dg12u3 = b3 + sumua*b3 + sumub*a3 
  dg12v1 = a1 + sumva*b1 + sumvb*a1  
  dg12v2 = a2 + sumva*b2 + sumvb*a2  
  dg12v3 = a3 + sumva*b3 + sumvb*a3

  !     These tractions, t0, t1, t2, t3, are used to calculate
  !     the derivatives of the stretch ratios with respect to
  !     the displacements. I think that the tractions correspond to
  !     the invariants of the right Cauchy strain tensor, but need to verify
  !     this.
  !     d Lambda1 / d u
  t0 = dsqrt((g11-g22)**2 + 4.d0*g12*g12)
  t1 = ((g11 + g22 + t0))
  jl1 = dsqrt(0.5d0*t1)
  jl2 = dsqrt(0.5d0*(g11+g22-t0))

  t2 = dg11u1 - dg22u1 
  if(dabs(t0) .gt. 1.0d-03) then 
     t3 = 0.5d0/t0*(2.d0*(g11-g22)*(dg11u1 - dg22u1) + &
          8.d0*g12*dg12u1)
  else       
     t3 = 0.0d0
  endif
  dl1du1 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1)*(t2+t3)
  dl2du1 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1-2.d0*t0)*(t2-t3)

  t2 = dg11u2 - dg22u2
  if(dabs(t0) .gt. 1.0d-03) then
     t3 = 0.5d0/t0*(2.d0*(g11-g22)*(dg11u2 - dg22u2) + &
          8.d0*g12*dg12u2)
  else
     t3 =0.0d0
  endif
  dl1du2 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1)*(t2+t3)
  dl2du2 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1-2.d0*t0)*(t2-t3)

  t2 = dg11u3 + dg22u3
  if(dabs(t0) .gt. 1.0d-03) then
     t3  = + 0.5d0/t0*(2.d0*(g11-g22)*(dg11u3 - dg22u3) &
          + 8.d0*g12*dg12u3 )
  else
     t3 = 0.0d0
  endif
  dl1du3 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1)*(t2+t3)
  dl2du3 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1-2.d0*t0)*(t2-t3)

  t2 = dg11v1 + dg22v1
  if(dabs(t0) .gt. 1.0d-03) then
     t3  = + 0.5d0/t0*(2.d0*(g11-g22)*(dg11v1 - dg22v1) &
          + 8.d0*g12*dg12v1 )
  else
     t3 = 0.0d0
  endif
  dl1dv1 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1)*(t2+t3)
  dl2dv1 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1-2.d0*t0)*(t2-t3)

  t2 = dg11v2 + dg22v2
  if(dabs(t0) .gt. 1.0d-03) then
     t3  = + 0.5d0/t0*(2.d0*(g11-g22)*(dg11v2 - dg22v2) &
          + 8.d0*g12*dg12v2 )
  else
     t3 = 0.0d0
  endif
  dl1dv2 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1)*(t2+t3)
  dl2dv2 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1-2.d0*t0)*(t2-t3)

  t2 = dg11v3 + dg22v3
  if(dabs(t0) .gt. 1.0d-03) then
     t3  = + 0.5d0/t0*(2.d0*(g11-g22)*(dg11v3 - dg22v3) &
          + 8.d0*g12*dg12v3 )
  else
     t3 = 0.0d0
  endif
  dl1dv3 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1)*(t2+t3)
  dl2dv3 = 0.5d0*dsqrt(0.5d0)/dsqrt(t1-2.d0*t0)*(t2-t3)

  Eh = $Eh$
  ji1 = jl1*jl1 + jl2*jl2 - 2.d0
  ji2 = jl1*jl1*jl2*jl2 - 1.d0
  !     strain energy for skalak red blood cell membrane
  !     Skalak, et al. J. o Biophys. 1973
  !     
  !     dwdl1 = (b/2.d0)*(ji1*jl1 + jl1 - jl1*jl2*jl2)
  !     &     + (c/2.d0)*ji2*jl1*jl2**2
  !     
  !     dwdl2 = (b/2.d0)*(ji1*jl2 + jl2 - jl1*jl1*jl2)
  !     &     + (c/2.d0)*ji2*jl2*jl1**2
  !     
  !     strain energy for Barthes-Besiel elastic membrane under
  !     small deformations JFM 1981
  !     
  !      dwdl1 = -(Eh)/(3.d0*jl1) + (Eh*jl1)/3.d0
  !     &     + (Eh)/(3.d0*jl1)*dlog(jl1**2*jl2**2)
  !     
  !      dwdl2 = -(Eh)/(3.d0*jl2) + (Eh*jl2)/3.d0
  !     &     + (Eh)/(3.d0*jl2)*dlog(jl1**2*jl2**2)
  !     
  !     
  !     strain energy function from evans and skalak
  !     mech and therm of biomembranes 1980 pg 98 (4.8.1)
  !     
  !     dwdl1 = cK*(jl1*jl2 - 1.d0)*jl2
  !     & + 0.5d0*bmu*(jl1*jl1 - jl2*jl2)/(jl1*jl1*jl2)
  !     
  !     dwdl2 = cK*(jl1*jl2 - 1.d0)*jl1 
  !     & - 0.5d0*bmu*(jl1*jl1 - jl2*jl2)/(jl1*jl2*jl2)
  !     
  !     
  !     linearized Evans and Skalak strain energy for Barthes-Besiel
  !     invariants under small deformations 
  !     
  !      dwdl1 = -(bmu)/(jl1) + (bmu*jl1)
  !     &     + (cK-bmu)/(jl1)*log(jl1*jl2)
  !     
  !      dwdl2 = -(bmu)/(jl2) + (bmu*jl2)
  !     &           + (cK-bmu)/(jl2)*log(jl1*jl2)
  !     
  !     
  !     Neo-Hookean -- added 5-24-07 by Alex
  dwdl1 = (Eh/3.d0)*(jl1-jl1**(-3)*jl2**(-2))
  dwdl2 = (Eh/3.d0)*(jl2-jl1**(-2)*jl2**(-3))
  !
  !     write(3,*) dwdl1,dwdl2
  !     
  !     Calculate nodal forces in cgs units.
  !     
  ! a0 is twice the area of the element, so the factor of 0.5 appears here
  ! to give the area of the triangular element. See comment in 
  fx1 =(dwdl1*dl1du1 + dwdl2*dl2du1)*0.5d0*a0
  fx2 =(dwdl1*dl1du2 + dwdl2*dl2du2)*0.5d0*a0
  fx3 =(dwdl1*dl1du3 + dwdl2*dl2du3)*0.5d0*a0
  fy1 =(dwdl1*dl1dv1 + dwdl2*dl2dv1)*0.5d0*a0
  fy2 =(dwdl1*dl1dv2 + dwdl2*dl2dv2)*0.5d0*a0
  fy3 =(dwdl1*dl1dv3 + dwdl2*dl2dv3)*0.5d0*a0
  fz1 = 0.0d0
  fz2 = 0.0d0
  fz3 = 0.0d0
  return
end subroutine elmfrc
!**********************************************************************
subroutine capsuletable(fineness, nnode, nelm)
  integer fineness, nnode, nelm
  nnode = 10*4**fineness+2
  nelm = 20*4**fineness
end subroutine capsuletable
!**********************************************************************
subroutine rectangle_table(my_rect_n1, my_rect_n2, my_rect_nnode)
  integer my_rect_n1, my_rect_n2, my_rect_nnode
  my_rect_nnode = my_rect_n1*my_rect_n2
end subroutine rectangle_table
!**********************************************************************
subroutine make_index_table_start_end(n, start, end)
  implicit none
  integer i
  integer n(:), start(:), end(:)

  end(1)=n(1) + start(1) - 1

  if (size(n) > 1) then
     do i=2,size(n)
        start(i) = end(i-1)+1
        end(i) = start(i) - 1 + n(i)
     end do
  end if
end subroutine make_index_table_start_end

