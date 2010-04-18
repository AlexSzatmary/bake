SUBROUTINE RESTART(LCUBE,NU,RHO,TD,UR,VR,WR, pr, &
     XFN,xpi,FIRSTN,NUMBER,NEXTN,elmnew,shpint,shpfs, xcenterold, ycenterold, &
     zcenterold, nsph, nellip, ellipa, ellipb, &
          ellipc, a_prestress, nvec_i)
  IMPLICIT NONE
  INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
  INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
  double precision :: FLNGX,FLNGY,FLNGZ
  parameter(lxng=$lngx$,lyng=lxng,lzng=lxng)
  PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
  PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
  PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
  PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
  double precision :: lcube,nu,td,rho
  double complex :: ur(0:,0:,0:)
  double complex :: vr(0:,0:,0:)
  double complex :: wr(0:,0:,0:)
  double complex :: pr(0:,0:,0:)
  integer firstn(:,:),number(:,:),nextn(:)
  double precision :: xfn(:,:)
  double precision :: xpi(:,:)
  integer elmnew(:,:)
  double precision :: shpint(:,:),shpfs(:,:)
  double precision :: xcenterold(:), ycenterold(:), zcenterold(:)
  integer nsph, nellip
  double precision :: ellipa(:), ellipb(:), ellipc(:), a_prestress(:)
  integer :: nvec_i(:,:)

  READ(101)LCUBE,NU,RHO,TD
  REWIND(101)
  READ(102) UR
  REWIND(102)
  READ(103) VR
  REWIND(103)
  READ(104) WR
  REWIND(104)
  READ(105) FIRSTN 
  REWIND(105)
  READ(106) NUMBER
  REWIND(106)
  READ(107) NEXTN
  REWIND(107)
  READ(109) elmnew
  REWIND(109)
  READ(111) shpfs
  REWIND(111)
  READ(110) shpint
  REWIND(110)
  READ(108) XFN
  REWIND(108)
  read(112) xpi
  rewind(112)
  read(113) xcenterold, ycenterold, zcenterold
  rewind(113)
  read(114) pr
  rewind(114)
  read(115) a_prestress
  rewind(115)
  read(116) nsph, nellip
  rewind(116)
  read(117) ellipa, ellipb, ellipc
  rewind(117)
  read(118) nvec_i
  rewind(118)
END SUBROUTINE RESTART
!**********************************************************************
SUBROUTINE WRSTART(LCUBE,NU,RHO,TD,KLOK,UR,VR,WR, pr, &
     XFN,xpi,FIRSTN,NUMBER,NEXTN,elmnew,shpint,shpfs, xcenterold, ycenterold, &
     zcenterold, nsph, nellip, ellipa, ellipb, &
          ellipc, a_prestress, nvec_i)
  IMPLICIT NONE
!todo Get rid of these magic parameters, they're not used.
  INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
  INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
  double precision :: FLNGX,FLNGY,FLNGZ
  parameter(lxng=$lngx$,lyng=lxng,lzng=lxng)
  PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
  PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
  PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
  PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
  double precision :: LCUBE,NU,TD,RHO
  INTEGER KLOK
  double COMPLEX :: UR(0:,0:,0:)
  double complex :: VR(0:,0:,0:)
  double COMPLEX :: WR(0:,0:,0:)
  double complex :: pr(0:,0:,0:)
  integer firstn(:,:),number(:,:),nextn(:)
  double precision :: xpi(:,:)
  double precision :: xfn(:,:)
  integer elmnew(:,:)
  double precision :: shpint(:,:),shpfs(:,:)
  double precision :: xcenterold(:), ycenterold(:), zcenterold(:)
  integer nsph, nellip
  double precision :: ellipa(:), ellipb(:), ellipc(:), a_prestress(:)
  integer :: nvec_i(:,:)

  open(100, position='rewind', form='unformatted')
  write(100) klok
  rewind(100)
  close(100)
  WRITE(101)LCUBE,NU,RHO,TD
  REWIND(101)
  WRITE(102) UR
  REWIND(102)
  WRITE(103) VR
  REWIND(103)
  WRITE(104) WR
  REWIND(104)
  WRITE(105) FIRSTN 
  REWIND(105)
  WRITE(106) NUMBER
  REWIND(106)
  WRITE(107) NEXTN
  REWIND(107)
  WRITE(108) XFN
  REWIND(108)
  WRITE(109) elmnew
  REWIND(109)
  WRITE(110) shpint
  REWIND(110)
  WRITE(111) shpfs
  REWIND(111)
  write(112) xpi
  rewind(112)
  write(113) xcenterold, ycenterold, zcenterold
  rewind(113)
  write(114) pr
  rewind(114)
  write(115) a_prestress
  rewind(115)
  write(116) nsph, nellip
  rewind(116)
  write(117) ellipa, ellipb, ellipc
  rewind(117)
  write(118) nvec_i
  rewind(118)
END SUBROUTINE WRSTART
!**********************************************************************

