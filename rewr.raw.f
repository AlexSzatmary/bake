      SUBROUTINE RESTART(LCUBE,NU,RHO,TD,UR,VR,WR,
     & XFN,xpi,FIRSTN,NUMBER,NEXTN,elmnew,shpint,shpfs)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ,NFSIZE,NFSIZE2
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      double precision :: FLNGX,FLNGY,FLNGZ
      parameter(lxng=$lngx$,lyng=lxng,lzng=lxng)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      parameter(nfsize=$nsnode$,nfsize2=$nselm$)
      double precision :: LCUBE,NU,TD,RHO
      double COMPLEX :: UR(0:NBX,0:NBY,0:NGZM1)
      double complex :: VR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: WR(0:NBX,0:NBY,0:NGZM1)
      INTEGER FIRSTN(1:NGX,1:NGY),NUMBER(1:NGX,1:NGY),NEXTN(1:NFSIZE)
      double precision :: XFN(1:3,1:NFSIZE)
      integer, parameter :: npl=$npl$
      double precision :: xpi(1:3,1:npl)
      INTEGER elmnew(1:3,1:NFSIZE2)
      double precision :: shpint(1:3,1:NFSIZE2),shpfs(1:7,1:NFSIZE2)
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
      END SUBROUTINE RESTART
!**********************************************************************
      SUBROUTINE WRSTART(LCUBE,NU,RHO,TD,KLOK,UR,VR,WR,
     & XFN,xpi,FIRSTN,NUMBER,NEXTN,elmnew,shpint,shpfs)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ,NFSIZE,NFSIZE2
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      double precision :: FLNGX,FLNGY,FLNGZ
      parameter(lxng=$lngx$,lyng=lxng,lzng=lxng)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      parameter(nfsize=$nsnode$,nfsize2=$nselm$)
      double precision :: LCUBE,NU,TD,RHO
      INTEGER KLOK
      double COMPLEX :: UR(0:NBX,0:NBY,0:NGZM1)
      double complex :: VR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: WR(0:NBX,0:NBY,0:NGZM1)
      INTEGER FIRSTN(1:NGX,1:NGY),NUMBER(1:NGX,1:NGY),NEXTN(1:NFSIZE)
      integer, parameter :: npl=$npl$
      double precision :: xpi(1:3,1:npl)
      double precision :: XFN(1:3,1:NFSIZE)
      INTEGER elmnew(1:3,1:NFSIZE2)
      double precision :: shpint(1:3,1:NFSIZE2),shpfs(1:7,1:NFSIZE2)
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
      END SUBROUTINE WRSTART
!**********************************************************************

