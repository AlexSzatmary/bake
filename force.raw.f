!*********************************************************************
      subroutine bodyfs(bfs,ur, vr, wr, vsc)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      integer i,j,k
      double precision :: bfs(3,3), vsc
      double COMPLEX :: UR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: VR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: WR(0:NBX,0:NBY,0:NGZM1)

      do j = 0, nby
         do k = 0, ngz - 1
            ur(ngx,j,k) = ur(ngx,j,k) + bfs(1,1)*vsc*flngx
            ur(1,j,k) = ur(1,j,k) - bfs(1,1)*vsc*flngx
            vr(ngx,j,k) = vr(ngx,j,k) + bfs(2,1)*vsc*flngx
            vr(1,j,k) = vr(1,j,k) - bfs(2,1)*vsc*flngx
            wr(ngx,j,k) = wr(ngx,j,k) + bfs(3,1)*vsc*flngx
            wr(1,j,k) = wr(1,j,k) - bfs(3,1)*vsc*flngx
         end do
      end do

      do i=0,nbx
         do k=0,ngzm1
            ur(i,ngy,k) = ur(i,ngy,k) + bfs(1,2)*vsc*flngy
            ur(i,1,k) = ur(i,1,k) - bfs(1,2)*vsc*flngy
            vr(i,ngy,k) = vr(i,ngy,k) + bfs(2,2)*vsc*flngy
            vr(i,1,k) = vr(i,1,k) - bfs(2,2)*vsc*flngy
            wr(i,ngy,k) = wr(i,ngy,k) + bfs(3,2)*vsc*flngy
            wr(i,1,k) = wr(i,1,k) - bfs(3,2)*vsc*flngy
         end do
      end do

      do i=0,nbx
         do j=0,nby
            ur(i,j,ngz-1) = ur(i,j,ngz-1) + bfs(1,3)*vsc*flngz
            ur(i,j,0) = ur(i,j,0) - bfs(1,3)*vsc*flngz
            vr(i,j,ngz-1) = vr(i,j,ngz-1) + bfs(2,3)*vsc*flngz
            vr(i,j,0) = vr(i,j,0) - bfs(2,3)*vsc*flngz
            wr(i,j,ngz-1) = wr(i,j,ngz-1) + bfs(3,3)*vsc*flngz
            wr(i,j,0) = wr(i,j,0) - bfs(3,3)*vsc*flngz
         end do
      end do

      RETURN
      END SUBROUTINE BODYFS
!*********************************************************************
      SUBROUTINE wrap(UR,VR,WR)
!     This subroutine wraps the velocities in the x and y directions.
!     The velocity arrays are dimensioned from 0 to ng+2, instead of
!     from 1 to ng. This is because the discrete Dirac delta has a
!     width of 4 nodes (two in each direction). When the discrete Dirac
!     delta is evaluated, it may need points outside of the flow field.
!     Rather than doing complicated and fragile mod arithmetic, we
!     simply say that u(0) = u(ng), u(ng+1)=u(1), and u(ng+2)=u(2)
!     by periodicity. This doesn't actually make sense, I suppose,
!     for cases in which the flow field is nonperiodic, but it keeps
!     the code running.
!     Commented by Alex Szatmary 25-08-2009
      IMPLICIT NONE
!     This subroutine contains microtasking directives.
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      INTEGER NGXP1,NGYP1,NGZP1
      INTEGER NGXP2,NGYP2,NGZP2
      INTEGER NGXB4,NGYB4,NGZB4
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(NGXP1=NGX+1,NGYP1=NGY+1,NGZP1=NGZ+1)
      PARAMETER(NGXP2=NGX+2,NGYP2=NGY+2,NGZP2=NGZ+2)
      PARAMETER(NGXB4=NGX/4,NGYB4=NGY/4,NGZB4=NGZ/4)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      INTEGER I,J,K
      double COMPLEX :: UR(0:NBX,0:NBY,0:NGZM1)
      double complex :: VR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: WR(0:NBX,0:NBY,0:NGZM1)

!$OMP   PARALLEL DO SHARED(UR,VR,WR) PRIVATE(J,K)
!     Wrap velocities for periodicity in the x direction
      do k=0,ngzm1
         do j=1,ngy
            ur(0, j,k) = ur(ngx,j,k)
            ur(ngx+1, j, k) = ur(1, j, k)
            ur(ngx+2, j, k) = ur(2, j, k)
            vr(0, j, k) = vr(ngx, j, k)
            vr(ngx+1, j, k) = vr(1, j, k)
            vr(ngx+2, j, k) = vr(2, j, k)
            wr(0, j, k) = wr(ngx, j, k)
            wr(ngx+1, j, k) = wr(1, j, k)
            wr(ngx+2, j, k) = wr(2, j, k)
         end do
      end do

!     Wrap velocities for periodicity in the y direction
      do k = 0, ngz - 1
         do i = 1, ngx
            ur(i, 0, k) = ur(i, ngy, k)
            ur(i, ngy+1, k) = ur(i, 1, k)
            ur(i, ngy+2, k) = ur(i, 2, k)
            vr(i, 0, k) = vr(i, ngy, k)
            vr(i, ngy+1, k) = vr(i, 1, k)
            vr(i, ngy+2, k) = vr(i, 2, k)
            wr(i, 0, k) = wr(i, ngy, k)
            wr(i, ngy+1, k) = wr(i, 1, k)
            wr(i, ngy+2, k) = wr(i, 2, k)
         end do
      end do

      RETURN
      END SUBROUTINE wrap
!*********************************************************************
      SUBROUTINE INPLANE(XPI,xfn)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      double precision :: FNGX,FNGY,FNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(FNGX=NGX,FNGY=NGY,FNGZ=NGZ)
      integer i,j
      double precision :: xfn(:,:)
      integer, parameter :: sqrtnpl=$sqrtnpl$

      double precision :: XPI(:,:)
      do i=0,sqrtnpl-1
         do j=0,sqrtnpl-1
            xpi(1,i*sqrtnpl+j+1)=(fngx-2*fngx/dble(sqrtnpl))*i/
     &           (dble(sqrtnpl)-1)+fngx/dble(sqrtnpl)+1.d0
            xpi(2,i*sqrtnpl+j+1)=$planey$
            xpi(3,i*sqrtnpl+j+1)=(fngz-2*fngz/dble(sqrtnpl))*j/
     &           (dble(sqrtnpl)-1)+fngz/dble(sqrtnpl)
            xfn(1,i*sqrtnpl+j+1)=xpi(1,i*sqrtnpl+j+1)
            xfn(2,i*sqrtnpl+j+1)=xpi(2,i*sqrtnpl+j+1)
            xfn(3,i*sqrtnpl+j+1)=xpi(3,i*sqrtnpl+j+1)
         end do
      end do
      return
      end subroutine inplane
!*********************************************************************
      SUBROUTINE PMHIST(XPI,XFN,FRC,FIRSTN,NEXTN,NUMBER,CONST,
     &     fp_start, fp_end, nfsize)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      INTEGER NGXP1,NGYP1,NGZP1
      INTEGER NGXP2,NGYP2,NGZP2
      INTEGER NGXB4,NGYB4,NGZB4
      INTEGER NFSIZE
      integer fp_start, fp_end
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG+1)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(NGXP1=NGX+1,NGYP1=NGY+1,NGZP1=NGZ+1)
      PARAMETER(NGXP2=NGX+2,NGYP2=NGY+2,NGZP2=NGZ+2)
      PARAMETER(NGXB4=NGX/4,NGYB4=NGY/4,NGZB4=NGZ/4)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)

      integer, parameter :: npl=$npl$
      double precision :: XPI(:,:)
      double precision :: XFN(:,:),FRC(:,:)
      integer firstn(:,:),number(:,:)
      integer nextn(:)
      integer i,ii,jj,n
      double precision :: CONST

      DO 1 I=fp_start,fp_end
      frc(1,I)=-CONST*(xfn(1,I)-XPI(1,I-fp_start+1))
      frc(2,I)=-CONST*(xfn(2,I)-XPI(2,I-fp_start+1))
      frc(3,I)=-CONST*(xfn(3,I)-XPI(3,I-fp_start+1))
 1    CONTINUE

!     Everything that follows in this sub might not need doing.
!     This part will be reassessed when performance profiling is taken.
!     Zero out all nextn?
      DO 16 N=1,nfsize
      NEXTN(N)=0
  16  CONTINUE

!     Zero out all firstn, nextn?
      DO 18 JJ=1,NGY
      DO 18 II=1,NGX
      FIRSTN(II,JJ)=0
      NUMBER(II,JJ)=0
  18  CONTINUE

!     Completely re-do firstn, nextn, number?
C SORT THE XFN DATA BY X-COORDINATE AND Y-COORDINATE USING LINKED LISTS
      DO 20 N=1,nfsize
         JJ=MOD(INT(xfn(2,N)+FLNGY)-1+NGY,NGY) + 1
         II=MOD(INT(xfn(1,N)+FLNGX)-1+NGX,NGX) + 1
         NEXTN(N)=FIRSTN(II,JJ)
         FIRSTN(II,JJ)=N
         NUMBER(II,JJ)=NUMBER(II,JJ)+1
  20  CONTINUE

      RETURN
      END SUBROUTINE PMHIST
!********************************************************************
      SUBROUTINE pushup(UR,VR,WR,XFN,FRC,FIRSTN,NUMBER,NEXTN)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      INTEGER NGXP1,NGYP1,NGZP1
      INTEGER NGXP2,NGYP2,NGZP2
      INTEGER NGXB4,NGYB4,NGZB4
      INTEGER NPTMAX
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(NGXP1=NGX+1,NGYP1=NGY+1,NGZP1=NGZ+1)
      PARAMETER(NGXP2=NGX+2,NGYP2=NGY+2,NGZP2=NGZ+2)
      PARAMETER(NGXB4=NGX/4,NGYB4=NGY/4,NGZB4=NGZ/4)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      PARAMETER(NPTMAX=256)
      double COMPLEX UR(0:,0:,0:),VR(0:,0:,0:), WR(0:,0:,0:)
      double precision :: XFN(:,:),FRC(:,:)
      INTEGER FIRSTN(:,:),NUMBER(:,:),NEXTN(:)
      double precision :: FLIZP1,FLJZP1,RAD1,RAD2,RAD3,ARG1,ARG2,ARG3
      INTEGER I,J,K,L,M,N,IZERO,JZERO,MZERO,NPT,NPOINTS,NUMREM
      INTEGER II,JJ,IJ,I3D,J3D,K3D,NEXTNOLD,NPREV,IR
      INTEGER MSHIFT,I0,IZ,JZ,KZ,IN,JN
      double precision,ALLOCATABLE :: XFN1OLD(:),XFN2OLD(:),XFN3OLD(:)
      double precision,ALLOCATABLE  :: FORCE1(:), FORCE2(:), FORCE3(:)
      double precision,ALLOCATABLE  ::  D1(:,:),   D2(:,:),   D3(:,:)
      double precision,ALLOCATABLE ::  D12(:,:,:)
      double precision,ALLOCATABLE ::DELTA(:,:)
      double precision,ALLOCATABLE  :: ULIN(:),VLIN(:),WLIN(:)
      INTEGER,ALLOCATABLE ::LINDX(:),ID(:,:)
      double precision,ALLOCATABLE :: FLKZP1(:)
C
C     THIS ROUTINE APPLIES THE FIBER FORCES TO THE FLUID.
C
C     FIBER POINT L INFLUENCES FLUID POINT (I,J,K) IFF
C     ALL OF THE FOLLOWING CONDITIONS ARE SATISFIED: 
C
C     ABS(I-XFN(1,L)) < 2.
C     ABS(J-XFN(2,L)) < 2.
C     ABS(K-XFN(3,L)) < 2.
C
C     USING FORTRAN THESE VALUES OF (I,J,K) ARE FOUND AS FOLLOWS: 
C
C     IZ=   XFN(1,L) -1.
C     JZ=   XFN(2,L) -1.
C     KZ=   XFN(3,L) -1.
C
C     THEN
C
C     I = IZ,IZ+3
C     J = JZ,JZ+3
C     K = KZ,KZ+3
C
C     IF THESE CONDITIONS ARE SATISFIED, THE COEFFICIENT THAT LINKS
C     FIBER POINT L TO THE FLUID POINT (I,J,K) IS: 
C
C     D = DEL(I-XFN(1,L)) * DEL(J-XFN(2,L)) * DEL(K-XFN(3,L))
C
C     WHERE
C     DEL(R) = (1. + COS((PI/2.)*R)) / 4.
C
C     ALGORITHM IS AS FOLLOWS: 
C
C     LET L=1,2, ... ,NP
C     FOR EACH L, LET (I,J,K) COVER RANGE OF VALUES DEFINED ABOVE.
C
C     FOR EACH VALUE OF (L,I,J,K)
C
C     U(I,J,K) = U(I,J,K) + D * FRC(1,L)
C     V(I,J,K) = V(I,J,K) + D * FRC(2,L)
C     W(I,J,K) = W(I,J,K) + D * FRC(3,L)
C
C     WHERE D WAS DEFINED ABOVE.
C

!      write(*,*) 'force l143'

      DO 89 IZERO=1,4
      DO 89 JZERO=1,4
!$OMP  PARALLEL DO &
!$OMP& SHARED(FIRSTN, NEXTN, NUMBER, IZERO, JZERO, XFN) &
!$OMP& PRIVATE(II, IJ, IN, JJ, JN, N, NEXTNOLD, NPREV)
      DO 85 IJ=0,NGXB4*NGYB4-1
      JJ   = JZERO + (IJ/NGYB4)*4
      II   = IZERO + (MOD(IJ,NGXB4))*4
      IF (NUMBER(II,JJ) .EQ. 0) GO TO 85
      NPREV = 0
      N     = FIRSTN(II,JJ)
   82 JN   = MOD((INT(XFN(2,N)+FLNGY)-1+NGY),NGY) + 1
      IN   = MOD((INT(XFN(1,N)+FLNGX)-1+NGX),NGX) + 1
      IF ((IN.NE.II) .OR. (JN.NE.JJ)) THEN
c       point N is in the wrong linked list
c       remember the pointer to the next one in the present linked list
        NEXTNOLD      = NEXTN(N)
c       add point N to the correct linked list
        NEXTN(N)      = FIRSTN(IN,JN)
        FIRSTN(IN,JN) = N
        NUMBER(IN,JN) = NUMBER(IN,JN) + 1
c       remove point N from the present linked list
        IF (NPREV .EQ. 0) THEN
          FIRSTN(II,JJ) = NEXTNOLD
        ELSE
          NEXTN(NPREV)  = NEXTNOLD
        END IF
        NUMBER(II,JJ) = NUMBER(II,JJ) - 1
c       advance to the remembered next point in the present linked list
c       the previous-point pointer is unchanged
        N             = NEXTNOLD
      ELSE
c       point N is in the correct linked list
c       advance the previous-point pointer and then the current-point pointer
        NPREV = N
        N     = NEXTN(N)
      END IF
      IF (N .NE. 0) GO TO 82
   85 CONTINUE
!$OMP END PARALLEL DO
   89 CONTINUE
C FIRSTN(I,J) CONTAINS THE INDEX OF THE FIRST POINT WHICH SIMULTANEOUSLY IS
C BETWEEN PLANES I AND I+1 AND IS BETWEEN PLANES J AND J+1.
C THIS POINT HAS COORDINATES:
C (XFN(1,FIRSTN(I,J)),XFN(2,FIRSTN(I,J)),XFN(3,FIRSTN(I,J))).
C NEXTN(FIRSTN(I,J)) CONTAINS THE INDEX OF THE SECOND SUCH POINT
C NEXTN(NEXTN(FIRSTN(I,J))) CONTAINS THE INDEX OF THE THIRD SUCH POINT
C ETC.
C IF FIRSTN(I,J) CONTAINS THE VALUE 0, THERE ARE NO SUCH POINTS.

!      write(*,*) 'force l179'

      MSHIFT = 16*NGZ
      DO 25 JZERO=1,4
      DO 25 IZERO=1,4
!$OMP  PARALLEL DO &
!$OMP& SHARED(FIRSTN,FRC,MSHIFT,NEXTN,NUMBER,UR,VR,WR,XFN,IZERO,JZERO)&
!$OMP& PRIVATE(ARG1,ARG2,ARG3,D1,D12,D2,D3,DELTA,&
!$OMP&        FLIZP1,FLJZP1,FLKZP1,FORCE1,FORCE2,FORCE3,&
!$OMP&        IJ,I,I0,I3D,II,IR,ID,IZ,J,J3D,JJ,JZ,K,K3D,KZ,L,&
!$OMP&        LINDX,M,MZERO,NPOINTS,NPT,NUMREM,RAD1,RAD2,RAD3,&
!$OMP&        ULIN,VLIN,WLIN,XFN1OLD,XFN2OLD,XFN3OLD)
      DO 20 IJ=0,NGXB4*NGYB4-1
      JJ = JZERO + (IJ/NGYB4)*4
      II = IZERO + (MOD(IJ,NGXB4))*4
      L = FIRSTN(II,JJ)
      IF (L.EQ.0) GO TO 20
      ALLOCATE  (ULIN(-15:16*NGZP2))
      ALLOCATE  (VLIN(-15:16*NGZP2))
      ALLOCATE  (WLIN(-15:16*NGZP2))
      ALLOCATE  (FORCE1(NPTMAX))
      ALLOCATE  (FORCE2(NPTMAX))
      ALLOCATE  (FORCE3(NPTMAX))
      ALLOCATE (XFN1OLD(NPTMAX),XFN2OLD(NPTMAX),XFN3OLD(NPTMAX))
      ALLOCATE   ( D1(NPTMAX,0:3),   D2(NPTMAX,0:3),   D3(NPTMAX,0:3))
      ALLOCATE   (D12(NPTMAX,0:3,0:3))
      ALLOCATE    (DELTA(0:64,NPTMAX))
      ALLOCATE (FLKZP1(NPTMAX))
      ALLOCATE  (lindx(NPTMAX))
      ALLOCATE  (ID(1:3,-15:16*NGZP2))
C CONSTRUCT THE INDEX ARRAY FOR SUBSEQUENT GATHERING AND SCATTERING OF THE
C VELOCITY DATA FOR THE APPROPRIATE 4 X 4 X NG COLUMN OF THE 3D GRID INTO
C AND OUT OF THE APPROPRIATE LINEAR ARRAYS, ONE FOR EACH DIRECTION.
C ZERO-FILL THE LINEAR ARRAYS. NOTICE THAT THE LINEAR DATA ARRAYS ARE LARGER
C THAN THE INDEX ARRAY; THIS FACILITATES SATISFYING THE PERIODIC B.C.
      DO 4 K=0,NGZM1
      DO 4 J=-1,2
      J3D   = MOD((JJ+J-1+NGY),NGY)+1
      I0    = K*16+(J+1)*4 + 2
      DO 4 I=-1,2
      I3D   = MOD((II+I-1+NGX),NGX)+1
      ID(1,I0+I)=I3D
      ID(2,I0+I)=J3D
      ID(3,I0+I)=K
    4 CONTINUE
!      write(*,*) 'force l219'
      DO 41 M=-15,0
      ULIN(M) = 0.0d0
      VLIN(M) = 0.0d0
      WLIN(M) = 0.0d0
   41 CONTINUE
!      write(*,*) 'force l225', izero, jzero
      DO 42 IR=1,16*NGZ
!         write(*,*) 'force l227', ir, id(1,ir), id(2,ir), id(3,ir)
      I3D=ID(1,IR)
      J3D=ID(2,IR)
      K3D=ID(3,IR)
!      write(*,*) 'force l231', ur(i3d, j3d, k3d), ur(i3d, j3d, k3d),
!     &     wr(i3d, j3d, k3d)
      ULIN(IR)=UR(I3D,J3D,K3D)
      VLIN(IR)=VR(I3D,J3D,K3D)
      WLIN(IR)=WR(I3D,J3D,K3D)
   42 CONTINUE
!      write(*,*) 'force l234'
      DO 43 M=16*NGZ+1,16*NGZP2
      ULIN(M) = 0.0d0
      VLIN(M) = 0.0d0
      WLIN(M) = 0.0d0
   43 CONTINUE
!      write(*,*) 'force l240'
      NUMREM = NUMBER(II,JJ)
    5 IF (NUMREM .GE. NPTMAX) THEN
        NPOINTS       = NPTMAX
        NUMREM        = NUMREM - NPTMAX
      ELSE
        NPOINTS       = NUMREM
        NUMREM        = 0
      END IF

!      write(*,*) 'force l250'
      DO 601 NPT=1,NPOINTS
      LINDX(NPT) = L
      L          = NEXTN(L)
  601 CONTINUE
!      write(*,*) 'force l255'
      DO 6 NPT=1,NPOINTS
      XFN1OLD(NPT)=XFN(1,LINDX(NPT))
      XFN2OLD(NPT)=XFN(2,LINDX(NPT))
      XFN3OLD(NPT)=XFN(3,LINDX(NPT))
      FORCE1(NPT)=FRC(1,LINDX(NPT))
      FORCE2(NPT)=FRC(2,LINDX(NPT))
      FORCE3(NPT)=FRC(3,LINDX(NPT))
    6 CONTINUE
      IZ     = INT(XFN1OLD(  1) - 1.d0 + FLNGX) - NGX
      JZ     = INT(XFN2OLD(  1) - 1.d0 + FLNGY) - NGY
      FLIZP1 = IZ+1
      FLJZP1 = JZ+1
      DO   7  NPT=1,NPOINTS
      KZ     = INT(XFN3OLD(NPT) - 1.d0 + FLNGZ) - NGZ
      FLKZP1(NPT) = KZ+1
    7 CONTINUE
      DO   8  NPT=1,NPOINTS
      ARG3      = XFN3OLD(NPT) - FLKZP1(NPT)
      ARG2      = XFN2OLD(NPT) - FLJZP1
      ARG1      = XFN1OLD(NPT) - FLIZP1
      RAD3      = dSQRT(1.d0+4.d0*ARG3*(1.d0-ARG3))
      RAD2      = dSQRT(1.d0+4.d0*ARG2*(1.d0-ARG2))
      RAD1      = dSQRT(1.d0+4.d0*ARG1*(1.d0-ARG1))
      D3(NPT,3) = (1.d0+2.d0*ARG3-RAD3)/8.d0
      D2(NPT,3) = (1.d0+2.d0*ARG2-RAD2)/8.d0
      D1(NPT,3) = (1.d0+2.d0*ARG1-RAD1)/8.d0
      D3(NPT,2) = (1.d0+2.d0*ARG3+RAD3)/8.d0
      D2(NPT,2) = (1.d0+2.d0*ARG2+RAD2)/8.d0
      D1(NPT,2) = (1.d0+2.d0*ARG1+RAD1)/8.d0
      D3(NPT,1) = (3.d0-2.d0*ARG3+RAD3)/8.d0
      D2(NPT,1) = (3.d0-2.d0*ARG2+RAD2)/8.d0
      D1(NPT,1) = (3.d0-2.d0*ARG1+RAD1)/8.d0
      D3(NPT,0) = (3.d0-2.d0*ARG3-RAD3)/8.d0
      D2(NPT,0) = (3.d0-2.d0*ARG2-RAD2)/8.d0
      D1(NPT,0) = (3.d0-2.d0*ARG1-RAD1)/8.d0
    8 CONTINUE
!      write(*,*) 'force l292'
      DO   9   J=0,3
      DO   9   I=0,3
      DO   9 NPT=1,NPOINTS
      D12(NPT,I,J) = D1(NPT,I)*D2(NPT,J)
    9 CONTINUE
      DO  10  K=0,3
      DO  10  J=0,3
      DO  10  I=0,3
      M        = K*16 + J*4 + I + 1
      DO  10  NPT=1,NPOINTS
      DELTA(M,NPT) = D12(NPT,I,J) * D3(NPT,K)
   10 CONTINUE
      DO 712 NPT=1,NPOINTS
      MZERO = 16*(INT(XFN3OLD(NPT) - 1.d0 + FLNGZ) - NGZ)
      DO 711   M=1,64
      ULIN(M+MZERO) = ULIN(M+MZERO) + DELTA(M,NPT)*FORCE1(NPT)
      VLIN(M+MZERO) = VLIN(M+MZERO) + DELTA(M,NPT)*FORCE2(NPT)
      WLIN(M+MZERO) = WLIN(M+MZERO) + DELTA(M,NPT)*FORCE3(NPT)
  711 CONTINUE
  712 CONTINUE
!      write(*,*) 'force l313'
      IF (NUMREM .NE. 0) GO TO 5
      DO 14 M=1,32
      ULIN(M) = ULIN(M) + ULIN(M+MSHIFT)
      VLIN(M) = VLIN(M) + VLIN(M+MSHIFT)
      WLIN(M) = WLIN(M) + WLIN(M+MSHIFT)
   14 CONTINUE
      DO 15 M=16*(NGZ-1)+1,16*NGZ
      ULIN(M) = ULIN(M) + ULIN(M-MSHIFT)
      VLIN(M) = VLIN(M) + VLIN(M-MSHIFT)
      WLIN(M) = WLIN(M) + WLIN(M-MSHIFT)
   15 CONTINUE
!      write(*,*) 'force l325'
      DO 16 IR=1,16*NGZ
      I3D=ID(1,IR)
      J3D=ID(2,IR)
      K3D=ID(3,IR)
      UR(I3D,J3D,K3D)=ULIN(IR)
      VR(I3D,J3D,K3D)=VLIN(IR)
      WR(I3D,J3D,K3D)=WLIN(IR)
   16 CONTINUE
!      write(*,*) 'force l334'
      DEALLOCATE (XFN1OLD,XFN2OLD,XFN3OLD)
      DEALLOCATE (D1,D2,D3)
      DEALLOCATE (D12,DELTA)
      DEALLOCATE (FLKZP1,LINDX,ID)
      DEALLOCATE  (ULIN,VLIN,WLIN)
      DEALLOCATE  (FORCE1)
      DEALLOCATE  (FORCE2)
      DEALLOCATE  (FORCE3)
20    CONTINUE
!$OMP END PARALLEL DO
25    CONTINUE
!      write(*,*) 'force l336'
      RETURN
      END SUBROUTINE pushup
!********************************************************************
      SUBROUTINE MOVE(UR,VR,WR,XFN,FIRSTN,NUMBER,NEXTN)
      IMPLICIT NONE
!     This subroutine contains microtasking directives.
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      INTEGER NGXP1,NGYP1,NGZP1
      INTEGER NGXP2,NGYP2,NGZP2
      INTEGER NGXB4,NGYB4,NGZB4
      INTEGER NPTMAX
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(NGXP1=NGX+1,NGYP1=NGY+1,NGZP1=NGZ+1)
      PARAMETER(NGXP2=NGX+2,NGYP2=NGY+2,NGZP2=NGZ+2)
      PARAMETER(NGXB4=NGX/4,NGYB4=NGY/4,NGZB4=NGZ/4)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      PARAMETER(NPTMAX=256)
      double complex :: ur(0:,0:,0:), vr(0:,0:,0:), wr(0:,0:,0:)
      double precision :: XFN(:,:)
      integer firstn(:,:),number(:,:),nextn(:)
      double precision :: FLIZP1,FLJZP1,RAD1,RAD2,RAD3,ARG1,ARG2,ARG3
      INTEGER I,J,K,L,M,MZERO,NPT,NPOINTS,NUMREM
      INTEGER II,JJ,IJ,I3D,J3D,K3D
      INTEGER I0,IZ,JZ,KZ
      double precision, ALLOCATABLE :: XFN1OLD(:),XFN2OLD(:),
     &     XFN3OLD(:)
      double precision :: xfn3oldold(nptmax)
      double precision,ALLOCATABLE :: D1(:,:),D2(:,:),D3(:,:)
      double precision,ALLOCATABLE :: D12(:,:,:)
      double precision,ALLOCATABLE :: DELTA(:,:)
      double precision,ALLOCATABLE :: ULIN(:),VLIN(:),WLIN(:)
      double precision,ALLOCATABLE :: UINT(:,:),VINT(:,:),WINT(:,:)
      double precision,ALLOCATABLE ::FLKZP1(:)
      INTEGER,ALLOCATABLE :: lindx(:)
C     THIS ROUTINE MOVES THE FIBER POINTS AT THE LOCAL FLUID VELOCITY.
C
C     FIBER POINT L IS INFLUENCED BY FLUID POINT (I,J,K) IFF
C     ALL THE FOLLOWING CONDITIONS ARE SATISFIED: 
C
C     ABS(I - XFN(1,L)) < 2.
C     ABS(J - XFN(2,L)) < 2.
C     ABS(K - XFN(3,L)) < 2.
C
C     USING FORTRAN, THESE VALUES OF (I,J,K) ARE FOUND AS FOLLOWS: 
C
C     IZ = XFN(1,L) -1.
C     JZ = XFN(2,L) -1.
C     KZ = XFN(3,L) -1.
C
C     THEN
C     I = IZ,IZ+3
C     J = JZ,JZ+3
C     K = KZ,KZ+3
C
C     IF THESE CONDITIONS ARE SATISFIED, THE COEFFICIENT THAT LINKS
C     FIBER POINT L TO FLUID POINT (I,J,K) IS: 
C
C     D = DEL(I-XFN(1,L)) * DEL(J-XFN(2,L)) * DEL(K-XFN(3,L))
C
C     WHERE
C
C     DEL(R) = (1. + COS((PI/2.)*R))/4. 
C
C     ALGORITHM IS AS FOLLOWS: 
C
C     LET L=1,2, ... ,NP
C     FOR EACH L, LET (I,J,K) COVER RANGE OF VALUES DEFINED ABOVE.
C
C     FOR EACH VALUE OF (L,I,J,K)
C
C     XFN(1,L) = XFN(1,L) + D * U(I,J,K)
C     XFN(2,L) = XFN(2,L) + D * V(I,J,K)
C     XFN(3,L) = XFN(3,L) + D * W(I,J,K)
C
C     WHERE D WAS DEFINED ABOVE.

!$OMP  PARALLEL DO &
!$OMP& SHARED(FIRSTN,NEXTN,NUMBER,UR,VR,WR,XFN)&
!$OMP& PRIVATE(ARG1,ARG2,ARG3,D1,D12,D2,D3,DELTA,&
!$OMP&        FLIZP1,FLJZP1,FLKZP1, &
!$OMP&        IJ,I,I0,I3D,II,IZ,J,J3D,JJ,JZ,K,K3D,KZ,L,&
!$OMP&        LINDX,M,MZERO,NPOINTS,NUMREM,NPT,RAD1,RAD2,RAD3,&
!$OMP&        UINT,ULIN,VLIN,VINT,WLIN,WINT,XFN1OLD,XFN2OLD,XFN3OLD)
      DO 20 IJ=0,NGX*NGY-1
      JJ = 1 + IJ/NGY
      II = 1 + MOD(IJ,NGX)
      L = FIRSTN(II,JJ)
      IF (L.EQ.0) GO TO 20
      ALLOCATE (XFN1OLD(NPTMAX),XFN2OLD(NPTMAX),XFN3OLD(NPTMAX))
      ALLOCATE   ( D1(NPTMAX,0:3),   D2(NPTMAX,0:3),   D3(NPTMAX,0:3))
      ALLOCATE   (D12(NPTMAX,0:3,0:3))
      ALLOCATE    (DELTA(0:64,NPTMAX))
      ALLOCATE (FLKZP1(NPTMAX))
      ALLOCATE  (lindx(NPTMAX))
      ALLOCATE (ULIN(-15:16*NGZP2),VLIN(-15:16*NGZP2),
     &     WLIN(-15:16*NGZP2))
      ALLOCATE (UINT(0:64,NPTMAX),VINT(0:64,NPTMAX),WINT(0:64,NPTMAX))
      DO 2 K=-1,NGZ+1
      K3D   = MOD(K+NGZ,NGZ)
      DO 2 J=-1,2
      J3D   = MOD(JJ+J-1+NGY,NGY)+1
      I0    = K*16+(J+1)*4 + 2
      DO 2 I=-1,2
      I3D   = MOD(II+I-1+NGX,NGX)+1
      ULIN(I0+I) = UR(I3D,J3D,K3D)
      VLIN(I0+I) = VR(I3D,J3D,K3D)
      WLIN(I0+I) = WR(I3D,J3D,K3D)
    2 CONTINUE
      NUMREM = NUMBER(II,JJ)
    5 IF(NUMREM .GE. NPTMAX) THEN
        NPOINTS       = NPTMAX
        NUMREM        = NUMREM - NPTMAX
      ELSE
        NPOINTS       = NUMREM
        NUMREM        = 0
      END IF
      DO 601 NPT=1,NPOINTS
      LINDX(NPT) = L
      L          = NEXTN(L)
  601 CONTINUE
      DO 6 NPT=1,NPOINTS
      XFN1OLD(NPT)=XFN(1,LINDX(NPT))
      XFN2OLD(NPT)=XFN(2,LINDX(NPT))
      XFN3OLD(NPT)=XFN(3,LINDX(NPT))
      xfn3oldold(npt) = xfn(3, lindx(npt))
    6 CONTINUE
      IZ     = INT(XFN1OLD(  1) - 1.d0 + FLNGX) - NGX
      JZ     = INT(XFN2OLD(  1) - 1.d0 + FLNGY) - NGY
      FLIZP1 = IZ+1
      FLJZP1 = JZ+1
      DO   7  NPT=1,NPOINTS
      KZ     = INT(XFN3OLD(NPT) - 1.d0 + FLNGZ) - NGZ
      FLKZP1(NPT) = KZ+1
    7 CONTINUE
      DO   8  NPT=1,NPOINTS
      ARG3      = XFN3OLD(NPT) - FLKZP1(NPT)
      ARG2      = XFN2OLD(NPT) - FLJZP1
      ARG1      = XFN1OLD(NPT) - FLIZP1
      RAD3      = dSQRT(1.d0+4.d0*ARG3*(1.d0-ARG3))
      RAD2      = dSQRT(1.d0+4.d0*ARG2*(1.d0-ARG2))
      RAD1      = dSQRT(1.d0+4.d0*ARG1*(1.d0-ARG1))
      D3(NPT,3) = (1.d0+2.d0*ARG3-RAD3)/8.d0
      D2(NPT,3) = (1.d0+2.d0*ARG2-RAD2)/8.d0
      D1(NPT,3) = (1.d0+2.d0*ARG1-RAD1)/8.d0
      D3(NPT,2) = (1.d0+2.d0*ARG3+RAD3)/8.d0
      D2(NPT,2) = (1.d0+2.d0*ARG2+RAD2)/8.d0
      D1(NPT,2) = (1.d0+2.d0*ARG1+RAD1)/8.d0
      D3(NPT,1) = (3.d0-2.d0*ARG3+RAD3)/8.d0
      D2(NPT,1) = (3.d0-2.d0*ARG2+RAD2)/8.d0
      D1(NPT,1) = (3.d0-2.d0*ARG1+RAD1)/8.d0
      D3(NPT,0) = (3.d0-2.d0*ARG3-RAD3)/8.d0
      D2(NPT,0) = (3.d0-2.d0*ARG2-RAD2)/8.d0
      D1(NPT,0) = (3.d0-2.d0*ARG1-RAD1)/8.d0
    8 CONTINUE
      DO   9   J=0,3
      DO   9   I=0,3
      DO   9 NPT=1,NPOINTS
      D12(NPT,I,J) = D1(NPT,I)*D2(NPT,J)
    9 CONTINUE
      DO  10  K=0,3
      DO  10  J=0,3
      DO  10  I=0,3
      M = K*16 + J*4 + I + 1
      DO  10 NPT=1,NPOINTS
      DELTA(M,NPT) = D12(NPT,I,J) * D3(NPT,K)
   10 CONTINUE
      DO 122 NPT=1,NPOINTS
      MZERO = 16*(INT(XFN3OLD(NPT) - 1.d0 + FLNGZ) - NGZ)
      DO 121 M=1,64
      UINT(M,NPT) = ULIN(MZERO+M)*DELTA(M,NPT)
      VINT(M,NPT) = VLIN(MZERO+M)*DELTA(M,NPT)
      WINT(M,NPT) = WLIN(MZERO+M)*DELTA(M,NPT)
  121 CONTINUE
  122 CONTINUE
      DO 124 M=1,64
      DO 123 NPT=1,NPOINTS
      XFN1OLD(NPT) = XFN1OLD(NPT) + UINT(M,NPT)
      XFN2OLD(NPT) = XFN2OLD(NPT) + VINT(M,NPT)
      XFN3OLD(NPT) = XFN3OLD(NPT) + WINT(M,NPT)
!      if (wint(m, npt) /= 0.0d0) write(*,*) 'force l461', m, npt, 
!     &     wint(m, npt), xfn3oldold(npt), xfn3old(npt)
  123 CONTINUE
  124 CONTINUE
      DO 17 NPT=1,NPOINTS
      XFN(1,LINDX(NPT))=XFN1OLD(NPT)
      XFN(2,LINDX(NPT))=XFN2OLD(NPT)
      XFN(3,LINDX(NPT))=XFN3OLD(NPT)
!      if (dabs(xfn3old(npt) - xfn3oldold(npt)) > 1.d-13) then
!         write(*,*) 'force l468',npt, xfn3old(npt), xfn3oldold(npt), 
!     &        xfn(3,lindx(npt)), xfn3fullold(lindx(npt))
!      end if
17    CONTINUE
      IF (NUMREM .NE. 0) GO TO 5
      DEALLOCATE (XFN1OLD,XFN2OLD,XFN3OLD)
      DEALLOCATE  (UINT,VINT,WINT)
      DEALLOCATE  (ULIN,VLIN,WLIN)
      DEALLOCATE (D1,D2,D3)
      DEALLOCATE (D12,DELTA)
      DEALLOCATE (FLKZP1,lindx)
20    CONTINUE
      RETURN
      END SUBROUTINE MOVE
!*********************************************************************

