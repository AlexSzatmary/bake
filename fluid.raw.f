      SUBROUTINE FLUIDUP(klok,UR,VR,WR, pr, QRFACT,DSQ, DX, DY, DZ)
      IMPLICIT NONE
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
      COMMON/FFT/FSCALE
      COMMON/FFTX/TABX(2*NGX,-1:1)
      COMMON/FFTX/TABY(2*NGY,-1:1)
      COMMON/FFTX/TABZ(2*NGZ,-1:1)
      double precision :: FSCALE
      INTEGER I,J,K,KLOK
      double COMPLEX :: TABX,TABY,TABZ
      double COMPLEX :: UR(0:NBX,0:NBY,0:NGZM1),VR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: WR(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: PR(0:NBX,0:NBY,0:NGZM1)
      double precision ::  QRFACT(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: DX(0:NBX,0:NBY,0:NGZM1),DY(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: DZ(0:NBX,0:NBY,0:NGZM1)
      double precision :: DSQ(0:NBX,0:NBY,0:NGZM1)

!     TAKE FORWARD TRANFORM OF UR, VR, WR, QR
      CALL CFFT3D(-1,FSCALE,UR,TABX,TABY,TABZ)
      CALL CFFT3D(-1,FSCALE,VR,TABX,TABY,TABZ)
      CALL CFFT3D(-1,FSCALE,WR,TABX,TABY,TABZ)
!     COMPUTE TRANSFORMED PR, UR, VR, WR

      DO 100 K=0,NGZM1
      DO 100 J=1,NGY
      DO 100 I=1,NGX
      UR(I,J,K) = UR(I,J,K)/QRFACT(I,J,K) 
      VR(I,J,K) = VR(I,J,K)/QRFACT(I,J,K)
      WR(I,J,K) = WR(I,J,K)/QRFACT(I,J,K)
      PR(I,J,K) = (UR(I,J,K)*DX(I,J,K)+VR(I,J,K)*DY(I,J,K)
     &  + WR(I,J,K)*DZ(I,J,K))*DSQ(I,J,K) 
      UR(I,J,K) = UR(I,J,K) + PR(I,J,K)*DX(I,J,K) 
      VR(I,J,K) = VR(I,J,K) + PR(I,J,K)*DY(I,J,K)
      WR(I,J,K) = WR(I,J,K) + PR(I,J,K)*DZ(I,J,K)
  100 CONTINUE
!$OMP END PARALLEL DO

!     TAKE BACKWARD TRANFORM OF UR, VR, WR, pressure
      CALL CFFT3D(+1,1.0d0,UR,TABX,TABY,TABZ)
      CALL CFFT3D(+1,1.0d0,VR,TABX,TABY,TABZ)
      CALL CFFT3D(+1,1.0d0,WR,TABX,TABY,TABZ)
      call cfft3d(+1,1.0d0,pr,tabx,taby,tabz)
      call dumpstatus(klok, 'fluid l60', 'status.txt')
      RETURN
      END SUBROUTINE FLUIDUP
!*********************************************************************
      SUBROUTINE INFLUIDU(VSC,QRFACT, DSQ,DX,DY,DZ)
      IMPLICIT NONE
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
      COMMON/FFT/FSCALE
      COMMON/FFTX/TABX(2*NGX,-1:1)
      COMMON/FFTX/TABY(2*NGY,-1:1)
      COMMON/FFTX/TABZ(2*NGZ,-1:1)
      double COMPLEX TABX,TABY,TABZ
      double precision :: FSCALE,VSC,FOURNU,PI,TWOPI
      INTEGER I,J,K,KS,MIDPTX,MIDPTY,MIDPTZ
      double precision :: SINSQX(-1:NGXP1),SINSQY(-1:NGYP1),
     &     SINSQZ(-1:NGZP1)
      double precision :: QRFACT(-1:NGXP1,-1:NGYP1,0:NGZM1) 
      double precision :: DXG(-1:NGXP1),DYG(-1:NGYP1),DZG(-1:NGZP1)
      double precision :: DXO(-1:NGXP1),DYO(-1:NGYP1),DZO(-1:NGZP1)
      double COMPLEX :: DX(-1:NGXP1,-1:NGYP1,0:NGZM1)
      double complex :: DY(-1:NGXP1,-1:NGYP1,0:NGZM1)
      double COMPLEX :: DZ(-1:NGXP1,-1:NGYP1,0:NGZM1)
      double precision :: DSQ(-1:NGXP1,-1:NGYP1,0:NGZM1)

      pi = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits
      TWOPI  = 2.d0*PI
      FOURNU = 4.d0*VSC
      MIDPTX = NGX/2
      MIDPTY = NGY/2
      MIDPTZ = NGZ/2
      DXG=dCMPLX(0.d0,0.d0)
      DYG=dCMPLX(0.d0,0.d0)
      DZG=dCMPLX(0.d0,0.d0)
      DXO=dCMPLX(0.d0,0.d0)
      DYO=dCMPLX(0.d0,0.d0)
      DZO=dCMPLX(0.d0,0.d0)
      DSQ=0.d0
      DO 11 KS=0,NGXM1
!      VXFACT(KS) =  dCMPLX(0.d0,-dSIN(TWOPI*KS/FLNGX))
      DXG(KS) =  dSIN(TWOPI*KS/FLNGX)*(dsqrt(2.d0)/2.d0+
     & (1.d0-dsqrt(2.d0)/2.d0)*dCOS(TWOPI*KS/FLNGX))
      DXO(KS) =  ((dCOS(PI*KS/FLNGX))**2)*(1.d0-2.d0* 
     & (1.d0-2.d0*dsqrt(2.d0)/PI)*((dSIN(PI*KS/FLNGX))**2))
   11 CONTINUE
      DO 12 KS=0,NGYM1
!      VYFACT(KS) =  dCMPLX(0.d0,-dSIN(TWOPI*KS/FLNGY))
      DYG(KS) =  dSIN(TWOPI*KS/FLNGY)*(dsqrt(2.d0)/2.d0+ 
     & (1.d0-dsqrt(2.d0)/2.d0)*dCOS(TWOPI*KS/FLNGY))
      DYO(KS) =  ((dCOS(PI*KS/FLNGY))**2)*(1.d0-2.d0* 
     & (1.d0-2.d0*dsqrt(2.d0)/PI)*((dSIN(PI*KS/FLNGY))**2))
   12 CONTINUE
      DO 13 KS=0,NGZM1
!      VZFACT(KS) =  dCMPLX(0.d0,-dSIN(TWOPI*KS/FLNGZ))
      DZG(KS) =  dSIN(TWOPI*KS/FLNGZ)*(dsqrt(2.d0)/2.d0+ 
     & (1.d0-dsqrt(2.d0)/2.d0)*dCOS(TWOPI*KS/FLNGZ))
      DZO(KS) =  ((dCOS(PI*KS/FLNGZ))**2)*(1.d0-2.d0* 
     & (1.d0-2.d0*dsqrt(2.d0)/PI)*((dSIN(PI*KS/FLNGZ))**2))
   13 CONTINUE
C$$$      VXFACT(-1) = (0.0d0,0.0d0)
C$$$      VXFACT(MIDPTX) = (0.0d0,0.0d0)
C$$$      VYFACT(-1) = (0.0d0,0.0d0)
C$$$      VYFACT(MIDPTY) = (0.0d0,0.0d0)
C$$$      VZFACT(-1) = (0.0d0,0.0d0)
C$$$      VZFACT(MIDPTZ) = (0.0d0,0.0d0)
      DXG(-1)=0.0d0
      DYG(-1)=0.0d0
      DZG(-1)=0.0d0
      DXG(MIDPTX)=0.0d0
      DYG(MIDPTY)=0.0d0
      DZG(MIDPTZ)=0.0d0
      DXO(MIDPTX)=0.0d0
      DYO(MIDPTY)=0.0d0
      DZO(MIDPTZ)=0.0d0
      DO 21 KS=0,NGXM1
      SINSQX(KS) = (dSIN(PI*KS/FLNGX))**2
   21 CONTINUE
      DO 22 KS=0,NGYM1
      SINSQY(KS) = (dSIN(PI*KS/FLNGY))**2
   22 CONTINUE
      DO 23 KS=0,NGZM1
      SINSQZ(KS) = (dSIN(PI*KS/FLNGZ))**2
   23 CONTINUE
      DO 25 K=0,NGZM1
      DO 25 J=0,NGYM1
      DO 25 I=0,NGXM1
      QRFACT(I,J,K) = 1.d0+FOURNU*(SINSQX(I)+SINSQY(J)+SINSQZ(K))
      DX(I,J,K)=dCMPLX(0.d0,-DXG(I)*DYO(J)*DZO(K))
      DY(I,J,K)=dCMPLX(0.d0,-DXO(I)*DYG(J)*DZO(K))
      DZ(I,J,K)=dCMPLX(0.d0,-DXO(I)*DYO(J)*DZG(K))
   25 CONTINUE
      DO 31 KS=0,NGXM1
      SINSQX(KS) = (dSIN(TWOPI*KS/FLNGX))**2
   31 CONTINUE
      DO 32 KS=0,NGYM1
      SINSQY(KS) = (dSIN(TWOPI*KS/FLNGY))**2
   32 CONTINUE
      DO 33 KS=0,NGZM1
      SINSQZ(KS) = (dSIN(TWOPI*KS/FLNGZ))**2
   33 CONTINUE
      SINSQX(MIDPTX) = 0.0d0
      SINSQY(MIDPTY) = 0.0d0
      SINSQZ(MIDPTZ) = 0.0d0
      DO 35 K=0,NGZM1
      DO 35 J=0,NGYM1
      DO 35 I=0,NGXM1
!      PRDENO(I,J,K) = (SINSQX(I)+SINSQY(J)+SINSQZ(K))
      IF(((I.EQ.0).AND.(J.EQ.0).AND.(K.EQ.0)).OR.
     & (I.EQ.MIDPTX).OR.(J.EQ.MIDPTY).OR.(K.EQ.MIDPTZ)) THEN
      DSQ(I,J,K)=0.d0
      ELSE
!     DSQ(I,J,K)=1.d0/((DX(I,J,K))**2+(DY(I,J,K))**2+(DZ(I,J,K))**2)
      DSQ(I,J,K)=1.d0/((DXG(I)*DYO(J)*DZO(K))**2+
     & (DXO(I)*DYG(J)*DZO(K))**2+(DXO(I)*DYO(J)*DZG(K))**2)
      ENDIF
   35 CONTINUE
C$$$      PRDENO(     0,     0,     0) = 1.0d0
C$$$      PRDENO(MIDPTX,     0,     0) = 1.0d0
C$$$      PRDENO(     0,MIDPTY,     0) = 1.0d0
C$$$      PRDENO(MIDPTX,MIDPTY,     0) = 1.0d0
C$$$      PRDENO(     0,     0,MIDPTZ) = 1.0d0
C$$$      PRDENO(MIDPTX,     0,MIDPTZ) = 1.0d0
C$$$      PRDENO(     0,MIDPTY,MIDPTZ) = 1.0d0
C$$$      PRDENO(MIDPTX,MIDPTY,MIDPTZ) = 1.0d0
      FSCALE = (1.d0/(FLNGX*FLNGY*FLNGZ))
!     DO 38 K=0,NGZ
!     DO 38 J=0,NGY
!     DO 38 I=0,NGX
!     write(231,36)i,j,k,DX(i,j,k)
!     write(232,36)i,j,k,DY(i,j,k)
!     write(233,36)i,j,k,DZ(i,j,k)
!     write(234,37)i,j,k,dsq(i,j,k)
!36     format(I3,I3,I3,X,e13.5,X,e13.5)
!37     format(I3,I3,I3,X,e13.5)
!38   CONTINUE
!     INITIALIZE THE ARRAYS TABLE AND WORK
      CALL CFFT3D(0,1.d0,1.d0,TABX,TABY,TABZ)
      RETURN
      END SUBROUTINE INFLUIDU
!*********************************************************************
      SUBROUTINE CFFT3D(ISIGN,FSCALE,A,UX,UY,UZ)
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      INTEGER NGXP1,NGYP1,NGZP1
      INTEGER NGXP2,NGYP2,NGZP2
      INTEGER NGXB4,NGYB4,NGZB4
      INTEGER NX,NY,NZ
      INTEGER MX,MY,MZ
      INTEGER NBXP1,NBYP1,NBZP1
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(NGXP1=NGX+1,NGYP1=NGY+1,NGZP1=NGZ+1)
      PARAMETER(NGXP2=NGX+2,NGYP2=NGY+2,NGZP2=NGZ+2)
      PARAMETER(NGXB4=NGX/4,NGYB4=NGY/4,NGZB4=NGZ/4)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      PARAMETER(MX=LXNG,MY=LYNG,MZ=LZNG)
      PARAMETER(NX=2**MX,NY=2**MY,NZ=2**MZ)
      PARAMETER(NBXP1=NBX+1,NBYP1=NBY+1,NBZP1=NBZ+1)
      double COMPLEX A(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX B(1:NBXP1,1:NBYP1,1:NGZ)
      double COMPLEX UX(2*NX,-1:1),UY(2*NY,-1:1),UZ(2*NZ,-1:1)
      double COMPLEX TX(1:NX),TY(1:NY),TZ(1:NZ),W
      double precision :: PI,SIGN,ANG,FSCALE
      INTEGER I,J,K,IS,ISIGN,L,LE,LE1,IW,IV,IP
      INTEGER NMX1,NVX2,NMY1,NVY2,NMZ1,NVZ2

      IF (ISIGN .EQ. 0) THEN
        PI = 3.141592653589793238462643d0
        DO 3 IS=-1,1,2
        SIGN = IS
        DO 2 L=1,MX
        LE        = 2**L
        LE1       = LE/2
        UX(LE1,IS) = (1.d0,0.d0)
        ANG       = SIGN*PI/LE1
        W         = dCMPLX(dCOS(ANG),dSIN(ANG))
        IF (LE1 .EQ. 1) THEN
          W = dcmplx(-1.d0, 0.d0)
        ELSEIF (LE1 .EQ. 2) THEN
          W = dCMPLX(0.d0, SIGN)
        END IF
        DO 2 J=1,LE1-1
        UX(LE1+J,IS) = UX(LE1+J-1,IS)*W
    2   CONTINUE
    3   CONTINUE
        DO 5 IS=-1,1,2
        SIGN = IS
        DO 4 L=1,MY
        LE        = 2**L
        LE1       = LE/2
        UY(LE1,IS) = dcmplx(1.d0,0.d0)
        ANG       = SIGN*PI/LE1
        W         = dCMPLX(dCOS(ANG),dSIN(ANG))
        IF (LE1 .EQ. 1) THEN
          W = dcmplx(-1.d0, 0.d0)
        ELSEIF (LE1 .EQ. 2) THEN
          W = dCMPLX(0.d0, SIGN)
        END IF
        DO 4 J=1,LE1-1
        UY(LE1+J,IS) = UY(LE1+J-1,IS)*W
    4   CONTINUE
    5   CONTINUE
        DO 7 IS=-1,1,2
        SIGN = IS
        DO 6 L=1,MZ
        LE        = 2**L
        LE1       = LE/2
        UZ(LE1,IS) = dcmplx(1.d0,0.d0)
        ANG       = SIGN*PI/LE1
        W         = dCMPLX(dCOS(ANG),dSIN(ANG))
        IF (LE1 .EQ. 1) THEN
          W = dcmplx(-1.d0, 0.d0)
        ELSEIF (LE1 .EQ. 2) THEN
          W = dCMPLX(0.d0, SIGN)
        END IF
        DO 6 J=1,LE1-1
        UZ(LE1+J,IS) = UZ(LE1+J-1,IS)*W
    6   CONTINUE
    7   CONTINUE
        RETURN
      END IF
      SIGN = ISIGN/IABS(ISIGN)
      IS   =  SIGN
!$OMP PARALLEL DO SHARED(A,B,FSCALE) PRIVATE(IW,IV,I)
      DO 11 IW=1,NZ 
      DO 11 IV=1,NY 
      DO 11  I=1,NX 
      B(I,IV,IW) = FSCALE*A(I,IV,IW-1)
   11 CONTINUE
!$OMP END PARALLEL DO
      NVX2 = NX/2
      NMX1 = NX-1
!DIR$ IVDEP
!$OMP PARALLEL DO SHARED(B,IS,NMX1,NVX2,UX) PRIVATE(IW,IV,IP,I,J,K,L,LE,LE1,TY)
      DO 30 IW=1,NZ
      J = 1
      DO 17 I=1,NMX1
      IF (I.LT.J) THEN
        DO 91 IV=1,NY
        TY(IV)      = B(J,IV,IW)
   91   CONTINUE
        DO 92 IV=1,NY
        B(J,IV,IW) = B(I,IV,IW)
        B(I,IV,IW) = TY(IV)
   92   CONTINUE
      END IF
      K = NVX2
   16 IF (K.LT.J) THEN
        J = J-K
        K = K/2
        GO TO 16
      END IF
      J = J+K
   17 CONTINUE
      DO 20 L=1,MX
      LE  = 2**L
      LE1 = LE/2
      DO 20 J=1,LE1
      DO 10 I=J,NX,LE
      IP    = I+LE1
      DO 95 IV=1,NY
      TY(IV)= B(IP,IV,IW)*UX(LE1+J-1,IS)
   95 CONTINUE
      DO 96 IV=1,NY
      B(IP,IV,IW) = B(I,IV,IW)-TY(IV)
      B( I,IV,IW) = B(I,IV,IW)+TY(IV)
   96 CONTINUE
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
!$OMP END PARALLEL DO
      NVY2 = NY/2
      NMY1 = NY-1
!DIR$ IVDEP
!$OMP PARALLEL DO SHARED(B,IS,NMY1,NVY2,UY) PRIVATE(IW,IV,IP,I,J,K,L,LE,LE1,TZ)
      DO 130 IW=1,NX
      J = 1
      DO 117 I=1,NMY1
      IF (I.LT.J) THEN
        DO 191 IV=1,NZ
        TZ(IV) = B(IW,J,IV)
  191   CONTINUE
        DO 192 IV=1,NZ
        B(IW,J,IV) = B(IW,I,IV)
        B(IW,I,IV) = TZ(IV)
  192   CONTINUE
      END IF
      K = NVY2
  116 IF (K.LT.J) THEN
        J = J-K
        K = K/2
        GO TO 116
      END IF
      J = J+K
  117 CONTINUE
      DO 120 L=1,MY
      LE  = 2**L
      LE1 = LE/2
      DO 120 J=1,LE1
      DO 110 I=J,NY,LE
      IP    = I+LE1
      DO 195 IV=1,NZ
      TZ(IV)= B(IW,IP,IV)*UY(LE1+J-1,IS)
  195 CONTINUE
      DO 196 IV=1,NZ
      B(IW,IP,IV) = B(IW,I ,IV)-TZ(IV)
      B(IW,I ,IV) = B(IW,I ,IV)+TZ(IV)
  196 CONTINUE
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
!$OMP END PARALLEL DO
      NVZ2 = NZ/2
      NMZ1 = NZ-1
!DIR$ IVDEP
!$OMP PARALLEL DO SHARED(A,B,IS,NMZ1,NVZ2,UZ) PRIVATE(IW,IV,IP,I,J,K,L,LE,LE1,TX)
      DO 230 IW=1,NY
      J = 1
      DO 217 I=1,NMZ1
      IF (I.LT.J) THEN
        DO 291 IV=1,NX
        TX(IV) = B(IV,IW,J)
  291   CONTINUE
        DO 292 IV=1,NX
        B(IV,IW,J) = B(IV,IW,I)
        B(IV,IW,I) = TX(IV)
  292   CONTINUE
      END IF
      K = NVZ2
  216 IF (K.LT.J) THEN
        J = J-K
        K = K/2
        GO TO 216
      END IF
      J = J+K
  217 CONTINUE
      DO 220 L=1,MZ
      LE  = 2**L
      LE1 = LE/2
      DO 220 J=1,LE1
      DO 210 I=J,NZ,LE
      IP    = I+LE1
      DO 295 IV=1,NX
      TX(IV)= B(IV,IW,IP)*UZ(LE1+J-1,IS)
  295 CONTINUE
      DO 296 IV=1,NX
      B(IV,IW,IP) = B(IV,IW,I )-TX(IV)
      B(IV,IW,I ) = B(IV,IW,I )+TX(IV)
  296 CONTINUE
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SHARED(A,B) PRIVATE(IW,IV,I)
      DO 401 IW=1,NZ
      DO 401 IV=1,NY
      DO 401  I=1,NX
      A(I,IV,IW-1) = B(I,IV,IW)
  401 CONTINUE
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE CFFT3D
!*************************************************************************

