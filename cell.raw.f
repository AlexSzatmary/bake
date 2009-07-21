C *********************************************************************
C *     Copyright 1994 Charles S. Peskin and David M. McQueen         *
C *     All rights reserved. No portion of this program may be        *
C *     transmitted to any persons without the prior permission       *
C *     of the Copyright holders.                                     *
C *     Contact: Charles S. Peskin   email:  peskin@cims.nyu.edu      *
C *              David M. McQueen    email: mcqueen@cims.nyu.edu      *
C *********************************************************************
C
C**********************************************************************
C  (FIBER ACTIVATION: AN INTEGER ARRAY, LAFLAG, FLAGS ACTIVE LINKS, [A AND V])
c   fixed bug in subroutine fluidup:
c   changed:
c     WR(I,J,K) = (WR(I,J,K) + PR(I,J,K)*PRFACT(K))/QRFACT(I,J,K)
c   to:
c     WR(I,J,K) = (WR(I,J,K) + PR(I,J,K)*PRFACT(K+1))/QRFACT(I,J,K)
c  (Combined upwind difference [hrtxp_61u6.f] and microtasking [hrtxp_61m3e.f])
C
C  OUTPUT FILES:
C      FORT.7  VELOCITY AND PRESSURE ON A COARSE GRID
C      FORT.8  VELOCITY AND PRESSURE ON A FINE   GRID
C      FORT.9  VELOCITY AND LOCATION ON THE VALVE FLOWMETER WEBS
C      FORT.10 FIBER  DATA
C      FORT.11 MARKER DATA
C
C  INPUT PARAMETERS: 
C      LCUBE=LENGTH OF SIDE OF CUBE (CM)
C      NU=KINEMATIC VISCOSITY (CM**2/SEC)
C      RHO=FLUID DENSITY (GM/CM**3)
C                (MU=NU*RHO=DYNAMIC VISCOSITY)
C      NG=NUMBER OF GRID POINTS IN EACH DIRECTION 
C      TD=TIME STEP (SEC)
C
C  PARAMETERS: 
C     PLANES ARE NG X NG
C     BUT ARE STORED IN ARRAYS WITH DIMENSIONS
C     (0:NB,1:NG)   WHERE NB=NG+2
C
C     A FIBER IS COMPOSED OF POINTS
C     A GROUP IS COMPOSED OF FIBERS HAVING THE SAME NUMBER OF POINTS
C     A BUNCH IS COMPOSED OF GROUPS
C
C     NBUNCH    = NUMBER OF BUNCHES IN THE ENTIRE STRUCTURE 
C     NGROUPS(J)= NUMBER OF GROUPS IN BUNCH J, J=1,...,NBUNCH   (NGROUPS(0)=0)
C     NFG(I)    = NUMBER OF FIBERS IN GROUP I, I=NGROUPS(J-1)+1,...,NGROUPS(J)
C     NPF(I)    = NUMBER OF POINTS IN A FIBER IN GROUP I
C     IMAX      = SUM(J=1,NBUNCH):NGROUPS(J)
C                 IMAX IS THE TOTAL NUMBER OF GROUPS IN THE ENTIRE STRUCTURE.
C
C     NFSIZE    = MAX(J=1,NBUNCH):SUM(I=ISTART(J),ISTOP(J)):NFG(I)*NPF(I)
C     ISTART(J) = NGROUPS(J-1)+1
C     ISTOP (J) = NGROUPS(J)
C     NGROUPS(0)=0
C
C     ASSUMED VALUE  MAX(I=1,IMAX):NFG(I)        =    64
C     ASSUMED VALUE  MAX(I=1,IMAX):NPF(I)        =   530
C     ASSUMED VALUE  MAX(I=1,IMAX):NFG(I)*NPF(I) = 33920
C
C     WARNING: THE ASSUMED VALUES GIVEN ABOVE
C     AND THE CONSTANT NFSIZE=2020 GIVEN
C     BELOW WERE DETERMINED EMPIRICALLY FOR A
C     PARTICULAR HEART ANATOMY. ANY CHANGE TO
C     THIS ANATOMY MAY REQUIRE ALTERATION OF
C     OF THESE VALUES.
C     A SIMILAR WARNING APPLIES TO NMSIZE.
C
c     nsrcs = the number of sources in the heart, reasonably 5:
c             (1) superior vena cava
c             (2) inferior vena cava
c             (3) pulmonary veins (left  atrium)
c             (4) pulmonary artery (normally thought of as a sink)
c             (5) aorta            (normally thought of as a sink)
c               
C *********************************************************************
C
c     NOTES ON CHANGES MADE TO PROGRAM PULSE3D TO CONVERT IT TO
c     RBC3D: SIMULATION CODE FOR THE STUDY OF ERYTHROCYTES IN 3D FLOW
c     NOTES STARTED ON 1/19/95
C
c  1/19/95   removed FIBERX routines that calculate heart fiber forces
c  replaced with MEMBNX and related routines that calculate forces
c  exreted by an erythrocyte membrane
c   added common block membrn
C**********************************************************************
      PROGRAM cell
      IMPLICIT NONE
      INTEGER LXNG,LYNG,LZNG,NGX,NGY,NGZ,NFSIZE,NFSIZE2
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      parameter(nfsize=$nsnode$,nfsize2=$nselm$)
      integer, parameter :: npl=$npl$
      integer, parameter :: m_start=$m_start$, m_end=$m_end$
      INTEGER KLOK,KLOK1,KLOK0,KLOKEND,NSTEP
      double precision :: T,H,h64,TD,VSC,TIME,RHO,PI,RADX,FOSTAR
      double precision :: xcenter, ycenter, zcenter
      double precision :: xcenterold, ycenterold, zcenterold
      double precision ::  LCUBE,NU,MU,MASS,LENGTH
!     Velocities are always expressed in program units
      double COMPLEX,ALLOCATABLE :: UR(:,:,:),VR(:,:,:),WR(:,:,:)
      double complex :: meanstablev
!     Handy look-up tables for fluid mechanics calculations.
      double COMPLEX, ALLOCATABLE :: VXFACT(:),VYFACT(:),VZFACT(:)
      double precision,ALLOCATABLE :: PRDENO(:,:,:),QRFACT(:,:,:)
      double COMPLEX :: DX(0:NBX,0:NBY,0:NGZM1)
      double complex :: DY(0:NBX,0:NBY,0:NGZM1)
      double COMPLEX :: DZ(0:NBX,0:NBY,0:NGZM1)
      double precision :: DSQ(0:NBX,0:NBY,0:NGZM1)
      double precision, parameter :: Eh = $Eh$
      double precision, parameter :: capillary_no = $capillary_no$
      integer, parameter :: FVS = $FVS$
      double precision :: planey = $planey$
      character*80 message

!     Array of linked lists representing the solid nodes for a given 
!     hxh tile, spanning (x,y) to (x+h, y+h)
      INTEGER,ALLOCATABLE :: FIRSTN(:,:),NUMBER(:,:),NEXTN(:)
!     xfn is the x coordinates of the solid nodes, in program units
!     frc is the force exerted at each node on the fluid; this is in program
!     units.
!     shpint and shpfs are arrays of finite element shape factor parameters.
!     These are in real dimensions.
      double precision,ALLOCATABLE :: XFN(:,:),FRC(:,:),shpint(:,:),
     &     shpfs(:,:)
      double precision xpi(1:3,1:npl)
!     This array associates 3 nodes with a numbered element; three corners
!     on a triangle.
      INTEGER,ALLOCATABLE :: elmnew(:,:)

!     Changes made 10-19-06 to implement BC Homog
      integer i, j, k
      double precision :: gamma_dot, gamma_dot_p
      double precision :: bfs(3,3), umean(3)

!     Handy string array for making file names.
      character*19 strfname

      integer ios
      logical lex

      double precision :: b = $b$

      ALLOCATE (UR(0:NBX,0:NBY,0:NGZM1),VR(0:NBX,0:NBY,0:NGZM1))
      ALLOCATE (WR(0:NBX,0:NBY,0:NGZM1))
      ALLOCATE (QRFACT(0:NBX,0:NBY,0:NGZM1),PRDENO(0:NBX,0:NBY,0:NGZM1))
      ALLOCATE (VXFACT(0:NBX),VYFACT(0:NBY),VZFACT(0:NBZ))
      ALLOCATE (FIRSTN(1:NGX,1:NGY),NUMBER(1:NGX,1:NGY))
      ALLOCATE (NEXTN(1:NFSIZE),XFN(1:3,1:NFSIZE),FRC(1:3,1:NFSIZE))
      ALLOCATE (elmnew(1:3,1:NFSIZE2))
      ALLOCATE (shpint(1:3,1:NFSIZE2),shpfs(1:7,1:NFSIZE2))

      pi = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits
!     Physical parameters -- using cgs system
      nstep = $nstep$ ! Number of timesteps
      radx = 3.338d-4 ! cell radius (cm)

      if (b < 0) b = 2*(1+2*$gapratio$)*ngy/(ngy-3)
      lcube = radx*b ! Length of one edge of the fluid domain cube (cm)
      mu = 1.20d-2 ! Viscosity of plasma (poise), Eggleton & Popel 1998
      rho = 1.025d0 ! Density of plasma (g/(cm^3)), Neofytou 2004
      nu = mu/rho ! Kinematic density

      
      gamma_dot = Eh*capillary_no/(radx*mu)
!     Shear rate in program units
      gamma_dot_p = $gamma_dot_p$

      td = gamma_dot_p/gamma_dot ! Timestep (s)
      h = lcube/flngx ! Fluid node spacing

!     Establish conversion factors from cgs units to program units
!     Running the program in these units simplifies several calculations;
!     specifically, with these scalings, 
!     delta_u = F (not F/rho)
!     delta_x = u (not u*delta_t)
      length = h
      time = td
!     This characteristic mass is the mass of fluid in a cube with a
!     side length of h.
      mass = rho*h**3

!     
      vsc = nu/((h**2)/td)
!     h64 is used in converting from length in program units to
!     non-dimensional units; currently, only used in shape.
      h64 = h/radx
!     Characteristic force; divide by this to get force in program units.
      fostar = (mass*h/td**2)

!     Added to implement BC homogenization
!     These are in the program units for 1/T
!     Plane Shear
      if (1 == $flow$) then
         bfs(1,:) = (/0.d0, 0.d0, 0.d0/)
         bfs(2,:) = (/0.d0, 0.d0, gamma_dot_p/2.d0/)
         bfs(3,:) = (/0.d0, gamma_dot_p/2.d0, 0.d0/)
      end if
!     Axisymmetric Extension
      if (2 == $flow$) then
         bfs(1,:) = (/-gamma_dot_p/2.d0, 0.d0, 0.d0/)
         bfs(2,:) = (/0.d0, -gamma_dot_p/2.d0, 0.d0/)
         bfs(3,:) = (/0.d0, 0.d0, gamma_dot_p/)
      end if
!     Wall
      if (3 == $flow$) then
         bfs(1,:) = (/0.d0, 0.d0, 0.d0/)
         bfs(2,:) = (/0.d0, 0.d0, 0.d0/)
         bfs(3,:) = (/0.d0, gamma_dot_p, 0.d0/)
      end if
      if (4 == $flow$) then
         mix = $mix$
         bfs(1,:) = (/-mix*gamma_dot_p/2.d0, gamma_dot_p*(1-mix), 
     &        0.d0/)
         bfs(2,:) = (/gamma_dot_p*(1-mix), -mix*gamma_dot_p/2.d0, 
     &        0.d0/)
         bfs(3,:) = (/0.d0, 0.d0, gamma_dot_p*mix/)
      end if
!     In program units for velocity (unitless)
      if (3 /= $flow$) then
         umean(:)= (/0.d0,0.d0,0.d0/)
      end if

!     This makes the flow zero at the wall
      if (3 == $flow$) then
         umean(:) = (/0.d0, 0.d0, bfs(3,2)*((flngy+1.d0)/2.d0-planey)/)
      end if

      klok = 0
      inquire(100, exist=lex, iostat=ios, recl=i)
      open(100,iostat=ios, form='unformatted')
      read(100, err=37, end=37) klok
      close(100)
      go to 38
 37   continue
 38   continue

      if (klok == 0) then
!     Initialize the solid arrays
      call inspher(lcube,radx,h,xfn,elmnew,shpint,shpfs)
      if ($npls$ > 0) then
         call inplane(xpi, xfn)
      end if
      call cellcenter(klok, xfn, xcenter, ycenter, zcenter)
      xcenterold = xcenter
      ycenterold = ycenter
      zcenterold = zcenter
!  Initialize the activation variables --
      write(206,*) nstep  ,' = nstep'
      write(206,*) lcube  ,' cm = lcube'
      write(206,*) nu     ,' cm**2/sec = nu'
      write(206,*) rho    ,' gm/cm**3 = rho'
      write(206,*) ngx,lxng,' = ngx l2ngx'
      write(206,*) td     ,' sec = td'
      write(206,*) pi,' = pi'
      write(206,*) h ,' cm = h '
      write(206,*) mu,' (gm/cm**3)*(cm**2/sec) = mu'
      write(206,*)'conversion factors --'
      write(206,*) mass  ,' =mass'
      write(206,*) length,' =length'
      write(206,*) time  ,' =time'
      write(206,*) vsc   ,' =vsc=kinematic viscosity in program units'

!     Initialize velocity
!     Throughout, u is in program units (normalized by h/td)
      if (fvs /= 0) then
         call fvssub(ur, vr, wr, -bfs, -umean)
      end if

      call wprofile(wr, 0)
!      call uvwdump(ur, vr, wr, 0)
      call makefilename('solidnodes', 0,'.txt',strfname)
      call saveallsolid(XFN,strfname)
      call makefilename('solidforce', 0,'.txt',strfname)
      call saveallsolid(frc,strfname)
      call shape(lcube,h64,klok,td,m_start,m_end,1,nfsize2,xfn,elmnew)
      call calculateDF(klok, xfn, m_start, m_end)
      else
         write(*,*) 'cell l184 Restarting'
         call restart(lcube, nu, rho,td,klok,ur,vr,wr,
     &        xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs)
      end if
      call inhist(xfn,firstn,number,nextn)
!     Initialize the fluid solver
      call influidu(vsc,vxfact,vyfact,vzfact,prdeno,qrfact, dsq,dx,
     &     dy,dz)

      klok1 = klok + 1
      klokend = nstep + klok/nstep*nstep

!     MAIN LOOP --
      DO 5 KLOK=KLOK1,KLOKEND
         message = 'cell l216'
         call dumpstatus(klok, message)
      T = T+TD
         CALL MEMBNX(KLOK,XFN,elmnew,shpint,shpfs,FRC,h,FOSTAR,RADX)
         message = 'cell l220'
         call dumpstatus(klok, message)
         if ($npls$ > 0) then
            call pmhist(XPI,XFN,FRC,FIRSTN,NEXTN,NUMBER,
     &           10240.d0/dble($npl$))
         end if
         message = 'cell l225'
         call dumpstatus(klok, message)
         call meanforce(klok, frc)
         message = 'cell l228'
         call dumpstatus(klok, message)
         call cellcenter(klok, xfn, xcenter, ycenter, zcenter)
         message = 'cell l232'
         call dumpstatus(klok, message)
         open(402, access='append')
         write(402,*) klok, xcenter - xcenterold, 
     &        ycenter - ycenterold, zcenter - zcenterold
         close(402)
         xcenterold = xcenter
         ycenterold = ycenter
         zcenterold = zcenter
         message = 'cell l240'
         call dumpstatus(klok, message)
         CALL pushup(KLOK,UR,VR,WR,XFN,FRC,FIRSTN,NUMBER,NEXTN)
         message = 'cell l243'
         call dumpstatus(klok, message)
         if (FVS == 0) then
            call bodyfs(klok, bfs, ur, vr, wr, vsc)
         end if
         message = 'cell l248'
         call dumpstatus(klok, message)
         if (fvs /= 0) then
            call fvssub(ur, vr, wr, bfs, umean)
         end if
         CALL FLUIDUP(KLOK,UR,VR,WR,VXFACT,VYFACT,VZFACT,PRDENO,QRFACT,
     &     DSQ, DX, DY, DZ)
         if (fvs /= 0) then
            call fvssub(ur, vr, wr, -bfs, -umean)
         end if
         message = 'cell l252'
         call dumpstatus(klok, message)
         call wrap(klok, ur, vr, wr)
         message = 'cell l255'
         call dumpstatus(klok, message)
         CALL MOVE(KLOK,UR,VR,WR,XFN,FIRSTN,NUMBER,NEXTN)
         message = 'cell l258'
         call dumpstatus(klok, message)
         CALL SHAPE(LCUBE,h64,KLOK,TD,m_start,m_end,1,nfsize2,XFN,
     &     elmnew)
         message = 'cell l262'
         call dumpstatus(klok, message)
         call calculateDF(klok, xfn, m_start, m_end)
      WRITE(206,*)' KLOK: ',KLOK,  ' ; TIME: ',T
         message = 'cell l266'
         call dumpstatus(klok, message)
         call wrstart(lcube, nu, rho,td,klok,ur,vr,wr,
     &        xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs)
         message = 'cell l270'
         call dumpstatus(klok, message)
         if ((klok/100)*100 == klok) then
         message = 'cell l273'
         call dumpstatus(klok, message)
         call wprofile(wr, klok)
         message = 'cell l276'
         call dumpstatus(klok, message)
         end if
         message = 'cell l279'
         call dumpstatus(klok, message)

      if ((klok/500)*500 == klok) then
!         call uvwdump(ur, vr, wr, KLOK)
         message = 'cell l284'
         call dumpstatus(klok, message)
         call makefilename('solidnodes',KLOK,'.txt',strfname)
         message = 'cell l287'
         call dumpstatus(klok, message)
         call saveallsolid(XFN,strfname)
         message = 'cell l290'
         call dumpstatus(klok, message)
         call makefilename('solidforce',KLOK,'.txt',strfname)
         message = 'cell l293'
         call dumpstatus(klok, message)
         call saveallsolid(frc,strfname)
         message = 'cell l296'
         call dumpstatus(klok, message)
      end if
         message = 'cell l299'
         call dumpstatus(klok, message)

         if ($vstable$ == 1) then
            if ($flow$ /= 3) then
               call meanfluidvelocity(ur, meanstablev)
               ur = ur - meanstablev
               call meanfluidvelocity(vr, meanstablev)
               vr = vr - meanstablev
               call meanfluidvelocity(wr, meanstablev)
               wr = wr - meanstablev
            end if
         end if
   5  CONTINUE
   6  CONTINUE

      message = 'cell l304'
      call dumpstatus(klok, message)

      DEALLOCATE (UR,VR,WR,FIRSTN,NUMBER,NEXTN,XFN,FRC)
      DEALLOCATE(VXFACT,VYFACT,VZFACT)
      DEALLOCATE(PRDENO,QRFACT,elmnew,shpint,shpfs)
  
      END PROGRAM cell
!**********************************************************************
