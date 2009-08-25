!     Some modifications made by Alex Szatmary. He claims copyright on
!     these changes, which, unfortunately, are not, at the moment,
!     clearly marked, but he's a nice guy and is generally happy to 
!     share.
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
C  PARAMETERS: 
C     PLANES ARE NG X NG
C     BUT ARE STORED IN ARRAYS WITH DIMENSIONS
C     (0:NB,1:NG)   WHERE NB=NG+2
C
C *********************************************************************
C
c     NOTES ON CHANGES MADE TO PROGRAM PULSE3D TO CONVERT IT TO
c     RBC3D: SIMULATION CODE FOR THE STUDY OF ERYTHROCYTES IN 3D FLOW
c     NOTES STARTED ON 1/19/95
C
C**********************************************************************
      PROGRAM cell
      IMPLICIT NONE

      interface
         subroutine importmesh(LCUBE,RAD,H,XFN, elmnew,shpint,shpfs,
     &        my_nnode, my_nelm, my_cap_center)
         implicit none
         double precision lcube, h, pi, rad
         double precision :: XFN(:,:)
         INTEGER elmnew(:,:)
         double precision :: shpint(:,:),shpfs(:,:)
         integer my_nnode, my_nelm
         double precision my_cap_center(3)
         end subroutine

         subroutine cellcenter(klok, xfn, my_nnode, cap_i, xcenter, 
     &     ycenter, zcenter)
         implicit none
         integer klok, my_nnode, cap_i
         double precision :: xfn(:,:)
         double precision :: xcenter, ycenter, zcenter
         end subroutine cellcenter

         subroutine inplane(xpi, xfn)
         implicit none
         double precision :: xpi(:,:), xfn(:,:)
         integer nnodes
         end subroutine inplane

         subroutine pmhist(xpi,xfn,frc,firstn,nextn,number,const,
     &     fp_start, fp_end, nfsize)
         implicit none
         double precision :: xpi(:,:)
         double precision :: xfn(:,:),frc(:,:)
         integer firstn(:,:),number(:,:)
         integer nextn(:)
         double precision :: const
         integer fp_start, fp_end
         integer nfsize
         end subroutine pmhist

         subroutine pushup(klok,ur,vr,wr,xfn,frc,firstn,number,nextn)
         implicit none
         integer klok
         double complex ur(:,:,:), vr(:,:,:), wr(:,:,:)
         double precision :: xfn(:,:),frc(:,:)
         integer firstn(:,:),number(:,:),nextn(:)
         end subroutine pushup

         subroutine move(klok,ur,vr,wr,xfn,firstn,number,nextn)
         implicit none
         integer klok
         double complex :: ur(:,:,:), vr(:,:,:), wr(:,:,:)
         double precision :: XFN(:,:)
         integer firstn(:,:),number(:,:),nextn(:)
         end subroutine move
      end interface

!**********************************************************************
!     Start variable declaration
!**********************************************************************

      integer lxng,lyng,lzng,ngx,ngy,ngz,nfsize,nfsize2
      INTEGER NGXM1,NGYM1,NGZM1,NBX,NBY,NBZ
      double precision :: FLNGX,FLNGY,FLNGZ
      PARAMETER(LXNG=$lngx$,LYNG=LXNG,LZNG=LXNG)
      PARAMETER(NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG)
      PARAMETER(NBX=NGX+2,NBY=NGY+2,NBZ=NGZ+2)
      PARAMETER(NGXM1=NGX-1,NGYM1=NGY-1,NGZM1=NGZ-1)
      PARAMETER(FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ)
      integer, parameter :: npl=$npl$
      double precision fnpl
      INTEGER KLOK,KLOK1,KLOK0,KLOKEND,NSTEP
      double precision :: T,H,h64,TD,VSC,TIME,RHO,PI,RADX,FOSTAR
      double precision, allocatable :: rad(:)
      double precision, allocatable :: xcenter(:), ycenter(:), 
     &     zcenter(:)
      double precision, allocatable :: xcenterold(:), ycenterold(:), 
     &     zcenterold(:)
      integer ncap
      double precision ::  LCUBE,NU,MU,MASS,LENGTH
!     Velocities are always expressed in program units
      double COMPLEX,ALLOCATABLE :: UR(:,:,:),VR(:,:,:),WR(:,:,:)
      double complex :: pr(0:NBX,0:NBY,0:NGZM1)
      double complex :: meanstablev
      double complex :: meanu, meanv, meanw
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

!     cap_n_start is an array of the indices to xfn, the *n*ode indices
!     for each capsule. cap_n_end is the indices for the end of each
!     capsule's node index count. cap_e_start is similar, but for
!     element indexing, and cap_e_end likewise.
      integer, allocatable :: cap_n_start(:), cap_n_end(:)
      integer, allocatable :: cap_e_start(:), cap_e_end(:)
      character*80 message

!     Array of linked lists representing the solid nodes for a given 
!     hxh tile, spanning (x,y) to (x+h, y+h)
      INTEGER,ALLOCATABLE :: FIRSTN(:,:),NUMBER(:,:),NEXTN(:)
!     xfn is the x coordinates of the solid nodes, in program units
!     frc is the force exerted at each node on the fluid; this is in program
!     units.
!     shpint and shpfs are arrays of finite element shape factor parameters.
!     These are in real dimensions.
      double precision, ALLOCATABLE :: XFN(:,:),FRC(:,:),shpint(:,:),
     &     shpfs(:,:)
      double precision xpi(1:3,1:npl)
!     This array associates 3 nodes with a numbered element; three corners
!     on a triangle.
      INTEGER, ALLOCATABLE :: elmnew(:,:)

!     Changes made 10-19-06 to implement BC Homog
      integer i, j, k
      double precision :: bfs(3,3), umean(3)

!     Handy string array for making file names.
      character*19 strfname

      integer ios
      logical lex

      double precision :: b = $b$
      double precision :: mix

      integer, allocatable :: fineness(:), nnode(:), nelm(:)
      double precision :: cap_center(3,$ncap$)

      integer fp_start, fp_end

!**********************************************************************
!     End variable declaration, start real code
!**********************************************************************

      ncap=$ncap$
      allocate(rad(ncap), xcenter(ncap), ycenter(ncap), zcenter(ncap),
     &     xcenterold(ncap), ycenterold(ncap), zcenterold(ncap))
      allocate(cap_n_start(ncap), cap_n_end(ncap), cap_e_start(ncap),
     &     cap_e_end(ncap))
      allocate(fineness(ncap), nnode(ncap), nelm(ncap))

      fnpl = $npl$
      message = '               '

      fineness = $fineness$
      do i=1,$ncap$
         call capsuletable(fineness(i),nnode(i),nelm(i))
      end do
      call make_cap_start_and_end(nnode, cap_n_start, cap_n_end,
     &     nelm, cap_e_start, cap_e_end)

      nfsize=cap_n_end($ncap$)+$npl$
      nfsize2=cap_e_end($ncap$)

      fp_start=cap_n_end($ncap$)+1
      fp_end=nfsize

      ALLOCATE(UR(0:NBX,0:NBY,0:NGZM1),VR(0:NBX,0:NBY,0:NGZM1))
      ALLOCATE(WR(0:NBX,0:NBY,0:NGZM1))
      ALLOCATE(QRFACT(0:NBX,0:NBY,0:NGZM1),PRDENO(0:NBX,0:NBY,0:NGZM1))
      ALLOCATE(VXFACT(0:NBX),VYFACT(0:NBY),VZFACT(0:NBZ))
      ALLOCATE(FIRSTN(1:NGX,1:NGY),NUMBER(1:NGX,1:NGY))
      ALLOCATE(NEXTN(1:NFSIZE),XFN(1:3,1:NFSIZE),FRC(1:3,1:NFSIZE))
      ALLOCATE(elmnew(1:3,1:NFSIZE2))
      ALLOCATE(shpint(1:3,1:NFSIZE2),shpfs(1:7,1:NFSIZE2))


      pi = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits
!     Physical parameters -- using cgs system
      nstep = $nstep$ ! Number of timesteps
      radx = $radx$ ! cell radius (cm)

      lcube = $lcube$ ! Length of one edge of the fluid domain cube (cm)
      mu = $mu$
      rho = $rho$
      nu = $nu$ ! Kinematic viscosity
      
      td = $td$ ! Timestep (s)
      h = lcube/flngx ! Fluid node spacing (cm)

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

      vsc = nu/((h**2)/td)
!     h64 is used in converting from length in program units to
!     non-dimensional units; currently, only used in shape.
      h64 = h/radx
!     Characteristic force; divide by this to get force in program units.
      fostar = (mass*h/td**2)

      rad = $rad$
      cap_center(1,:)=$xc_cap$
      cap_center(2,:)=$yc_cap$
      cap_center(3,:)=$zc_cap$

!     Added to implement BC homogenization
!     These are in the program units for 1/T
!     bfs is the velocity gradient matrix, \nabla u (no dot)
!     This can be applied to FVS or to the bodyfs subroutine.
!     FVS allows for direct imposition of a linear flow field.
!     bodyfs imposes stresses on the faces of the domain; these stresses are
!     picked to generate the same flow field that FVS would, if there were
!     no immersed body.
      mix = $mix$
      bfs(1,:) = $bfs1$
      bfs(2,:) = $bfs2$
      bfs(3,:) = $bfs3$

!     This is the mean fluid velocity vector. For unbounded flows, this is
!     typically zero, so that a capsule in the center is immobilized.
!     In program units for velocity (unitless)
      umean(:)= $umean$

!     This automatically restarts an aborted run.
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
         do i = 1,$ncap$
            write(*,*) 'l290', cap_n_start(i), cap_n_end(i), 
     &           cap_e_start(i),cap_e_end(i), cap_center(:,i)
            call importmesh(lcube,rad(i),h,
     &           xfn(1:3,cap_n_start(i):cap_n_end(i)),
     &           elmnew(1:3,cap_e_start(i):cap_e_end(i)),
     &           shpint(1:3,cap_e_start(i):cap_e_end(i)),
     &           shpfs(1:7,cap_e_start(i):cap_e_end(i)), 
     &           nnode(i), nelm(i), cap_center(:,i))
            write(*,*) 'cell l309'
         end do
!     $npls$ is the number of planes. If there is one, it should be
!     initialized.
         if ($npls$ > 0) then
            call inplane(xpi, xfn(1:3,fp_start:fp_end))
         end if
         write(*,*) 'cell l316'
!     Get an initial measurement of the center of the capsule
         do i=1,$ncap$
            call cellcenter(klok,xfn(1:3,cap_n_start(i):cap_n_end(i)), 
     &           nnode(i), i, xcenter(i), ycenter(i), zcenter(i))
            write(*,*) 'cell l319'
!     Use these values to measure velocity. It's a crappy measure, it's
!     a backward difference.
            xcenterold(i) = xcenter(i)
            ycenterold(i) = ycenter(i)
            zcenterold(i) = zcenter(i)
         end do

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
         write(206,*) vsc   ,' =vsc=kinematic viscosity, program units'

!     Initialize velocity
!     Throughout, u is in program units (normalized by h/td)
         if (fvs /= 0) then
            call fvssub(ur, vr, wr, -bfs, -umean)
            call poiseuille(ur, vr, wr, pr, 
     &           $dpdz$, vsc)
         end if
         write(*,*) 'cell l352'

         call wprofile(wr, 0)
         write(*,*) 'cell l355'
         call uvwpdump(ur, vr, wr, pr, 0)
         call makefilename('solidnodes', 0,'.txt',strfname)
         call saveallsolid(XFN,strfname)
         call makefilename('solidforce', 0,'.txt',strfname)
         call saveallsolid(frc,strfname)
         do i = 1,$ncap$
            call shape(lcube,h64,klok,td,cap_n_start(i),cap_n_end(i),
     &           cap_e_start(i), cap_e_end(i),xfn,elmnew)
            call calculateDF(klok, i, xfn, cap_n_start(i), 
     &           cap_n_end(i))
         end do
         t=0.d0
      else
         write(*,*) 'cell l367 Restarting'
         call restart(lcube, nu, rho,td,klok,ur,vr,wr,
     &        xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs)
         T=klok*time
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
         call dumpstatus(klok, message, 'status.txt')
         T = T+TD
         do i=1,$ncap$
            CALL MEMBNX(KLOK,XFN,elmnew,shpint,shpfs,FRC,h,FOSTAR,
     &           RAD(i), cap_n_start(i), cap_n_end(i))
            write(message, *) 'frc(1,1)',frc(1,1)
            call dumpstatus(klok, message, 'thumbprint.txt')
            write(message, *) 'frc(3,nfsize)',frc(1,nfsize)
            call dumpstatus(klok, message, 'thumbprint.txt')
            write(message, *) 'frc(2,nfsize)',frc(1,nfsize/2)
            call dumpstatus(klok, message, 'thumbprint.txt')
         end do

         message = 'cell l220'
         call dumpstatus(klok, message, 'status.txt')
         if ($npls$ > 0) then
            call pmhist(xpi,xfn,frc,firstn,nextn,number,
     &           10240.d0/fnpl, fp_start, fp_end, nfsize)
         end if
         message = 'cell l225'
         call dumpstatus(klok, message, 'status.txt')
         call meanforce(klok, frc)
         message = 'cell l228'
         call dumpstatus(klok, message, 'status.txt')
         do i=1,$ncap$
            call cellcenter(klok,xfn(1:3,cap_n_start(i):cap_n_end(i)), 
     &           nnode(i), i, xcenter(i), ycenter(i), zcenter(i))
            message = 'cell l232'
            call dumpstatus(klok, message, 'status.txt')
!     todo Make filename change with capsule index
            call makefilename('capsulev__', i,'.txt',strfname)
            open(402,file=strfname, access='append')
            write(402,*) klok, xcenter(i) - xcenterold(i),
     &           ycenter(i) - ycenterold(i), zcenter(i) - zcenterold(i)
            close(402)
            xcenterold(i) = xcenter(i)
            ycenterold(i) = ycenter(i)
            zcenterold(i) = zcenter(i)
         end do

         call meanfluidvelocity(ur, meanu)
         call meanfluidvelocity(vr, meanv)
         call meanfluidvelocity(wr, meanw)
         open(403, access='append')
         write(403,*) klok, dreal(meanu), dreal(meanv), dreal(meanw)
         close(403)
         message = 'cell l240'
         call dumpstatus(klok, message, 'status.txt')
         CALL pushup(KLOK,UR,VR,WR,XFN,FRC,FIRSTN,NUMBER,NEXTN)
         message = 'cell l243'
         call dumpstatus(klok, message, 'status.txt')
         write(message, *) 'ur(3*ngx/4,ngy/2,ngz/4)', 
     &        ur(3*ngx/4,ngy/2,ngz/4)
         call dumpstatus(klok, message, 'thumbprint.txt')
         write(message, *) 'vr(3*ngx/4,ngy/2,ngz/4)',
     &        vr(3*ngx/4,ngy/2,ngz/4)
         call dumpstatus(klok, message, 'thumbprint.txt')
         write(message, *) 'wr(3*ngx/4,ngy/2,ngz/4)',
     &        wr(3*ngx/4,ngy/2,ngz/4)
         call dumpstatus(klok, message, 'thumbprint.txt')
         if (FVS == 0) then
            call bodyfs(klok, bfs, ur, vr, wr, vsc)
         end if
         message = 'cell l248'
         call dumpstatus(klok, message, 'status.txt')
         if (fvs /= 0) then
            call fvssub(ur, vr, wr, bfs, umean)
            call poiseuille(ur, vr, wr, pr,
     &           -$dpdz$, vsc)
         end if
         CALL FLUIDUP(KLOK,UR,VR,WR, pr, VXFACT,VYFACT,VZFACT,
     &        PRDENO,QRFACT, DSQ, DX, DY, DZ)
         if (fvs /= 0) then
            call fvssub(ur, vr, wr, -bfs, -umean)
            call poiseuille(ur, vr, wr, pr, 
     &           $dpdz$, vsc)
         end if
         message = 'cell l252'
         call dumpstatus(klok, message, 'status.txt')
         call wrap(klok, ur, vr, wr)
         write(message, *) 'ur(3*ngx/4,ngy/2,ngz/4)', 
     &        ur(3*ngx/4,ngy/2,ngz/4)
         call dumpstatus(klok, message, 'thumbprint.txt')
         write(message, *) 'vr(3*ngx/4,ngy/2,ngz/4)',
     &        vr(3*ngx/4,ngy/2,ngz/4)
         call dumpstatus(klok, message, 'thumbprint.txt')
         write(message, *) 'wr(3*ngx/4,ngy/2,ngz/4)',
     &        wr(3*ngx/4,ngy/2,ngz/4)
         call dumpstatus(klok, message, 'thumbprint.txt')
         message = 'cell l255'
         call dumpstatus(klok, message, 'status.txt')
         CALL MOVE(KLOK,UR,VR,WR,XFN,FIRSTN,NUMBER,NEXTN)
         message = 'cell l258'
         call dumpstatus(klok, message, 'status.txt')
         write(message, *) 'xfn(1,1)',xfn(1,1)
         call dumpstatus(klok, message, 'thumbprint.txt')
         write(message, *) 'xfn(3,nfsize)',xfn(3,nfsize)
         call dumpstatus(klok, message, 'thumbprint.txt')
         write(message, *) 'xfn(2,nfsize/2)',xfn(2,nfsize/2)
         call dumpstatus(klok, message, 'thumbprint.txt')
!     End of the loop, do post-processing stuff
         do i=1,$ncap$
            CALL SHAPE(LCUBE,h64,KLOK,TD,cap_n_start(i),cap_n_end(i),1,
     &           nfsize2,XFN, elmnew)
            message = 'cell l262'
            call dumpstatus(klok, message, 'status.txt')
            call calculateDF(klok, i, xfn,cap_n_start(i), cap_n_end(i))
         end do
         WRITE(206,*)' KLOK: ',KLOK,  ' ; TIME: ',T
         message = 'cell l266'
         call dumpstatus(klok, message, 'status.txt')
         call wrstart(lcube, nu, rho,td,klok,ur,vr,wr,
     &        xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs)
         message = 'cell l270'
         call dumpstatus(klok, message, 'status.txt')
         if ((klok/$smalldumpint$)*$smalldumpint$ == klok) then
            message = 'cell l273'
            call dumpstatus(klok, message, 'status.txt')
            call wprofile(wr, klok)
            message = 'cell l276'
            call dumpstatus(klok, message, 'status.txt')
         end if
         message = 'cell l279'
         call dumpstatus(klok, message, 'status.txt')

         if ((klok/$bigdumpint$)*$bigdumpint$ == klok) then
            call uvwpdump(ur, vr, wr, pr, KLOK)
            message = 'cell l284'
            call dumpstatus(klok, message, 'status.txt')
            call makefilename('solidnodes',KLOK,'.txt',strfname)
            message = 'cell l287'
            call dumpstatus(klok, message, 'status.txt')
            call saveallsolid(XFN,strfname)
            message = 'cell l290'
            call dumpstatus(klok, message, 'status.txt')
            call makefilename('solidforce',KLOK,'.txt',strfname)
            message = 'cell l293'
            call dumpstatus(klok, message, 'status.txt')
            call saveallsolid(frc,strfname)
            message = 'cell l296'
            call dumpstatus(klok, message, 'status.txt')
      end if
         message = 'cell l299'
         call dumpstatus(klok, message, 'status.txt')

C          if ($vstable$ == 1) then
C             if ($flow$ /= 3) then
C                call meanfluidvelocity(ur, meanstablev)
C                ur = ur - meanstablev
C                call meanfluidvelocity(vr, meanstablev)
C                vr = vr - meanstablev
C                call meanfluidvelocity(wr, meanstablev)
C                wr = wr - meanstablev
C             end if
C          end if
   5  CONTINUE
   6  CONTINUE

      message = 'cell l304'
      call dumpstatus(klok, message, 'status.txt')

      DEALLOCATE (UR,VR,WR,FIRSTN,NUMBER,NEXTN,XFN,FRC)
      DEALLOCATE(VXFACT,VYFACT,VZFACT)
      DEALLOCATE(PRDENO,QRFACT,elmnew,shpint,shpfs)
      deallocate(fineness, nnode, nelm)
      deallocate(cap_n_start, cap_n_end, cap_e_start,
     &     cap_e_end)
      deallocate(rad, xcenter, ycenter, zcenter,
     &     xcenterold, ycenterold, zcenterold)
      END PROGRAM cell
!**********************************************************************
