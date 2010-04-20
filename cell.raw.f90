!     Some modifications made by Alex Szatmary. He claims copyright on
!     these changes, which, unfortunately, are not, at the moment,
!     clearly marked, but he's a nice guy and is generally happy to 
!     share.
! *********************************************************************
! *     Copyright 1994 Charles S. Peskin and David M. McQueen         *
! *     All rights reserved. No portion of this program may be        *
! *     transmitted to any persons without the prior permission       *
! *     of the Copyright holders.                                     *
! *     Contact: Charles S. Peskin   email:  peskin@cims.nyu.edu      *
! *              David M. McQueen    email: mcqueen@cims.nyu.edu      *
! *********************************************************************
!
!**********************************************************************
!  PARAMETERS: 
!     PLANES ARE NG X NG
!     BUT ARE STORED IN ARRAYS WITH DIMENSIONS
!     (0:NB,1:NG)   WHERE NB=NG+2
!
! *********************************************************************
!
!     NOTES ON CHANGES MADE TO PROGRAM PULSE3D TO CONVERT IT TO
!     RBC3D: SIMULATION CODE FOR THE STUDY OF ERYTHROCYTES IN 3D FLOW
!     NOTES STARTED ON 1/19/95
!
!**********************************************************************
PROGRAM cell
  IMPLICIT NONE

  ! This interface block is needed to pass dynamically allocated arrays
  ! Useful comments will be found near the actual subroutines.
  interface
     subroutine inhist(xfn, firstn, number, nextn)
       double precision :: xfn(:,:)
       integer firstn(:,:),number(:,:),nextn(:)       
     end subroutine inhist

     subroutine generatecapsule(my_ellipa, my_ellipb, my_ellipc, H, XFN, &
          elmnew, shpint, shpfs, &
          my_cap_center, my_fineness, cap_i, a_prestress, my_nvec_i)
       implicit none
       double precision :: h
       double precision :: xfn(:,:), shpint(:,:), shpfs(:,:)
       integer :: elmnew(:,:)
       double precision :: my_ellipa(:), my_ellipb(:), my_ellipc(:)
       double precision :: my_cap_center(:), a_prestress
       integer :: my_fineness, cap_i
       integer :: my_nvec_i(:)
     end subroutine generatecapsule

     subroutine cellcenter(klok, xfn, my_nnode, cap_i, xcenter, &
          ycenter, zcenter) 
       implicit none 
       integer klok, my_nnode, cap_i 
       double precision :: xfn(:,:)
       double precision :: xcenter, ycenter, zcenter 
     end subroutine cellcenter

     subroutine inrect(xpi, xfn, my_rect_n1, my_rect_n2, my_recty)
       implicit none
       double precision :: xpi(:,:), xfn(:,:)
       integer my_rect_n1, my_rect_n2
       double precision :: my_recty
     end subroutine inrect

     subroutine pmhist(xpi,xfn,frc,const)
       implicit none
       double precision :: xpi(:,:)
       double precision :: xfn(:,:),frc(:,:)
       double precision :: const
     end subroutine pmhist

     subroutine pushup(ur,vr,wr,xfn,frc,firstn,number,nextn)
       implicit none
       double complex ur(:,:,:), vr(:,:,:), wr(:,:,:)
       double precision :: xfn(:,:),frc(:,:)
       integer firstn(:,:),number(:,:),nextn(:)
     end subroutine pushup

     subroutine move(ur,vr,wr,xfn,firstn,number,nextn)
       implicit none
       double complex :: ur(:,:,:), vr(:,:,:), wr(:,:,:)
       double precision :: XFN(:,:)
       integer firstn(:,:),number(:,:),nextn(:)
     end subroutine move

     subroutine shape(h64, klok, cap_i, xfn, nnode, elmnew, nelm)
       implicit none
       double precision :: h64, xfn(:,:)
       integer :: klok, nnode, elmnew(:,:), nelm, cap_i
     end subroutine shape

     subroutine calculateDF(clock, cap_i, xfn, my_nnode, lambda1, lambda2, &
          xcenter, ycenter, zcenter)
       implicit none
       integer clock, cap_i, my_nnode
       double precision :: xfn(:,:)
       double precision :: lambda1(:), lambda2(:)
       double precision :: xcenter, ycenter, zcenter
     end subroutine calculateDF

     subroutine profile(cap_i, xfn, clock)
       implicit none
       integer cap_i, clock
       double precision :: xfn(:,:)
     end subroutine profile

     subroutine membnx(xfn, elmnew, shpint, shpfs, FRC,h,FOSTAR, &
          lambda1, lambda2)
       implicit none
       double precision :: xfn(:,:), shpint(:,:),&
            shpfs(:,:), frc(:,:), h, fostar
       double precision :: lambda1(:), lambda2(:)
       integer :: elmnew(:,:)
     end subroutine membnx

     subroutine restart(lcube,nu,rho,td,ur,vr,wr,pr, &
          xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs, &
          xcenterold, ycenterold, zcenterold, nsph, nellip, ellipa, ellipb, &
          ellipc, a_prestress, nvec_i)
       implicit none
       double precision :: lcube,nu,td,rho
       double COMPLEX :: UR(:,:,:)
       double complex :: VR(:,:,:)
       double COMPLEX :: WR(:,:,:)
       double COMPLEX :: pr(:,:,:)
       double precision :: xfn(:,:)
       double precision :: xpi(:,:)
       integer firstn(:,:),number(:,:),nextn(:)
       integer elmnew(:,:)
       double precision :: shpint(:,:),shpfs(:,:)
       double precision :: xcenterold(:), ycenterold(:), zcenterold(:)
       integer nsph, nellip
       double precision :: ellipa(:), ellipb(:), ellipc(:), a_prestress(:)
       integer :: nvec_i(:,:)
     end subroutine restart

     subroutine wrstart(lcube,nu,rho,td,klok,ur,vr,wr,pr, &
          xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs, &
          xcenterold, ycenterold, zcenterold, nsph, nellip, ellipa, ellipb, &
          ellipc, a_prestress, nvec_i)
       implicit none
       double precision :: lcube,nu,td,rho
       integer klok
       double COMPLEX :: UR(:,:,:)
       double complex :: VR(:,:,:)
       double COMPLEX :: WR(:,:,:)
       double COMPLEX :: pr(:,:,:)
       double precision :: xfn(:,:)
       double precision :: xpi(:,:)
       integer firstn(:,:),number(:,:),nextn(:)
       integer elmnew(:,:)
       double precision :: shpint(:,:),shpfs(:,:)       
       double precision :: xcenterold(:), ycenterold(:), zcenterold(:)
       integer nsph, nellip
       double precision :: ellipa(:), ellipb(:), ellipc(:), a_prestress(:)
       integer :: nvec_i(:,:)
     end subroutine wrstart

     subroutine saveallsolid(xfn, strfname)
       character(len=*) strfname
       double precision :: xfn(:,:)
     end subroutine saveallsolid

     subroutine meanforce(clock, frc)
       integer clock
       double precision :: frc(:,:)
     end subroutine meanforce

     subroutine make_index_table_start_end(n, start, end)
       implicit none
       integer :: n(:), start(:), end(:)
     end subroutine make_index_table_start_end

     subroutine minmaxxfn(clock, xfn, cap_i)
       implicit none
       double precision :: xfn(:,:)
       integer :: cap_i, clock
     end subroutine minmaxxfn

     subroutine volume_area(klok,xfn,elmnew, cap_i)
       integer klok, cap_i
       double precision :: xfn(:,:)
       integer elmnew(:,:)
     end subroutine volume_area


     subroutine dumplambdas(lambda1, lambda2, strfname)
       double precision :: lambda1(:), lambda2(:)
       character(len=*) strfname
     end subroutine dumplambdas

     subroutine dumpnvec(clock, xfn, xcenter, ycenter, zcenter, my_nvec_i, &
          cap_i)
       integer clock
       double precision :: xfn(:,:), xcenter, ycenter, zcenter
       integer my_nvec_i(:)
       integer cap_i
     end subroutine dumpnvec
  end interface

  !**********************************************************************
  !     Start variable declaration
  !**********************************************************************


  ! ngx_im1 is the number of grid nodes in the i-direction, minus one.
  ! This is only really useful for the z-direction.
  ! nbx_i is ngx_i+2, which is only really useful in the x and y directions.
  ! lx_ing is the log_2(ngx_i)
  integer, parameter :: lxng=$lngx$, lyng=lxng, lzng=lxng
  ! ngx_i is the number of fluid grid nodes in the i-direction
  integer, parameter :: NGX=2**LXNG,NGY=2**LYNG,NGZ=2**LZNG
  integer, parameter :: nbx=ngx+2,nby=ngy+2,nbz=ngz+2
  integer, parameter :: ngxm1=ngx-1,ngym1=ngy-1,ngzm1=ngz-1

  ! These are just the ngx_i in double precision
  double precision, PARAMETER :: FLNGX=NGX,FLNGY=NGY,FLNGZ=NGZ

  ! nfsize is the number of solid nodes
  ! nfsize2 is the number of elements
  integer nfsize, nfsize2

  integer nrects
  integer nrectnodes
  double precision fnrectnodes
  INTEGER KLOK,KLOK1,KLOKEND,NSTEP
  double precision :: T,H,h64,TD,VSC,TIME,RHO,PI,RADX,FOSTAR
  double precision, allocatable :: rad(:), a_prestress(:)
  ! The lengths of the semi-axes of all the capsules
  double precision, allocatable :: ellipa(:), ellipb(:), ellipc(:)
  ! The orientations of the semi-axes of all the capsules
  double precision, allocatable :: xcenter(:), ycenter(:), &
       zcenter(:)
  double precision, allocatable :: xcenterold(:), ycenterold(:), & 
       zcenterold(:)
  ! This is a 6xncap array of integers. The integers are indices for xfn, they
  ! represent nodes, the six nodes that are at the tips of the semi-major
  ! axes of an ellipse. For a sphere, these indices correspond to the points
  ! at the +/- extremes in the x, y, z directions of the sphere before the flow
  ! is applied.
  integer, allocatable :: nvec_i(:,:)
  ! ncap is the total number of elastic capsules: all elastic capsules are
  ! modeled as ellipsoids
  integer ncap
  ! nsph is the number of spheres input from the bp files, nellip is the number
  ! of ellipsoids input from the bp files. ncap=nsph+nellip, and all of
  ! the arrays that talk about the parameters of the capsules have entries for
  ! the spheres, then for the ellipsoids.
  integer nsph, nellip
  double precision ::  LCUBE,NU,MU,MASS,LENGTH
  !     Velocities are always expressed in program units
  double COMPLEX,ALLOCATABLE :: UR(:,:,:),VR(:,:,:),WR(:,:,:)
  double complex :: pr(0:NBX,0:NBY,0:NGZM1)
  double complex :: meanu, meanv, meanw
  !     Handy look-up tables for fluid mechanics calculations.
  double precision,ALLOCATABLE :: QRFACT(:,:,:)
  double COMPLEX :: DX(0:NBX,0:NBY,0:NGZM1)
  double complex :: DY(0:NBX,0:NBY,0:NGZM1)
  double COMPLEX :: DZ(0:NBX,0:NBY,0:NGZM1)
  double precision :: DSQ(0:NBX,0:NBY,0:NGZM1)
  double precision, parameter :: Eh = $Eh$

  integer, parameter :: FVS = $FVS$

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
  double precision, ALLOCATABLE :: XFN(:,:),FRC(:,:),Foptical(:,:),shpint(:,:), &
       shpfs(:,:), rays(:,:)
  double precision, allocatable :: lambda1(:), lambda2(:)
  double precision, allocatable :: xpi(:,:)
  !     This array associates 3 nodes with a numbered element; three corners
  !     on a triangle.
  INTEGER, ALLOCATABLE :: elmnew(:,:)

  !     Changes made 10-19-06 to implement BC Homog
  integer i
  double precision :: bfs(3,3), umean(3)

  !     Handy string array for making file names.
  character*19 strfname

  integer ios
  logical lex

  integer, allocatable :: fineness(:), nnode(:), nelm(:)
  double precision, allocatable :: cap_center(:,:)

  integer, allocatable :: rect_n1(:), rect_n2(:), rect_nnode(:), &
       rect_n_start(:), rect_n_end(:)
  double precision, allocatable :: recty(:)

  !**********************************************************************
  !     Optical variables declaration
  !**********************************************************************
  double precision, parameter :: nm = $nm$ , np = $np$ , omega_zero = $omega_zero$, epsilono = $epsilono$ 
  double precision, parameter :: lambda= $lambda$ , opticalPower = $opticalPower$, lightSpeed = $lightSpeed$ 
  integer, parameter ::  numberOfMaxReflections = $numberOfMaxReflections$, index = $index$
  double precision :: disp  , z0
  integer numberOfRays
  !  integer sysclockstart, sysclockend, rate

  !**********************************************************************
  !     End variable declaration, start real code
  !**********************************************************************

  nfsize = 0
  nfsize2 = 0

  ncap=$ncap$
  nsph = $nsph$
  nellip = $nellip$
  if (ncap > 0) then
     if (nsph > 0) allocate(rad(nsph))
     allocate(ellipa(3*ncap), ellipb(3*ncap), ellipc(3*ncap))
     allocate(xcenter(ncap), ycenter(ncap), zcenter(ncap), &
          xcenterold(ncap), ycenterold(ncap), zcenterold(ncap))
     allocate(cap_n_start(ncap), cap_n_end(ncap), cap_e_start(ncap), &
          cap_e_end(ncap))
     allocate(fineness(ncap), nnode(ncap), nelm(ncap), cap_center(3,ncap))
     allocate(a_prestress(ncap), nvec_i(6, ncap))
     fineness = $fineness$
     do i=1,ncap
        call capsuletable(fineness(i),nnode(i),nelm(i))
     end do
     cap_n_start(1) = 1
     call make_index_table_start_end(nnode, cap_n_start, cap_n_end)

     cap_e_start(1) = 1
     call make_index_table_start_end(nelm, cap_e_start, cap_e_end)
     nfsize = nfsize + cap_n_end(ncap)
     nfsize2 = nfsize2 + cap_e_end(ncap)
  end if
  print *, 'cell l278', nfsize, nfsize2
  nrects = $npls$
  nrectnodes = 0
  if (nrects > 0) then 
     allocate(rect_nnode(nrects), rect_n1(nrects), rect_n2(nrects), &
          recty(nrects), rect_n_start(nrects), rect_n_end(nrects))
     rect_n1 = (/$pl_n1$/)
     rect_n2 = (/$pl_n2$/)
     recty = (/$recty$/)
     do i=1,nrects
        call rectangle_table(rect_n1(i), rect_n2(i), rect_nnode(i))
     end do
     if (ncap > 0) then
        rect_n_start(1) = cap_n_end(ncap)+1
     else
        rect_n_start(1) = 0   
     end if

     call make_index_table_start_end(rect_nnode, rect_n_start, rect_n_end)

     nrectnodes = rect_n_end(nrects) - rect_n_start(1) + 1
     fnrectnodes = nrectnodes
     nfsize = rect_n_end(nrects)
  end if

  print *, 'cell l298', nfsize, nfsize2
  do i = 1, ncap
     print *, i, nnode(i), cap_n_start(i), cap_n_end(i)
     print *, i, nelm(i), cap_e_start(i), cap_e_end(i)
  end do

  do i = 1, nrects
     print *, 'rectangle', i, rect_nnode(i), rect_n_start(i), rect_n_end(i)
  end do
  message = '               '


  ALLOCATE(UR(0:NBX,0:NBY,0:NGZM1),VR(0:NBX,0:NBY,0:NGZM1))
  ALLOCATE(WR(0:NBX,0:NBY,0:NGZM1))
  ALLOCATE(QRFACT(0:NBX,0:NBY,0:NGZM1))
  ALLOCATE(FIRSTN(1:NGX,1:NGY),NUMBER(1:NGX,1:NGY))
  ALLOCATE(NEXTN(1:NFSIZE),XFN(1:3,1:NFSIZE),FRC(1:3,1:NFSIZE),Foptical(1:3,1:NFSIZE))
  ALLOCATE(elmnew(1:3,1:NFSIZE2))
  ALLOCATE(shpint(1:3,1:NFSIZE2),shpfs(1:7,1:NFSIZE2))
  allocate(lambda1(nfsize2),lambda2(nfsize2))
  allocate(xpi(3,nrectnodes))

  pi = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits
  !     Physical parameters -- using cgs system
  z0 = PI*omega_zero*omega_zero/lambda
  nstep = $nstep$ ! Number of timesteps
  radx = $radx$ ! cell radius (cm)
  disp = index * RADX*1e-2/30

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

  if (ncap > 0) then
     ! All objects described as spheres are re-posed as ellipsoids.
     if (nsph > 0) then
        rad = $rad$
        do i = 1,nsph
           ellipa(1+(i-1)*3:3+(i-1)*3) = rad(i)*(/1.d0, 0.d0, 0.d0/)
           ellipb(1+(i-1)*3:3+(i-1)*3) = rad(i)*(/0.d0, 1.d0, 0.d0/)
           ellipc(1+(i-1)*3:3+(i-1)*3) = rad(i)*(/0.d0, 0.d0, 1.d0/)
        end do
     end if

     if (nellip > 0) then
        ellipa(nsph*3+1:ncap*3) = $ellipa$
        ellipb(nsph*3+1:ncap*3) = $ellipb$
        ellipc(nsph*3+1:ncap*3) = $ellipc$
     end if

     cap_center(1,:)=$xc_cap$
     cap_center(2,:)=$yc_cap$
     cap_center(3,:)=$zc_cap$
     a_prestress = $a_prestress$
  end if

  !todo take these lines out once ellipsoid code is validated
  print *, "fineness", fineness
  print *, "ellipa", ellipa
  print *, "ellipb", ellipb
  print *, "ellipc", ellipc
  print *, "cap_center(1,:)", cap_center(1,:)
  print *, "cap_center(2,:)", cap_center(2,:)
  print *, "cap_center(3,:)", cap_center(3,:)

  !     Added to implement BC homogenization
  !     These are in the program units for 1/T
  !     bfs is the velocity gradient matrix, \nabla u (no dot)
  !     This can be applied to FVS or to the bodyfs subroutine.
  !     FVS allows for direct imposition of a linear flow field.
  !     bodyfs imposes stresses on the faces of the domain; these stresses are
  !     picked to generate the same flow field that FVS would, if there were
  !     no immersed body.
  bfs(1,:) = $bfs1$*td
  bfs(2,:) = $bfs2$*td
  bfs(3,:) = $bfs3$*td

  !     This is the mean fluid velocity vector. For unbounded flows, this is
  !     typically zero, so that a capsule in the center is immobilized.
  !     In program units for velocity (unitless)
  umean(:)= $umean$*td/h

  !     This automatically restarts an aborted run.
  klok = 0
  inquire(100, exist=lex, iostat=ios, recl=i)
  open(100,iostat=ios, form='unformatted')
  read(100, err=37, end=37) klok
  close(100)
  go to 38
37 continue
38 continue

  klok1 = klok + 1
  klokend = min(nstep + klok/nstep*nstep, $nend$)

  if (klok >= klokend) then
     print *, 'Maximum number of runs done; not restarting. Terminating job.'
     stop
  end if

  if (klok == 0) then
     !     Initialize the solid arrays
     do i = 1,ncap
        write(*,*) 'l290', cap_n_start(i), cap_n_end(i), &
             cap_e_start(i),cap_e_end(i), cap_center(:,i)
        call generatecapsule(ellipa((i-1)*3+1:(i-1)*3+3), &
             ellipb((i-1)*3+1:(i-1)*3+3), &
             ellipc((i-1)*3+1:(i-1)*3+3), &
             h, &
             xfn(1:3,cap_n_start(i):cap_n_end(i)), &
             elmnew(1:3,cap_e_start(i):cap_e_end(i)), &
             shpint(1:3,cap_e_start(i):cap_e_end(i)), &
             shpfs(1:7,cap_e_start(i):cap_e_end(i)), &
             cap_center(:,i), fineness(i), i, a_prestress(i), nvec_i(:, i))
        ! Get an initial measurement of the center of the capsule
        ! This is so that the first capsule center velocity measurement has
        ! some physical meaning.
        call cellcenter(klok,xfn(1:3,cap_n_start(i):cap_n_end(i)), &
             nnode(i), i, xcenter(i), ycenter(i), zcenter(i))
        !     Use these values to measure velocity. It's a crappy measure, it's
        !     a backward difference.
        xcenterold(i) = xcenter(i)
        ycenterold(i) = ycenter(i)
        zcenterold(i) = zcenter(i)
        ! Save the n vectors, the six vectors from the center of the capsule
        ! to the six tracer points.
        call dumpnvec(klok, xfn(:,cap_n_start(i):cap_n_end(i)), xcenter(i), &
             ycenter(i), zcenter(i), nvec_i(:, i), i)
        call shape(h64,klok, i, xfn(1:3, cap_n_start(i): &
             cap_n_end(i)), &
             nnode(i), elmnew(1:3, cap_e_start(i):cap_e_end(i)), &
             nelm(i))
        call calculateDF(klok, i, xfn(1:3, cap_n_start(i): &
             cap_n_end(i)), nnode(i), &
             lambda1(cap_e_start(i):cap_e_end(i)), &
             lambda2(cap_e_start(i):cap_e_end(i)), &
             xcenter(i), ycenter(i), zcenter(i))
        call profile(i, xfn(1:3,cap_n_start(i):cap_n_end(i)), &
             klok)
     end do
     !     nrects is the number of rectangles. If there is one, it should be
     !     initialized.
     do i= 1, nrects
        call inrect(xpi, xfn(1:3,rect_n_start(i):rect_n_end(i)), rect_n1(i), &
             rect_n2(i), recty(i))
     end do


     !  Initialize the activation variables --
     open(206,file='checkinit.txt',access='append')
     write(206,*) $nscale$  ,' = nscale'
     write(206,*) $nend$  ,' = nend'
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
     close(206)
     !     Initialize velocity
     !     Throughout, u is in program units (normalized by h/td)
     if (fvs /= 0) then
        call fvssub(ur, vr, wr, -bfs, -umean)
        call poiseuille(wr, pr, $dpdz$, vsc)
     end if

     call wprofile(wr, 0)
     call uvwpdump(ur, vr, wr, pr, 0)
     call makefilename('solidnodes', 0,'.txt',strfname)
     call saveallsolid(XFN,strfname)
     call makefilename('solidforce', 0,'.txt',strfname)
     call saveallsolid(frc,strfname)

     t=0.d0

     if ($optical$ /= 0) then
        !        CALL findRays (XFN , elmnew ,zcenter(1), numberOfrays )
        ALLOCATE (rays(1:9,1:numberOfRays))
        !        CALL initializeRays(rays, numberOfRays)
        !        call capsuleForce(XFN , Foptical, shpfs , elmnew , rays , FOSTAR, &
        !        zcenter(1), &
        !             RADX, H,cap_center(:,1),z0,disp,numberOfrays)
     end if

     call inhist(xfn,firstn,number,nextn)
  else
     write(*,*) 'cell l367 Restarting'
     call restart(lcube, nu, rho,td,ur,vr,wr, pr, &
          xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs, xcenterold, ycenterold, zcenterold, nsph, nellip, ellipa, ellipb, &
          ellipc, a_prestress, nvec_i)
     T=klok*time
  end if


  !     Initialize the fluid solver
  call influidu(vsc,qrfact, dsq,dx,dy,dz)


  !     MAIN LOOP --
  DO KLOK=KLOK1,KLOKEND
     message = 'cell l216'
     call dumpstatus(klok, message, 'status.txt')
     T = T+TD
     do i=1,ncap
        CALL MEMBNX(XFN(:,cap_n_start(i):cap_n_end(i)), &
             elmnew(:,cap_e_start(i):cap_e_end(i)), & 
             shpint(:,cap_e_start(i):cap_e_end(i)), &
             shpfs(:,cap_e_start(i):cap_e_end(i)), &
             FRC(:,cap_n_start(i):cap_n_end(i)),h,FOSTAR,&
             lambda1(cap_e_start(i):cap_e_end(i)), &
             lambda2(cap_e_start(i):cap_e_end(i)))
     end do

     message = 'cell l220'
     call dumpstatus(klok, message, 'status.txt')
     do i=1,nrects
        call pmhist(xpi,xfn(:,rect_n_start(i):rect_n_end(i)), &
             frc(:,rect_n_start(i):rect_n_end(i)),10240.d0/fnrectnodes)
     end do
     if (nrects > 0) call inhist(xfn, firstn, number, nextn)

     write(message, *) 'frc(1,1)',frc(1,1)
     call dumpstatus(klok, message, 'thumbprint.txt')
     write(message, *) 'frc(3,nfsize)',frc(1,nfsize)
     call dumpstatus(klok, message, 'thumbprint.txt')
     write(message, *) 'frc(2,nfsize)',frc(1,nfsize/2)
     call dumpstatus(klok, message, 'thumbprint.txt')

     message = 'cell l225'
     call dumpstatus(klok, message, 'status.txt')
     call meanforce(klok, frc)
     message = 'cell l228'
     call dumpstatus(klok, message, 'status.txt')

     call meanfluidvelocity(ur, meanu)
     call meanfluidvelocity(vr, meanv)
     call meanfluidvelocity(wr, meanw)
     open(403, file='meanfluidv.txt',  access='append')
     write(403,'(i6,3(x,es24.17))') klok, dreal(meanu), dreal(meanv), &
          dreal(meanw)
     close(403)

     message = 'cell l240'
     call dumpstatus(klok, message, 'status.txt')

     CALL pushup(UR,VR,WR,XFN,FRC,FIRSTN,NUMBER,NEXTN)
     message = 'cell l243'

     call dumpstatus(klok, message, 'status.txt')
     write(message, *) 'ur(3*ngx/4,ngy/2,ngz/4)', &
          ur(3*ngx/4,ngy/2,ngz/4)
     call dumpstatus(klok, message, 'thumbprint.txt')

     write(message, *) 'vr(3*ngx/4,ngy/2,ngz/4)', &
          vr(3*ngx/4,ngy/2,ngz/4)
     call dumpstatus(klok, message, 'thumbprint.txt')

     write(message, *) 'wr(3*ngx/4,ngy/2,ngz/4)', &
          wr(3*ngx/4,ngy/2,ngz/4)
     call dumpstatus(klok, message, 'thumbprint.txt')

     if (FVS == 0) then
        call bodyfs(bfs, ur, vr, wr, vsc)
     end if

     message = 'cell l248'
     call dumpstatus(klok, message, 'status.txt')

     if (fvs /= 0) then
        call fvssub(ur, vr, wr, bfs, umean)
        call poiseuille(wr, pr, -$dpdz$, vsc)
     end if

     !     call system_clock(sysclockstart)
     CALL FLUIDUP(KLOK,UR,VR,WR, pr, QRFACT, DSQ, DX, DY, DZ)
     if (fvs /= 0) then
        call fvssub(ur, vr, wr, -bfs, -umean)
        call poiseuille(wr, pr, $dpdz$, vsc)
     end if
     !     call system_clock(sysclockend, rate)
     !     print *, 'clock', klok, 'fluidup time', &
     !          float(sysclockend - sysclockstart)/rate
     !     print *, 'start', sysclockstart, 'end', sysclockend
     message = 'cell l252'
     call dumpstatus(klok, message, 'status.txt')

     call wrap(ur, vr, wr)

     write(message, *) 'ur(3*ngx/4,ngy/2,ngz/4)', &
          ur(3*ngx/4,ngy/2,ngz/4)
     call dumpstatus(klok, message, 'thumbprint.txt') 
     write(message, *) 'vr(3*ngx/4,ngy/2,ngz/4)', &
          vr(3*ngx/4,ngy/2,ngz/4)
     call dumpstatus(klok, message, 'thumbprint.txt')
     write(message, *) 'wr(3*ngx/4,ngy/2,ngz/4)', &
          wr(3*ngx/4,ngy/2,ngz/4)
     call dumpstatus(klok, message, 'thumbprint.txt')

     message = 'cell l255'
     call dumpstatus(klok, message, 'status.txt')
     CALL MOVE(UR,VR,WR,XFN,FIRSTN,NUMBER,NEXTN)
     if ($centerstable$ == 1) then
        xfn(1,:) = xfn(1,:) - sum(xfn(1,:))/size(xfn,2) + cap_center(1,1)
        xfn(2,:) = xfn(2,:) - sum(xfn(2,:))/size(xfn,2) + cap_center(2,1)
        xfn(3,:) = xfn(3,:) - sum(xfn(3,:))/size(xfn,2) + cap_center(3,1)
     end if
     message = 'cell l258'
     call dumpstatus(klok, message, 'status.txt')
     write(message, *) 'xfn(1,1)',xfn(1,1)
     call dumpstatus(klok, message, 'thumbprint.txt')
     write(message, *) 'xfn(3,nfsize)',xfn(3,nfsize)
     call dumpstatus(klok, message, 'thumbprint.txt')
     write(message, *) 'xfn(2,nfsize/2)',xfn(2,nfsize/2)
     call dumpstatus(klok, message, 'thumbprint.txt')

     !     End of the loop, do post-processing stuff
     do i=1,ncap
        call cellcenter(klok,xfn(1:3,cap_n_start(i):cap_n_end(i)), &
             nnode(i), i, xcenter(i), ycenter(i), zcenter(i))
        message = 'cell l232'
        call dumpstatus(klok, message, 'status.txt')
        call makefilename('capsulev__', i,'.txt',strfname)
        open(402,file=strfname, access='append')
        write(402,'(i6,3(x,es24.17))') klok, xcenter(i) - xcenterold(i), &
             ycenter(i) - ycenterold(i), zcenter(i) - zcenterold(i)
        close(402)
        xcenterold(i) = xcenter(i)
        ycenterold(i) = ycenter(i)
        zcenterold(i) = zcenter(i)
        ! Save the n vectors, the six vectors from the center of the capsule
        ! to the six tracer points.
        call dumpnvec(klok, xfn(:,cap_n_start(i):cap_n_end(i)), xcenter(i), &
             ycenter(i), zcenter(i), nvec_i(:, i), i)
        call shape(h64,klok, i, xfn(1:3, cap_n_start(i): &
             cap_n_end(i)), &
             nnode(i), elmnew(1:3, cap_e_start(i):cap_e_end(i)), &
             nelm(i))
        message = 'cell l262'
        call dumpstatus(klok, message, 'status.txt')
        call calculateDF(klok, i, xfn(1:3, cap_n_start(i): &
             cap_n_end(i)), nnode(i), &
             lambda1(cap_e_start(i):cap_e_end(i)), &
             lambda2(cap_e_start(i):cap_e_end(i)), &
             xcenter(i), ycenter(i), zcenter(i))

        call minmaxxfn(klok, xfn(:,cap_n_start(i):cap_n_end(i)), i)
        call volume_area(klok, xfn(:,cap_n_start(i):cap_n_end(i)), &
             elmnew(:, cap_e_start(i):cap_e_end(i)), i)

        call makefilename('lambda1___', i, '.txt', strfname)
        open(400, file=strfname, access='append')
        write(400,'(i6,x,es24.17,x,es24.17)') klok, minval(lambda1),&
             maxval(lambda1)
        close(400)

        call makefilename('lambda2___', i, '.txt', strfname)
        open(400, file=strfname, access='append')
        write(400,'(i6,x,es24.17,x,es24.17)') klok, minval(lambda2),&
             maxval(lambda2)
        close(400)

        if ((klok/$smalldumpint$)*$smalldumpint$==klok) then
           call profile(i, xfn(1:3,cap_n_start(i):cap_n_end(i)), &
                klok)
        end if
     end do

     open(206, file='checkinit.txt', access='append')
     WRITE(206,*)' KLOK: ',KLOK,  ' ; TIME: ',T
     close(206)

     message = 'cell l266'
     call dumpstatus(klok, message, 'status.txt')

     call wrstart(lcube, nu, rho,td,klok,ur,vr,wr, pr, &
          xfn,xpi,firstn,number,nextn,elmnew,shpint,shpfs, xcenterold, ycenterold, zcenterold, nsph, nellip, ellipa, ellipb, &
          ellipc, a_prestress, nvec_i)
     message = 'cell l270'

     call dumpstatus(klok, message, 'status.txt')
     if ((klok/$smalldumpint$)*$smalldumpint$ == klok) then
        message = 'cell l273'
        call dumpstatus(klok, message, 'status.txt')
        call wprofile(wr, klok)
        message = 'cell l276'
        call dumpstatus(klok, message, 'status.txt')

        message = 'cell l279'
        call dumpstatus(klok, message, 'status.txt')

        ! Files that are saved every smalldumpint, but only kept every bigdumpint
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
        call makefilename('stretches_',KLOK,'.txt',strfname)
        call dumplambdas(lambda1, lambda2, strfname)

        ! This doesn't kill the files that are zero, one, or two smalldumpints old
        ! or any files that have a klok count that is an even multiple of bigdumpint
        if ((((klok-2*$smalldumpint$)/$bigdumpint$)*$bigdumpint$) /= klok - &
             2*$smalldumpint$ .and. (klok-2*$smalldumpint$) > 0) then

           ! Delete old copies of these files
           ! Delete old velocity/pressure information
           call makefilename('uvwpdump__', klok-2*$smalldumpint$,'.txt',&
                strfname)
           open(300, file=strfname, form = 'unformatted')
           close(300, status='delete')

           ! Delete the total configuration of the solid nodes
           call makefilename('solidnodes',klok-2*$smalldumpint$,'.txt', &
                strfname)
           open(300, file=strfname, form = 'unformatted')
           close(300, status='delete')

           ! Delete the total configuration of the solid forces
           call makefilename('solidforce',klok-2*$smalldumpint$,'.txt', &
                strfname)
           open(300, file=strfname, form = 'unformatted')
           close(300, status='delete')

           ! Delete the file storing stretch ratios for each solid node
           call makefilename('stretches_',klok-2*$smalldumpint$,'.txt', &
                strfname)
           open(300, file=strfname, form = 'unformatted')
           close(300, status='delete')

        end if
     end if

     message = 'cell l299'
     call dumpstatus(klok, message, 'status.txt')

     if (minval(xfn(1,:)) < 0.5d0) then
        message = 'Nodes exited box in -x direction, terminating job.'
        call dumpstatus(klok, message, 'status.txt')
        stop
     end if
     if (maxval(xfn(1,:)) > flngx + 0.5d0) then
        message = 'Nodes exited box in +x direction, terminating job.'
        call dumpstatus(klok, message, 'status.txt')
        stop
     end if
     if (minval(xfn(2,:)) < 0.5d0) then
        message = 'Nodes exited box in -y direction, terminating job.'
        call dumpstatus(klok, message, 'status.txt')
        stop
     end if
     if (maxval(xfn(2,:)) > flngy + 0.5d0) then
        message = 'Nodes exited box in +y direction, terminating job.'
        call dumpstatus(klok, message, 'status.txt')
        stop
     end if
     if (minval(xfn(3,:)) < -0.5d0) then
        message = 'Nodes exited box in -z direction, terminating job.'
        call dumpstatus(klok, message, 'status.txt')
        stop
     end if
     if (maxval(xfn(3,:)) > flngz - 0.5d0) then
        message = 'Nodes exited box in +z direction, terminating job.'
        call dumpstatus(klok, message, 'status.txt')
        stop
     end if
  end do

  message = 'cell l304'
  call dumpstatus(klok, message, 'status.txt')

  deallocate(xpi)
  deallocate(lambda1,lambda2)
  DEALLOCATE (UR,VR,WR,FIRSTN,NUMBER,NEXTN,XFN,FRC)
  DEALLOCATE(QRFACT,elmnew,shpint,shpfs)
  if (nrects > 0) then
     deallocate(rect_nnode, rect_n1, rect_n2, recty,rect_n_start, &
          rect_n_end)
  end if
  if (ncap > 0) then 
     deallocate(a_prestress, nvec_i)
     deallocate(fineness, nnode, nelm)
     deallocate(cap_n_start, cap_n_end, cap_e_start, cap_e_end)
     deallocate(xcenter, ycenter, zcenter, &
          xcenterold, ycenterold, zcenterold)
     deallocate(ellipa, ellipb, ellipc)
     if (nsph > 0) deallocate(rad) 
  end if
END PROGRAM cell
!**********************************************************************
