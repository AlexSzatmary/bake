!**********************************************************************
! Makes a profile out of the capsule by looking
! to see which nodes are within one fluid node of the assumed
! centerline.
!m This subroutine is obsolescent. Instead, a subroutine
! that plots zy in shape should be used.
subroutine profile(cap_i, xfn, clock)
  implicit none
  integer, parameter :: lngz = $lngz$
  integer, parameter :: ngz = 2**lngz - 1
  double precision, parameter :: fngz = ngz
  character*21 strfname
  integer cap_i, clock, i
  double precision :: xfn(:, :)
  integer nnode

  nnode = size(xfn,2)

  write(strfname,'(A6,I4,A1,I6,A4)') 'cappro', cap_i, '_', clock, &
       '.txt'
  call padzeros(strfname)
  open(25, file=strfname, access='append')
  DO i=1, nnode
     IF((XFN(3,i) > 0.5d0*(fngz-1.d0) - 0.5d0) .AND.(XFN(3,i) &
          < 0.5d0*(fngz-1.d0) + 0.5d0)) THEN
        write(25,'(es24.17,x,es24.17)') XFN(1,i),XFN(2,i)
     ENDIF
  ENDDO
  close(25)
end subroutine profile
!**********************************************************
SUBROUTINE SHAPE(KLOK, cap_i, xfn,nnode,elmnew,nelm)
  implicit none
  integer, parameter :: lxng=$lngx$,lyng=$lngy$,lzng=$lngz$
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  integer, parameter :: fngx=ngx,fngy=ngy,fngz=ngz
  double precision :: XFN(:,:)
  INTEGER elmnew(:,:)
  double precision :: ra, cgx, cgy, cgz, xo, xa, ya, za
  double precision ::  xb, yb, zb, xc, yc, zc, check2, check3
  integer :: klok, nnode, nelm, i, icount, k, j1, j2, j3
  double precision, allocatable :: xfp(:, :), zy(:, :)
  integer cap_i

  ALLOCATE (XFP(3,nnode))
  ALLOCATE (zy(2,8000))

  !     define GROUP and Point parameters
  ! Make xfp. If xfn represents a set of nodes, xfp is those nodes, but
  ! translated so that their center is at the origin.
  ra = 1.0d0
  cgx = 0.0d0
  cgy = 0.0d0
  cgz = 0.0d0
  do i = 1, nnode
     cgx = cgx + XFN(1,i)
     cgy = cgy + XFN(2,i)
     cgz = cgz + XFN(3,i)
  enddo
  cgx = cgx/dble(nnode)
  cgy = cgy/dble(nnode)
  cgz = cgz/dble(nnode)
  do i = 1, nnode
     XFP(1,i) = XFN(1,i) - cgx
     XFP(2,i) = XFN(2,i) - cgy
     XFP(3,i) = XFN(3,i) - cgz
  enddo

  xo = 0.0d0
  icount = 0

  ! This loop picks out each element and looks at each edge; if the edge
  ! crosses the y-z plane (cutting through the center of the nodes) then
  ! the point on the edge that has x=0 has its coordinates recorded.
  ! This way, rather than getting just a ring of nodes near the center, the
  ! precise nodes that pass through x=0 are selected.
  do k = 1, nelm
     j1 = elmnew(1,k)
     j2 = elmnew(2,k)
     j3 = elmnew(3,k)
     xa = XFP(1,j1)
     ya = XFP(2,j1)
     za = XFP(3,j1)
     xb = XFP(1,j2)
     yb = XFP(2,j2)
     zb = XFP(3,j2)
     xc = XFP(1,j3)
     yc = XFP(2,j3)
     zc = XFP(3,j3)
     if((xa-xo)*(xb-xo) < 0.d0)then
        icount = icount + 1
        zy(1,icount) = za + (zb-za)*(xo-xa)/(xb-xa)
        zy(2,icount) = ya + (yb-ya)*(xo-xa)/(xb-xa)
     endif
     if((xa-xo) .eq. 0.d0) then
        icount = icount + 1
        zy(1,icount) = za
        zy(2,icount) = ya
     endif
     if((xb-xo) .eq. 0.d0) then
        icount = icount + 1
        zy(1,icount) = zb
        zy(2,icount) = yb
     endif
     check2= (xa-xo)*(xc-xo)
     if(check2.lt. dabs(check2)) then
        icount = icount + 1
        zy(1,icount) = za + (zc-za)*(xo-xa)/(xc-xa)
        zy(2,icount) = ya + (yc-ya)*(xo-xa)/(xc-xa)
     endif
     if((xa-xo) .eq. 0.d0) then
        icount = icount + 1
        zy(1,icount) = za
        zy(2,icount) = ya
     endif
     if((xc-xo) .eq. 0.d0) then
        icount = icount + 1
        zy(1,icount) = zc
        zy(2,icount) = yc
     endif
     check3= (xb-xo)*(xc-xo)
     if(check3 .lt. dabs(check3)) then
        icount = icount + 1
        zy(1,icount) = zb + (zc-zb)*(xo-xb)/(xc-xb)
        zy(2,icount) = yb + (yc-yb)*(xo-xb)/(xc-xb)
     endif
     if((xb-xo) .eq. 0.d0) then
        icount = icount + 1
        zy(1,icount) = zb
        zy(2,icount) = yb
     endif
     if((xc-xo) .eq. 0.d0) then
        icount = icount + 1
        zy(1,icount) = zc
        zy(2,icount) = yc
     endif
  enddo
  !     printout theta and DF based on dzy
  ! The points on x=0 are passed to dzy, which does some output on the DF
  ! seen in the y-z plane, and the angle of inclination from the z-axis
  call dzy(zy,icount,klok, cap_i)
  DEALLOCATE (zy,XFP)
  return 
end subroutine SHAPE
!********************************************************
subroutine dzy(zy,icount,klok, cap_i)
  implicit none
  ! This subroutine outputs the DF
  ! seen in the y-z plane, and the angle of inclination from the z-axis
  ! Points in dzy are from shape, which picks nodes in the y-z plane.
  double precision :: zy(2,8000)
  double precision :: rmax, rmin, pi, r, r2, theta, ddzy
  integer icount, klok
  integer k, cap_i
  character*15 strfname

  !     !m Change these numbers so that they're a lot bigger
  rmax = -100.d0
  rmin =  100.d0
  pi = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits

  do k = 1,icount
     r2 = zy(1,k)**2 + zy(2,k)**2
     r = dsqrt(r2)
     if(r .gt. rmax) then
        rmax = r
        theta = datan((zy(2,k)/zy(1,k)))/pi
     endif
     if(r .lt. rmin) rmin = r
  end do

  ddzy = (rmax - rmin)/(rmax+rmin)

  write(strfname,'(A7,I4,A4)') 'theta__', cap_i, '.txt'
  call padzeros(strfname)
  open(202, file=strfname, access='append')
  write(202,'(i6,es24.17)')klok,theta
  close(202)

  write(strfname,'(A7,I4,A4)') 'planedf', cap_i, '.txt'
  call padzeros(strfname)
  open(203, file=strfname, access='append')
  write(203,'(i6,es24.17)')klok, ddzy
  close(203)
  return
end subroutine Dzy
!**********************************************************************
subroutine calculateDF(clock, cap_i, xfn, my_nnode, lambda1, lambda2, &
     xcenter, ycenter, zcenter)
  implicit none
  ! This subroutine calculates the DF based on the points nearest to and
  ! farthest from the centroid.
  ! It also does other interesting things, like giving the stretch ratios
  ! at those points, and recording the directions to those points as mvec
  integer my_nnode
  double precision :: xfn(:, :)
  integer cap_i, i
  integer clock
  double precision :: rmin, rmax, r
  double precision :: cgx, cgy, cgz
  double precision :: DF
  character*19 strfname
  integer imax, imin
  double precision :: lambda1(:), lambda2(:)
  double precision :: xcenter, ycenter, zcenter

  rmin = 1000d0
  rmax = -1000d0

  cgx = 0.0d0
  cgy = 0.0d0
  cgz = 0.0d0
  do i = 1, my_nnode
     cgx = cgx + xfn(1,i)
     cgy = cgy + xfn(2,i)
     cgz = cgz + xfn(3,i)
  enddo
  cgx = cgx/dble(my_nnode)
  cgy = cgy/dble(my_nnode)
  cgz = cgz/dble(my_nnode)

  do i = 1, my_nnode
     r = dsqrt((xfn(1,i)-cgx)**2 + (xfn(2,i)-cgy)**2 + &
          (xfn(3,i)-cgz)**2)
     if(r > rmax) then
        rmax = r
        imax = i
     end if
     if(r < rmin) then 
        rmin = r
        imin = i
     end if
  end do

  DF = (rmax - rmin)/(rmax+rmin)

  call makefilename('TaylorDF__', cap_i,'.txt',strfname)
  open(204, file=strfname, access='append')
  write(204,'(i6,3(x,es24.17))') clock, DF, rmin, rmax
  close(204)

  call makefilename('lambdarmax', cap_i,'.txt',strfname)
  open(204, file=strfname, access='append')
  write(204,'(i6,3(x,es24.17))') clock, lambda1(imax), lambda2(imax)
  close(204)

  call makefilename('lambdarmin', cap_i,'.txt',strfname)
  open(204, file=strfname, access='append')
  write(204,'(i6,3(x,es24.17))') clock, lambda1(imin), lambda2(imin)
  close(204)

  call makefilename('mvecrmax__', cap_i,'.txt',strfname)
  open(204, file=strfname, access='append')
  write(204,'(i6,3(x,es24.17))') clock, xfn(1, imax) - xcenter, &
          xfn(2, imax) - ycenter, xfn(3, imax) - zcenter
  close(204)

  call makefilename('mvecrmin__', cap_i,'.txt',strfname)
  open(204, file=strfname, access='append')
  write(204,'(i6,3(x,es24.17))') clock, xfn(1, imin) - xcenter, &
          xfn(2, imin) - ycenter, xfn(3, imin) - zcenter
  close(204)
  close(204)

  return
end subroutine calculateDF
!**********************************************************************
subroutine saveallsolid(xfn,strfname)
  !     saveallsolid
  !     Alex Szatmary
  !     8-24-06
  !     Records the coordinates of each solid node
  ! It's since been stretched to output the force on each solid node.
  !     Format: x, y, z

  implicit none
  character(len=*) strfname
  double precision :: xfn(:,:)
  integer i
  open(25,file=strfname,status='unknown')
  do i = 1,size(xfn,2)
     write(25,'(es24.17,2(x,es24.17))') xfn(1,i),xfn(2,i),xfn(3,i)
  end do
  close(25)
end subroutine saveallsolid
!**********************************************************************
subroutine makefilename(str1, clock, str2, strfname)
  ! This subroutine should be deemed obsolescent; idioms seen in, say,
  ! dumpnvec, using padzeros, are simpler.
  !     makefilename
  !     Alex Szatmary
  !     8-24-06
  !     Makes an appropriate filename for a file.
  !     str1 should be a description of the file (i.e. 'sideview__')
  !     It must be 10 characters long (use underscores or similar for padding)
  !     clock is the current timestep
  !     str2 is the extension, which must be 4 characters long. I recommend '.txt'
  !     strfname is the filename. Pass this into a function that uses a filename
  !     as a parameter.
  !     Warnings:
  !     * The character limits are hard; this function will have to be modified
  !     for longer file names, shorter file names should use underscores for
  !     padding.
  !     * clock values greater than 99999 are not permitted
  !     !m
  !     No error checking is done on the value of clock.

  implicit none
  character*10 str1
  character*4 str2
  character*19 strfname
  integer clock

  select case (clock)
  case(0:9)
     write(strfname,'(A10,A4,I1,A4)') str1,'0000',clock,str2
  case(10:99)
     write(strfname,'(A10,A3,I2,A4)') str1,'000',clock,str2
  case(100:999)
     write(strfname,'(A10,A2,I3,A4)') str1,'00',clock,str2
  case(1000:9999)
     write(strfname,'(A10,A1,I4,A4)') str1,'0',clock,str2
  case(10000:99999)
     write(strfname,'(A10,I5,A4)') str1,clock,str2
  end select

end subroutine makefilename
!**********************************************************************
subroutine wprofile(wr,clock)
  implicit none
  ! Gives the w velocity profile along i=ngx/2, y = ngz/2-1
  integer, parameter :: lxng=$lngx$,lyng=$lngy$,lzng=$lngz$
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  integer i, j, k, clock
  double complex :: wr(0:ngx+2,0:ngy+2,0:ngz-1)
  character*19 wprofname

  call makefilename('wprofile__', clock, '.txt', wprofname)

  open(25,file=wprofname,status='unknown')
  i = ngx/2
  !      do i = 1, ngx
  !      j = ngy/2
  do j = 1, ngy
     k = ngz/2-1
     !            do k = 0, ngz-1
     write(25,'(3(i5x),x,e24.17,x,e24.17)') i,j,k, dble(wr(i,j,k)), dimag(wr(i,j,k))
     !            end do
  end do
  !      end do
  close(25)

end subroutine wprofile
!**********************************************************************
subroutine wdump(wr,clock)
  implicit none
  ! Saves the w velocity at all nodes; deprecated and unused.
  integer, parameter :: lxng=$lngx$,lyng=$lngy$,lzng=$lngz$
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  integer i, j, k, clock
  double complex :: wr(0:ngx+2,0:ngy+2,0:ngz-1)
  character*11 wdumpname

  if((clock >= 0).and.(clock .lt. 10)) then
     write(wdumpname,601)'wdu000',clock,'.txt' 
  end if
601 format(a6,i1,a4)
  if(clock  >=  10 .and. clock .lt. 100) then
     write(wdumpname,602)'wdu00',clock,'.txt'
  end if
602 format(a5,i2,a4)
  if(clock >= 100.and.clock.lt.1000) then
     write(wdumpname,603)'wdu0',clock,'.txt'
  end if
603 format(a4,i3,a4)
  if(clock >= 1000.and.clock.lt.10000) then
     write(wdumpname,604)'wdu',clock,'.txt'
  end if
604 format(a3,i4,a4)

  open(25,file=wdumpname,status='unknown')
  do i =1,ngx
     do j = 1,ngy
        do k = 0,ngz-1
           write(25,'(3(i4),2(es24.17))') i,j,k, dreal(wr(i,j,k)), &
                dimag(wr(i,j,k))
        end do
     end do
  end do
  close(25)

end subroutine wdump
!**********************************************************************
subroutine uvwpdump(ur, vr, wr, pr, clock)
  implicit none
  ! Saves the velocity and pressure at every fluid node
  integer, parameter :: lxng=$lngx$,lyng=$lngy$,lzng=$lngz$
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  integer i, j, k, clock
  double complex :: ur(0:ngx+2,0:ngy+2,0:ngz-1)
  double complex :: vr(0:ngx+2,0:ngy+2,0:ngz-1)
  double complex :: wr(0:ngx+2,0:ngy+2,0:ngz-1)
  double complex :: pr(0:ngx+2,0:ngy+2,0:ngz-1)

  character*19 uvwpdumpname

  call makefilename('uvwpdump__', clock,'.txt', uvwpdumpname)

  open(25,file=uvwpdumpname,status='unknown')
  do i =1,ngx
     do j = 1,ngy
        do k = 0,ngz-1
           write(25,'(3(i4),4(x,es24.17))') i,j,k, dreal(ur(i,j,k)), &
                dreal(vr(i,j,k)), &
                dreal(wr(i,j,k)), dreal(pr(i,j,k))
        end do
     end do
  end do
  close(25)

end subroutine uvwpdump
!**********************************************************************
subroutine wlindump(wlin,clock)
  implicit none
  ! Saves wlin velocities, deprecated and unused
  integer, parameter :: lxng=$lngx$,lyng=$lngy$,lzng=$lngz$
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  integer i, clock
  double precision :: wlin(-15:16*(ngz+2))
  character*11 wlindumpname

  if((clock >= 0).and.(clock .lt. 10)) then
     write(wlindumpname,601)'wli000',clock,'.txt' 
  end if
601 format(a6,i1,a4)
  if(clock  >=  10 .and. clock .lt. 100) then
     write(wlindumpname,602)'wli00',clock,'.txt'
  end if
602 format(a5,i2,a4)
  if(clock >= 100.and.clock.lt.1000) then
     write(wlindumpname,603)'wli0',clock,'.txt'
  end if
603 format(a4,i3,a4)
  if(clock >= 1000.and.clock.lt.10000) then
     write(wlindumpname,604)'wli',clock,'.txt'
  end if
604 format(a3,i4,a4)

  open(25,file=wlindumpname,status='unknown')
  do i =-15,16*(ngz+2)
     write(25,'(i6,es24.17)') i, wlin(i)
  end do
  close(25)

end subroutine wlindump
!**********************************************************
subroutine meanforce(clock, frc)
  implicit none
  !     !! Watch out; this gives mean force on all nodes, including in walls
  integer i, clock
  double precision :: frc(:,:)
  double precision :: forcex, forcey, forcez
  forcex = 0.d0
  forcey = 0.d0
  forcez = 0.d0
  do i=1,size(frc,2)
     forcex=forcex + frc(1,i)
     forcey=forcey + frc(2,i)
     forcez=forcez + frc(3,i)
  end do
  forcex = forcex/size(frc,2)
  forcey = forcex/size(frc,2)
  forcez = forcex/size(frc,2)

  open(451, file='meanforce.txt', access='append')
  write(451,'(i6,3(xes24.17))') clock, forcex, forcey, forcez
  close(451)
end subroutine meanforce
!**********************************************************
subroutine cellcenter(klok, xfn, my_nnode, cap_i, xcenter, &
     ycenter, zcenter)
  implicit none
  ! Calculates the center of a given capsule and saves these coordinates
  integer cap_i
  double precision :: xcenter, ycenter, zcenter
  double precision :: xfn(:,:)
  integer i, klok, my_nnode
  character*19 strfname

  xcenter = 0.d0
  ycenter = 0.d0
  zcenter = 0.d0

  do i = 1, my_nnode
     xcenter = xcenter + XFN(1,i)
     ycenter = ycenter + XFN(2,i)
     zcenter = zcenter + XFN(3,i)
  end do

  xcenter = xcenter/dble(my_nnode)
  ycenter = ycenter/dble(my_nnode)
  zcenter = zcenter/dble(my_nnode)

  call makefilename('capsulex__', cap_i,'.txt',strfname)
  open(401, file=strfname, access='append')
  write(401,'(i6,3(x,es24.17))') klok, xcenter,ycenter,zcenter
  close(401)
  return
end subroutine cellcenter
!**********************************************************
subroutine dumpstatus(clock, message, fname)
  integer clock
  ! Quick and dirty way to dump a string to disk, normally used for debugging.
  character(len=*) fname
  character(len=*) message
  open(500,file=fname,status='unknown', access='append')
  write(500, *) clock
  write(500, *) message
  close(500)
end subroutine dumpstatus
!**********************************************************
subroutine meanfluidvelocity(ur, meanstablev)
  ! Calculates the mean fluid velocity; this is applied through separate calls
  ! on u, v, and w. 
  integer, parameter :: lxng=$lngx$,lyng=$lngy$,lzng=$lngz$
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  double complex :: ur(0:ngx+2,0:ngy+2,0:ngz-1)
  integer i, j, k
  double complex :: meanstablev
  meanstablev = 0.
  do i = 1,ngx
     do j = 1,ngy
        do k=0,ngz-1
           meanstablev = meanstablev + ur(i,j,k)
        end do
     end do
  end do
  meanstablev = meanstablev/(ngx*ngy*ngz)
  return
end subroutine meanfluidvelocity
!**********************************************************************
subroutine padzeros(str)
  implicit none
  ! Takes a string and replaces any space characters in it with zeros.
  integer i
  character(len=*) str
  do i=1,len(str)
     if (str(i:i) == ' ') str(i:i) = '0'
  end do
end subroutine padzeros
!**********************************************************************
subroutine minmaxxfn(clock, xfn, cap_i)
  implicit none
  ! Outputs the coordinates for a rectanguloid that would bound a given 
  ! capsule 
  double precision :: xfn(:,:)
  integer :: cap_i
  character*19 strfname
  integer clock

  call makefilename('minmaxx___', cap_i,'.txt',strfname)
  open(400, file=strfname, access='append')
  write(400,'(i6,x,es24.17,x,es24.17)') clock, minval(xfn(1,:)),&
       maxval(xfn(1,:))
  close(400)
  
  call makefilename('minmaxy___', cap_i,'.txt',strfname)
  open(400, file=strfname, access='append')
  write(400,'(i6,x,es24.17,x,es24.17)') clock, minval(xfn(2,:)),&
       maxval(xfn(2,:))
  close(400)
  
  call makefilename('minmaxz___', cap_i,'.txt',strfname)
  open(400, file=strfname, access='append')
  write(400,'(i6,x,es24.17,x,es24.17)') clock, minval(xfn(3,:)),&
       maxval(xfn(3,:))
  close(400)
  
end subroutine minmaxxfn
!**********************************************************************
SUBROUTINE volume_area(KLOK,XFN,elmnew, cap_i)
  IMPLICIT NONE
  ! Calculates and outputs the volume and surface area of a given capsule.
  INTEGER NFSIZE,NFSIZE2
  INTEGER i,KLOK, cap_i
  integer j1,j2,j3
  double precision :: x1,y1,z1,x2,y2,z2,x3,y3,z3,xc1,yc1,zc1
  double precision :: a1,a2,a3,b1,b2,b3,c1,c2,c3,vol,area
  double precision :: XFN(:,:)
  INTEGER ELMNEW(:,:)
  character*19 strfname

  nfsize = size(xfn, 2)
  nfsize2 = size(elmnew, 2)

  area=0.d0
  vol=0.d0
  xc1 = 0.0d0
  yc1 = 0.0d0
  zc1 = 0.0d0
  do i = 1, size(xfn,2)
     xc1 = xc1 + XFN(1,i)
     yc1 = yc1 + XFN(2,i)
     zc1 = zc1 + XFN(3,i)
  end do
  xc1 = xc1/dble(nfsize)
  yc1 = yc1/dble(nfsize)
  zc1 = zc1/dble(nfsize)

!!$OMP  PARALLEL DO &
!!$OMP& SHARED(XFN,ELMNEW,xc1,yc1,zc1,vol,vol2,area) &
!!$OMP& PRIVATE(a1,a2,a3,b1,b2,b3,c1,c2,c3,i,j1,j2,j3, &
!!$OMP& X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3)
  do i = 1, nfsize2
     j1 = ELMNEW(1,i)
     j2 = ELMNEW(2,i)
     j3 = ELMNEW(3,i)
     X1 = XFN(1,j1)
     Y1 = XFN(2,j1)
     Z1 = XFN(3,j1)
     X2 = XFN(1,j2)
     Y2 = XFN(2,j2)
     Z2 = XFN(3,j2)
     X3 = XFN(1,j3)
     Y3 = XFN(2,j3)
     Z3 = XFN(3,j3)
     a1 = X2 - X1
     a1 = X2 - X1
     a2 = Y2 - Y1
     a3 = Z2 - Z1
     b1 = X3 - X1
     b2 = Y3 - Y1
     b3 = Z3 - Z1
     c1 = a2*b3-b2*a3
     c2 = b1*a3-a1*b3
     c3 = a1*b2-b1*a2
     vol=vol+(1.d0/6.d0)*dABS(c1*(XC1-X1)+c2*(YC1-Y1)+c3*(ZC1-Z1))
     area=area+0.5d0*dabs(dsqrt(c1**2+c2**2+c3**2))
  end do
!!$OMP END PARALLEL DO

  call makefilename('volumearea', cap_i,'.txt',strfname)
  open(400, file=strfname, access='append')
  write(400,'(i6,x,es24.17,x,es24.17)') klok, vol, area
  close(400)

  return
end subroutine volume_area
!**********************************************************************
subroutine dumplambdas(lambda1, lambda2, strfname)
  !     Records the stretch ratio, lambda, on each element.
  implicit none
  character(len=*) strfname
  double precision :: lambda1(:), lambda2(:)
  integer i
  open(25,file=strfname,status='unknown')
  do i = 1,size(lambda1)
     write(25,'(i8,2(x,es24.17))') i, lambda1(i), lambda2(i)
  end do
  close(25)
end subroutine dumplambdas
!**********************************************************************
subroutine dumpnvec(clock, xfn, xcenter, ycenter, zcenter, my_nvec_i, cap_i)
  implicit none
  ! Saves the orientation, the n vectors, of the capsule
  integer clock
  double precision :: xfn(:,:), xcenter, ycenter, zcenter
  integer my_nvec_i(:)
  integer j
  character*14 strfname
  integer cap_i
  do j = 1, 6
     write(strfname,'(A4,I1,A1,I4,A4)') 'nvec', j, '_', cap_i, '.txt'
     call padzeros(strfname)
     open(401, file=strfname, access='append')
     write(401,'(i6,3(x,es24.17))') clock, xfn(1, my_nvec_i(j)) - xcenter, &
          xfn(2, my_nvec_i(j)) - ycenter, &
          xfn(3, my_nvec_i(j)) - zcenter
     close(401)
  end do
  return
end subroutine dumpnvec
