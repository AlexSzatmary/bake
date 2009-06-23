!**********************************************************
      SUBROUTINE SHAPE(LCUBE,H64,KLOK,TD,NBEG,NENDS,NELM1,NELM2,XFN,
     &     elmnew)
      implicit none
      integer, parameter :: nfsize=$nsnode$,nfsize2=$nselm$
      integer, parameter :: m_start=$m_start$, m_end=$m_end$

      integer, parameter :: lxng=$lngx$,lyng=lxng,lzng=lxng
      integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
      integer, parameter :: fngx=ngx,fngy=ngy,fngz=ngz
      double precision :: LCUBE
      CHARACTER*15 filepro
      double precision :: XFN(1:3,1:NFSIZE)
      INTEGER elmnew(1:3,1:NFSIZE2)
      double precision :: h64, td, ra, cgx, cgy, cgz, xo, ro, 
     &     t, xa, ya, za
      double precision ::  xb, yb, zb, xc, yc, zc, check2, check3
      integer klok, nbeg, nends, nelm1, nelm2, i, icount, k, j1, j2, 
     &     j3
      double precision, allocatable :: xfp(:, :), zy(:, :)

      ALLOCATE (XFP(3,NFSIZE))
      ALLOCATE (zy(2,8000))

      IF(NBEG.EQ.1) THEN
         IF ((KLOK/400)*400.EQ.KLOK) THEN
            if(klok == 0) 
     &           write(filepro,305)'rbcshpk0000.pro'
 305        format(a15)
            if((KLOK.ge.1).and.(KLOK .lt. 10)) 
     &           write(filepro,301)'rbcshpk000',KLOK,'.pro'
 301        format(a10,i1,a4)
            if(KLOK .ge. 10 .and. KLOK .lt. 100) 
     &           write(filepro,302)'rbcshpk00',KLOK,'.pro'
 302        format(a9,i2,a4)
            if(KLOK.ge.100.and.KLOK.lt.1000) 
     &           write(filepro,303)'rbcshpk0',KLOK,'.pro'
 303        format(a8,i3,a4)
            if(KLOK.ge.1000.and.KLOK.lt.10000) 
     &           write(filepro,304)'rbcshpk',KLOK,'.pro'
 304        format(a7,i4,a4)
            open(25,file=filepro,status='unknown')
            DO i=m_start,m_end
               IF((XFN(1,i) >fngx/2.d0-1.d0).AND.(XFN(1,i) 
     &              <fngx/2.d0+1.d0)) THEN
                  write(25,300) XFN(2,i),XFN(3,i)
 300              format(e12.5,3x,e12.5)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF((KLOK/1)*1.EQ.KLOK) THEN
!     define GROUP and Point parameters
         ra = 1.0d0
         cgx = 0.0d0
         cgy = 0.0d0
         cgz = 0.0d0
         do i = NBEG,NENDS
            cgx = cgx + XFN(1,i)
            cgy = cgy + XFN(2,i)
            cgz = cgz + XFN(3,i)
         enddo
         cgx = cgx/dble(m_end-m_start+1)
         cgy = cgy/dble(m_end-m_start+1)
         cgz = cgz/dble(m_end-m_start+1)
         do i = NBEG,NENDS
            XFP(1,i) = XFN(1,i) - cgx
            XFP(2,i) = XFN(2,i) - cgy
            XFP(3,i) = XFN(3,i) - cgz
         enddo
         write(200,*) klok,cgx,cgy,cgz
         xo = 0.0d0
         ro =  ra/H64
         icount = 0
         t = dble(KLOK)*TD
         do k = NELM1,NELM2
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
 110        format(i5,x,i4,x,i4,x,6(f7.4,x))
 111        format(i5,x,3(e11.4,x))
         enddo
!     printout profile in the zy plane and value of Dxy
         call Dzy(zy,t,ro,icount,KLOK,NBEG)
      ENDIF
      DEALLOCATE (zy,XFP)
      return 
      end subroutine SHAPE
!********************************************************
      subroutine Dzy(zy,t,ro,icount,KLOK,NBEG)
      implicit none
      double precision :: zy(2,8000)
      INTEGER NBEG
      double precision :: rmax, rmin, pi, r, r2, ro, t, theta, ddzy
      integer icount, klok
      integer k
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

      write(202,*)klok,theta
      open(203, access='append')
      write(203,*)klok,ddzy
      close(203)
      return
      end subroutine Dzy
!**********************************************************************
      subroutine calculateDF(clock, xfn, nbeg, nend)
      implicit none
      integer nbeg, nend
      integer, parameter :: nfsize = $nsnode$
      integer, parameter :: nmem = $nmem$
      double precision :: xfn(3, nfsize)
      integer i
      integer clock
      double precision :: rmin, rmax, r
      double precision :: cgx, cgy, cgz
      double precision :: DF

      rmin = 1000d0
      rmax = -1000d0

!      if (clock == 1) then
!         write(*,*) 'nbeg', nbeg, 'nbeg', nend, 'nmem', nmem
!         write(*,*) 'rmin', rmin, 'rmax', rmax
!      end if

      cgx = 0.0d0
      cgy = 0.0d0
      cgz = 0.0d0
      do i = nbeg,nend
         cgx = cgx + xfn(1,i)
         cgy = cgy + xfn(2,i)
         cgz = cgz + xfn(3,i)
      enddo
      cgx = cgx/dble(nmem)
      cgy = cgy/dble(nmem)
      cgz = cgz/dble(nmem)
C       if (clock == 1) then
C          write(*,*) 'cgx', cgx, 'cgy', cgy, 'cgz', cgz
C       end if

      do i = nbeg, nend
         r = dsqrt((xfn(1,i)-cgx)**2 + (xfn(2,i)-cgy)**2 + 
     &        (xfn(3,i)-cgz)**2)
C          if (clock == 1) then
C             write(*,*) 'r', r
C          end if
         if(r > rmax) rmax = r
         if(r < rmin) rmin = r
C          if (clock == 1) then
C             write(*,*) 'rmin', rmin, 'rmax', rmax
C             write(*,*) 'i', i, 'x', xfn(1,i), 'y', xfn(2,i), 
C      &           'z', xfn(3,i)
C          end if
      end do

      DF = (rmax - rmin)/(rmax+rmin)

      open(204, access='append')
      write(204,*) clock, ',', DF
      close(204)

      return
      end subroutine calculateDF
!**********************************************************************
      subroutine saveallsolid(xfn,strfname)
!     saveallsolid
!     Alex Szatmary
!     8-24-06
!     Records the coordinates of each solid node
!     !? Is this the actual output format?
!     Format: x, y, z

      implicit none
      character*19 strfname
      integer, parameter :: nfsize=$nsnode$
      double precision :: xfn(1:3,1:nfsize)
      integer i
      open(25,file=strfname,status='unknown')
      do i = 1,nfsize
         write(25,*) xfn(1,i),xfn(2,i),xfn(3,i)
      end do
      close(25)
      end subroutine saveallsolid
!**********************************************************************
      subroutine makefilename(str1, clock, str2, strfname)
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
      integer, parameter :: lngx=$lngx$,lngy=lngx,lngz=lngx
      integer, parameter :: ngx=2**lngx,ngy=2**lngy,ngz=2**lngz
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
               write(25,*) i,j,k, dble(wr(i,j,k)), dimag(wr(i,j,k))
!            end do
         end do
!      end do
      close(25)

      end subroutine wprofile
!**********************************************************************
      subroutine wdump(wr,clock)
      implicit none
      integer, parameter :: lngx=$lngx$,lngy=lngx,lngz=lngx
      integer, parameter :: ngx=2**lngx,ngy=2**lngy,ngz=2**lngz
      integer i, j, k, clock
      double complex :: wr(0:ngx+2,0:ngy+2,0:ngz-1)
      character*11 wdumpname

      if((clock >= 0).and.(clock .lt. 10)) then
         write(wdumpname,601)'wdu000',clock,'.txt' 
      end if
 601  format(a6,i1,a4)
      if(clock  >=  10 .and. clock .lt. 100) then
         write(wdumpname,602)'wdu00',clock,'.txt'
      end if
 602  format(a5,i2,a4)
      if(clock >= 100.and.clock.lt.1000) then
         write(wdumpname,603)'wdu0',clock,'.txt'
      end if
 603  format(a4,i3,a4)
      if(clock >= 1000.and.clock.lt.10000) then
         write(wdumpname,604)'wdu',clock,'.txt'
      end if
 604  format(a3,i4,a4)
      
      open(25,file=wdumpname,status='unknown')
      do i =1,ngx
         do j = 1,ngy
            do k = 0,ngz-1
               write(25,*) i,j,k, dreal(wr(i,j,k)), dimag(wr(i,j,k))
            end do
         end do
      end do
      close(25)

      end subroutine wdump
!**********************************************************************
      subroutine uvwdump(ur, vr, wr,clock)
      implicit none
      integer, parameter :: lngx=$lngx$,lngy=lngx,lngz=lngx
      integer, parameter :: ngx=2**lngx,ngy=2**lngy,ngz=2**lngz
      integer i, j, k, clock
      double complex :: ur(0:ngx+2,0:ngy+2,0:ngz-1)
      double complex :: vr(0:ngx+2,0:ngy+2,0:ngz-1)
      double complex :: wr(0:ngx+2,0:ngy+2,0:ngz-1)
      character*19 uvwdumpname

      call makefilename('uvwdump___', clock,'.txt', uvwdumpname)
      
      open(25,file=uvwdumpname,status='unknown')
      do i =1,ngx
         do j = 1,ngy
            do k = 0,ngz-1
               write(25,*) i,j,k, dreal(ur(i,j,k)), dreal(vr(i,j,k)), 
     &              dreal(wr(i,j,k))
            end do
         end do
      end do
      close(25)

      end subroutine uvwdump
!**********************************************************************
      subroutine wlindump(wlin,clock)
      implicit none
      integer, parameter :: lngx=$lngx$,lngy=lngx,lngz=lngx
      integer, parameter :: ngx=2**lngx,ngy=2**lngy,ngz=2**lngz
      integer i, clock
      double precision :: wlin(-15:16*(ngz+2))
      character*11 wlindumpname

      if((clock >= 0).and.(clock .lt. 10)) then
         write(wlindumpname,601)'wli000',clock,'.txt' 
      end if
 601  format(a6,i1,a4)
      if(clock  >=  10 .and. clock .lt. 100) then
         write(wlindumpname,602)'wli00',clock,'.txt'
      end if
 602  format(a5,i2,a4)
      if(clock >= 100.and.clock.lt.1000) then
         write(wlindumpname,603)'wli0',clock,'.txt'
      end if
 603  format(a4,i3,a4)
      if(clock >= 1000.and.clock.lt.10000) then
         write(wlindumpname,604)'wli',clock,'.txt'
      end if
 604  format(a3,i4,a4)
      
      open(25,file=wlindumpname,status='unknown')
      do i =-15,16*(ngz+2)
         write(25,*) i, wlin(i)
      end do
      close(25)

      end subroutine wlindump
!**********************************************************
      subroutine meanforce(clock, frc)
      implicit none
      integer, parameter :: nfsize=$nsnode$
!     !! Watch out; this gives mean force including the plane.
      integer i, clock
      double precision :: frc(3,nfsize)
      double precision :: forcex, forcey, forcez
      forcex = 0.d0
      forcey = 0.d0
      forcez = 0.d0
      do i=1,nfsize
         forcex=forcex + frc(1,i)
         forcey=forcey + frc(2,i)
         forcez=forcez + frc(3,i)
      end do

      open(451, access='append')
      write(451,*) clock, forcex, forcey, forcez
      close(451)
      end subroutine meanforce
!**********************************************************
      subroutine cellcenter(klok, xfn, xcenter, ycenter, zcenter)
      implicit none
      integer, parameter :: nfsize=$nsnode$
      integer, parameter :: m_start=$m_start$, m_end=$m_end$
      double precision :: xcenter, ycenter, zcenter
      double precision :: xfn(3, nfsize)
      integer :: nmem=$nmem$
      integer i, klok

      do i = m_start, m_end
         xcenter = xcenter + XFN(1,i)
         ycenter = ycenter + XFN(2,i)
         zcenter = zcenter + XFN(3,i)
      end do

      xcenter = xcenter/dble(nmem)
      ycenter = ycenter/dble(nmem)
      zcenter = zcenter/dble(nmem)

      open(401, access='append')
      write(401,*) klok, xcenter,ycenter,zcenter
      close(401)
      return
      end subroutine cellcenter
!**********************************************************
      subroutine dumpstatus(clock, message)
      integer clock
      character*11 fname
      character*80 message
      fname='status.txt'
      open(500,file=fname,status='unknown', access='append')
      write(500, *) clock
      write(500, *) message
      close(500)
      end subroutine dumpstatus
!**********************************************************
      subroutine meanfluidvelocity(ur, meanstablev)
      integer, parameter :: lngx=$lngx$,lngy=lngx,lngz=lngx
      integer, parameter :: ngx=2**lngx,ngy=2**lngy,ngz=2**lngz
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
