      subroutine fvssub(ur, vr, wr, bfs, umean)
      integer, parameter :: lxng=$lngx$,lyng=lxng,lzng=lxng
      integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
      double precision, parameter :: flngx=ngx, flngy=ngy, flngz=ngz
      double complex :: ur(0:ngx+2, 0:ngy+2,0:ngz-1), 
     &     vr(0:ngx+2, 0:ngy+2, 0:ngz-1), wr(0:ngx+2, 0:ngy+2, 0:ngz-1)
      double precision :: bfs(3,3), umean(3)
      integer i, j, k

      do i = 1,ngx
         do j = 1,ngy
            do k = 0,ngz-1
               ur(i,j,k) = ur(i,j,k) - (bfs(1,1)*(i - (flngx+1)/2) + 
     &              bfs(1,2)*(j - (flngy+1)/2) + 
     &              bfs(1,3)*(k - (flngz-1)/2) + umean(1))
               vr(i,j,k) = vr(i,j,k) - (bfs(2,1)*(i - (flngx+1)/2) + 
     &              bfs(2,2)*(j - (flngy+1)/2) + 
     &              bfs(2,3)*(k - (flngz-1)/2) + umean(2))
               wr(i,j,k) = wr(i,j,k) - (bfs(3,1)*(i - (flngx+1)/2) + 
     &              bfs(3,2)*(j - (flngy+1)/2) + 
     &              bfs(3,3)*(k - (flngz-1)/2) + umean(3))
            end do
         end do
      end do

      return
      end subroutine fvssub
