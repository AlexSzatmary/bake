program rbc3d
  implicit none

  integer, parameter :: lxng=6,lyng=lxng,lzng=lxng+1
  integer, parameter :: ngx=2**lxng,ngy=2**lyng,ngz=2**lzng
  real, parameter :: flngx=ngx,flngy=ngy,flngz=ngz
  integer, parameter :: nfsize=10242, nfsize2=20480
  integer, parameter :: mic=252,mrec=50

  integer clock,clock1,clock0,clockend,nstep,n1,n2,i, j
  real t,h,h64,dt,vsc,time,rho,radx,fostar,bfs,const
  real, parameter :: pi = 3.141592654
  real lcube,nu,mu,mass,length
  integer :: seed(2)
  complex, allocatable :: ur(:,:,:),vr(:,:,:),wr(:,:,:)
  real, allocatable :: prdeno(:,:,:),qrfact(:,:,:)
  real, allocatable :: dxg(:),dyg(:),dzg(:)
  real, allocatable :: dxo(:),dyo(:),dzo(:)
  complex, allocatable :: dx(:,:,:),dy(:,:,:),dz(:,:,:)
  real, allocatable :: dsq(:,:,:)
  integer, allocatable :: firstn(:,:),number(:,:),nextn(:)
  real, allocatable :: xpl(:,:),fpl(:,:),xpi(:,:),xtn(:,:),ftn(:,:)
  real, allocatable :: xfi(:,:),xfn(:,:),frc(:,:),shpint(:,:),shpfs(:,:)
  integer, allocatable :: imic(:),irec(:,:),elmv(:,:),ihist(:)
  real,allocatable :: rlbs(:,:,:),rlbf(:,:)
  integer,allocatable :: elmnew(:,:),edr(:),ltb(:,:,:)
  logical denova
  integer k
  allocate(ur(0:ngx+2,0:ngy+2,0:ngz-1),vr(0:ngx+2,0:ngy+2,0:ngz-1))
  allocate(wr(0:ngx+2,0:ngy+2,0:ngz-1))
  allocate(qrfact(0:ngx+2,0:ngy+2,0:ngz-1),prdeno(0:ngx+2,0:ngy+2,0:ngz-1))
  allocate(dx(0:ngx+2,0:ngy+2,0:ngz-1),dy(0:ngx+2,0:ngy+2,0:ngz-1))
  allocate(dz(0:ngx+2,0:ngy+2,0:ngz-1),dsq(0:ngx+2,0:ngy+2,0:ngz-1))
  allocate(dxg(0:ngx+2),dyg(0:ngy+2),dzg(0:ngz+2))
  allocate(dxo(0:ngx+2),dyo(0:ngy+2),dzo(0:ngz+2))
  allocate(firstn(1:ngx,1:ngy),number(1:ngx,1:ngy))
  allocate(nextn(1:2*(ngx-2)*(ngz-2)+nfsize))
  allocate(xfi(1:3,1:nfsize))
  allocate(xfn(1:3,1:nfsize),frc(1:3,1:nfsize))
  allocate(xtn(1:3,1:2*(ngx-2)*(ngz-2)+nfsize))
  allocate(ftn(1:3,1:2*(ngx-2)*(ngz-2)+nfsize))
  allocate(xpi(1:3,1:2*(ngx-2)*(ngz-2)))
  allocate(xpl(1:3,1:2*(ngx-2)*(ngz-2)),fpl(1:3,1:2*(ngx-2)*(ngz-2)))
  allocate(elmnew(1:3,1:nfsize2),edr(1:nfsize2),ltb(1:2,1:mrec,1:mic))
  allocate(shpint(1:3,1:nfsize2),shpfs(1:7,1:nfsize2))
  allocate(imic(1:mic),irec(1:mrec,1:mic),elmv(1:12,1:mic),ihist(0:50))
  allocate(rlbs(1:3,1:mrec,1:mic),rlbf(1:mrec,1:mic))
  k = 2
  seed(1)=260
  seed(2)=215
  call random_seed
  call random_seed(size=k)
  call random_seed(put=seed(1:k))

  nstep = 250000
  lcube = 15.
  radx = 3.75
  nu = 8.0e+05
  rho = 1.0e-06
  dt = 1.0e-06
  clock = 0
  bfs = 1.*880.*(dt/1.e-4)**2
  const = 10.0
  ur=0.
  vr=0.
  wr=0.
  mu = nu*rho
  h  = lcube/flngx
  mass = rho*h**3 
  length = h
  time = dt
  vsc = nu/((length**2)/time)
  h64 = h/radx
  fostar = (mass*length/time**2)
  call inspher(lcube,radx,h,n1,n2,xfi,xfn,elmnew,shpint,shpfs)
  call inmv(lcube,radx,h,xfn,elmnew,n1,n2,elmv,imic,irec,rlbs,ltb)
  do i = 1,12
     print *, elmv(i,1)
  end do
  do i = 1,12
     print *, elmnew(:, elmv(i,1))
  end do
  do i = 1,12
        print *, xfn(:, elmv(i,1))
  end do

  call fmv(clock,mass,length,time,xfn,frc,elmv, &
       imic,irec,rlbs,rlbf,ihist,ltb)
  print *, mic
  deallocate(dx,dy,dz,dxo,dyo,dzo,dxg,dyg,dzg,dsq)
  deallocate(ur,vr,wr,edr,firstn,number,nextn,xfi,xfn,frc)
  deallocate(prdeno,qrfact,elmnew,shpint,shpfs)
  deallocate(xpi,xpl,fpl,xtn,ftn)
  deallocate(irec,imic,rlbs,ltb,ihist,rlbf)

end program rbc3d
