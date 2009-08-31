SUBROUTINE rayIntensity(ray , x , y , z , N , nm, omega_zero, epsilono ,z0, RADX, disp,opticalPower) 

IMPLICIT NONE
double precision nm, omega_zero, epsilono ,z0, RADX, disp
double precision ray(1:9) , z ,x , y ,  N(1:3), iRay(1:3)

double precision intensity  , opticalPower

double precision cos_theta ,PowerDist , PCoeff , PowerD , PowerDis

double precision omega

!parameter(omega_zero=$omega_zero$, epsilono=$epsilono$ ,opticalPower=$opticalPower$,nm=$nm$)
iRay = ray(4:6)
cos_theta = iRay(1)*N(1)+iRay(2)*N(2)+iRay(3)*N(3)

!	 x = (P(1) - 0.5000d0*LCUBE/H)*H*1e-6 !
!	 y = (P(2) - 0.5000d0*LCUBE/H)*H*1e-6 ! Subtract original center position
!	 z = (P(3) - 0.5000d0*LCUBE/H)*H*1e-6

	omega = omega_zero*(1+(z/z0)**2)**0.5
	PowerDist = nm**2*epsilono/2*(omega_zero/omega*exp(-(x-disp)**2/omega**2))**2
	PCoeff = 6.79246d20*opticalPower/0.1*5e-4/RADX  ! how to find this value when totalPower = 0.1
	PowerDis = PowerDist*PCoeff;
	PowerD =  PowerDis;
	intensity = PowerD * abs(cos_theta);
!	call calcReflectance(cos_theta, transmittance, n1 , n2 )   ! new transmittance function of incidence angle
   ray(7) = intensity
	ray(9) = acos(abs(cos_theta))


RETURN
END SUBROUTINE rayIntensity
!********************************
SUBROUTINE fresnelEquations(cos_theta, transmittance, n1 , n2 )

IMPLICIT NONE
double precision cos_theta , cos_theta_out, sin_theta
double precision n1 
double precision n2  
double precision Rperp , Rpara , Reflectance , incidentAngle , refractedAngle
double precision transmittance 

cos_theta = abs(cos_theta) ! Making sure WE ARE USING THE SMALLEST ANGEL BETWEEN 0 and PI/2
sin_theta = sqrt(1-cos_theta**2)
cos_theta_out = sqrt(1-((n1/n2)*sin_theta)**2) ! n1sin1=n1sin2 -> sin2 = n1/n2sin1 but cos2 = sqrt(1-sin2^2) = sqrt(1-(n1/n2sin1)^2)

If((n1.le.n2) ) then ! less to more , refraction
	Rperp = ((n1*cos_theta-n2*cos_theta_out)/((n1*cos_theta + n2*cos_theta_out)))**2
	Rpara = ((n1*cos_theta_out-n2*cos_theta)/((n1*cos_theta_out + n2*cos_theta)))**2

else !! more to less , possibility of reflection
	if ((sin_theta.ge.(n2/n1))) then ! total internal reflection
		Rperp = 1
		Rpara = 1
	else	! refraction
		Rperp = ((n1*cos_theta-n2*cos_theta_out)/((n1*cos_theta + n2*cos_theta_out)))**2
		Rpara = ((n1*cos_theta_out-n2*cos_theta)/((n1*cos_theta_out + n2*cos_theta)))**2

	end if
end if

Reflectance= (Rperp + Rpara)/2; 
Transmittance = 1 - Reflectance;

RETURN
END SUBROUTINE fresnelEquations

!**********************************************************************
SUBROUTINE rayStress(ray, Selement, x,y,z, N , n1 , n2, lightSpeed) 

IMPLICIT NONE
double precision ray(1:9) , N(1:3)
double precision Transmittance , trans, n1,n2 , lightSpeed
double precision Selement(1:4)  ! para. perpen then theta
double precision iRay(1:3), tRay(1:3),rRay(1:3)
!
double precision   intensity
double precision  iStress, rStress, tStress , z ,x , y
double precision Fparallel ,Fperpendicular,F
double precision  cos_theta , cos_theta_out , sin_theta, THETA

double precision cosiAngle , cosrAngle , siniAngle , sinrAngle ,costAngle , sintAngle

!ray(1) = x ; ! Update rays origins
!ray(2) = y ; 
!ray(3) = z ; 
iRay = ray(4:6)


cos_theta = iRay(1)*N(1)+iRay(2)*N(2)+iRay(3)*N(3)
sin_theta = sqrt(1-cos_theta**2)

cos_theta_out = sqrt(1-((n1/n2)*sin_theta)**2) ! n1sin1=n1sin2 -> sin2 = n1/n2sin1 but cos2 = sqrt(1-sin2^2) = sqrt(1-(n1/n2sin1)^2)
rRay = iRay - (2*(cos_theta))*N

if ( cos_theta.lt.-2e-2) then
tRay  = (n1/n2)*iRay-(cos_theta_out+n1/n2*cos_theta)*N
else
tRay  = (n1/n2)*iRay+(cos_theta_out-n1/n2*cos_theta)*N
end if

CALL normalize(tRay);

CALL normalize(rRay);

If((n1.gt.n2) ) then
! Possibility of Total internal refraction
if ((sin_theta.ge.(n2/n1))) then
	ray(4:6) = rRay	! Update rays direction vector
	tRay = rRay
else 
	ray(4:6) = tRay	! Update rays direction vector
end if
else
	ray(4:6) = tRay	! Update rays direction vector
end if


! Finding angles with laser axis 
cosiAngle  = iRay(3)
costAngle  = tRay(3)
cosrAngle  = (rRay(3))
!
siniAngle  = sqrt(1- cosiAngle**2)
sintAngle  = sqrt(1- costAngle**2)
sinrAngle  = sqrt(1- cosrAngle**2)

!	 x = (P(1) - 0.5000d0*LCUBE/H)*H*1e-6 !
!	 y = (P(2) - 0.5000d0*LCUBE/H)*H*1e-6
!	 z = (P(3) - 0.5000d0*LCUBE/H)*H*1e-6

  intensity = ray(7)
  cos_theta = cos(ray(9))	
  Transmittance = ray(8)
! We dont trace relfected rays therefore transmittance is zero and no effect
	iStress = intensity*n1*Transmittance/lightSpeed ;  ! transmittance upong hitting
! check this later
	call fresnelEquations(cos_theta, trans, n1 , n2 )   ! new transmittance function of incidence angle
	
	tStress = intensity*n2*Trans*Transmittance/lightSpeed; !
	
	rStress = intensity*n1*(1-trans)*Transmittance/lightSpeed;

	Transmittance = Transmittance  * trans					! updating transmittance Product of before
	
	ray(8) = Transmittance  
 	 
   Theta = atan2( y , x) 
 
	Fparallel     =       (iStress*cosiAngle - rStress*cosrAngle - tStress*costAngle)  ! possible minus value

	Fperpendicular =  (iStress*siniAngle - rStress*sinrAngle +  tStress*sintAngle) 
	
	F = (Fparallel**2 + Fperpendicular**2)**0.5;

	Selement(4) = Selement(4)  + F
	Selement(3) =  Theta 
	Selement(2) =  Selement(2) + Fperpendicular   
	Selement(1) =  Selement(1) + Fparallel  
RETURN
END SUBROUTINE rayStress


!*******************************************
SUBROUTINE rayTrace(XFN , Selements, elmnew , rays, numberOfRays, &
   nm, np,nbReflections,omega_zero, epsilono ,z0, RADX, H, disp,my_cap_center, opticalPower,lightSpeed)
IMPLICIT none

integer nbNodes, nbElements, numberOfRays, nbReflections
parameter(nbNodes=$nsnode$,nbElements=$nselm$)

double precision nm, np , stress(1:4)
double precision  XFN(1:3,1:nbNodes) , Selements(1:4,1:nbElements)
double precision A(1:3)  , B(1:3),C(1:3),FA(1:4)   , mag
double precision rays(1:9,1:numberOfRays), P(1:3),N1(1:3), nearestPoint(1:3), normal(1:3) , nearestNormal(1:3)
double precision e(1:3) , d(1:3),r(1:3), hither,  yon, t , rAngle , iAngle
logical result 
integer elmnew(1:3,1:nbElements), i  , iRays , iReflection, nearestElement , s , k
double precision nearestIntersection
double precision omega_zero, epsilono ,z0, RADX, disp,opticalPower, lightSpeed
double precision x , y , z, my_cap_center(1:3), H
hither = 1e-5
yon	 = 100.

do iReflection = 1 , nbReflections
!print *, iReflection, numberOfRays
	! Loop over rays
	do iRays = 1, numberOfRays
		e   =  rays(1:3,iRays);   ! point
		d   =  rays(4:6,iRays);	! direction vector

		result = .false.
		t = 0. 
		nearestIntersection = 100.
		! Loop over elements to check for intersection
		do i = 1 , nbElements  
			A   = (/ XFN(1,elmnew(1,i)) , XFN(2,elmnew(1,i)) ,  XFN(3,elmnew(1,i)) /)
			B   = (/ XFN(1,elmnew(2,i)) , XFN(2,elmnew(2,i)) ,  XFN(3,elmnew(2,i)) /)
			C   = (/ XFN(1,elmnew(3,i)) , XFN(2,elmnew(3,i)) ,  XFN(3,elmnew(3,i)) /)
			CALL raytri_intersect(A,B,C, e, d,  hither,  yon, result , P ,normal,t)
			if (result.eqv.(.true.))then
						if ((t < nearestIntersection).and.(t.gt.1e-3)) then
						nearestIntersection = t;
						nearestElement = i;
						nearestPoint = P;
						nearestNormal = normal ;
						end if
			end if
			result = .false.
		end do
! 	End  Loop over elements , nearest intersection found

	! Update ray vector and point
	! print *, ' intersecion '
    mag = (P(1)**2+P(2)**2+P(3)**2)
	if ((nearestIntersection.gt.0).and.(mag.ne.0).and.(nearestIntersection.lt.100.)) then
			stress = Selements(1:4,nearestElement)
			rays(1,iRays) = nearestPoint(1)
			rays(2,iRays) = nearestPoint(2)
			rays(3,iRays) = nearestPoint(3)

			x = (nearestPoint(1) - my_cap_center(1))*H*1d-2 ! center = 0.5000d0*LCUBE/H
			y = (nearestPoint(2) - my_cap_center(2))*H*1d-2 !
			z = (nearestPoint(3) - my_cap_center(3))*H*1d-2 !
			
			if (iReflection.eq.1) then
				call rayIntensity(rays(1:9,iRays) , x,y,z , nearestNormal , nm, omega_zero, epsilono ,z0, RADX, disp,opticalPower) 
			end if
			if ( mod(iReflection,2).eq.1) then
			call rayStress(rays(1:9,iRays), stress,x,y,z, nearestNormal,nm,np,lightSpeed) !n one and n two has to swap
			else 
			call rayStress(rays(1:9,iRays), stress,x,y,z, nearestNormal,np,nm,lightSpeed) !n one and n two has to swap
			end if
			Selements(1:4,nearestElement) = stress 
	end if

	end do
	end do
	
!print *,'done rayTrace'


end SUBROUTINE rayTrace
!***************************************
SUBROUTINE raytri_intersect(A,B,C, e, d,  hither,  yon, result, P, normal , t)

IMPLICIT NONE

    double precision aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, t, beta, gamma, denom;
    double precision eihf, gfdi, dheg, akjb, jcal, blkc;
	 double precision P(1:3), N1(1:3), normal(1:3)
	 double precision A(1:3) , B(1:3) , C(1:3)
	 double precision e(1:3), d(1:3)
	 double precision hither,  yon
	 logical result
	 result  = .true.
	 P = 10000.
    aa = A(1) - B(1);
    bb = A(2) - B(2);
    cc = A(3) - B(3);
    dd = A(1) - C(1);
    ee = A(2) - C(2);
    ff = A(3) - C(3);
    gg = d(1);
    hh = d(2);    ! V  = (gg, hh, ii)
    ii = d(3);
    jj = A(1) - e(1);
    kk = A(2) - e(2);
    ll = A(3) - e(3);
    eihf = ee*ii-hh*ff;
    gfdi = gg*ff-dd*ii;
    dheg = dd*hh-ee*gg;
    akjb = aa*kk-jj*bb;
    jcal = jj*cc-aa*ll;
    blkc = bb*ll-kk*cc;
    
	 denom = 1. / (aa*eihf+bb*gfdi+cc*dheg + 1e-20);
    
	 t = -(ff*akjb+ee*jcal+dd*blkc) * denom;  ! t = -(P0.N+d)/(V.N) ; P = P0+tV
  
	 if ((t.lt.hither).or.(t.gt.yon)) then ! early termination condintion 
	 result = .false.
	 RETURN
	 end if
	 
	 gamma = (ii*akjb+hh*jcal+gg*blkc)*denom;
	
    if ((gamma .lt. (0.)) .or. (gamma.gt. 1.0)) then
	 result = .false.
	 return
	 end if

    beta = (jj*eihf+kk*gfdi+ll*dheg) * denom;

    if ((beta .lt.(0.)) .or. (beta .gt.(1.0 - gamma))) then
     result = .false.
	 return
	 end if

	  P = e + d * t
	 CALL findNormal(A,B,C,normal)  
     result = .true.
return
end SUBROUTINE raytri_intersect
!************************************
SUBROUTINE initializeRays(rays , numberOfRays )
IMPLICIT NONE
! no need for first hit 
double precision  rays(1:9,1:numberOfRays) 
integer  k , numberOfRays
double precision  x , y  , z
open(232,file='rays.dat',status='unknown')
!print *, numberOfRays
do k = 1 , numberOfRays
	read(232,*) x, y ,z
	rays(1,k) =  x ! P0_x   
	rays(2,k) =  y ! P0_y
	rays(3,k) =  z + 1. ! P0_z 
	rays(4,k) =  0.;    ! V_x
	rays(5,k) =  0.;	  ! V_y 
	rays(6,k) =  -1.0;  ! V_z
	rays(7,k)  = 0.0 ;  ! intensity = 0 
 	rays(8,k)  = 1.0    ! transmitanece  = 1  
	rays(9,k)  = 0.0    ! cos(theta) forgot y 
end do
!print *,'done initializeRays'
close(232)
RETURN
END SUBROUTINE initializeRays
!**********************************************************************
SUBROUTINE findRays (XFN , elmnew  , zc , numberOfrays)
IMPLICIT NONE
integer iElement, numberOfrays, nbElements, nbNodes

parameter(nbNodes=$nsnode$,nbElements=$nselm$)

double precision N(1:3,1:nbElements) , zc, PI

double precision A(1:3),B(1:3),C(1:3),center(1:3),tail(1:3) , t
double precision XFN(1:3,1:nbNodes)
integer elmnew(1:3,1:nbElements)

open(232,file='rays.dat')
PI = 3.14159265358979323846d0 ! Taken from Wikipedia; 20 digits

numberOfrays = 0

do iElement = 1, nbElements! Loop over rays
			A   = (/ XFN(1,elmnew(1,iElement)) , XFN(2,elmnew(1,iElement)) ,  XFN(3,elmnew(1,iElement)) /)
			B   = (/ XFN(1,elmnew(2,iElement)) , XFN(2,elmnew(2,iElement)) ,  XFN(3,elmnew(2,iElement)) /)
			C   = (/ XFN(1,elmnew(3,iElement)) , XFN(2,elmnew(3,iElement)) ,  XFN(3,elmnew(3,iElement)) /)
			CALL findNormal(A, B, C , N(1:3,iElement))
			if((acos(N(3,iElement))*180/PI.le.90)) then
			center = (A+B+C)/3
			if (center(3).ge.zc) then
			numberOfrays = numberOfrays + 1
			write(232, *) center(1), center(2),center(3)
			end if
			end if
end do
!print *, nbElements, numberOfRays
!print *,'done initializeRays'

close (232)

END SUBROUTINE findRays

!*****************************************
      SUBROUTINE capsuleForce(XFN , Fnodes, shpfs , elmnew ,rays , FOSTAR, zc, &
	   RADX, H,my_cap_center,z0,disp,numberOfrays)	
		IMPLICIT NONE
		
		integer nbNodes, nbElements, numberOfrays
		double precision  nm, np, omega_zero, epsilono ,z0, RADX, H, disp,my_cap_center(1:3), opticalPower,lightSpeed
		integer nbReflections
		parameter(nbNodes=$nsnode$,nbElements=$nselm$)
		parameter(opticalPower=$opticalPower$,lightSpeed=$lightSpeed$)
		parameter(nm=$nm$,np=$np$,nbReflections=$numberOfMaxReflections$,omega_zero=$omega_zero$, epsilono=$epsilono$)

		integer elmnew(1:3,1:nbElements), ii
		double precision  XFN(1:3,1:nbNodes) , Fnodes(1:3,1:nbNodes)
		
		double precision  shpfs(1:7,1:nbElements), FX ,FY ,FZ

		double precision Felements(1:3,1:nbElements),Selements(1:4,1:nbElements)
		double precision x,y,z, phi,theta
		double precision	F1i ,	F1t ,	F1r ,	F2i,	F2t ,	F2r ,Fparallel ,Fperpendicular   
		double precision total, Fphi, FOSTAR
		REAL (8) FX_Back	,	FZ_Back	, 	FX_FRONT	,	FZ_FRONT, xc,yc,zc
		double precision rays(1:9,1:numberOfRays)
		integer i,j ,k
		

		Felements = 0d0
		Fnodes = 0d0
		Selements = 0d0
!print *,'rayTrace'
    CALL rayTrace(XFN , Selements, elmnew , rays, numberOfRays, &
   nm, np,nbReflections,omega_zero, epsilono ,z0, RADX, H, disp,my_cap_center, opticalPower,lightSpeed)
!print *,'End rayTrace'

			
		total = 0
		!RADIUS = RADX*1e-6
		FZ_FRONT = 0D0
		FX_FRONT = 0D0
		FZ_Back  = 0D0
		FX_Back  = 0D0
		Fphi = 0d0
		
		! Calculating force from stress vector and nondimensiolaize
		do i=1, nbElements
			!stressMag = Selements(4,i) !sqrt( Selements(1,i)**2 + Selements(2,i)**2) 
			Felements(1,i) = Selements(4,i)*cos(Selements(3,i))*1d1/FOSTAR*shpfs(7,i)/2 
			Fphi = atan2(Selements(2,i),Selements(1,i))
			Felements(3,i) = Selements(4,i)*cos(Fphi)*1d1/FOSTAR*shpfs(7,i)/2 
		end do
	
		! Interpolating element to nodes forces
		do i=1,nbNodes
				Fnodes(1,i) = 0d0
				Fnodes(2,i) = 0d0
				Fnodes(3,i) = 0d0
				ii=0
				do j=1,nbElements
					do k=1,3
						if (elmnew(k,j) == i) then
							ii=ii+1
							Fnodes(1,i)=Fnodes(1,i)+Felements(1,j)/3
							Fnodes(2,i)=Fnodes(2,i)+Felements(2,j)/3 
							Fnodes(3,i)=Fnodes(3,i)+Felements(3,j)/3 
						end if					
					end do
				end do
		end do
	
	FZ_FRONT = 0D0
	FX_FRONT= 0D0
	FZ_Back = 0D0
	FX_Back = 0D0
	
	! Calculating total back and front forces
	do i=1, nbNodes
      
	 z = XFN(3,i)
	if (z.le.zc) then
		FX_Back	=  FX_Back	+  Fnodes(1,I)
		FZ_Back	=  FZ_Back	+  Fnodes(3,I)	
	else
		FX_Front =  FX_Front +  Fnodes(1,I)
		FZ_Front =  FZ_Front +  Fnodes(3,I)
	end if
	end do

!	print *, 'Direction		: ' , 'X',							'					Z'
!	print *, 'Using ray tracer front part : ' , FX_Front *FOSTAR/1000 , FZ_Front *FOSTAR/1000
!	print *, 'Using ray tracer back part : '  , FX_Back  *FOSTAR/1000  , FZ_Back *FOSTAR/1000
!	print *, 'Using ray tracer Net forces : '  , (FX_Front + FX_Back) *FOSTAR /1000 , (FZ_Front+FZ_Back)*FOSTAR/1000

!  saving back and front total forces 

   write(1981,*)  FX_Front *FOSTAR/1d-7 , FZ_Front *FOSTAR/1d-7
   write(1982,*)  FX_Back  *FOSTAR/1d-7  , FZ_Back *FOSTAR/1d-7

   write(*,*)  FX_Front *FOSTAR/1d-7 , FZ_Front *FOSTAR/1d-7
   write(*,*)  FX_Back  *FOSTAR/1d-7  , FZ_Back *FOSTAR/1d-7

END SUBROUTINE capsuleForce

