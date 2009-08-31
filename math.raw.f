!**********************************************************************
SUBROUTINE findNormal(A, B, C, N )
IMPLICIT NONE

double precision A(1:3), B(1:3) , C(1:3), N(1:3)
double precision AC(1:3) , AB(1:3),BC(1:3), dot2,dot1,dot3
double precision    testi, testj, testk, dotprod

	   AC  = C - A ;	 
	   BC  = C - B
	   AB  = B - A ;	 

	  CALL findCross(AB , AC, N);
          CALL normalize(N);

RETURN 
END SUBROUTINE findNormal

!**********************************************************************

SUBROUTINE findCross(V1, V2, cross)

IMPLICIT NONE
	
double precision V1(1:3), V2(1:3), cross(1:3)

	cross(1) = V1(2) * V2(3) - V1(3) * V2(2) ;
	cross(2) = V1(3) * V2(1) - V1(1) * V2(3) ;
	cross(3) = V1(1) * V2(2)-  V1(2) * V2(1) ;

END SUBROUTINE findCross

!**********************************************************************
SUBROUTINE findDot(V1, V2, dot)

IMPLICIT NONE
	
double precision V1(1:3), V2(1:3), dot

	dot = V1(1) * V2(1) + V1(2) * V2(2) + V1(3) * V2(3) ;
	

END SUBROUTINE findDot
!**************************
SUBROUTINE findMagnitude(vector, mag)

IMPLICIT NONE

double precision vector(1:3), mag

	mag = sqrt(vector(1) *vector(1) + vector(2) * vector(2) + vector(3) * vector(3));

END SUBROUTINE findMagnitude

!**************************
SUBROUTINE normalize(vector)

IMPLICIT NONE

double precision vector(1:3)
double precision mag;

	CALL findMagnitude(vector, mag);
	vector(1) = vector(1)/(mag+1e-20) ;
	vector(2) = vector(2)/(mag+1e-20) ;
	vector(3) = vector(3)/(mag+1e-20) ;

END SUBROUTINE normalize
!**************************