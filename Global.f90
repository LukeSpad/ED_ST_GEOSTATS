MODULE global
	IMPLICIT NONE
	SAVE

	!Derived types
	TYPE point
		INTEGER :: oid,sid
		REAL :: x,y,t,value,secondary,h,dt,theta
	END TYPE point

	TYPE pairs
		REAL :: h,dt,theta
	END TYPE pairs

	!Global parameter declarations
	REAL,PARAMETER :: pi=3.141592654

	!Global variable declarations
	INTEGER :: method
	REAL :: sill_h,sill_t,sill_gl
	REAL,DIMENSION(2) :: tau
	REAL,DIMENSION(3) :: k
	CHARACTER(20),DIMENSION(7) :: names
	LOGICAL :: correction=.TRUE.

	!Global array declarations
	TYPE(point),DIMENSION(:),ALLOCATABLE :: observations
	REAL,DIMENSION(:),ALLOCATABLE :: theta1,theta2,theta3,theta4,rotation
	CHARACTER(20),DIMENSION(:),ALLOCATABLE :: structflag

END MODULE global