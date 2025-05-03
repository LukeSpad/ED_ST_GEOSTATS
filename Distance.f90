SUBROUTINE distance(output,input,n)
	USE global
	USE interfaces
	IMPLICIT NONE

	!Dummy argument definitions
	TYPE(point),DIMENSION(:),INTENT(in) :: input
	TYPE(pairs),DIMENSION(:,:),INTENT(out) :: output
	INTEGER,INTENT(in) :: n

	!Local variable declarations
	INTEGER :: i,j

	!Subroutine to find the Euclidean distances, between two vectors of points.
	
	DO i=1,n
		DO j=i,n
	
		!Find spatial distance between points using pythagoras theorem
		output(i,j)%h = SQRT((input(i)%x-input(j)%x)**2+(input(i)%y-input(j)%y)**2)
		output(j,i)%h = output(i,j)%h 

		!Find the temporal separation
		output(i,j)%dt = input(i)%t-input(j)%t
		output(j,i)%dt = output(i,j)%dt

		!Calculate angle in radians between points if anisotropy specified in any 
		!of the semivariogram functions.
		IF (ANY(theta2-theta3 .NE. 0)) THEN
			IF ((input(i)%x-input(j)%x) .NE. 0) THEN
				output(i,j)%theta = ATAN((input(i)%y-input(j)%y)/(input(i)%x-input(j)%x))
				output(j,i)%theta = output(i,j)%theta 
			
			ELSE
				!if delta x =0, points on north/south vector: theta = 90 degrees
				output(i,j)%theta = 1.570796
				output(j,i)%theta = output(i,j)%theta
				
			END IF
		END IF

		END DO
	END DO

END SUBROUTINE distance