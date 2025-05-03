SUBROUTINE separation(point1,point2,out,db)
	USE global
	IMPLICIT NONE

	!Dummy variable declaration
	INTEGER,INTENT(in) :: db
	TYPE(point),INTENT(in) :: point2
	TYPE(point),DIMENSION(:),INTENT(in) :: point1
	TYPE(point),DIMENSION(:),INTENT(out) :: out

	!Subroutine to find the Euclidean distances between a point in space and
	!a vector of other points. Also returns angles if anisotropy present, 
	!or in debug mode.

	!Pass on values and object id
	out = point1

	!Find distance between points using pythagoras theorem
	out%h = SQRT((point1%x-point2%x)**2+(point1%y-point2%y)**2)
	
	!Find the temporal separation
	out%dt = point1%t-point2%t

	!Return the angle between points1 and point2 if anisotropy present in 
	!specified semivariance functions, or in debug mode.
	IF (ANY(theta2-theta3 .NE. 0) .OR. db==0) THEN
		WHERE ((point1%x-point2%x) .NE. 0)
			out%theta = ATAN((point1%y-point2%y)/(point1%x-point2%x))
		ELSEWHERE
			!if delta x =0, points on north/south vector: theta = 90 degrees
			out%theta = 1.570796
		END WHERE
	END IF

END SUBROUTINE separation