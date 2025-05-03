SUBROUTINE writeout(prediction,localmu,est,var,alpha,beta,controlpts,p,db)
	USE global
	IMPLICIT NONE

	!Parameter declarations
	INTEGER,PARAMETER :: out=40

	!Dummy arguments
	INTEGER,INTENT(in) :: p,db
	REAL,INTENT(in) :: est,var,localmu,alpha,beta
	TYPE(point),INTENT(in) :: prediction
	TYPE(point),DIMENSION(:),INTENT(in) :: controlpts

	!Write results to output file
	IF (db .NE. 2) THEN
		WRITE(out,FMT=203) prediction%x,prediction%y,prediction%secondary,&
		prediction%t,localmu,est,var,SQRT(var),alpha,beta,controlpts(1)%h,&
		SUM(controlpts%h)/p,SUM(controlpts%secondary)/p,controlpts(1)%SID
	ELSE
		WRITE(out,FMT=204) prediction%sid,prediction%x,prediction%y,&
						   prediction%secondary,prediction%t,&
						   prediction%value,localmu,est,var,SQRT(var),alpha,beta,&
						   prediction%value-est,controlpts(1)%h,&
						   SUM(controlpts%h)/p,controlpts(1)%SID
	END IF

	!Format definitions
	203 FORMAT (2(F15.4,1X),F10.4,1X,F8.2,1X,6(F10.4,1X),3(F8.2,1X),I6)
	204 FORMAT (I5,1X,2(F15.4,1X),F10.4,1X,F8.2,1X,7(F10.4,1X),3(F8.2,1X),I5)

END SUBROUTINE writeout