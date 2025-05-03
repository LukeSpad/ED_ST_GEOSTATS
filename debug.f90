SUBROUTINE debug(weights,controlpts,prediction,est,var,n,cv,localmu,alpha,beta)
	USE global
	IMPLICIT NONE

	!Dummy argument declarations
	INTEGER,INTENT(in) :: n
	REAL,INTENT(in) :: est,var,localmu,alpha,beta
	REAL,DIMENSION(:),INTENT(in) :: weights,cv
	TYPE(point),INTENT(in) :: prediction
	TYPE(point),DIMENSION(:),INTENT(in) :: controlpts

	!Local variable declarations
	INTEGER :: i

	!Summary of internal data for kriging. Results are displayed for all prediction
	!locations. Provided for diagnostic reasons.

	!Display observations
	PRINT '(/"Observations:"/1X,&
			"---------------------------------------------------------------------------------------------")'
	PRINT '(2X,"OID",8X,"SID",6X,"Value",6X,"Drift",4X,"Distance",5X,"Angle",&
			6X,"Time",6X,"Gamma",6X,"Weight",/1X,&
			"=============================================================================================")'
	DO i=1,n
		PRINT '(2X,I6,5X,I3,1X,5(F10.2,1X),F10.4,1X,F10.4)',&
		controlpts(i)%oid,controlpts(i)%sid,controlpts(i)%value,controlpts(i)%secondary,controlpts(i)%h,&
		(controlpts(i)%theta)/pi*180,controlpts(i)%dt,sill_gl-cv(i),weights(i)
	END DO
	PRINT '(1X,"---------------------------------------------------------------------------------------------")'

	PRINT '(11X,"Mean:",3(1X,F10.2),29X,"Sum:",1X,F10.4)',SUM(controlpts%value)/n,SUM(controlpts%secondary)/n,&
			SUM(controlpts%h)/n,SUM(weights(1:n))

	!Display prediction
	PRINT '(/"Estimate:"/1X,&
			"--------------------------------------------------------------------------------------------------------------")'
	PRINT '(3X,"East",7X,"North",7X,"Drift",7X,"Time",6X,"Local MU",3X,&
			"Estimate",4X,"Variance",3X,"St.Dev",5X"Alpha",5X,"Beta"/1X,&
			"==============================================================================================================")'

	SELECT CASE(method)
	CASE(1)
		PRINT '(1X,8(F10.4,1X),4X,"NA",8X,"NA")',prediction%x,prediction%y,prediction%secondary,prediction%t,localmu,est,var,SQRT(var)

	CASE DEFAULT
		PRINT '(1X,10(F10.4,1X))',prediction%x,prediction%y,prediction%secondary,prediction%t,localmu,est,var,SQRT(var),alpha,beta

	END SELECT
	PRINT '(1X,"--------------------------------------------------------------------------------------------------------------"/)'

	PRINT*,"Type Exit to quit, or hit enter to continue."
	PAUSE

END SUBROUTINE debug