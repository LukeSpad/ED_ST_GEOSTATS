SUBROUTINE dataframe(datafile,n)
	USE global
	USE interfaces
	IMPLICIT NONE

	!Dummy argument definitions
	INTEGER,INTENT(out) :: n
	CHARACTER(100),INTENT(in) :: datafile

	!Local parameter declarations
	INTEGER,PARAMETER :: in=35

	!Local variable declarations
	INTEGER :: ios,i
	CHARACTER(100) :: dataname

	!Initialises the data frame: All observations are read into a global vector

	!Open observations file
	OPEN(UNIT=in,FILE=datafile,IOSTAT=ios)
	IF(ios .NE. 0) THEN
		PRINT '("Invalid data file"/)'
!		PAUSE
		STOP
	END IF

	!Skip header
	READ(in,*),dataname
	READ(in,*),names

	PRINT '(/"Data file:")'
	PRINT*, TRIM(datafile),dataname

	!Find total number of data points
	n=-1
	DO WHILE (ios==0)
		READ(in,*,IOSTAT=ios)
		n=n+1
	END DO

	!Go back to the beginning and skip header
	REWIND(in)
	DO i=1,2
		READ(in,*)
	END DO

	!Build observations vector
	ALLOCATE(observations(n))
	DO i=1,n
		READ(in,*) observations(i)%oid,observations(i)%sid,observations(i)%x,observations(i)%y,&
				   observations(i)%t,observations(i)%value,observations(i)%secondary
	END DO

	!Close the file
	CLOSE(in)

	!Produce summary statistics
	CALL summary(observations)

END SUBROUTINE dataframe