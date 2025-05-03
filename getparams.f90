SUBROUTINE getparams(datafile,outfile,db,cv_flag,ns_h,ns_t,neighbors,window,n)
	USE global
	IMPLICIT NONE

	!Dummy argument declarations
	INTEGER,INTENT(out) :: ns_h,ns_t,neighbors,window,db,cv_flag,n
	CHARACTER(100),INTENT(out) :: datafile,outfile

	!Local parameter declarations
	INTEGER,PARAMETER :: in=30,out=40

	!Local variable declarations
	INTEGER :: i,ios
	CHARACTER(100) :: parameterfile,gridfile

	!A subroutine to retrieve and report the kriging parameters from a specified
	!parameter file. Semivariance parameters are held in the global module for
	!ease of access by other subroutines.

	!Read parameters
	PRINT '("Enter parameter file name")'
	DO
		READ '(A)', parameterfile
		IF(parameterfile=="QUIT") STOP
		OPEN(UNIT=in,FILE=parameterfile,STATUS="OLD",IOSTAT=ios)
		IF(ios==0) EXIT
		PRINT '(/"Unable to read parameter file: Please try again, or type QUIT to exit."/)'
	END DO
	PRINT*, ""

	!Read file names and parameters
	READ(UNIT=in,FMT=101) datafile,gridfile,outfile
	READ(in,'(I//)')method
	READ(in,'(I//)')db
	READ(in,*)neighbors,window
	READ(in,'(/)')
	READ(in,*)ns_h,ns_t
	READ(in,'(/)')
	READ(in,*)tau(1),tau(2)
	READ(in,'(//F//)')sill_gl
	READ(in,'(I///)')cv_flag

	!Prepare the observations vector
	CALL dataframe(datafile,n)

	!Allocate semivariogram parameter arrays
	ALLOCATE(structflag(ns_h+ns_t))
	ALLOCATE(theta1(ns_h+ns_t))
	ALLOCATE(theta2(ns_h+ns_t))
	ALLOCATE(theta3(ns_h+ns_t))
	ALLOCATE(theta4(ns_h+ns_t))
	ALLOCATE(rotation(ns_h))

	!Get spatial semivariogram parameters
	DO i=1,ns_h
		READ(in,*,IOSTAT=ios) structflag(i),theta1(i),theta2(i),theta3(i),theta4(i),rotation(i)
 		IF(ios .NE. 0) THEN
			PRINT '("Error in parameters file: Improper spatial semivariance model specification"/)'
!			PAUSE
			STOP
		END IF
	END DO

	READ(in,'(//)')

	!Get temporal semivariogram parameters
	DO i=ns_h+1,ns_h+ns_t
		READ(in,*,IOSTAT=ios) structflag(i),theta1(i),theta2(i),theta4(i)
		theta3(i)=theta2(i)
 		IF(ios .NE. 0) THEN
			PRINT '("Error in parameters file: Improper temporal semivariance model specification"/)'
!			PAUSE
			STOP
		END IF
	END DO

	!Close the file
	CLOSE(UNIT=in)

	!Open grid file
	IF(db .NE. 2) THEN
		OPEN(UNIT=in,FILE=gridfile,STATUS="OLD",IOSTAT=ios)
		IF(ios .NE. 0) THEN
			PRINT '("Invalid prediction locations file!"/)'
!			PAUSE
			STOP
		END IF
	ELSE
		OPEN(UNIT=in,FILE=datafile,STATUS="OLD",IOSTAT=ios)
		IF(ios .NE. 0) THEN
			PRINT '("Invalid data file!"/)'
!			PAUSE
			STOP
		END IF
		
		!Skip header	
		READ(in,*)
	END IF

	!Skip header
	READ(in,*)

	!Open output file
	OPEN(UNIT=out,FILE=outfile,IOSTAT=ios)
	IF(ios .NE. 0) THEN
		PRINT '("Invalid output file name!"/)'
!		PAUSE
		STOP
	END IF

	!Report interpolation Method
	SELECT CASE(method)
	CASE DEFAULT
		PRINT '(/"Interpolation via Simple Kriging.")'
	CASE(1)
		PRINT '(/"Interpolation via Ordinary Kriging.")'
	CASE(2)
		PRINT '(/"Interpolation via Kriging with a Trend.")'
	END SELECT

	!Report operating mode
	IF (db == 0) THEN
		PRINT '("Operating in debug mode."/)'
	ELSE IF (db == 1) THEN
		PRINT '("Operating in default mode."/)'
	ELSE
		PRINT '("Operating in jack knife mode."/)'
	END IF

	!Summarise semivariance model specification
	PRINT '(/"Spatial semivariogram model with",1X,I1,1X,"Structures:"/)',ns_h 
	PRINT '("---------------------------------------------------------------------------------")'
	PRINT '("Structure",18X,"Theta1",5X,"Theta2",6X,"Theta3",6X,"Theta4",3X,"Rotations")'
	PRINT '("=================================================================================")'
	PRINT '("Nugget",17X,F10.3)',tau(1)
	DO i=1,ns_h
		SELECT CASE(structflag(i))
			CASE ("Dampened","dampened","DAMPENED")
				PRINT '(A20,5(2X,F10.2))',structflag(i),theta1(i),theta2(i),theta3(i),theta4(i),rotation(i)

			CASE ("Linear","linear","LINEAR")
				PRINT '(A20,2X,F10.2,4(10X,"NA"))',structflag(i),theta1(i)

			CASE DEFAULT
				PRINT '(A20,3(2X,F10.2),10X,"NA",2X,F10.2)',structflag(i),theta1(i),theta2(i),theta3(i),rotation(i)

		END SELECT
	END DO
	PRINT '("---------------------------------------------------------------------------------")'

	PRINT '(/"Temporal semivariogram model with",1X,I1,1X,"Structures:"/)',ns_t 
	PRINT '("---------------------------------------------------------")'
	PRINT '("Structure",18X,"Theta1",5X,"Theta2",6X,"Theta4")'
	PRINT '("=========================================================")'
	PRINT '("Nugget",17X,F10.3)',tau(2)
	DO i=ns_h+1,ns_h+ns_t
		SELECT CASE(structflag(i))
			CASE ("Dampened","dampened","DAMPENED")
				PRINT '(A20,3(2X,F10.2))',structflag(i),theta1(i),theta2(i),theta4(i)

			CASE ("Linear","linear","LINEAR")
				PRINT '(A20,2X,F10.2,2(10X,"NA"))',structflag(i),theta1(i)

			CASE DEFAULT
				PRINT '(A20,2(2X,F10.2),10X,"NA")',structflag(i),theta1(i),theta2(i)

		END SELECT
	END DO
	PRINT '("---------------------------------------------------------"/)'
	
	!Calculate the sill variances
	sill_h = SUM(theta1(1:ns_h))+tau(1)
	sill_t = SUM(theta1(ns_h+1 : ns_h+ns_t))+tau(2)
	
	!Screen prompt for Deutsch non-zero corrections; convex estimator
	!Deutsch, C.V. (1996). Computers & Geosciences 22(7):765-773
	IF(correction==.TRUE.) THEN
		PRINT '("Using Deutsch convex estimator corrections: Estimator constrained to data range.")'
		PRINT '("Deutsch, C.V. (1996). Computers & Geosciences 22(7):765-773"/)'
	END IF

	!Calculate constants for covariance models
	SELECT CASE (cv_flag)
		CASE DEFAULT
			!Covariance model defaults to the product method
			PRINT '("Using the Product covariance model."/)'
			k(1) = sill_gl/(sill_h*sill_t)
		CASE (1)
			!Product-sum covariance parameters
			PRINT '("Using the Product-Sum covariance model:")'
			IF (db==0) THEN
				PRINT '("-------------------------")'
				PRINT '("Parameters")'
				PRINT '("=========================")'
				PRINT '("Spatial sill  ="1X,F8.3)',sill_h
				PRINT '("Temporal sill ="1X,F8.3)',sill_t
			END IF
			PRINT '("Global sill   ="1X,F8.3/)',sill_gl
			
			k(1) = (sill_h+sill_t-sill_gl)/(sill_h*sill_t)
			k(2) = (sill_gl-sill_t)/sill_h
			k(3) = (sill_gl-sill_h)/sill_t
			
			IF (db==0) THEN
				PRINT '("k1 ="1X,F8.3)',k(1)
				PRINT '("k2 ="1X,F8.3)',k(2)
				PRINT '("k3 ="1X,F8.3)',k(3)
				PRINT '("-------------------------")'
			END IF

			IF(k(1)<=0 .OR. k(2)<0 .OR. k(3)<0) THEN
				PRINT '("Error in parameters file: Improper semivariance model specification")'
				PRINT '("Check sills produce permissable k values:"/)'
				PRINT '("k1 must be strictly (non-zero) positive,"/"whilst k2 and k3 must be non-negative."/)'
!				PAUSE
				STOP
			END IF
	END SELECT

	!Convert rotations to radians
	rotation=rotation/180*pi
	
	!Write header on output file
	IF(db .NE. 2) THEN
		WRITE(out,FMT=201) TRIM(names(7))
	ELSE
		WRITE(out,FMT=202) TRIM(names(7)),TRIM(names(6))
	END IF

	!Format definitions
	!Inputs
	101 FORMAT (/A100/A100/A100//)

	!Outputs
	201 FORMAT (3X,"East",3X,"North",3X,A20,3X,"Time",3X,"Local_Mu",3X,"Estimate",3X,&
				"Variance",3X,"St.Dev."3X,"Alpha",3X,"Beta",3X,"Nearest",3X,"Mu_Dist",3X,&
				"Mu_Drift",3X,"Nearest.ID")
	202 FORMAT (3X,"SID",3X,"East",3X,"North",3X,A20,3X,"Time",3X,A20,3X,&
				"Local_Mu",3X,"Estimate",3X,"Variance",3X,"St.Dev."3X,"Alpha",3X,&
				"Beta",3X,"Error",3X,"Nearest",3x,"Mu_Dist",3X,"Nearest.ID")

END SUBROUTINE getparams