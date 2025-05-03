SUBROUTINE covariance(dist,angle,dt,out,ns_h,ns_t,cv_flag)
	USE global
	IMPLICIT NONE

	!Dummy variable declarations
	REAL,DIMENSION(:),INTENT(in) :: dist,angle,dt
	REAL,DIMENSION(:),INTENT(out) :: out
	INTEGER,INTENT(in) :: ns_h,ns_t,cv_flag

	!Local variable declarations
	REAL,DIMENSION(SIZE(dist)) :: sv,h_prime,c_h,c_t
	INTEGER :: i

	!Subroutine to calculate a covariance vector from a semivariogram model
	!and a vector of separation distances. Calculate semivariances function as
	!the sum of all the semivariogram contribution functions. This is purely a fn
	!of separation distance. Called before the 2d version, so contains error checks,
	!and extra documentation/comments. Semivariogram parameters are held in global module. 

	sv = tau(1)
	DO i=1,ns_h
		!Check for anisotropy
		IF (theta2(i) == theta3(i)) THEN
			h_prime = dist
		ELSE
			!Check minor semiaxis is less than major semiaxis
			IF (theta3(i)<theta2(i)) THEN
				PRINT '(/"Major semiaxis is less than minor semiaxis:")'
				PRINT '("Check parameter file!"/)'
!				PAUSE
				STOP
			END IF
			
			!Transform distance vector by applying anisotropy correction factor:
			!The minor axis of variation reaches sill before major axis of variation;
			!this is achieved by adding extra distance to the point separations 
			!to 'pull back' the range towards the origin. Separation distances 
			!are reprojected into a transformed coordinate system h'. 
			!First angle is rotated to the major axis of semivariance lies on the horizontal.
			!The anisotropy correction factor is calculated by multiplying the SIN of the 
			!angle between the major axis of variation and the point separation angle by the 
			!eccentricity of the anisotropy ellipsoid (major semiaxis/minor semiaxis). 
			!The transformed coordinate system h' is then the product of the separation 
			!and the correction factor.
			h_prime = dist+dist*ABS(SIN(angle-rotation(i))*(theta3(i)/theta2(i)))
		END IF

		SELECT CASE(structflag(i))
			CASE("Linear","linear","LINEAR")
				!Linear model: Only requires contribution parameter theta1
				sv = sv+theta1(i)*h_prime

			CASE("Power","power","POWER")
				!Power law model: Theta4 (the power parameter usually referred to as omega)
				!must be > 0 and < 2
				IF (theta4(i)>0 .AND. theta4(i)<2) THEN
					sv = sv+theta1(i)*h_prime**theta4(i)
				ELSE
					PRINT '("Power model must satisfy 0 < w < 2:")'
					PRINT '("Check theta2 in parameter file!"/)'
!					PAUSE
					STOP
				END IF

			CASE("Spherical","spherical","SPHERICAL")
				!Spherical model
				WHERE(h_prime < theta3(i))
					sv = sv+theta1(i)*((1.5*(h_prime/theta3(i)))-(0.5*(h_prime/theta3(i))**3))
				ELSEWHERE
					sv = sv+theta1(i)
				END WHERE

			CASE("Exponential","exponential","EXPONENTIAL")
				!Exponential model
				IF(theta4(i)<=0 .OR. theta4(i)>2) THEN
					PRINT '("Powered Exponential model requires 0 < w <= 2")'
					PRINT '("Check parameter file!"/)'
				END IF

				sv = sv+theta1(i)*(1-EXP(-3*(h_prime/theta3(i))**theta4(i)))
				
			CASE("Gaussian","gaussian","GAUSSIAN")
				!Gaussian model: MUST have a nugget effect to function without instability
				IF (tau(1) == 0) THEN
					PRINT '("Gaussian model requires a small nugget "/&
							"effect to avoid numerical instability:")'
					PRINT '("Check parameter file!"/)'
!					PAUSE
					STOP
				END IF

				sv = sv+theta1(i)*(1-EXP(-3*h_prime**2/theta3(i)**2))

			CASE("Quadratic","quadratic","QUADRATIC")
				!Rational quadratic model
				sv = sv+theta1(i)*h_prime**2/(1+h_prime**2/theta3(i))

			CASE("Hole","hole","HOLE")
				!Hole effect model
				WHERE(h_prime .NE. 0)
					sv = sv+theta1(i)*(1-COS((h_prime/theta3(i))*pi))
				END WHERE

			CASE("Dampened","dampened","DAMPENED")
				!Dampened hole effect model
				IF(theta4(i) < theta2(i)) THEN
					PRINT '("Dampening distance is less than spatial periodicity:")'
					PRINT '("Check theta4 in parameter file!"/)'
!					PAUSE
					STOP
				END IF

				WHERE(h_prime .NE. 0)
					sv = sv+theta1(i)*(1-EXP(-3*h_prime/theta4(i))*COS((h_prime/theta3(i))*pi))
				END WHERE

			CASE DEFAULT
				!Error message
				PRINT '("Unrecognised semivariogram structure: ",A20)',structflag(i)
				PRINT '("Check parameter file!"/)'
!				PAUSE
				STOP	

		END SELECT 
	END DO
	
	!Convert spatial sv to covariance by subtracting sv from the sill
	c_h=sill_h-sv

	!Temporal covariance
	sv = tau(2)
	DO i=ns_h+1,ns_h+ns_t		
		SELECT CASE(structflag(i))
			CASE("Linear","linear","LINEAR")
				!Linear model: Only requires contribution parameter theta1
				sv = sv+theta1(i)*dt

			CASE("Power","power","POWER")
				!Power law model: Theta4 (the power parameter usually referred to as omega)
				!must be > 0 and < 2
				IF (theta4(i)>0 .AND. theta4(i)<2) THEN
					sv = sv+theta1(i)*dt**theta4(i)
				ELSE
					PRINT '("Power model must satisfy 0 < w < 2:")'
					PRINT '("Check theta2 in parameter file!"/)'
!					PAUSE
					STOP
				END IF

			CASE("Spherical","spherical","SPHERICAL")
				!Spherical model
				WHERE(dt < theta2(i))
					sv = sv+theta1(i)* ((1.5*(dt/theta2(i))) - (0.5*(dt/theta2(i))**3))
				ELSEWHERE
					sv = sv+theta1(i)
				END WHERE

			CASE("Exponential","exponential","EXPONENTIAL")
				!Exponential model
				IF(theta4(i)<=0 .OR. theta4(i)>2) THEN
					PRINT '("Powered Exponential model requires 0 < w <= 2")'
					PRINT '("Check parameter file!"/)'
				END IF
				sv = sv+theta1(i)*(1-EXP(-3*(dt/theta2(i)))**theta4(i))

			CASE("Gaussian","gaussian","GAUSSIAN")
				!Gaussian model: MUST have a nugget effect to function without instability
				IF (tau(2) == 0) THEN
					PRINT '("Gaussian model requires a small nugget "/&
							"effect to avoid numerical instability:")'
					PRINT '("Check parameter file!"/)'
!					PAUSE
					STOP
				END IF

				sv = sv+theta1(i)*(1-EXP(-3*dt**2/theta2(i)**2))

			CASE("Quadratic","quadratic","QUADRATIC")
				!Rational quadratic model
				sv = sv+theta1(i)*dt**2/(1+dt**2/theta2(i))

			CASE("Hole","hole","HOLE")
				!Hole effect model
				WHERE(dt .NE. 0)
					sv = sv+theta1(i)*(1-COS((dt/theta2(i))*pi))
				END WHERE

			CASE("Dampened","dampened","DAMPENED")
				!Dampened hole effect model
				IF(theta4(i) < theta2(i)) THEN
					PRINT '("Dampening distance is less than spatial periodicity:")'
					PRINT '("Check theta4 in parameter file!"/)'
!					PAUSE
					STOP
				END IF
				
				WHERE(dt .NE. 0)
					sv = sv+theta1(i)*(1-EXP(-3*dt/theta4(i))*COS((dt/theta2(i))*pi))
				END WHERE

			CASE DEFAULT
				!Error message
				PRINT '("Unrecognised semivariogram structure: ",A20)',structflag(i)
				PRINT '("Check parameter file!"/)'
!				PAUSE
				STOP	

		END SELECT 
	END DO

	!Convert spatial sv to covariance by subtracting sv from the sill
	c_t=sill_t-sv

	SELECT CASE (cv_flag)
		CASE DEFAULT
			!Covariance model defaults to the product method
			out = k(1)*c_h*c_t
		CASE (1)
			!If cv_flag = 1, use the product-sum method
			out = (k(1)*c_h*c_t)+(k(2)*c_h)+(k(3)*c_t)
	END SELECT

END SUBROUTINE covariance