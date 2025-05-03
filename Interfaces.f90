MODULE interfaces
	IMPLICIT NONE
	SAVE

	!Explicit interfaces
	INTERFACE
		SUBROUTINE dataframe(datafile,n)
		USE global
		IMPLICIT NONE
			INTEGER,INTENT(out) :: n
			CHARACTER(100),INTENT(in) :: datafile
		END SUBROUTINE dataframe
	END INTERFACE

	INTERFACE
		SUBROUTINE getparams(datafile,outfile,db,cv_flag,ns_h,ns_t,neighbors,window,n)	
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(out) :: ns_h,ns_t,neighbors,window,db,cv_flag,n
			CHARACTER(100),INTENT(out) :: datafile,outfile
		END SUBROUTINE getparams
	END INTERFACE

	INTERFACE
		SUBROUTINE covariance(dist,angle,dt,out,ns_h,ns_t,cv_flag)
			USE global
			IMPLICIT NONE
			REAL,DIMENSION(:),INTENT(in) :: dist,angle,dt
			REAL,DIMENSION(:),INTENT(out) :: out
			INTEGER,INTENT(in) :: ns_h,ns_t,cv_flag
		END SUBROUTINE covariance
	END INTERFACE

	INTERFACE
		SUBROUTINE covariance2d(dist,angle,dt,out,ns_h,ns_t,cv_flag)
			USE global
			IMPLICIT NONE
			REAL,DIMENSION(:,:),INTENT(in) :: dist,angle,dt
			REAL,DIMENSION(:,:),INTENT(out) :: out
			INTEGER,INTENT(in) :: ns_h,ns_t,cv_flag
		END SUBROUTINE covariance2d
	END INTERFACE

	INTERFACE
		SUBROUTINE separation(point1,point2,out,db)
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(in) :: db
			TYPE(point),INTENT(in) :: point2
			TYPE(point),DIMENSION(:),INTENT(in) :: point1
			TYPE(point),DIMENSION(:),INTENT(out) :: out
		END SUBROUTINE separation
	END INTERFACE

	INTERFACE	
		SUBROUTINE sorted(input,ns,nt,p,controlpts,hobs,last)
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(in) :: ns,nt,p
			INTEGER,DIMENSION(:),INTENT(in) :: last
			TYPE(point),DIMENSION(:),INTENT(inout) :: input
			TYPE(point),DIMENSION(:),INTENT(out) :: controlpts
			TYPE(pairs),DIMENSION(:,:),INTENT(out) :: hobs
		END SUBROUTINE sorted
	END INTERFACE

	INTERFACE
		SUBROUTINE sort(input,output,flag)
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(in) :: flag
			TYPE(point),DIMENSION(:),INTENT(in) :: input
			TYPE(point),DIMENSION(:),INTENT(out) :: output
		END SUBROUTINE sort
	END INTERFACE

	INTERFACE
		SUBROUTINE distance(obsdist,input,n)
			USE global
			IMPLICIT NONE
			TYPE(point),DIMENSION(:),INTENT(in) :: input
			TYPE(pairs),DIMENSION(:,:),INTENT(out) :: obsdist
			INTEGER,INTENT(in) :: n
		END SUBROUTINE distance
	END INTERFACE

	INTERFACE
		SUBROUTINE matinv(input,output)
			USE global
			IMPLICIT NONE
			REAL,DIMENSION(:,:),INTENT(in) :: input
			REAL,DIMENSION(:,:),INTENT(out) :: output
		END SUBROUTINE matinv
	END INTERFACE

	INTERFACE
		SUBROUTINE krige(cvobs,cvest,weights,controlpts,est,var,ck,n,localmu,alpha,beta,db)
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(in) :: n,db
			REAL,INTENT(out) :: est,var,ck,localmu,alpha,beta
			REAL,DIMENSION(:),INTENT(in) :: cvest
			REAL,DIMENSION(:),INTENT(out) :: weights
			REAL,DIMENSION(:,:),INTENT(in) :: cvobs
			REAL,DIMENSION(:),INTENT(in) :: controlpts
		END SUBROUTINE krige
	END INTERFACE

	INTERFACE
		SUBROUTINE debug(weights,controlpts,prediction,est,var,n,cv,localmu,alpha,beta)
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(in) :: n
			REAL,INTENT(in) :: est,var,localmu,alpha,beta
			REAL,DIMENSION(:),INTENT(in) :: weights,cv
			TYPE(point),INTENT(in) :: prediction
			TYPE(point),DIMENSION(:),INTENT(in) :: controlpts
		END SUBROUTINE debug
	END INTERFACE

	INTERFACE
		SUBROUTINE summary(input)
			USE global
			IMPLICIT NONE
			TYPE(point),DIMENSION(:),INTENT(in) :: input
		END SUBROUTINE summary
	END INTERFACE

	INTERFACE
		SUBROUTINE writeout(prediction,localmu,est,var,alpha,beta,controlpts,p,db)
			USE global
			IMPLICIT NONE
			INTEGER,INTENT(in) :: p,db
			REAL,INTENT(in) :: est,var,localmu,alpha,beta
			TYPE(point),INTENT(in) :: prediction
			TYPE(point),DIMENSION(:),INTENT(in) :: controlpts
		END SUBROUTINE writeout
	END INTERFACE

END MODULE interfaces