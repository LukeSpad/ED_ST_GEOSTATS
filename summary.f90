SUBROUTINE summary(input)
	USE global
	USE interfaces
	IMPLICIT NONE

	!Dummy arguments
	TYPE(point),DIMENSION(:),INTENT(in) :: input

	!Local variable declarations
	REAL :: n
	REAL,DIMENSION(2) :: mu,var,sd,skew,cvar,kurt,median,q1,q3,iqr,ps,jb
	TYPE(point),DIMENSION(:),POINTER :: obs

	!A subroutine to summarise a vector of input data points. A fairly complete
	!set of descriptive statistics is produced and printed to screen.

	!Generate and report global summary stats
	n = SIZE(input)
	ALLOCATE(obs(NINT(n)))
	CALL sort(input,obs,1)
	
	!Parametric stats of normal distribution:
	mu(1)   = SUM(obs%value)/n
	mu(2)   = SUM(obs%secondary)/n
	var(1)  = SUM((obs%value - mu(1))**2)/(n-1)
	var(2)  = SUM((obs%secondary - mu(2))**2)/(n-1)
	sd(1)   = SQRT(var(1))
	sd(2)   = SQRT(var(2))
	
	!Robust stats:
	IF (MOD(n,2.0)==0.0) THEN
		median(1) = (obs(n/2)%value+obs(n/2+1)%value)/2
		median(2) = (obs(n/2)%secondary+obs(n/2+1)%secondary)/2
	ELSE
		median(1) = obs(NINT(n/2))%value
		median(2) = obs(NINT(n/2))%secondary
	END IF
	
	IF (MOD(n,2.0)==0.0 .OR. MOD(n,2.0)==0.5) THEN
		q1(1)     = (obs(FLOOR(n/4))%value+obs(FLOOR(n/4+1))%value)/2
		q1(2)     = (obs(FLOOR(n/4))%secondary+obs(FLOOR(n/4+1))%secondary)/2
		q3(1)     = (obs(FLOOR((3*n)/4))%value+obs(FLOOR((3*n)/4+1))%value)/2
		q3(2)     = (obs(FLOOR((3*n)/4))%secondary+obs(FLOOR((3*n)/4+1))%secondary)/2
	ELSE
		q1(1)     = obs(CEILING(n/4))%value
		q1(2)     = obs(CEILING(n/4))%secondary
		q3(1)     = obs(CEILING((3*n)/4))%value
		q3(2)     = obs(CEILING((3*n)/4))%secondary

	END IF
	
	iqr(1)    = ABS(q3(1)-q1(1))
	iqr(2)    = ABS(q3(2)-q1(2))	

	!Higher order stats:
	skew(1) = (SUM((obs%value - mu(1))**3)/n)/sd(1)**3
	skew(2) = (SUM((obs%secondary - mu(2))**3)/n)/sd(2)**3
	ps(1)   = (3*(mu(1)-median(1)))/sd(1)
	ps(2)   = (3*(mu(2)-median(2)))/sd(2)
	kurt(1) = ((SUM((obs%value - mu(1))**4)/n)/var(1)**2)-3
	kurt(2) = ((SUM((obs%secondary - mu(2))**4)/n)/var(2)**2)-3
	cvar(1) = sd(1)/mu(1)
	cvar(2) = sd(2)/mu(2)
	
	!Jarque-Bera test
	jb(1)   = (n/6)*(skew(1)**2+(kurt(1)**2/4))
	jb(2)   = (n/6)*(skew(2)**2+(kurt(2)**2/4))

	PRINT '(/"----------------------------------------------------"/)'
	PRINT*,  "Statistic                       ",TRIM(ADJUSTL(names(6))),"     ",TRIM(ADJUSTL(names(7)))
	PRINT '( "===================================================="/&
			 "n                           =",5X,I6//&
	         "Min                         =",1X,F10.4,1X,F10.4/&
			 "First Quartile              =",1X,F10.4,1X,F10.4/&
			 "Median                      =",1X,F10.4,1X,F10.4/&
			 "Third Quartile              =",1X,F10.4,1X,F10.4/&
			 "Max                         =",1X,F10.4,1X,F10.4//&
			 "Interquartile range         =",1X,F10.4,1X,F10.4//&
			 "Mean                        =",1X,F10.4,1X,F10.4/&
			 "Variance                    =",1X,F10.4,1X,F10.4/&
			 "Standard Deviation          =",1X,F10.4,1X,F10.4//&
			 "Skew                        =",1X,F10.4,1X,F10.4/&
			 "Pearsons skewness           =",1X,F10.4,1X,F10.4/&
			 "Sample Kurtosis             =",1X,F10.4,1X,F10.4/&
			 "Coefficient of Variation    =",1X,F10.4,1X,F10.4//&
			 "Jarque-Bera test statistic  =",1X,F10.4,1X,F10.4/&
			 "----------------------------------------------------")',&
			 NINT(n),MINVAL(input%value),MINVAL(input%secondary),&
			 q1,median,q3,MAXVAL(input%value),MAXVAL(input%secondary),&
			 iqr,mu,var,sd,skew,ps,kurt,cvar,jb

END SUBROUTINE