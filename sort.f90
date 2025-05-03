SUBROUTINE sort(input,output,flag)
	USE global
	IMPLICIT NONE
	
	!Dummy argument declaration
	INTEGER,INTENT(in) :: flag
	TYPE(point),DIMENSION(:),INTENT(in) :: input
	TYPE(point),DIMENSION(:),INTENT(out) :: output
	
	!Local variable declarations
	INTEGER :: i,n
	TYPE(point),DIMENSION(:),POINTER :: work
	TYPE(point) :: dummy
	
	!Subroutine to sort a vector of points using the Heap-sort algorithm. 
	!Data are sorted by value if flag =1, or by separation distance if flag = 0.
	!Actually, only a pointer is sorted, as this is more efficient. Subroutine
	!Adapted from Numerical Recipes in Fortran 90: Press, Teukolsky, Vetterling and
	!Flannery (1996) pg. 1171, ISBN: 0-201-54446-6 

	n=SIZE(input)
	ALLOCATE(work(n))
	work=input

	DO i=n/2,1,-1
		!Loop down the left range of the sift-down: The element to be sifted
		!is decremented to n/2 to 1 during "Hiring" in heap creation
		CALL sift_down(i,n,flag)
	END DO

	DO i=n,2,-1
		!Loop down the right range of the sift-down: Decremented from n-1 to
		!1 during the "Retirement and promotion" stage of heap creation
		
		!Clear space at the top of the array and retire the top of the heap into it
		dummy=work(1)
		work(1)=work(i)
		work(i)=dummy

		CALL sift_down(1,i-1,flag)
	END DO

	!Return sorted array and tidy up
	output=work
	DEALLOCATE(work)

CONTAINS

	SUBROUTINE sift_down(l,r,flag)
		IMPLICIT NONE
		!Dummy argument declaration
		INTEGER,INTENT(in) :: l,r,flag
	
		!Local variable declarations
		INTEGER :: j,old
		Type(point) :: a
	
		!Carry out sift-down on element array(l) to maintain heap structure
		
		!Get element to sift
		a=work(l)
		old=l
		j=l+l
		
		SELECT CASE(flag)

		CASE(0)
			!If flag = 0, sort by h
			!Do while j<=r
			DO
				IF(j>r) EXIT
				IF(j<r) THEN
					!Compare to the better underling
					IF(work(j)%h<work(j+1)%h) j=j+1
				END IF

				!If found a's level, terminate the sift-down, else demote and continue
				IF(a%h >= work(j)%h) EXIT
				work(old) = work(j)
				old=j
				j=j+j
			END DO
		
		CASE(1)
			!If flag = 1, sort by value
			!Do while j<=r
			DO
				IF(j>r) EXIT
				IF(j<r) THEN
					!Compare to the better underling
					IF(work(j)%value<work(j+1)%value) j=j+1
				END IF

				!If found a's level, terminate the sift-down, else demote and continue
				IF(a%value >= work(j)%value) EXIT
				work(old) = work(j)
				old=j
				j=j+j
			END DO
		
		END SELECT

		!Put into its slot
		work(old)=a

	END SUBROUTINE sift_down

END SUBROUTINE sort