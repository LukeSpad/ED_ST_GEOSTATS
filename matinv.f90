SUBROUTINE matinv(input,output)
	USE global
	USE interfaces
	IMPLICIT NONE

	!Dummy variable declarations
	REAL,DIMENSION(:,:),INTENT(in) :: input
	REAL,DIMENSION(:,:),INTENT(out) :: output

	!Local variable declarations
	INTEGER :: i,row,n,column
	INTEGER,DIMENSION(SIZE(input,1),SIZE(input,2)) :: identity
	
	REAL(8) :: a,b,pivot
	REAL(8),DIMENSION(SIZE(input,1)*2,SIZE(input,2)) :: work
	REAL(8),DIMENSION(SIZE(input,1)*2) :: tmp

	!Matinv V1.1: A generic subroutine to invert a matrix via the Gauss-Jordan method. 
	!Copyright (C) 2006  Luke Spadavecchia
	!This program is free software; you can redistribute it and/or modify
    !it under the terms of the GNU General Public License as published by
    !the Free Software Foundation; either version 2 of the License, or
    !(at your option) any later version.
    !This program is distributed in the hope that it will be useful,
    !but WITHOUT ANY WARRANTY; without even the implied warranty of
    !MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !GNU General Public License for more details.
    !You should have received a copy of the GNU General Public License along
    !with this program; if not, write to the Free Software Foundation, Inc.,
	!51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

	!A generic subroutine to invert a square 2D matrix by pivoting ("Gauss-Jordan" method).
	!An n*2,n work array is assembled, with the array of interest on the left half, and the
	!identity matrix on the right half. A check is made to ensure the matrix is invertable,
	!and the matrix inverse is returned, providing the matrix is not singular. The matrix
	!is invertable if, after pivoting and row reduction, the identity matrix shifts to the
	!left half of the work array. Initially the code was written in integer form to avoid
	!fractions in the work array, but numbers rapidly become inconcevably large. The 
	!implemented solution therefore works on fractional pivot rows, with pivots factored
	!to leading ones. Double precision arrays have been used to maintain numerical stability.

	!Get size of input
	n=SIZE(input,1)

	!Build identity matrix
	identity=0
	
	!Make diagonals = 1
	DO i=1,n
		identity(i,i)=1
	END DO

	!Build work array:	
	!Input array on the left half. 
	work(1:n,:)=input

	!Identity matrix on the right half
	work(n+1:n*2,:) = identity

	!Step through the work array rows
	DO row=1,n
		!Step through row j's columns
		DO column=1,n
			!Find first non-zero element of row and record it as the pivot column
			IF (work(column,row) .NE. 0) THEN				
				pivot = work(column,row)
				EXIT
			END IF
		END DO

		!Divide pivot row by pivot to rescale for leading 1
		work(:,row)=work(:,row)/pivot

		!Loop over rows to 'Clear' the pivot column using the pivot value
		DO i=1,n
			!Skip over if i is the pivot row
			IF (i == row) CYCLE

			!Get clearance row scaling factor
			b = work(column,i)
				
			!Replace row with the following:
			!Difference between clearance row and (b * pivot row)
			work(:,i) = work(:,i) - b*work(:,row)
		END DO 
	END DO

	!Validate that the inversion is successful, by checking left side of work array
	!is equal to the identity matrix.
	IF (ANY(work(1:n,:) .NE. identity)) THEN
		PRINT '("Inversion not possible: Matrix is degenerate!"/)'
!		PAUSE
		STOP
	END IF

	!Output the inverse of input array (right hand side of work array).
	output = work(n+1:2*n,:) 

END SUBROUTINE matinv