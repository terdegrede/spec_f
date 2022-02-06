module utils
	contains 

        function hamming(P, M, N) 
        implicit none
        ! Allows to find Hamming distances vector-wise
        ! Inputs:
        ! 		P: Matrix  
        ! 		M: Number of matrix rows
        !		N: Number of matrix columns
        !
        ! Outputs:
        !		hamming: Hamming matrix
        !

			! Parameters 
			! **********************************************************************
			integer:: k, l, l_0                  ! Loop index rows/cols
			integer, intent(in):: M , N          ! Matrix dimensions
			integer, dimension(M, N), intent(in):: P       ! Declaring Input matrix
			integer, dimension(M, M):: hamming ! Declaring Output matrix

			! Loop to fill matrix
			! **********************************************************************
			hamming = 0
			l_0 = 1
			do k = 1, M       !rows
				do l = l_0, M ! cols

					! Hamming distances
					hamming(k, l) = sum(abs(P(k, :) - P(l, :)))
					hamming(l, k) = hamming(k, l) ! Full matrix (not only triangular)

					if (k == l) then
						hamming(k, l) = N + 1
					end if
				end do 
				l_0 = l_0 + 1
			
				!print*, k, "out of ", M
			end do
			end function hamming

end module utils