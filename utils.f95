module utils
	contains 

        function hamming_tr(P, M, N) result(H)
        ! Allows to find Hamming distances vector-wise
        ! Inputs:
        ! 		P: Matrix  
        ! 		M: Number of matrix rows
        !		N: Number of matrix columns
        !
        ! Outputs:
        !		H: Hamming matrix
        !

			! Parameters 
			! **********************************************************************
			integer:: k, l, l_0 = 1               ! Loop index rows/cols
			! integer, parameter:: M , N          ! Matrix dimensions
			integer, dimension(M, N), intent(in):: P       ! Declaring Input matrix
			integer, dimension(M, M):: H                   ! Declaring Output matrix

			! Initializing
			H = 0

			! Loop to fill matrix
			! **********************************************************************
			do k = 1, M       !rows
				do l = l_0, M ! cols

					! Hamming distances
					H(k, l) = sum(abs(P(k, :) - P(l, :)))
					H(l, k) = H(k, l) ! Full matrix (not only triangular)

					if (k == l) then
						H(k, l) = N + 1
					end if
				end do 
				l_0 = l_0 + 1

				print*, k, "out of ", M
			end do

			end function hamming_tr

end module utils