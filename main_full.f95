program main_full

	implicit none

	! ****************************************
	! 				Parameters
	! ****************************************	
	integer:: ix_al, ix_bt, l, q, r, t, j = 1, G = 50, tf = 500 !, B = 1000
	integer, parameter:: M = 100, N = 1000, B = 1000		! Population matrix size
	integer, dimension(M, N):: P, Pn!, Pb	! Allocation(?)
	integer, dimension(N):: al, bt, tmp, off	! Individuals
	integer, dimension(M, M):: H, Hb 			! Hamming distances matrix
	integer, dimension(M, B):: Pb
	integer, allocatable:: C(:)					! Index of compatible individuals

	real:: r_n, mu = 0.0025 						! Random number, mutation rate
	real, dimension(N):: r_v					! Random vector
	real, dimension(M, N):: test_val			! Random matrix

	! Creates Population matrix for testing	
	! call random_number(test_val)
	! P = nint(test_val)
	P = 0

	! ****************************************
	!              Main Loop
	! ****************************************
	! Opens .dat file for population and B
	open(1, file = "pop.dat") 
	open(2, file = "vmat.dat") 

	do t = 1, tf
		Pb = P(:, 1:B)				! Matrix of mating segments
		Hb = hamming(Pb, M, B)   ! H. distances between mating segments

		j = 1
		do while (j <= M)
			! Random index for 1st individual
			call random_number(r_n)
			ix_al = 1+ floor(M*r_n)

			! 1st individual
			al = P(ix_al, :)
			
			! Creates list of compatible individuals ----------
			allocate(C(0))
			do q = 1, M
				if(Hb(ix_al, q) <= G)	C = [C, q]
			end do 
			
			! -------------------------------------------------
			
			if (size(C) /= 0) then
				
				call random_number(r_n)
				ix_bt = C(1 + floor(r_n*size(C))) ! Random index for beta from the compatible individuals list
				

				! 2nd individual
				bt = P(ix_bt, :)

				! Recombination + mutation
				call random_number(r_v)
				tmp = floor(r_v*2)

				! do l = 1, N
				! 	call random_number(r_n)
				!	off(l) = (tmp(l)*al(l)+ (1 - tmp(l))*bt(l) - merge(1, 0, r_n < mu))**2 
				!	! if (merge(1, 0, r_n < mu) == 1) print*, "j: ", j, "l: ", l, "Mutó"
				! end do
				 do l = 1, N
				 	call random_number(r_n)
					
				 	! off(l) = (tmp(l)*al(l)+ (1 - tmp(l))*bt(l) - merge(1, 0, r_n < mu))**2 
				 	off(l) = tmp(l)*al(l)+ (1 - tmp(l))*bt(l)
				 	if (r_n < mu)	off(l) = 1 - off(l)
				 	
				 end do 
				
				! Fills new population
				Pn(j, :) = off
				j = j + 1
			end if

			deallocate(C)
			
		end do ! End while
		P = Pn
	
	print*, t, "Out of ", tf
	end do ! End generations loop

	! Save last population
	do r = 1, M
		write(1, *) P(r, :)		
	end do

	write(2, *) B


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

end program main_full

