program main_full
use utils

	implicit none

	! ****************************************
	! 				Parameters
	! ****************************************	
	integer:: ix_al, k, l, r, j = 0, G = 50, T = 500, B = 1000
	integer, parameter:: M = 100, N = 1000		! Population matrix size
	integer, dimension(M, N):: P, Pn, Pb, H			! Allocation(?)
	integer, dimension(N):: al, bt, tmp, off	! Individuals

	real:: r_n, mu = 0.001 						! Random number, mutation rate
	real, dimension(N):: r_v					! Random vector
	real, dimension(M, N):: test_val			! 

	! Creates Population matrix for testing	
	call random_number(test_val)
	P = nint(test_val)
	P = 0

	! ****************************************
	!              Main Loop
	! ****************************************
	open(1, file = "pop.dat") ! Opens .dat file 

	do t = 0, T

		Pb = P(:, 1:B)
		do while (j <= M)
			! Random index
			call random_number(r_n)
			ix_al = 1+ floor((M+ 1)*r_n)

			! 1st individual
			al = P(ix_al, :)
			
			do k = 1, M
				
				if (sum(abs(Pb(ix_al, :) - Pb(k, :))) <= G) then
					
					! 2nd individual
					bt = P(k, :)
					
					! Offspring + mutation
					call random_number(r_v)
					tmp = floor(r_v*2)

					do l = 1, N
						call random_number(r_n)
						
						off(l) = (tmp(l)*al(l)+ (1 - tmp(l))*bt(l) - merge(1, 0, r_n < mu))**2 

					end do 

					! Fills new population
					Pn(j, :) = off
					j = j + 1
				end if 
			end do
			print*, j 

		end do ! End while
		P = Pn
		
		!print*, t, "out of ", 500
	end do

	! H = hamming_tr(P, M, N)
	do r = 1, M
		write(1, *) P(r, :)
	end do

end program main_full

