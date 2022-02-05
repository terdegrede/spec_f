program main_full
use utils

	implicit none

	! ****************************************
	! 				Parameters
	! ****************************************	
	integer:: ix_al, ix_bt, k, l, q, r, t, j = 0, G = 50, tf = 500, B = 1000
	integer, parameter:: M = 100, N = 1000		! Population matrix size
	integer, dimension(M, N):: P, Pn, Pb	! Allocation(?)
	integer, dimension(N):: al, bt, tmp, off	! Individuals
	integer, dimension(M, M):: H 				! Hamming distances matrix
	integer, allocatable:: C(:)					! Index of compatible individuals
	allocate(C(0))

	real:: r_n, mu = 0.1 						! Random number, mutation rate
	real, dimension(N):: r_v					! Random vector
	real, dimension(M, N):: test_val			! Random matrix

	! Creates Population matrix for testing	
	call random_number(test_val)
	P = nint(test_val)
	P = 0

	! ****************************************
	!              Main Loop
	! ****************************************
	open(1, file = "pop.dat") ! Opens .dat file 

	do t = 0, tf

		Pb = P(:, 1:B)				! Matrix of mating segments
		Hb = hamming_tr(Pb, M, B)   ! H. distances between mating segments

		do while (j <= M)

			! Random index for 1st individual
			call random_number(r_n)
			ix_al = 1+ floor((M+ 1)*r_n)

			! 1st individual
			al = P(ix_al, :)
			
			! Creates list of compatible individuals ----------
			do q = 1, M
				if(Hb(ix_al, q) <= G)	C = [C, q]
			end do 
			! -------------------------------------------------
			
			if (size(A) /= 0) then
				call random_number(r_n)
				ix_bt = C(1+ floor(r_n*size(C))) ! Random index for beta from the compatible individuals list

				! 2nd individual
				bt = P(ix_bt, :)

				! Offspring + mutation
				call random_number(r_v)
				tmp = floor(r_v*2)

				!print*, tmp
				!print*, off
				do l = 1, N
					call random_number(r_n)
					
					! off(l) = (tmp(l)*al(l)+ (1 - tmp(l))*bt(l) - merge(1, 0, r_n < mu))**2 
					off(l) = tmp(l)*al(l)+ (1 - tmp(l))*bt(l)
					if (r_n < mu)	off(l) = 1 - off(l)

				end do 
				!print*, off

			end if 
			! Fills new population
			Pn(j, :) = off
			j = j + 1
			
			! print*, j 

		end do ! End while
		P = Pn
	
	print*, t, "Out of ", tf
	end do

	! H = hamming_tr(P, M, N)
	do r = 1, M
		write(1, *) P(r, :)
	end do

end program main_full

