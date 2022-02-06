program main_full
use utils

	implicit none

	! ****************************************
	! 				Parameters
	! ****************************************	
	integer:: ix_al, ix_bt, k, l, q, r, t, j = 1, G = 50, tf = 500 !, B = 1000
	integer, parameter:: M = 100, N = 1000, B = 1000		! Population matrix size
	integer, dimension(M, N):: P, Pn!, Pb	! Allocation(?)
	integer, dimension(N):: al, bt, tmp, off	! Individuals
	integer, dimension(M, M):: H, Hb 			! Hamming distances matrix
	integer, dimension(M, B):: Pb
	integer, allocatable:: C(:)					! Index of compatible individuals

	real:: r_n, mu = 0.001 						! Random number, mutation rate
	real, dimension(N):: r_v					! Random vector
	real, dimension(M, N):: test_val			! Random matrix

	!	allocate(C(0))


	! Creates Population matrix for testing	
	! call random_number(test_val)
	! P = nint(test_val)
	P = 0

	! ****************************************
	!              Main Loop
	! ****************************************
	open(11, file = "pop.dat") ! Opens .dat file 
	!open(12, file = "hamm.dat")	! TEST FILE
	do t = 1, tf
		Pb = P(:, 1:B)				! Matrix of mating segments
		Hb = hamming_tr(Pb, M, B)   ! H. distances between mating segments

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
				
				! print*, t, size(off)
				! Fills new population
				Pn(j, :) = off
				j = j + 1
			end if

			deallocate(C)
			
		end do ! End while
		P = Pn
	
	! print*, t, "Out of ", tf
	end do ! End generations loop

	!H = hamming_tr(P, M, N)
	do r = 1, M
		write(11, *) P(r, :)
	!	write(12, *) H(r, :)
	end do

end program main_full

