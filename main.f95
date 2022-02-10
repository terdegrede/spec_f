program main
use utils
implicit none

	! *************************************************************************************
	! 							Declaration of variables
	! *************************************************************************************	
	integer:: ix_al, ix_bt								! Indexes to choose individuals
	integer:: l, q, r, t, j = 1							! Loop indexes
	integer:: G, tf = 500 						! Fixed parameters
	integer, parameter:: M = 500, N = 10000, B = 4000	! Population/mating matrix size parameters
	integer, dimension(M, N):: P, Pn					! Population matrix declaration
	integer, dimension(N):: al, bt, tmp, off			! Individuals, recombination and offspring
	integer, dimension(M, M):: H, Hb 					! Hamming distances matrix
	integer, dimension(M, B):: Pb						! Mating matrix
	integer, allocatable:: C(:)							! List of compatible individuals indexes

	real:: r_n, mu = 0.001 						! Random number, mutation rate
	real, dimension(N):: r_v					! Random vector
	real, dimension(M, N):: test_val			! Random matrix

    
	! Initialization
	P = 0
	G = 0.05*B 									! Mating distance (dB)
	print*, G
	! Opens .dat file for population and B
	open(1, file = "pop.dat") 
	open(2, file = "vmat.dat") 
	

	! *************************************************************************************
	! 									Main Loop
	! *************************************************************************************	

	do t = 1, tf

		Pb = P(:, 1:B)							! Matrix of mating segments
		Hb = hamming(Pb, M, B)   				! H. distances between mating segments

		j = 1									! Loop index refresh
		do while(j <= M)						! Start of mating loop
			
			! 1st individual
			call random_number(r_n)
			ix_al = 1+ floor(M*r_n)

			al = P(ix_al, :)

			! List of individuals compatible with 1st individual
			allocate(C(0))
			do q = 1, M
				if(Hb(ix_al, q) <= G)	C = [C, q]
			end do 

			if (size(C) /= 0) then				! Condition for mating
				
				! 2nd individual
				call random_number(r_n)
				ix_bt = C(1 + floor(r_n*size(C))) 

				bt = P(ix_bt, :)

				! Recombination + mutation
				call random_number(r_v)
				tmp = floor(r_v*2)

				do l = 1, N
					call random_number(r_n)
					off(l) = (tmp(l)*al(l)+ (1 - tmp(l))*bt(l) - merge(1, 0, r_n < mu))**2 
				end do

				! New population
				Pn(j, :) = off
				j = j + 1						! M-j individuals remaining in the new population
			end if 

			deallocate(C)
		end do 									! End of mating loop
		
		P = Pn 									! Redefining population for new generation
	
	! if (mod(t, 50) == 0) print*, t, "Out of ", tf
	print*, t, "Out of ", tf
	end do
	

	! Save last population
	
	write(1, *) P
	write(2, *) B, M, N
	

end program main