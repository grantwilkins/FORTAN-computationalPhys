!Grant Wilkins
!This is the FORTRAN version of the second difference approximation for the heat
!Diffusion equation in one dimension.
	program heatApprox	
		implicit none

		real(8), dimension(:,:), allocatable :: heat
		real(8), dimension(:), allocatable :: sums
		real(8) :: tau, k, l, h, const, time, i, sumHeat
		integer(4) :: n, tsteps, x, y

		l = 1.0
		k = 1.0
		tau = 0.0001
		n = 61
		time = 0.1
		tsteps = int(time/tau)
		h = real(l/n)
		const = (k*tau)/((h)**2)
		i = 1.0
		
		!For the 1D case all data can be stored in this one data file
		!If using gnuplot splot the data with the option to plot as matrix.
		open(unit = 100, file ='heatData.dat')
		open(unit = 110, file = 'heatSum.dat')
		!Stores the heat in the appropriate dimensions.
		allocate(heat(tsteps, n), sums(tsteps))

		heat(:,1) = 0
		heat(:,n) = 0
		heat(1,:) = 0
		do x = 1, n
			heat(1, x) = abs(10.0*sin(i))
			i = i + 0.25
		end do
		do x = 2, tsteps
			do y = 2, n - 1
				heat(x, y) = heat(x - 1, y) + const*(heat(x -1, y +1 ) - 2.0*heat(x - 1, y) + heat(x - 1, y - 1))

			end do
		end do
		
		do x = 1, tsteps
			do y = 1, n
				sumHeat = sumHeat + heat(x,y)
			end do
			sums(x) = sumHeat
			sumHeat = 0
		end do

		do x = 1, tsteps
			write(110, *) x, sums(x)
			write(100,*) heat(x,:)
		end do
	
					
	
	end program heatApprox
