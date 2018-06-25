	
	program gaussianUnitsAreForRealMen

		implicit none
		real(8), dimension(3) :: r, v, a, E, B, rtemp, vtemp, rtemp2, vtemp2, rtemp3, vtemp3, atemp, atemp2, atemp3
		real(8) :: m, q, tstep, t
		integer(4) :: i, n, x

		t = 0.00
		!Declaration of initial conditions
		!1st position is x, 2nd y, 3rd z
		m = 1.0
		q = 1.0
		n = 1000000
		tstep = 0.001
		E(1) = 0.01
		E(2) = 0.01
		E(3) = 0.1
		B(1) = 0.0
		B(2) = 0.01
		B(3) = 0.1
		r(1) = 0
		r(2) = 0
		r(3) = 0
		v(1) = 1.0
		v(2) = 0
		v(3) = 0
		
		!Open the files
		open(unit = 100, file = 'rElec.dat')
		open(unit = 110, file = 'vElec.dat')
		open(unit = 120, file = 'aElec.dat')

		do i = 1, n - 1
			call rk4(r,v)
               		write(100,*) r(1), r(2), r(3)
			write(110,*) v(1), v(2), v(3)
			write(120,*) a(1), a(2), a(3)

		end do
	contains
	function cross(a, b) result(c)
		real(8), dimension(3), intent(in) :: a, b
		real(8), dimension(3) :: c
		c(1) = a(2)*b(3) - a(3)*b(2)
		c(2) = -a(1)*b(3) + a(3)*b(1)
		c(3) = a(1)*b(2) - a(2)*b(1)
	end function cross

	subroutine rk4(r, v)
		real(8), dimension(3) :: r, v
		B = t*B
		a = (q/m)*(E + cross(v,B))
		rtemp = r + 0.5*tstep*v
		vtemp = v + 0.5*tstep*a
		
		atemp = (q/m)*(E+cross(v,B))
		rtemp2 = r + 0.5*tstep*vtemp
		vtemp2 = v + 0.5*tstep*atemp

		atemp2 = (q/m)*(E + cross(v,B))
		rtemp3 = r + tstep*vtemp2
		vtemp3 = (q/m)*(E + cross(v,B))

		r = r + tstep*(v + 2.0*vtemp + 2.0*vtemp2 + vtemp3)/6.0
		v = v + tstep*(a + 2.0*atemp + 2.0*atemp2 + atemp3)/6.0
		t = t + tstep

	end subroutine rk4
		
	end program gaussianUnitsAreForRealMen
