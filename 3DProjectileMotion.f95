!Grant Wilkins
!6/24/18
!This program will simulate a projectile motion problem with air drag using the Verlet
!integration method. The user input will be specified based on knowledge of the given system.
program projectile
        implicit none
	real(8), dimension(1 : 3) :: r, v, a, vmid, aD, ag	
	real(8) :: t, tstep, theta, pi, g, coeffD, maxh, Area, m, rho, d, Cd
        
        open(unit = 100, file = 'xyz.dat')
        open(unit = 110, file = 'vxvyvz.dat')
        tstep = 0.01
        pi = 4.0*atan(1.0)
        g = 9.81
        !X is the 1st position, Y the 2nd and Z the 3rd for all vectors
        theta = 25.0*(pi/180.0)
        r(1) = 0
        r(2) = 0
        r(3) = 1.0075

        v(1) = 6.9668*cos(theta)
        v(2) = 0
        v(3) = 6.9668*sin(theta)

        ag(1) = 0
        ag(2) = 0
        ag(3) = -g

	!Initialization of the air drag information.
	!d = diameter, m = mass, rho is the density of the medium,
	!Cd is the constant of drag for the specified object.
        d = 0.02554
	Area = 0.25*pi*(d)**2.0

        m = 0.00774
	rho = 1.225
	Cd = 0.47
        coeffD = -0.5*d*Area*rho/m
        a = ag + coeffD*v*norm2(v)
        t = 0.0

	!The posistions will be calculated until the ball hits the ground.
        do while (r(3) > 0)
                vmid = v
                v = v +tstep*a
                call windV(r, v)               
		aD = coeffD*v*norm2(v)
                a = ag + aD
                if (r(3) > maxh) then
                        maxh = r(3)
                end if
                r = r + tstep*0.5*(vmid+v)
                t = t + tstep
                
                write(100,*)  r(1), r(2), r(3)
                write(110,*)  v(1), v(2), v(3)
        end do

	!Writes to the terminal the maximum height reached and the time spent in the flight.
        write(*,*) t
        write(*,*) maxh

contains

!Functions as a path-dependent vector field for the wind.
!If you do not want this enabled remove the call statement in the do-loop.
subroutine  windV(pos, vel)                
	real(8), dimension(3), intent(in) :: pos
        real(8), dimension(1 : 3) :: res, vel
        res(1) = log(pos(1)**2.0 + 1)
        res(2) = exp(-pos(2))
       	res(3) = exp(sin(pos(3)**2 + pos(1)))
	vel = vel - res
end subroutine windV

end program
