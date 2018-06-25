        program orbitalMotion
                implicit none
                
                real(8), dimension(3) :: r, v, a, vtemp, vtemp2, vtemp3, rtemp, rtemp2, rtemp3, atemp, atemp2, atemp3               
                real(8) :: tstep, gm, t, mass, KE, UE, pi
                integer(4) :: i, n
                
                !Initialization of constants
                pi = 4.0*atan(1.0)
                gm = 4.0*pi**2.0
                t = 0.0
                tstep = 0.0001
                n = 50000
                mass = 1.0

                !File opening for the .dat plots
                open(unit = 100, file = 'rOrb.dat')
                open(unit = 110, file = 'vOrb.dat')
                open(unit = 120, file = 'KtOrb.dat')
                open(unit = 130, file = 'UtOrb.dat')
                open(unit = 140, file = 'EtOrb.dat')

                r(1) = 1.0
                r(2) = 0.0
                r(3) = 0.0

                v(1) = pi
                v(2) = 2.0*pi
                v(3) = 0.0
               
	        !runs the Runge-Kutta 4 method for the position of the planet
		!during this time period.
                do i = 1, n
                        !K2
                        a = (-gm/( norm2(r)**3.0))*r
                        rtemp = r + 0.5*tstep*v
                        vtemp = v + 0.5*tstep*a
                        atemp = (-gm/ (norm2(rtemp)**3.0))*rtemp
                        
                        !K3
                        rtemp2 = r + 0.5*tstep*vtemp
                        vtemp2 = v + 0.5*tstep*atemp
                        atemp2 = (-gm/(norm2(rtemp2)**3.0))*rtemp2

                        !K4
                        rtemp3 = r + tstep*vtemp2
                        vtemp3 = v + tstep*atemp2
                        atemp3 = (-gm/(norm2(rtemp3)**3.0))*rtemp3

                        !The updating of the vectors using the RK4
                        !method 1/6*(k1+k2+k3+k4)*tstep
                        r = r + tstep*(v + 2.0*vtemp + 2.0*vtemp2 + vtemp3)/6.0
                        v = v + tstep*(a + 2.0*atemp + 2.0*atemp2 + atemp3)/6.0
                        t = tstep + t
                        
                        KE = 0.5*mass*norm2(v)**2.0
                        UE = -gm/(norm2(r))

                        !Writing to datafiles
                        write(100,*) r(1), r(2), r(3)
                        write(110, *) v(1), v(2), v(3)
                        write(120,*) t, KE                     
                        write(130, *) t, UE
                        write(140,*) t, KE + UE

                end do
              
        end program orbitalMotion
