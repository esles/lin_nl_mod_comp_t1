
NOT TO BE RUN
FOR REFERENCE ONLY

!=====================================================

!This inputs cpar ie get parameters, so start there.

	subroutine initial_data(cpar,InD,Ni,Rin,Din)
	implicit none
	integer I,Ni
	double precision InD(10), Din(Ni), Rin(Ni)
  double precision cpar(30)
  double precision zo,zz,zf,cto,ctf,ri,di

! InD(1) = initial time instant
! InD(2) = background's density
! InD(3) = background's expansion rate
! InD(4) = final time instant
! InD(5) = final background's density
! InD(6) = cosmological constant
! InD(7) = virialisation
! InD(8) = time step


! initial values: redshift, time instant, density, and expansion rate (the LCDM model assumed)
	zo = 1090.0d0
	zz = (zo+1.0d0)
	call timelcdm(zo,cto)
	InD(1) = cto
	InD(2) = cpar(5)*(zz**3)
	InD(3) = 3.0d0*cpar(2)*dsqrt(cpar(3)*(zz**3) + cpar(4))

! final time instants
	zf = 0.0
        call timelcdm(zf,ctf)
	InD(4) = ctf
	zz = (zf+1.0d0)
	InD(5) = cpar(5)*(zz**3)

! initial vector with density contrasts
!Generate a simple example of initial conditions.
!Modify this to read in a more realistic set of initial conditions, e.g. from the Millenium simulation initial conditions as in arXiv:1708.09143

do I=1,Ni

	   ri = I*1d0
	   call initial_profile(ri,di)
	   Rin(I) = ri
	   Din(I) = di

enddo

! other parameters
	InD(6) = cpar(7)
	InD(7) = cpar(10)
	InD(8) = cpar(11)

	end

!=====================================================
	subroutine initial_profile(r,d)
	implicit none
	double precision r,d,dc,sc

!call get initial conditions?

	dc = -0.002
	sc = 20

	d = dc*exp(-(r/sc)**2)

	end
!=====================================================

  subroutine get_parameters(cpar)
	implicit none
	double precision cpar(30)
	double precision omega_matter,omega_lambda,H0,age
	double precision pi, mu,lu,tu,gcons,cs,kap,kapc2,Ho,gkr,lb

	pi = 4d0*datan(1.0d0)

! cosmological parameters / Planck 2015 (TT+lowP+lensing)
	omega_matter = 0.308
	omega_lambda = 1.0 - omega_matter
	H0 = 67.810d0


	  if(omega_lambda.ne.(1.0 - omega_matter)) then
	   print *, 'The code uses the LCDM model to set up '
	   print *, ' the initial conditions '
	   print *, 'if you want to use non-LCDM models '
	   print *, 'then change the subroutine *timelcdm* '
  	   print *, '---calculations are being aborted---'
	   stop
	  endif

! units: time in 10^6 years, length in kpc, mass in 10^15 M_{\odot}
	mu=1.989d45
	lu=3.085678d19
	tu=31557600*1d6
! and other constants
	gcons= 6.6742d-11*((mu*(tu**2))/(lu**3))
	cs=299792458*(tu/lu)
	kap=8d0*pi*gcons*(1d0/(cs**4))
	kapc2=8d0*pi*gcons*(1d0/(cs**2))
	Ho=(tu/(lu))*H0
	gkr=3d0*(((Ho)**2)/(8d0*pi*gcons))
	lb=3d0*omega_lambda*(((Ho)**2)/(cs*cs))
	gkr=kapc2*gkr*omega_matter

	cpar(1) = H0*1d-2
	cpar(2) = Ho/cs
	cpar(3) = omega_matter
	cpar(4) = omega_lambda
	cpar(5) = gkr
	cpar(6) = cs
	cpar(7) = lb
	cpar(8) = pi
	cpar(9) = kapc2

! virialisation type: 1=turnaround, 2=near singularity, 3=stable halo
	cpar(10) = 1.0d0
! for option 2 ("collapsed"), please use either a fixed step (option 2 below), or please decrease the time step -- with default setting the results may not be accurate

! fix vs dynamical time step: 1=dynamical, 2=fixed
	cpar(11) = 1.d0
! dynamical step does not work well for some extreme cases
! so always test if this choice works well with your system


	end
!=====================================================
