

program linear_void_model
!Initialisation
implicit none
integer :: i, j, n
real::OmegaM, OmegaL, H0, Hb, ai, zi, a0, dt, ax, tf, dx,xmax, h
real::rho_crit, rhobg_0, rhobg_i, t0, ti, ts, Hu, epsil,a1,a2,t1,t2,c,G,pi
real::k1,k2,k3,k4,p1,p2,p3,p4
real::rmax, rmin, rs, del0, sigma, Hi
integer :: sz, xn,tnum, rnum
real, dimension(:), allocatable ::tv              !v for vector
real, dimension(:), allocatable ::av
real, dimension(:), allocatable ::rv,del_iv,deld_iv,del_v,deld_v
!real,dimension(:,:),allocatable::del_rt         !matrix (rn x tn)
!real,dimension(:,:),allocatable::deld_rt         !matrix (rn x tn)

real,dimension(:,:),allocatable::dlt_rt,dltdt_rt
real,dimension(:,:),allocatable::ddlt_r
real,dimension(:),allocatable::ddlt_n,ddlt_np1,ddlt_tmp,ddlt_tmp_prime

real,Dimension(4)::rk4wt1, rk4wt2	!rk4 spacing and weighting
!real,dimension(:,:),allocatable::delvn,delvnp1,delvtmp,deldvtmp

real,Dimension(0:4)::kvec_1
real,Dimension(0:4,2)::kvec_2


real::a_tmp,t_tmp,a_n,t_n,a_np1,t_np1,a_tmp_prime
integer::rk4s
!real,dimension(:,:,:),allocatable::deltadrt       !matrix (rn x rn x tn)
!real,dimension(:,:,:),allocatable::deltadrt       !matrix (rn x rn x tn)
!real, dimension(1:2)::deltaw,deltaw1,deltaw2,q1,q2,q3,q4,l1,l2,l3,l4
logical::static_time_step=.false.

write(*,*) 'static_time_step = ',static_time_step

write(*,*) 'shape(kvec_1) = ',shape(kvec_1)
write(*,*) 'shape(kvec_2) = ',shape(kvec_2)


c = 2.97*(10**8)
G = 6.674/(10**11)
pi = 4.D0*DATAN(1.D0)

!testtesttest

write(*,*) 'pi value', pi
write(*,*) 'gravitational constant', G
write(*,*) 'light speed', c


!Initial Conditions

zi = 1000
ai = 1/(zi+1)
write(*,*) 'Initial scalefactor',ai

!
tnum = 30000
!if static then 50000 is ok
!if dynamic then need less (depending on epsil)


!=====

rnum = 250
!250

rmin = 0
rmax = 50

allocate(rv(1:rnum))

allocate(del_iv(1:rnum))
allocate(deld_iv(1:rnum))

allocate(del_v(1:rnum))
allocate(deld_v(1:rnum))

!

!allocate(delvn(2,1:rnum))
!allocate(delvnp1(2,1:rnum))

!allocate(del_rt(0:tnum,1:rnum))
!allocate(deld_rt(0:tnum,1:rnum))

!allocate(deltadrt(1:rnum,1:rnum,0:tnum))

!


allocate(dlt_rt(0:tnum,1:rnum))
allocate(dltdt_rt(0:tnum,1:rnum))

allocate(ddlt_r(2,1:rnum))

allocate(ddlt_n(2))
allocate(ddlt_np1(2))
allocate(ddlt_tmp(2))
allocate(ddlt_tmp_prime(2))

!


!
kvec_1(0) = 0

kvec_2(0,1) = 0
kvec_2(0,2) = 0
!kvec_2(0,:) = 0 ??

rs = real(rmax - rmin)/rnum

!====

del0 = 0.1
sigma = 10

!=====

epsil = 1.0/(10**3)
!epsil = 1.0/(10**4)

write(*,*) 'epsil',epsil

allocate(av(0:tnum))
allocate(tv(0:tnum))

write(*,*) 'shape of tv',shape(tv)
!
a0 = 1;

!Omega Factors
OmegaM = 0.27
!OmegaL = 0.3
OmegaL = 1 - OmegaM     !LCDM model constraint.
write(*,*) 'Dimensionless Matter density', OmegaM
write(*,*) 'Dimensionless DE density', OmegaL

!Hubble Parameter
H0 = 69                                     !Units: km/s/Mpc
Hb = (H0*1e3/(3.086e16*1e6))*31557600*1e9   !Units: 1/(Gyr)
write(*,*) 'Present day Hubble parameter in units 1/Gyr ', Hb

!The critical density (flat universe)
rho_crit = 3*(Hb**2)/(8*pi*G)
write(*,*) 'Critical Density', rho_crit
!Background density
rhobg_0 = rho_crit*OmegaM         !Present day value
rhobg_i = rhobg_0*((1+zi)**3)      !Initial value

call get_t_from_rhobg(rhobg_0,OmegaM,OmegaL,Hb,c,G,pi,t0)
!call get_t_from_z(0.0,OmegaM,OmegaL,H0,t0)
write(*,*) 'Age of the universe',t0

call get_t_from_rhobg(rhobg_i,OmegaM,OmegaL,Hb,c,G,pi,ti)
!call get_t_from_z(zi,OmegaM,OmegaL,H0,ti)
write(*,*) 'Time corresponding to initial redshift',ti

ts = ti - t0

if (static_time_step.eqv..true.) then
    dt = 0.0001
    !write(*,*) 'Static time step = TRUE w size = ', dt

    rk4wt1=(/0.,dt/2,dt/2,dt/)
    rk4wt2=(/dt/6.,dt/3.,dt/3.,dt/6./)	!
endif

tv(0) = ts
av(0) = ai

!====

!SUBROUTINE GET_INITIAL_DENSITY_PROFILE

!Initial profile for density perturbation

call get_Hubble(ti,ai,OmegaM,OmegaL,Hb,Hi)

do j = 0,rnum
  rv(j) = rmin + j*rs

  del_iv(j) = (-1)*del0*EXP((-1)*(rv(j)/sigma)**2)
  deld_iv(j) = Hi*del_iv(j)
enddo

write(*,*) 'rv', shape(rv)

open(1, file="data1.dat", status="new")
  do j=1,rnum
     write(1,*) rv(j), del_iv(j)
  end do
close(1)

do j =1,rnum
  dlt_rt(0,j) = del_iv(j)
  dltdt_rt(0,j) = deld_iv(j)
enddo


!============

!SUBROUTINE scale_factor ???

do n = 0,(tnum-1)   !time loop

!  write(*,*) 'debug', 0

  a_n = av(n)
  t_n = tv(n)

  !replace static_time_step with dt_opt
  if (static_time_step.eqv..false.) then
      call get_Hubble(t_n,a_n,OmegaM,OmegaL,Hb,Hu)
      dt = epsil*(1.0/Hu) !or
      !write(*,*) 'static time step = FALSE w size =', dt
      rk4wt1=(/0.0,dt/2.0,dt/2.0,dt/)
      rk4wt2=(/dt/6.0,dt/3.0,dt/3.0,dt/6.0/)	!
  endif

  !call get_Hubble(tn,an,OmegaM,OmegaL,Hb,Hu)

  !write(*,*) 'debug', 1

  do rk4s=1,4			              !call the 4 steps of the rk4 routine
    t_tmp = t_n + rk4wt1(rk4s)
    a_tmp = a_n + rk4wt1(rk4s)*kvec_1(rk4s-1)
    call get_aprime(t_tmp,a_tmp,OmegaM,OmegaL,Hb,a_tmp_prime)
    kvec_1(rk4s) = a_tmp_prime
    a_np1 = a_n + rk4wt2(rk4s)*kvec_1(rk4s)
  enddo                         !End RK4 loop

  t_np1 = t_n + dt

  av(n+1)=a_np1
  tv(n+1)=t_np1

enddo               !end time loop


!SUBROUTINE LINEAR_EVOLUTION

!av,tv,OmegaM,OmegaL,Hb,dlt_rt,ddlt_rt,tnum,rnum,static_time_step,epsil,
!Hu,

do j =1,rnum      !Start space loop

  !to copy kbs, the subroutine goes here

  do n = 0,(tnum-1)   !Start time loop

      t_n = tv(n)
      a_n = av(n)

      ddlt_r(1,j) = dlt_rt(n,j)
      ddlt_r(2,j) = dltdt_rt(n,j)

      ddlt_n(1) = ddlt_r(1,j)
      ddlt_n(2) = ddlt_r(2,j)


      if (static_time_step.eqv..false.) then
          call get_Hubble(t_n,a_n,OmegaM,OmegaL,Hb,Hu)
          dt = epsil*(1.0/Hu) !or
          !write(*,*) 'static time step = FALSE w size =', dt
          rk4wt1=(/0.0,dt/2.0,dt/2.0,dt/)
          rk4wt2=(/dt/6.0,dt/3.0,dt/3.0,dt/6.0/)
      endif

      do rk4s=1,4			            !call the 4 steps of the rk4 routine
        !write(*,*) 'debug', rk4s
        t_tmp = t_n + rk4wt1(rk4s)
        ddlt_tmp = ddlt_n + rk4wt1(rk4s)*kvec_2(rk4s-1,:)
        !HERE
        call get_deltaprime(t_tmp,ddlt_tmp,a_n,OmegaM,OmegaL,Hb,ddlt_tmp_prime)
        kvec_2(rk4s,:) = ddlt_tmp_prime
        ddlt_np1 = ddlt_n + rk4wt2(rk4s)*kvec_2(rk4s,:)
      enddo                       !end rk4 loop

      dlt_rt((n+1),j) = ddlt_np1(1)
      dltdt_rt((n+1),j) = ddlt_np1(2)

  enddo   !end time loop

  

enddo   !end space loop


write(*,*) 'shape of scale factor', SHAPE(av)
write(*,*) 'shape of dlt_rt', SHAPE(dlt_rt)


open(2, file="data2.dat", status="new")
  do n=0,(tnum-1)
     write(2,*) tv(n), av(n)
  end do
close(2)

open(3, file="data3.dat", status="new")
  do j=1,rnum
     write(3,*) rv(j), dlt_rt(1,j)
  end do
close(3)

open(4, file="data4.dat", status="new")
  do j=1,rnum
     write(4,*) rv(j), dlt_rt(int(tnum/2),j)
  end do
close(4)

open(5, file="data5.dat", status="new")
  do j=1,rnum
     write(5,*) rv(j), dlt_rt(tnum,j)
  end do
close(5)

!This, the scale factor, and the time and space domain vectors are crucial
open(6, file="data6.dat", status="replace")
  do j=1,rnum
     write(6,*)  (dlt_rt(n,j), n=0,(tnum-1))
  end do
close(6)


print*, 'end program'



end program linear_void_model

!=======================

subroutine get_t_from_rhobg(rhobg,OmegaM,OmegaL,H0,c,G,pi,t)
  !From the lambda-cdm model theory
  implicit none
  real,intent(in)::rhobg,OmegaM,OmegaL,H0
  real,intent(out)::t
  real::c,G,pi,lambda,xi,ash,td

!  write(*,*) 'pi value', pi
!  write(*,*) 'gravitational constant', G

  lambda = (3*(H0**2)*OmegaL)/(c**2)
  xi = sqrt((lambda/((8*pi*G/(c**2))*rhobg)))
  ash = log(xi + sqrt(xi*xi + 1))
  t = (1/c)*sqrt(4/(3*lambda))*ash

  return
end subroutine get_t_from_rhobg

!==========================

subroutine get_aprime(t,a,OmegaM,OmegaL,H0,ap)
  implicit none
  real,intent(in)::t,a,OmegaM,OmegaL,H0
  real,intent(out)::ap

  ap = H0*sqrt((OmegaM/a) + (OmegaL*(a**2)))

  return
end subroutine get_aprime

!=========================

subroutine get_Hubble(t,a,OmegaM,OmegaL,H0,Hu)
  implicit none
  real,intent(in)::t,a,OmegaM,OmegaL,H0
  real,intent(out)::Hu

  Hu = H0*sqrt((OmegaM/(a**3))+OmegaL);

  return
end subroutine get_Hubble

!==========================

subroutine get_deltaprime(t,deltaw,a,OmegaM,OmegaL,H0,deltawprime)
  implicit none
  real,intent(in)::t,a,OmegaM,OmegaL,H0
  real,intent(in),dimension(2)::deltaw        !vector
  real,intent(out),dimension(2)::deltawprime  !vector
  real::Hu,delta,deltadot,v,vdot

  delta = deltaw(1)
  deltadot = deltaw(2)

  !write(*,*) 'deltaw',deltaw

  call get_Hubble(t,a,OmegaM,OmegaL,H0,Hu)

  v = deltadot
  vdot = (3.0/2.0)*OmegaM*(H0**2.0)*(1.0/(a**3.0))*delta - (2.0*Hu)*v

  deltawprime(1) = v
  deltawprime(2) = vdot

  return
end subroutine get_deltaprime

!==============================
