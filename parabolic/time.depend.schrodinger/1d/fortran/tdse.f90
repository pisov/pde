program tdse
implicit none

double precision, parameter :: PI = 3.14159265359d0

double complex, allocatable, dimension(:) :: Psi, PsiPrime

! Hamiltonian matix definition since is hertian we 
! defien the matrix with diagonal - d(n+1)
! and first neighbour diagonal - e(n)
! the right side of the equation (1 + i dt H / 2) * PsiNew = (1 - i dt H / 2) Psi 
! is vector b(n+1) = (1 - i dt H / 2) Psi


double complex, allocatable, dimension(:) :: D, DL, DU

double complex, allocatable, dimension(:) :: Pot


double complex :: Im, One, Zero

integer :: i, k, n, steps, cnt, info

double precision xmin, xmax, l0, T0, x, dx, t, dt, tout, sig, norm, xpmin, xpmax, Factor

double precision :: Ekin, Epot, V0, k0

double precision, allocatable, dimension(:) :: TrapCoef

character (len=32) :: filename

!Set domain size
xmin = -25.d0
xmax =  25.d0

!Set the potential energy barrier position left xpmin point
!and xpmax right point
xpmin = 4.d0
xpmax = 6.d0

!Wave packet initial position offset
l0   = -5.d0
!Wave packet spread
sig  = 0.5d0
!Define constants
Im   = cmplx(0.d0, 1.d0)
One  = cmplx(1.d0, 0.d0)
Zero = cmplx(0.d0, 0.d0)

write(0,*)'Please eneter grid resolution dx = '
read(*,*)dx

norm = dx **2 / ( 1 + v0 * 2 * dx **2)

write(0,*)'Please enter time integration step (dt < ', norm  ,') dt = '
read(*,*)dt

!if ( dt.gt.norm) then
!  write(0,*)'Time integration step dt is out of range !!!'
!  stop
!end if

!write(0,*)'Please enter kinetic energy 1/8 > Ekin  = '
write(0,*)'Please enter kinetic energy 1/2 > Ekin  = '
read(*,*)Ekin

if ( Ekin.lt.0.5d0) then
   write(0,*)'Kinetic energy should be > 0.5'
  stop
end if

!Calculate the impulse k0
!k0 = 0.5d0*dsqrt(8*Ekin-1.d0)
k0 = dsqrt(2*Ekin-1.d0)

write(0,*)'Enter potential well deep V0 = '
read(*,*)V0

write(0,*)'Enter integration period T = '
read(*,*)T0

write(0,*)'Enter write out period (Tout < ',T0,') Tout = '
read(*,*)tout

n = nint((xmax - xmin) / dx) - 1

steps = tout / dt

!Calculate Factor multiplier
Factor = 0.5d0 * dt / dx ** 2


allocate(Psi(0:n+1))
allocate(PsiPrime(1:n))
allocate(Pot(0:n+1))
allocate(D(1:n))
allocate(DL(1:n-1))
allocate(DU(1:n-1))
allocate(TrapCoef(1:n))

!Initialise coefficients for Trapezoidal rule
TrapCoef(2:n-1) = 2.d0
TrapCoef(1) = 1.d0
TrapCoef(n) = 1.d0

do i = 0, n + 1
  x = xmin + i * dx
!Initialise wave function
  Psi(i) = (exp( -(x-l0) ** 2 / (4 * sig **2)  ) / ( dsqrt(dsqrt (2 * PI) * sig ))) * cmplx(cos(k0*(x-l0)), sin(k0*(x-l0)))
!Initialise the potential energy
  if (x.ge.xpmin.and.x.le.xpmax) then
    Pot(i) = cmplx(v0,   0.d0)
  else
    Pot(i) = cmplx(0.d0, 0.d0)
  end if
end do

! Apply boundary conditions
Psi(0)   = Zero
Psi(n+1) = Zero

norm = dsqrt(dx*real(sum(Psi(:) * conjg(Psi(:)))))
Psi(:) = Psi(:) / norm

t = 0.d0
cnt = 0



do while ( t < T0)
  !Solve the right side of Crank Nicolson equation
  ! PsiRight(t) = (1 - 0.5 * H * dt) * Psi(t)
  Psi(1:n) = (One - Im * Factor - Im * 0.5d0 * dt * Pot(1:n)) * Psi(1:n) + 0.5d0 * Im * Factor * (Psi(0:n-1) + Psi(2:n+1))

  D(1:n)    = One +         Im * Factor + Im * 0.5d0 * dt * Pot(1:n) 
  DL(1:n-1) =     - 0.5d0 * Im * Factor
  DU(1:n-1) =     - 0.5d0 * Im * Factor
  !Solve linear system (1 + 0.5 * H * dt) * Psi(t+dt) = PsiRight(t)
  !Crank Nicolson equation left side 
  call zgtsv(n, 1, DL, D, DU, Psi(1:n), n, info)

  !Calculate first derivative of Psi
  PsiPrime(1:n) = 0.5d0*(Psi(2:n+1) - Psi(0:n-1))/dx

  !Calculate the kinetic energy
  Ekin = 0.25d0*dx*(sum(TrapCoef(1:n)*PsiPrime(1:n)*conjg(PsiPrime(1:n))))

  !Calculate the potential energy
  Epot = 0.5d0*dx*sum(TrapCoef(1:n)*conjg(Psi(1:n))*Pot(1:n)*Psi(1:n))

  !Calculate the wave function norm
  norm = 0.5d0*dx*sum(TrapCoef(1:n)*conjg(Psi(1:n))*Psi(1:n))

  !Write data file every 'steps' count
  if (mod(cnt,steps).eq.0) then
    write(*,'(5F15.7)')t,Ekin,Epot,Ekin+Epot,norm
    write(filename,'(A3I0.6A4)')'wf-',cnt,'.dat'
    open(unit=10, file=filename)
    write(10,'(A1F15.7)')'#',t
    do i = 0, n + 1
      x = xmin + i * dx
      write(10,'(7F15.7)')x,real(Psi(i)), aimag(Psi(i)), cdabs(Psi(i)),real(Pot(i)),Ekin,Epot
    end do
    close(10)
  end if

  t = t + dt
  cnt = cnt + 1
end do

deallocate(Psi)
deallocate(PsiPrime)
deallocate(Pot)
deallocate(D)
deallocate(DL)
deallocate(DU)
deallocate(TrapCoef)

end program tdse
