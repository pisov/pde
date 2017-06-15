program tdse
use utils
implicit none

double precision, parameter :: PI = 3.14159265359d0

double complex, allocatable, dimension(:,:) :: Psi2D, Pot2D

double complex :: Im, One, Zero

integer :: i, j, k, nx, ny, steps, cnt, info

double precision :: xmin, xmax, ymin, ymax, T0, x, y, Delta, t, dt, tout, sig, Factor

double precision :: x0, y0, norm

double precision :: Ekin, V0, kx, ky

double precision :: potWidth, potGateWidth

character (len=32)  :: filename, filenamePPM

! Set computational box dimensions
xmin = -25.d0
xmax =  25.d0
ymin = -15.d0
ymax =  15.d0
! Set potential energy barrier width -> potWidth
! and the gates width -> potGateWidth
potWidth = 1.d0
potGateWidth = 1.d0

! Wave packet initial conditions x0, y0
x0 = -10.d0
y0 =   0.d0
! Wave packet spread SIGMA
sig  = 0.5d0
! Complex constants
Im   = cmplx(0.d0, 1.d0)
One  = cmplx(1.d0, 0.d0)
Zero = cmplx(0.d0, 0.d0)

! Read initial data from standard input
write(0,*)'Please eneter grid resolution Delta = '
read(*,*)Delta


write(0,*)'Please enter time integration step (dt < ', norm  ,') dt = '
read(*,*)dt


write(0,*)'Please enter Kx  = '
read(*,*)kx

write(0,*)'Please enter Ky  = '
read(*,*)ky


write(0,*)'Enter potential barriaer height V0 = '
read(*,*)V0

write(0,*)'Enter integration period T = '
read(*,*)T0

write(0,*)'Enter write out period (Tout < ',T0,') Tout = '
read(*,*)tout

! Calculate image dimension
nx = nint((xmax - xmin) / Delta) - 2
ny = nint((ymax - ymin) / Delta) - 2

steps = tout / dt

!Calculate Factor multiplier

Factor = dt 0.5d0 / Delta ** 2

! Allocate arrays Psi2D and Pot2D
allocate(Psi2D(0:nx+1, 0:ny+1))
allocate(Pot2D(0:nx+1, 0:ny+1))

Pot2D(0:nx+1, 0:ny+1) = zero;

!Initialise wave function and potential energy
do j = 0, ny + 1
  y = ymin + j * Delta
  do i = 0, nx + 1
    x = xmin + i * Delta
    Psi2D(i, j) = (exp( -((x-x0)**2 + (y-y0)**2) / (4 * sig **2)  ) / ( dsqrt(dsqrt (2 * PI) * sig ))) &
                * cmplx(cos(kx*(x-x0)), sin(kx*(x-x0)))&
                * cmplx(cos(ky*(y-y0)), sin(ky*(y-y0)))
    if ((dabs(x).lt.(0.5d0*potGateWidth)).and.((dabs(y).lt.(1.5d0*potGateWidth)).or.(dabs(y).gt.(2.5d0*potGateWidth)))) then
      Pot2D(i, j) = cmplx(v0,   0.d0)
    else
      Pot2D(i, j) = cmplx(0.d0, 0.d0)
    end if
  end do
end do

! create image of pottential energy potEnergy.ppm
filenamePPM='potEnergy.ppm'
call ppmwrite(filenamePPM, cdabs(Pot2D(:,:)),nx+2, ny+2)

! Apply boundary conditions
Psi2D(0,:)   = Zero
Psi2D(nx+1,:)= Zero
Psi2D(:,0)   = Zero
Psi2D(:,ny+1)= Zero

! norm the wave function
norm = Delta*dsqrt(real(sum(Psi2D(:,:) * conjg(Psi2D(:,:)))))
Psi2D(:,:) = Psi2D(:,:) / norm

t = 0.d0
cnt = 0

do while ( t .le. T0)

  if (mod(cnt,steps).eq.0) then
    write(0,'(A6F15.5)')'Time: ',t
    write(filenamePPM,'(A3I0.6A4)')'wf-',cnt,'.ppm'
    call ppmwrite(filenamePPM, cdabs(Psi2D(:,:)),nx+2, ny+2)
  end if
! Calculate excplicit Shrodinger timedependent PDE

  t = t + dt
  cnt = cnt + 1
end do


deallocate(Psi2D)
deallocate(Pot2D)

end program tdse
