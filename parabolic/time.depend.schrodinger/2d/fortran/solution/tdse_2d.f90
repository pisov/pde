module solver
implicit none
double complex, parameter  :: Im   = cmplx(0.d0, 1.d0)
double complex, parameter  :: One  = cmplx(1.d0, 0.d0)
double complex, parameter  :: Zero = cmplx(0.d0, 0.d0)

contains

! Subroutine jacobi solving left side of Crank - Nicolson equation
! using Jacobi iterative solver (1 + i dt H / 2) * PsiNew = Psi
! Psi2D - on input contain the right side of equaion
! as result returns PsiNew
! nx, ny - dimension 
subroutine jacobi(Psi2D, Pot2D, nx, ny, Factor, dt, eps)
implicit none
double complex, dimension(0:nx+1, 0:ny+1), intent(inout) :: Psi2D
double complex, dimension(0:nx+1, 0:ny+1), intent(in) :: Pot2D
integer, intent(in) :: nx, ny
double precision, intent(in) :: Factor
double precision, intent(in) :: dt
double precision, intent(in) :: eps

double complex, dimension(0:nx+1, 0:ny+1) :: PsiNew, Psi
double precision :: diff


Psi(:,:) = Psi2D(:,:)
diff = 2 * eps
do while (diff .gt. eps)
  PsiNew(1:nx, 1:ny) = (Psi2D(1:nx, 1:ny) + 0.5d0 * Im * Factor * (Psi(0:nx-1, 1:ny) + Psi(2:nx+1, 1:ny) + &
                       Psi(1:nx, 0:ny-1) + Psi(1:nx, 2:ny+1))) / &
                       (One + 2 * Im * Factor + Im * 0.5d0 * dt * Pot2D(1:nx, 1:ny))
  diff = cdabs(sum(PsiNew(1:nx, 1:ny) - Psi(1:nx,1:ny)))
  Psi(1:nx,1:ny) = PsiNew(1:nx, 1:ny)
end do

Psi2D(1:nx,1:ny) = Psi(1:nx,1:ny)

end subroutine jacobi

end module solver

program tdse
use solver
use utils
implicit none

double precision, parameter :: PI = 3.14159265359d0

double complex, allocatable, dimension(:)   :: Psi1D, Pot1D
double complex, allocatable, dimension(:,:) :: Psi2D, Pot2D

! Hamiltonian matix definition since is hertian we 
! defien the matrix with diagonal - d(n+1)
! and first neighbour diagonal - e(n)
! the right side of the equation (1 + i dt H / 2) * PsiNew = (1 - i dt H / 2) Psi 
! is vector b(n+1) = (1 - i dt H / 2) Psi

double complex, allocatable, dimension(:) :: D, DL, DU
integer :: i, j, k, nx, ny, steps, cnt, info
double precision :: xmin, xmax, ymin, ymax, T0, x, y, Delta, t, dt, tout, sig, Factor
double precision :: x0, y0, norm
double precision :: Ekin, V0, kx, ky
double precision :: potWidth, potGateWidth
character (len=32)  :: filename, filenamePPM

xmin = -25.d0
xmax =  25.d0
ymin = -15.d0
ymax =  15.d0
potWidth = 1.d0
potGateWidth = 1.d0

x0 = -10.d0
y0 =   0.d0
sig  = 0.5d0

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

nx = nint((xmax - xmin) / Delta) - 2
ny = nint((ymax - ymin) / Delta) - 2

steps = tout / dt

!Calculate Factor multiplier
Factor = 0.5d0 * dt / Delta ** 2

allocate(Psi2D(0:nx+1, 0:ny+1))
allocate(Pot2D(0:nx+1, 0:ny+1))

Pot2D(0:nx+1, 0:ny+1) = zero;

!Initialise wave function
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

filenamePPM='potEnergy.ppm'
call ppmwrite(filenamePPM, cdabs(Pot2D(:,:)),nx+2, ny+2)

! Apply boundary conditions
Psi2D(0,:)   = Zero
Psi2D(nx+1,:)= Zero
Psi2D(:,0)   = Zero
Psi2D(:,ny+1)= Zero

norm = Delta*dsqrt(real(sum(Psi2D(:,:) * conjg(Psi2D(:,:)))))
Psi2D(:,:) = Psi2D(:,:) / norm

t = 0.d0
cnt = 0

do while ( t .le. T0)
! Calculate right side of Crank - Nicolson equation
   Psi2D(1:nx,1:ny) = (One - 2 * Im * Factor - Im * 0.5d0 * dt * Pot2D(1:nx, 1:ny)) * Psi2D(1:nx, 1:ny) + &
                    0.5d0 * Im * Factor * (Psi2D(0:nx-1, 1:ny) + Psi2D(2:nx+1, 1:ny) + &
                    Psi2D(1:nx, 0:ny-1) + Psi2D(1:nx, 2:ny+1))
! Solve left side of Crank-Nicolson equation using Jacobi solver
  call jacobi(Psi2D, Pot2D, nx, ny, Factor, dt, 1.d-4)

  if (mod(cnt,steps).eq.0) then
    norm = Delta*dsqrt(real(sum(Psi2D(:,:) * conjg(Psi2D(:,:)))))
    write(0,'(A6F15.5A9F15.5)')'Time: ', t, 'Norm: ', norm
    write(filenamePPM,'(A3I0.6A4)')'wf-',cnt,'.ppm'
    call ppmwrite(filenamePPM, cdabs(Psi2D(:,:)),nx+2, ny+2)
  end if
  t = t + dt
  cnt = cnt + 1
end do

deallocate(Psi2D)
deallocate(Pot2D)

end program tdse
