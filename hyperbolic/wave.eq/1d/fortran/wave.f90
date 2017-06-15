program wave
implicit none

double precision, allocatable, dimension(:) :: u, unew, r, rnew, s, snew
double precision :: x, y, h, dt, t, tend, tout, v, Q
double precision :: xmin, xmax, alpha, sigma
integer :: i, j, k, nsteps, np, step, stepout
character (len=32) :: filename

xmin = 0.d0
xmax = 1.d0
sigma = 0.05d0

write(0,*)'Please enter velocity v = '
read(*,*)v

write(0,*)'Please enter number of points of discretization N = '
read(*,*)np

h = (xmax - xmin) / (np + 1)

write(0,'(A32F15.7A4)')'Please enter time step dt < (',h/v,') = '
read(*,*)dt

alpha = dt*v/h
Q = (1.d0 - alpha) / (1.d0 + alpha)

if (alpha.ge.1.d0) then
  write(0,*)'Wrong discretization parameters dt*v/h > 1 :',alpha
end if

write(0,*)'Please enter the total time period of integration Tend = '
read(*,*)Tend

write(0,*)'Please enter the writeout period Tout = '
read(*,*)Tout

stepout = nint(Tout / dt)

!Allocate memory arrays
allocate(   u(0:np+1))
allocate(unew(0:np+1))
allocate(   r(0:np+1))
allocate(rnew(0:np+1))
allocate(   s(0:np+1))
allocate(snew(0:np+1))

!Set the initial state of r, s and u

do i = 0, np+1
  x = i * h + xmin
  u(i) = exp(-(x-0.5d0)**2/(2*sigma**2))
!  r(i) = v * (0.5d0-x)*u(i)/sigma**2 
end do

s(:)=0.d0
r(1:np) = v * (u(2:np+1) - u(0:np-1)) / (2 * h)

!Main iteration loop
t = 0.d0
step = 0
do while ( t.le.tend )
  !Evolve solution with one time step
  ! evolve r
  rnew(1:np) = r(1:np) + alpha*0.5d0*(s(2:np+1)-s(0:np-1)+alpha*(r(2:np+1)-2*r(1:np)+r(0:np-1)))
  ! evolve s
  snew(1:np) = s(1:np) + alpha*0.5d0*(r(2:np+1)-r(0:np-1)+alpha*(s(2:np+1)-2*s(1:np)+s(0:np-1)))
  ! evolve u (leapfrog)
  unew(0:np+1) = u(0:np+1) + 0.5d0*dt*(snew(0:np+1)+s(0:np+1))

  rnew(np+1) = r(np) - rnew(np) * Q + r(np+1) * Q 
  rnew(0)    = r(1)  - rnew(1)  * Q + r(0)    * Q 
  snew(np+1) = s(np) - snew(np) * Q + s(np+1) * Q 
  snew(0)    = s(1)  - snew(1)  * Q + s(0)    * Q 

  if (mod(step,stepout).eq.0) then
    write(filename,'(A5I0.6A4)')'wave-',step,'.dat'
    open(unit=10, file=filename)
    write(10,'(A1F15.7)')'#',t
    do i = 0, np+1
      x = i * h + xmin
      write(10,'(7F15.7)')x,u(i)
    end do
    close(10)
  end if

  r(:) = rnew(:)
  s(:) = snew(:)
  u(:) = unew(:)
  
  t = t + dt
  step = step + 1
end do

!Deallocate memory arrays
deallocate(u)
deallocate(unew)
deallocate(r)
deallocate(rnew)
deallocate(s)
deallocate(snew)

end program wave
