program wave2d
implicit none

double precision, allocatable, dimension(:,:) :: u, unew, r, rnew, s, snew
double precision :: x, y, h, dt, t, tend, tout
double precision :: xmin, xmax, ymin, ymax, alpha, sigma
integer :: i, j, k, nsteps, np, step, stepout
character (len=32) :: filename

xmin = 0.d0
xmax = 1.d0
ymin = 0.d0
ymax = 1.d0

sigma = 0.01d0

write(0,*)'Please enter number of points of discretization N = '
read(*,*)np

h = (xmax - xmin) / (np + 1)

write(0,'(A32F15.7A4)')'Please enter time step dt < (',2*h,') = '
read(*,*)dt

alpha = dt/(2*h)

if (alpha.ge.1.d0) then
  write(0,*)'Wrong discretization parameters dt/(2*h) > 1 :',alpha
end if

write(0,*)'Please enter the total time period of integration Tend = '
read(*,*)Tend

write(0,*)'Please enter the writeout period Tout = '
read(*,*)Tout

stepout = nint(Tend / Tout)

!Allocate memory arrays
allocate(   u(0:np+1,0:np+1))
allocate(unew(0:np+1,0:np+1))
allocate(   r(0:np+1,0:np+1))
allocate(rnew(0:np+1,0:np+1))
allocate(   s(0:np+1,0:np+1))
allocate(snew(0:np+1,0:np+1))

!Set the initial state u(x,y,0)
do j = 1, np
  y = j * h + ymin
  do i = 1, np
    x = i * h + xmin
    u(i,j) = exp(-((x-0.5d0)**2+(y-0.5d0)**2)/(2*sigma))
  end do
end do

!Set zero boundary conditions
u(0,:) = 0.d0
u(np+1,:) = 0.d0
u(:,0) = 0.d0
u(:,np+1) = 0.d0

!Main iteration loop
t = 0.d0
step = 0
do while ( t.lt.tend )
  !Evolve solution with one time step
  ! evolve r
  rnew(1:np) = r(1:np) + alpha*0.5d0(s(2:np+1)-s(0:np-1)-alpha*(r(2:np+1)-2*r(1:np)+r(0:np-1)))
  ! evolve s
  snew(1:np) = s(1:np) + alpha*0.5d0(r(2:np+1)-r(0:np-1)-alpha*(s(2:np+1)-2*s(1:np)+s(0:np-1)))
  ! evolve u (leapfrog)
  

  unew(1:np,1:np) = u(1:np,1:np)
  if (mod(step,stepout).eq.0) then
    write(filename,'(A5I0.6A4)')'wave-',step,'.dat'
    open(unit=10, file=filename)
    write(10,'(A1F15.7)')'#',t
    do j = 0, np+1
      y = j * h + ymin
      do i = 0, np+1
        x = i * h + xmin
        write(10,'(7F15.7)')x,y,u(i,j)
      end do
    end do
    close(10)
  end if
  u(1:np,1:np) = unew(1:np,1:np)
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

end program wave2d
