program poisson
implicit none

integer :: i, j, n, m, iter
double precision :: hx, hy, x, y, eps, xmin, xmax, ymin, ymax, diff

double precision, allocatable, dimension(:, :) :: uold, unew, source, uexact


write(0,*)'Please enter n = '
read(*,*)n

write(0,*)'Please enter m = '
read(*,*)m


allocate(uold(0:n+1, 0:m+1))
allocate(unew(0:n+1, 0:m+1))
allocate(source(0:n+1, 0:m+1))
allocate(uexact(0:n+1, 0:m+1))

xmin = 0.d0
xmax = 1.d0
ymin = 0.d0
ymax = 1.d0

hx = (xmax - xmin) / (n+1)
hy = (ymax - ymin) / (m+1)

do j = 1, m
  y = j * hy + ymin
  do i = 1, n
    x = i * hx + xmin
    source(i, j) = 2 * x **3 - 6 * x * y * (1.d0 - y)
    uexact(i, j) = y * (1.d0 - y) * x ** 3
  end do
  uold(n+1, j) = y * (1.d0 - y)
end do

uold(1:n, 1:m) = source(1:n, 1:m)
uold(0, :) = 0.d0
uold(:, 0) = 0.d0
uold(:, n+1) = 0.d0
unew(:,:) = 0.d0

eps = 1.d-5
diff = sum(abs(unew(1:n, 1:m) - uold(1:n, 1:m)))

do while (diff .gt. eps)
! Apply here the iteartive formula
! for Poisson equation
! using Jakobi relaxation method
  diff = sum(abs(unew(1:n, 1:m) - uold(1:n, 1:m)))
  uold(1:n, 1:m) = unew(1:n, 1:m)
end do

do j = 1, m
  y = j * hy + ymin
  do i = 1, n
    x = i * hx + xmin
    write(*,'(4E20.10)')x,y,uold(i,j),uexact(i,j)
  end do
end do

deallocate(uold)
deallocate(unew)
deallocate(source)
deallocate(uexact)

end program poisson

