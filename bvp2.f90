program bvp2
implicit none
integer, parameter : : N = 20-1
real, parameter : : pi = 3.141592653589793
real, parameter : : c = 10., k = 900., a = 0., b = 0., omega = 10*pi
integer, parameter : : itmax = 200
real, parameter : : tol = 1e-6
real : : u(N),y(N),dx,e(0: itmax),x(N+2),xi(N),f(N)
type sparse
integer, allocatable, dimension(: ) : : col ind,row ind
real , allocatable, dimension(: ) : : val
end type
type(sparse) : : S
integer : : M,i,j(3)
dx = 1.0/(N+1)
x = (/(i*dx,i=0,N+1)/)
xi = x(2: N+1)
f = rhs(xi,omega)
M = 3*N-2 ! storage for tridiagonal matrix
call allocate sparse(M,S)
S%val((/1,2/)) = (/-2,1/)/dx**2+c*(/0,1/)/(2*dx)+k*(/1,0/)
S%row ind((/1,2/)) = (/1,1/)
S%col ind((/1,2/)) = (/1,2/)
f(1) = f(1)-a*(1/dx**2-c/(2*dx))
do i=2,N-1
j = 3*i-(/3,2,1/)
S%val(j) = (/1,-2,1/)/dx**2+c*(/-1,0,1/)/(2*dx)+k*(/0,1,0/)
S%row ind(j) = (/i,i,i/)
S%col ind(j) = (/i-1,i,i+1/)
end do
S%val(3*N-(/3,2/)) = (/1,-2/)/dx**2+c*(/-1,0/)/(2*dx)+k*(/0,1/)
10
S%row ind(3*N-(/3,2/)) = (/N,N/)
S%col ind(3*N-(/3,2/)) = (/N-1,N/)
f(N) = f(N)-b*(1/dx**2+c/(2*dx))
call print system(S,f)
u = 0.0
call gs(S,f,u,tol,itmax,e)
print *, ’solution u of au = b’
print ’(f12.7)’, u
print *, ’residual b-a*u’
print ’(f12.7)’, f-ax(S,u)
open(unit=1,file=”bvp sol f90.dat”,action=”write”,status=”replace”)
write(1,*) 0.0,a ! BC at x=0
do i=1,size(u,1)
write(1,*) xi(i),u(i)
end do
write(1,*) 1.0,b ! BC at x=1
close(1)
open(unit=1,file=”bvp gs cv.dat”,action=”write”,status=”replace”)
do i=1,size(pack(e,e.ne.huge(1.0)),1)
write(1,*) e(i)
end do
close(1)
call deallocate sparse(S)
contains
subroutine allocate sparse(m,a)
implicit none
integer, intent(in) : : m
type(sparse), intent(inout) : : a
allocate(a%val(m))
allocate(a%row ind(m))
allocate(a%col ind(m))
end subroutine allocate sparse
subroutine deallocate sparse(a)
implicit none
type(sparse), intent(inout) : : a
deallocate(a%row ind)
deallocate(a%col ind)
deallocate(a%val)
11
end subroutine deallocate sparse
subroutine print system(a,b)
implicit none
type(sparse), intent(in) : : a
real, intent(in) : : b(: )
print *, ’Matrix a’
print *, ’ i j a(i,j)’
print *, ’——————-’
do i=1,size(a%val,1)
print ’(2(i4),f12.3)’,a%row ind(i),a%col ind(i),a%val(i)
end do
print *
print *, ’rhs b’
print ’(f12.7)’, b
print *
end subroutine print system
real elemental function rhs(x,omega)
implicit none
real, intent(in) : : x, omega
rhs = sin(omega*x)
end function rhs
function ax(a,x) result(y)
implicit none
type(sparse), intent(in) : : a
real, intent(in) : : x(: )
real : : y(size(x,1))
integer : : i,j,k,m,n
n = size(x,1)
m = size(a%val,1)
y = 0
do k=1,m
i = a%row ind(k); j = a%col ind(k)
y(i) = y(i)+a%val(k)*x(j)
end do
end function ax
subroutine gs(a,b,x,tol,itmax,e)
implicit none
type(sparse), intent(in) : : a
real, intent(in) : : b(: )
12
real, intent(inout) : : x(: )
real, intent(in) : : tol
real, intent(out) : : e(0: itmax)
integer, intent(in) : : itmax
integer : : i,j,k,n,nn
real : : r(size(b,1)), x old(size(x,1))
n = 0
nn = size(a%val,1)
e = huge(1.0)
do while (e(n)>tol .and. n<itmax)
r = b
x old = x
do k=1,nn
i = a%row ind(k); j = a%col ind(k) ! a(k) = A(i,j)
if (i .lt. j) r(i) = r(i)-a%val(k)*x old(j)
end do
do k=1,nn
i = a%row ind(k); j = a%col ind(k)
if (i .gt. j) then
r(i) = r(i)-a%val(k)*x(j)
elseif (i .eq. j) then
x(i) = r(i)/a%val(k)
end if
end do
n = n+1
e(n) = norm(x-x old)/norm(x)
end do
if (n==itmax .and. e(n)>tol) print *, ’no convergence in gs’
end subroutine gs
real function norm(x)
implicit none
real, intent(in) : : x(: )
norm = sqrt(dot product(x,x))
end function norm
end program bvp2
