program linsolve
use working precision
implicit none
real(wp) : : a(8,8),b(8),x(8),a save(8,8),d
integer : : i,j,p(8)
a = reshape((/((1.0/(i+j+1.0),i=0,7),j=0,7)/),(/8,8/))
b = matmul(a,(/(17/7),(11/6),(10/7),(6/5),1,(8/9),(4/5),(5/7)/))
print *, ’Matrix a’
print ’(8(g22.8))’, transpose(a)
a save = a
call mylup(a,p,d)
print *, ’determinant(A) = ’,d
1
2
print *, ’LU decomposition with pivoting’
print ’(8(g22.8))’, transpose(a)
print ’(”p =” 8i2)’, p
call mysolve(a,b(p),x)
print *, ’Solution x of system’
print ’(g22.8)’, x
contains
subroutine mysolve v1(a,b,x)
implicit none
integer : : i,j,n
real(wp), intent(in) : : a(: ,: ),b(: )
real(wp), intent(out) : : x(size(b))
n = size(a,1)
x = b
do i = 2,n
do j = 1,i-1
x(i) = x(i)-a(i,j)*x(j)
enddo
enddo
do i = n,1,-1
do j = i+1,n
x(i) = x(i)-a(i,j)*x(j)
enddo
x(i) = x(i)/a(i,i)
enddo
end subroutine mysolve v1
subroutine mysolve v2(a,b,x)
implicit none
integer : : i,n
real(wp), intent(in) : : a(: ,: ),b(: )
real(wp), intent(out) : : x(size(b))
n = size(a,1)
forall(i=1: n) x(i) = b(i)-dot product(a(i,1: i-1),x(1: i-1))
forall(i=n: 1: -1) x(i) = (x(i)-dot product(a(i,i+1: n),x(i+1: n)))/a(i,i)
end subroutine mysolve v2
subroutine mylu v1(a)
implicit none
integer : : i,j,k,n
real(wp),intent(inout) : : a(: ,: )
n = size(a,1)
do k = 1,n-1
do i = k+1,n
3
a(i,k) = a(i,k)/a(k,k)
do j = k+1,n
a(i,j) = a(i,j)-a(i,k)*a(k,j)
enddo
enddo
enddo
end subroutine mylu v1
pure subroutine mylu v2(a)
implicit none
integer : : i,j,k,n
real(wp),intent(inout) : : a(: ,: )
n = size(a,1)
do k = 1,n-1
forall(i=k+1: n) a(i,k) = a(i,k)/a(k,k)
forall(i=k+1: n, j=k+1: n) a(i,j) = a(i,j)-a(i,k)*a(k,j)
enddo
end subroutine mylu v2
subroutine mylu v3(a)
use basic, only: outer product
implicit none
integer : : k,n
real(wp), intent(inout) : : a(: ,: )
n = size(a,1)
do k = 1,n-1
a(k+1: n,k) = a(k+1: n,k)/a(k,k)
a(k+1: n,k+1: n) = a(k+1: n,k+1: n)-outer product(a(k+1: n,k),a(k,k+1: n))
enddo
end subroutine mylu v3
subroutine mylup(a,p,d)
use basic
implicit none
real(wp), intent(inout) : : a(: ,: )
real(wp), intent(out) : : d
integer, intent(out) : : p(size(a,1))
integer : : n,k,imax
n = size(a,1)
p = (/(k,k=1,n)/)
d = 1. wp
do k = 1,n-1
imax = k-1+maxloc(abs(a(k: n,k)),1)
if (k/=imax) then
call swap(a(k,: ),a(imax,: )) ! overloaded rswap
4
call swap(p(k),p(imax)) ! overloaded iswap
d = -d
end if
a(k+1: n,k) = a(k+1: n,k)/a(k,k)
a(k+1: n,k+1: n) = a(k+1: n,k+1: n)-outer product(a(k+1: n,k),a(k,k+1: n))
d = d*a(k,k)
end do
d = d*a(n,n)
end subroutine mylup
subroutine mysolve(a,b,x)
implicit none
integer : : i,n,j
real(wp), intent(in) : : a(: ,: ),b(: )
real(wp), intent(out) : : x(size(a,2))
n = size(a,1)
x = b
do j=1,n
do i=j+1,n
x(i) = x(i)-a(i,j)*x(j)
end do
end do
do j=n,1,-1
x(j) = x(j)/a(j,j)
do i=1,j-1
x(i) = x(i)-a(i,j)*x(j)
end do
end do
!fe: do i=1,n
! x(i) = b(i)-dot product(a(i,1: i-1),x(1: i-1))
!end do fe
!bs: do i=n,1,-1
! x(i) = (x(i)-dot product(a(i,i+1: n),x(i+1: n)))/a(i,i)
!end do bs
end subroutine mysolve
end program linsolve
