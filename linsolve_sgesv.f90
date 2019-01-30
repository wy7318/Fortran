program linsolve_sgesv
implicit none
integer, parameter : : lda=50,ldb=lda,nrhs=1
real : : a(lda,lda),b(ldb,nrhs)
integer : : i,j,n,ipiv(lda),info
character(30) : : myformat
n = 8
myformat = ’(g22.8)’
a(1: n,1: n) = reshape((/((1.0/(i+j+1.0),i=0,7),j=0,7)/),(/8,8/))
b(1: n,1) = (/(17/7),(11/6),(10/7),(6/5),1,(8/9),(4/5),(5/7)/)
b(1: n,1) = matmul(a(1: n,1: n),b(1: n,1))
call sgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
if (info==0) then
print *, ’— LU is —————’
do i=1,n
print ’(50’//myformat//’)’,a(i,1: n)
end do
print *, ’— x is —————-’
print myformat, b(1: n,1)
print *, ’— ipiv is ————-’
print ’(a7,50(i8))’, ’ipiv = ’,ipiv(1: n-1)
else
print *, ’info = ’,info
end if
end program linsolve_sgesv
