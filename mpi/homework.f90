module Homework
implicit none
contains
subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
implicit none
include "mpif.h"
real(8), intent(in), dimension(:,:) :: A
integer(4), intent(out) :: x1,y1,x2,y2
integer(4) :: mpiErr, mpiSize, mpiRank
integer(4) :: n,m,l,i,j,k,c,b,num1,num2, coord(4)
real(8) m1,m2
real(8),allocatable, dimension(:) :: p, mm1
integer, allocatable, dimension(:) :: mx1, mx2, my1, my2
 n=size(A,1)
 m=size(A,2)
 allocate(p(n))
 m1=A(1,1)
 x1=1; x2=1; y1=1; y2=1
call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpiErr)
call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpiErr)
 i=mpiSize
 num1=i*i
 k=m/num1
 if (mpiRank+1.ne.mpiSize) then
  c=1+k*mpiRank*mpiRank
  b=(mpiRank+1)*(mpiRank+1)*k
 else 
  c=1+k*mpiRank*mpiRank
  b=m
 endif
 do i=c, b
  do j=1,n
   p(j)=0
  enddo
  do k=i, m 
   m2=0; l=1;
   do j=1,n
    p(j)=p(j)+a(j,k)
    m2=m2+p(j)
    if (m2>m1) then 
     m1=m2; x1=i; y1=l; x2=k; y2=j 
    endif
    if (m2<0) then 
      m2=0; l=j+1
    endif
   enddo
  enddo
 enddo
 deallocate(p)
 num1=mpiSize
 if (mpiRank==0) then
  allocate(mm1(num1)); allocate(mx1(num1)); allocate(mx2(num1)); allocate(my1(num1)); allocate(my2(num1))
  mm1=0
 endif
 call mpi_barrier(MPI_COMM_WORLD, mpiErr)
 call mpi_gather(m1, 1, MPI_REAL8, mm1, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
 call mpi_gather(x1, 1, MPI_INTEGER4, mx1, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
 call mpi_gather(x2, 1, MPI_INTEGER4, mx2, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
 call mpi_gather(y1, 1, MPI_INTEGER4, my1, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
 call mpi_gather(y2, 1, MPI_INTEGER4, my2, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
 call mpi_barrier(MPI_COMM_WORLD, mpiErr)  
 if (mpiRank==0) then
   m2=mm1(1)
   j=1
   do i=1,num1
    if(mm1(i)>m2) then
     m2=mm1(i); j=i
    endif
   enddo
   coord(1)=mx1(j); coord(2)=mx2(j); coord(3)=my1(j); coord(4)=my2(j)
   x1=coord(1); x2=coord(2); y1=coord(3); y2=coord(4) 
 endif
 call mpi_barrier(MPI_COMM_WORLD, mpiErr)
 call mpi_bcast(coord, 4, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
 x1=coord(1); x2=coord(2); y1=coord(3); y2=coord(4)  
 if (mpiRank==0) then
  deallocate(mm1); deallocate(mx1); deallocate(mx2); deallocate(my1); deallocate(my2)
 endif
end subroutine 
end module
  
