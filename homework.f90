 module Homework
 implicit none
 contains
 subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
 implicit none
 include "mpif.h"
 real(8), intent(in), dimension(:,:) :: A
 integer(4), intent(out) :: x1,y1,x2,y2
 integer(4) :: mpiErr, mpiSize, mpiRank
 integer(4) :: n,m,l,i,j,k,left_border,right_border, coord(4)
 real(8) candidate_for_max_summ,max_summ
 real(8),allocatable, dimension(:) :: array_of_continious_summ, array_of_candidate_for_max_summ
 integer, allocatable, dimension(:) :: array_of_candidate_for_x1, array_of_candidate_for_x2, array_of_candidate_for_y1, array_of_candidate_for_y2
   n=size(A,1)
   m=size(A,2)
   allocate(array_of_continious_summ(n))
   max_summ=A(1,1)
   x1=1; x2=1; y1=1; y2=1
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpiErr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpiErr)
   k=m/(mpiSize*mpiSize)
   if (mpiRank+1.ne.mpiSize) then
     left_border=1+k*mpiRank*mpiRank
     right_border=(mpiRank+1)*(mpiRank+1)*k
   else 
     left_border=1+k*mpiRank*mpiRank
     right_border=m
   endif
   do i=left_border, right_border
     do j=1,n
       array_of_continious_summ(j)=0
     enddo
     do k=i, m 
      candidate_for_max_summ=0; l=1;
      do j=1,n
        array_of_continious_summ(j)=array_of_continious_summ(j)+a(j,k)
        candidate_for_max_summ=candidate_for_max_summ+array_of_continious_summ(j)
        if (candidate_for_max_summ>max_summ) then 
         max_summ=candidate_for_max_summ; x1=i; y1=l; x2=k; y2=j 
        endif
        if (candidate_for_max_summ<0) then 
          candidate_for_max_summ=0; l=j+1
        endif
      enddo
     enddo
   enddo
   deallocate(array_of_continious_summ)
   if (mpiRank==0) then
     allocate(array_of_candidate_for_max_summ(mpiSize)); 
     allocate(array_of_candidate_for_x1(mpiSize))
     allocate(array_of_candidate_for_x2(mpiSize))
     allocate(array_of_candidate_for_y1(mpiSize))
     allocate(array_of_candidate_for_y2(mpiSize))
     array_of_candidate_for_max_summ=0
   endif
   call mpi_gather(max_summ, 1, MPI_REAL8, array_of_candidate_for_max_summ, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
   call mpi_gather(x1, 1, MPI_INTEGER4, array_of_candidate_for_x1, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
   call mpi_gather(x2, 1, MPI_INTEGER4, array_of_candidate_for_x2, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
   call mpi_gather(y1, 1, MPI_INTEGER4, array_of_candidate_for_y1, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
   call mpi_gather(y2, 1, MPI_INTEGER4, array_of_candidate_for_y2, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)  
   if (mpiRank==0) then
     candidate_for_max_summ=array_of_candidate_for_max_summ(1)
     j=1
     do i=1,mpiSize
       if(array_of_candidate_for_max_summ(i)>candidate_for_max_summ) then
        candidate_for_max_summ=array_of_candidate_for_max_summ(i); j=i
       endif
     enddo
     coord(1)=array_of_candidate_for_x1(j)
     coord(2)=array_of_candidate_for_x2(j)
     coord(3)=array_of_candidate_for_y1(j)
     coord(4)=array_of_candidate_for_y2(j)
     x1=coord(1); x2=coord(2); y1=coord(3); y2=coord(4) 
   endif
   call mpi_bcast(coord, 4, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
   x1=coord(1); x2=coord(2); y1=coord(3); y2=coord(4)  
   if (mpiRank==0) then
     deallocate(array_of_candidate_for_max_summ)
     deallocate(array_of_candidate_for_x1)
     deallocate(array_of_candidate_for_x2)
     deallocate(array_of_candidate_for_y1)
     deallocate(array_of_candidate_for_y2)
   endif
 end subroutine 
 end module
   
