

integer grid,gridz

real*4,allocatable ::   rho(:,:,:)

double precision ss,ss2,totval

character*80  filename

print*, ' Input densfield'

read (*,'(A)') filename

open(1,file=filename, form='unformatted',status='old')

read (1) grid,grid, gridz

allocate(rho(grid,grid,gridz))


 do kk=1,gridz
    read(1)  ((RHO(II,JJ,KK), ii=1,GRID),jj=1,GRID)
 enddo

close(1)


ss=0.

!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE(I,J,K) REDUCTION(+:ss)
do K=1,grid
   do J=1,grid
      do I=1,grid
         ss = ss + rho(i,j,k)
      enddo
   enddo
enddo

totval = float(grid)**3

ss = ss / totval

print*, ' Mean valuel=', ss

ss2=0.

!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE(I,J,K) REDUCTION(+:ss2)
do K=1,grid
   do J=1,grid
      do I=1,grid
         ss2 = (ss-rho(i,j,k))**2 + ss2
      enddo
   enddo
enddo

print*, ' Min-Max values ', Minval(rho), Maxval(rho)



ss2 = sqrt(ss2) / totval

print*, ' RMS=', ss2

stop
end



