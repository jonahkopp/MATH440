program main

  use omp_lib
  use mod_file

  implicit none

  !Declare variables and parameters
  integer, parameter :: my_kind = selected_real_kind(16,300)
  real(kind=my_kind), allocatable, dimension(:,:) :: K, K_inv, K_orig, E, E1
  integer, allocatable, dimension(:,:) :: M
  character,allocatable,dimension(:,:) :: C
  integer :: i,j,N
  real(kind=my_kind) :: timeseed,time_start,time_end
  character(len=20) :: file_name

  call cpu_time(timeseed)

  !Seed rng with the cpu time found above
  call srand(int(timeseed))

  !Prompt user for the file name containing the message to be encoded
  print *, "Enter the file name:"
  read(*,*) file_name

  !Set N to be the number of rows/cols of the key matrix, the number of rows of rows of the message matrix/encoded matrix
  N = 20

  !Generate the original message matrix using the following subroutine
  call gen_msg_mat(N,file_name,M)

  !Allocate matrices
  allocate(K(N,N),K_orig(N,N))
  allocate(E(size(M(:,1)),size(M(1,:))),E1(size(M(:,1)),size(M(1,:))))

  !Generate the key matrix using random numbers reals between 1 and 2
  do i=1,N
     do j=1,N
        K(i,j) = rand()+1
     end do
  end do

  !Calculate the encoded matrix E using mpi matrix portioning parallel algorithm
  E = matmul(K,M)

  K_orig = K

  !Sequential row reduction algorithm (uncomment to compare run time with omp version):
  !time_start = omp_get_wtime()
  !call row_red(K,K_inv)
  !time_end = omp_get_wtime()

  !print *,"seq time = ",time_end-time_start
  !print *,''

  !deallocate(K_inv)

  !The following line is necessary only when running sequential (above) AND parallel (below):
  !K = K_orig

  !Call the omp parallel row reduction subroutine to get K's inverse
  time_start = omp_get_wtime()
  call row_red_omp(K,K_inv)
  time_end = omp_get_wtime()
   
  print *,"parallel time = ",time_end-time_start
  print *,''

  !Compute the original message message again by multiplying the inverse key and the encoded message, E
  E1 = matmul(K_inv,E)

  !Convert the decoded matrix of reals into the nearest whole numbers
  M = nint(E1)

  !Convert the integer decoded matrix of integers into their ASCII equivalent characters
  C = char(M)

  print *, C(1,:)
  print *, C(2,:)
  print *, C(3,:)
  print *, C(4,:)
 
  !Deallocate all arrays used in main
  deallocate(K_inv)
  deallocate(K)
  deallocate(K_orig)
  deallocate(E)
  deallocate(M)

end program main
