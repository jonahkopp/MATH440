program main

  use mpi
  use omp_lib
  use mod_file

  implicit none

  !Declare variables and parameters
  integer, parameter :: my_kind = selected_real_kind(16,300)
  real(kind=my_kind), allocatable, dimension(:,:) :: K, K_inv, K_orig, K_core, K_core_inv, E, E1, M_core
  integer, allocatable, dimension(:,:) :: M
  character,allocatable,dimension(:,:) :: C
  integer,dimension(mpi_status_size) :: mpi_status
  integer :: i,j,N,ierror,my_rank,num_cores,master,div,rem,tag
  real(kind=my_kind) :: timeseed,time_start,time_end
  character(len=20) :: file_name
  integer,allocatable,dimension(:) :: num_rows,start_vec
  
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,my_rank,ierror)
  call mpi_comm_size(mpi_comm_world,num_cores,ierror)

  master = 0
  tag = 0
  
  if(my_rank == master) then
  
  call cpu_time(timeseed)

  !Seed rng with the cpu time found above
  call srand(int(timeseed))

  !Prompt user for the file name containing the message to be encoded
  print *, "Enter the file name:"
  
  read(*,*) file_name

  end if
   
  call mpi_bcast(file_name,len(file_name),mpi_double,master,mpi_comm_world,ierror)
  
  !Set N to be the number of rows/cols of the key matrix, the number of rows of rows of the message matrix/encoded matrix
  N = 1000
  
  !Generate the original message matrix using the following subroutine
  call gen_msg_mat(N,file_name,M)

  if (my_rank == master) then
  
  !Allocate matrices
  allocate(K(N,N),K_orig(N,N))
  !allocate(E(size(M(:,1)),size(M(1,:))),E1(size(M(:,1)),size(M(1,:))))

  !Generate the key matrix using random numbers reals between 1 and 2
  do i=1,N
     do j=1,N
        K(i,j) = rand()+1
     end do
  end do

  allocate(num_rows(num_cores),start_vec(num_cores))
  
  div = N/num_cores
  rem = mod(N,num_cores)

  do i = 1,num_cores
     num_rows(i) = div
  end do
  
  do i = 1,rem
     num_rows(i) = num_rows(i) + 1
  end do

  start_vec(1) = 1
  
  do i = 2,num_cores
     start_vec(i) = num_rows(i-1)+start_vec(i-1)
  end do
  
  
  do i = 2,num_cores
     call mpi_send(K(start_vec(i):(start_vec(i)+num_rows(i)),:),num_rows(i)*N,mpi_double,i-1,tag,mpi_comm_world,ierror)
  end do
  
  end if

  if (my_rank .ne. master) then
     allocate(K(num_rows(my_rank+1),N))
     call mpi_recv(K,num_rows(my_rank+1)*N,mpi_double,master,tag,mpi_comm_world,mpi_status,ierror)
  end if

  if (my_rank == master) then
     allocate(K_core(num_rows(1),N))
     !allocate(E(num_rows(1),size(M(1,:))))
     K_core = reshape(K(start_vec(1):(start_vec(2)-1),N),(/num_rows(1),N/))
     call par_mat_mul_int(K_core,M,E)
  else
     !allocate(E(num_rows(my_rank+1),size(M(1,:))))
     call par_mat_mul_int(K,M,E)
  end if

  !Calculate the encoded matrix E using mpi matrix portioning parallel algorithm
  !E = matmul(K,M)

  !K_orig = K

  !Sequential row reduction algorithm (uncomment to compare run time with omp version):
  !time_start = omp_get_wtime()
  !call row_red(K,K_inv)
  !time_end = omp_get_wtime()

  !print *,"seq time = ",time_end-time_start
  !print *,''

  !deallocate(K_inv)

  !The following line is necessary only when running sequential (above) AND parallel (below):
  !K = K_orig

  if (my_rank == master) then
     !Call the omp parallel row reduction subroutine to get K's inverse
     time_start = omp_get_wtime()
     call row_red_omp(K,K_inv)
     time_end = omp_get_wtime()
   
     print *,"parallel time = ",time_end-time_start
     print *,''
 

     do i = 2,num_cores
        call mpi_send(K_inv(start_vec(i):(start_vec(i)+num_rows(i)),:),num_rows(i)*N,mpi_double,i-1,tag,mpi_comm_world,ierror)
     end do

  end if

  if (my_rank .ne. master) then
     allocate(K_inv(num_rows(my_rank+1),N))
     call mpi_recv(K_inv,num_rows(my_rank+1)*N,mpi_double,master,tag,mpi_comm_world,mpi_status,ierror)
  end if

  if (my_rank == master) then
     allocate(K_core_inv(num_rows(1),N))
     K_core_inv = reshape(K(start_vec(1):(start_vec(2)-1),N),(/num_rows(1),N/))
     call par_mat_mul(K_core_inv,E,M_core)
  else
     call par_mat_mul(K_inv,E,M_core)
  end if

  
  !Compute the original message message again by multiplying the inverse key and the encoded message, E
  !E1 = matmul(K_inv,E)

  !Convert the decoded matrix of reals into the nearest whole numbers
  M = nint(M_core)

  !Convert the integer decoded matrix of integers into their ASCII equivalent characters
  C = char(M)
  
  do i = 0,num_cores-1
     if (my_rank == master) then
        open(unit=3,file='output.txt')
     end if
     
     if (my_rank == i) then
        write(3,*) M_core
     end if
     call mpi_barrier(mpi_comm_world,ierror)
  end do
  
 
  !Deallocate all arrays used in main
  deallocate(K_inv)
  deallocate(K)
  deallocate(K_orig)
  deallocate(E)
  deallocate(M)

end program main
