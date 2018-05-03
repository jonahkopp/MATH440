program main
  
  use mpi
  use omp_lib
  use mod_file

  implicit none

  !Declare variables and parameters
  integer, parameter :: my_kind = kind(0.0d0)
  real(kind=my_kind), allocatable, dimension(:,:) :: K, K_core,K_inv, K_orig, K_core_inv,M, E, E1, M_core, M_whole, E_whole
  real(kind=my_kind),allocatable,dimension(:) :: K_core_vec,K_vec
  character,allocatable,dimension(:,:) :: C
  integer,dimension(mpi_status_size) :: mpi_status
  integer :: i,j,N,ierror,my_rank,num_cores,master,div,rem,tag,msg_rows,msg_cols
  real(kind=my_kind) :: timeseed,time_start,time_end
  character(len=20) :: file_name
  integer,allocatable,dimension(:) :: num_rows,start_vec
  
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,my_rank,ierror)
  call mpi_comm_size(mpi_comm_world,num_cores,ierror)
  
  master = 0
  tag = 0

  !Set N to be the number of rows/cols of the key matrix, the number of rows of rows of the message matrix/encoded matrix
  N = 400

  call cpu_time(timeseed)

  call srand(int(exp(real(my_rank))))
  
  if (my_rank == master) then
  
     !Prompt user for the file name containing the message to be encoded
     print *, "Enter the file name:"
  
     read(*,*) file_name

     !Generate the original message matrix using the following subroutine
     call gen_msg_mat(N,file_name,M)

     msg_rows = size(M(:,1))
     msg_cols = size(M(1,:))

  end if

  !Broadcast the size of the message to each core
  call mpi_bcast(msg_rows,1,mpi_integer,master,mpi_comm_world,ierror)
  call mpi_bcast(msg_cols,1,mpi_integer,master,mpi_comm_world,ierror)
  
  call mpi_barrier(mpi_comm_world,ierror)

  !All other cores beside master allocate their own message matrix
  if (my_rank .ne. master) then
     allocate(M(msg_rows,msg_cols))
  end if
  
  !call mpi_barrier(mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)
  
  !Broadcast the message matrix to all cores
  call mpi_bcast(M,size(M(1,:))*size(M(:,1)),mpi_double,master,mpi_comm_world,ierror)
  
  allocate(num_rows(num_cores),start_vec(num_cores))

  !Set div to the initial number of rows of K that each core will have
  div = N/num_cores

  !Set rem to be the number of leftover rows
  rem = mod(N,num_cores)

  !Num_rows is a vector that holds the number of rows that each core will work with
  do i = 1,num_cores
     num_rows(i) = div
  end do

  !Keep distributing the remainder of rows to the cores until all rows of K are covered
  do i = 1,rem
     num_rows(i) = num_rows(i) + 1
  end do

  !Start vector contains the starting index of each core's block of rows of K
  start_vec(1) = 1
  do i = 2,num_cores
     start_vec(i) = num_rows(i-1)+start_vec(i-1)
  end do
  
  !Allocate matrices
  allocate(K_core(num_rows(my_rank+1),N))
  allocate(K_core_vec(num_rows(my_rank+1)*N))

  !Generate the key matrix using random numbers reals between 1 and 2
  do i=1,num_rows(my_rank+1)
     do j=1,N
        K_core(i,j) = rand()+1
     end do
  end do

  !Allocate the full key matrix
  if (my_rank == master) then
     allocate(K(N,N))
     allocate(K_vec(N*N))
  end if

  call par_mat_mul(K_core,M,E)

  do i = 1,num_rows(my_rank+1)
     do j = 1,N
        K_core_vec((i-1)*N+j) = K_core(i,j)
     end do
  end do

  print *,'sum of K = ',sum(K_core),'vector',sum(K_core_vec)  
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  call mpi_gather(K_core_vec,N*N/num_cores,mpi_double,K_vec,N*N/num_cores,mpi_double,master,mpi_comm_world,ierror)
  
  call mpi_barrier(mpi_comm_world,ierror)

  !Testing if the full K is correct (it is not as of now)
  if (my_rank == master) then
     do i = 1,N
        do j = 1,N
           K(i,j) = K_vec((i-1)*N+j)
        end do
     end do
     print *, K(1,1:5),K(101,1:5),K(201,1:5),K(301,1:5)
  end if

  if (my_rank == master) then
     print *,'sum of all K ',sum(K),'vector: ',sum(K_vec)
  end if

  !call mpi_barrier(mpi_comm_world,ierror)





  
  
  allocate(K_core_inv(num_rows(my_rank+1),N))
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  call mpi_scatterv(K_inv,num_rows*N,(start_vec-1)*N+1,mpi_double, &
       K_core_inv,(N/num_cores)*N,mpi_double,master,mpi_comm_world,ierror)
  
  allocate(E_whole(msg_rows,msg_cols))
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  !Trying to get all of E to all of the cores (CHANGE THIS TO ALLGATHERV WHEN WE FIGURE OUT HOW TO HAVE VARYING NUM_ROWS)
  call mpi_allgather(E,msg_rows*msg_cols,mpi_double,E_whole,msg_rows*msg_cols, &
       mpi_double,mpi_comm_world,ierror)
  
  M = real(M,my_kind)

  print *, E_whole(1,:)
  
  call par_mat_mul(K_core_inv,E_whole,M)

  !Compute the original message message again by multiplying the inverse key and the encoded message, E
  !E1 = matmul(K_inv,E)

  !Convert the decoded matrix of reals into the nearest whole numbers

  if (my_rank == master) then
     allocate(M_whole(msg_rows,msg_cols))
  end if

  print *, size(M)

  if (my_rank == master) then
     print *, size(M_whole)
  end if
  
  !gathering the decoded messaage in master
  call mpi_gather(M,size(M),mpi_double,M_whole, &
       size(M_whole),mpi_double,master,mpi_comm_world,ierror)

  !print *, ierror, my_rank
  
  print *, my_rank, "made it past last gather call"
  
  !Convert the integer decoded matrix of integers into their ASCII equivalent characters
  if (my_rank == master) then
     C = char(nint(M_whole))
     open(unit=3,file='output.txt')
     write(3,*) C
     close(3)
  end if
   
  !Deallocate all arrays used in master core
  if (my_rank == master) then
     deallocate(K_inv)
     deallocate(K)
     deallocate(C)
     deallocate(M_whole)

  end if
     

     !Deallocate the arrays that all cores use
     deallocate(start_vec)
     deallocate(num_rows)
     deallocate(E)
     deallocate(M)
     deallocate(K_core)
     deallocate(K_core_inv)
     deallocate(E_whole)
  
  call mpi_finalize(ierror)
  
end program main
