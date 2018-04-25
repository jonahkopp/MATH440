program main
  
  use mpi
  use omp_lib
  use mod_file

  implicit none

  !Declare variables and parameters
  integer, parameter :: my_kind = selected_real_kind(16,300)
  real(kind=my_kind), allocatable, dimension(:,:) :: K, K_inv, K_orig, K_core, K_core_inv,M, E, E1, M_core, M_whole, E_whole
  !integer, allocatable, dimension(:,:) :: M
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
  
  if (my_rank == master) then
  
     call cpu_time(timeseed)

     !Seed rng with the cpu time found above
     call srand(int(timeseed))

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

  !Each core finds its num_rows(my_rank+1) rows of the key matrix
  !if (my_rank == master) then
  
     !Allocate matrices
     allocate(K_core(num_rows(my_rank+1),N))
     !allocate(K_orig(N,N))
     !allocate(E(size(M(:,1)),size(M(1,:))),E1(size(M(:,1)),size(M(1,:))))

     !Generate the key matrix using random numbers reals between 1 and 2
     do i=1,num_rows(my_rank+1)
        do j=1,N
           K_core(i,j) = rand()+1
        end do
     end do
     
     !Master will send each block of K to the respective core with rank i-1
     !do i = 2,num_cores
      !  print *, start_vec(i)+num_rows(i)-1
      !  call mpi_send(K(start_vec(i):(start_vec(i)+num_rows(i)-1),:),num_rows(i)*N,mpi_double,i-1,tag,mpi_comm_world,ierror)
     !end do
     
     !K = transpose(K)
     
  !end if
  
  !allocate(K_core(num_rows(my_rank+1),N))
 
  call mpi_barrier(mpi_comm_world,ierror)
  
  !call mpi_scatterv(K,num_rows*N,(start_vec-1)*N+1,mpi_double,K_core &
   !    ,(N/num_cores)*N,mpi_double,master,mpi_comm_world,ierror)

  !call mpi_barrier(mpi_comm_world,ierror)

  !Each core calculates its section of the encoded msg matrix, E
  call par_mat_mul(K_core,M,E)

  call mpi_barrier(mpi_comm_world,ierror)
  
  !If not the master core, then allocate the block of K and call mpi recv so that each core obtains the block of K that was sent by master
  !if (my_rank .ne. master) then
   !  allocate(K(num_rows(my_rank+1),N))
    ! call mpi_recv(K,num_rows(my_rank+1)*N,mpi_double,master,tag,mpi_comm_world,mpi_status,ierror)
  !end if

  !If master core, then set its block of K (upper-most section), called K_core, to the correct set of rows of original K
  !if (my_rank == master) then
     !allocate(K_core(num_rows(1),N))
     !allocate(E(num_rows(1),size(M(1,:))))

     !K_core = reshape(K(start_vec(1):(start_vec(2)-1),N),(/num_rows(1),N/))

     !Call the parallel matrix multiplication subroutine with the top section of K and the message matrix M
   !  call par_mat_mul_int(K_core,M,E)
  !else
     !allocate(E(num_rows(my_rank+1),size(M(1,:))))

     !Parallel matrix mult. with the other cores' blocks of K and the msg matrix M
    ! call par_mat_mul_int(K_core,M,E)
  !end if

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

  !Allocate the full key matrix
  allocate(K(N,N))

  !THIS GATHER IS PRODUCING A FULL K MATRIX WITH VERY SMALL VALUES RATHER THAN VALUES BETW 1 AND 2:
  call mpi_gatherv(K_core,num_rows(my_rank+1)*N,mpi_double,K,num_rows*N, &
       (start_vec-1)*N+1,mpi_double,master,mpi_comm_world,ierror)

  !Testing if the full K is correct (it is not as of now)
  if (my_rank == master) then
     print *, K(1:5,1:5)
     print *, size(K(1,:)), size(K(:,1))
  end if

  !Obtain K inverse matrix
  if (my_rank == master) then
     !Call the omp parallel row reduction subroutine to get K's inverse
     time_start = omp_get_wtime()
     
     call row_red_omp(K,K_inv,N)
     
     time_end = omp_get_wtime()
   
     print *,"parallel time = ",time_end-time_start
     print *,''

     !Master sends each block of K inverse to the other cores
     !do i = 2,num_cores
        !call mpi_send(K_inv(start_vec(i):(start_vec(i)+num_rows(i)),:),num_rows(i)*N,mpi_double,i-1,tag,mpi_comm_world,ierror)
     !end do

  end if

  call mpi_barrier(mpi_comm_world,ierror)

  allocate(K_core_inv(num_rows(my_rank+1),N))
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  call mpi_scatterv(K_inv,num_rows*N,(start_vec-1)*N+1,mpi_double, &
       K_core_inv,(N/num_cores)*N,mpi_double,master,mpi_comm_world,ierror)

  !All other cores recv the blocks of K inverse, as above with the original blocks of K
  !if (my_rank .ne. master) then
     !allocate(K_inv(num_rows(my_rank+1),N))
     !call mpi_recv(K_inv,num_rows(my_rank+1)*N,mpi_double,master,tag,mpi_comm_world,mpi_status,ierror)
  !end if

  !Parallel matrix mult. using the subroutine again, but this time with the inverse of K
  !if (my_rank == master) then
     !allocate(K_core_inv(num_rows(1),N))
     !K_core_inv = reshape(K(start_vec(1):(start_vec(2)-1),N),(/num_rows(1),N/))
     !call par_mat_mul(K_core_inv,E,M_core)
  !else
     !call par_mat_mul(K_inv,E,M_core)
  !end if

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
  
  !Each core writes the decoded message to the output file in order (based on my_rank)
  !do i = 0,num_cores-1
     !Only the master opens the file for writing
     !if (my_rank == master) then
        !open(unit=3,file='output.txt')
     !end if
     
     !if (my_rank == i) then
        !write(3,*) M
     !end if
     !Barrier assures that the cores write in order so message makes sense to reader
     !call mpi_barrier(mpi_comm_world,ierror)
  !end do
   
  !Deallocate all arrays used in master core
  if (my_rank == master) then
     deallocate(K_inv)
     deallocate(K)
     deallocate(C)
     deallocate(M_whole)
  end if

  !Deallocate the arrays that all cores use
  !deallocate(K_orig)
  deallocate(start_vec)
  deallocate(num_rows)
  deallocate(E)
  deallocate(M)
  deallocate(K_core)
  deallocate(K_core_inv)
  deallocate(E_whole)
  
  call mpi_finalize(ierror)
  
end program main
