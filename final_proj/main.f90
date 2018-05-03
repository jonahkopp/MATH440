program main
  
  use mpi
  use omp_lib
  use mod_file

  implicit none

  !Declare variables and parameters
  integer, parameter :: my_kind = kind(0.0d0)
  real(kind=my_kind), allocatable, dimension(:,:) :: K, K_core,K_inv, K_orig, K_core_inv,M, E, E1, M_core, M_whole, E_whole
  real(kind=my_kind),allocatable,dimension(:) :: K_core_vec,K_vec,K_inv_vec,K_core_inv_vec,E_vec,E_whole_vec,M_vec,M_whole_vec
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
  N = 100

  !Not sure we need this line anymore
  call cpu_time(timeseed)

  !Seed rng with exp(my_rank+1) so each core has different seed
  call srand(int(exp(real(my_rank+1))))

  !Master core does the following:
  if (my_rank == master) then
  
     !Prompt user for the file name containing the message to be encoded
     print *, "Enter the file name:"
     read(*,*) file_name

     !Generate the original message matrix using the following subroutine
     call gen_msg_mat(N,file_name,M)

     !Save the size of the message matrix so it can be known to all other cores later
     msg_rows = size(M(:,1))
     msg_cols = size(M(1,:))

  end if

  !Broadcast the size of the message to each core
  call mpi_bcast(msg_rows,1,mpi_integer,master,mpi_comm_world,ierror)
  call mpi_bcast(msg_cols,1,mpi_integer,master,mpi_comm_world,ierror)

  !All other cores beside master allocate their own message matrix
  if (my_rank .ne. master) then
     allocate(M(msg_rows,msg_cols))
  end if

  !Make sure all cores are caught up before broadcast call below
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
  
  !Allocate each core's portion of the key matrix, K, in 1-d and 2-d form
  allocate(K_core(num_rows(my_rank+1),N))
  allocate(K_core_vec(num_rows(my_rank+1)*N))

  !Generate the key matrix using random numbers reals between 1 and 2
  do i=1,num_rows(my_rank+1)
     do j=1,N
        K_core(i,j) = rand()+1
     end do
  end do

  !Allocate the full key matrix in 2-d and 1-d form
  if (my_rank == master) then
     allocate(K(N,N))
     allocate(K_vec(N*N))
  end if

  !Get each core's portion of the encoded message, E, by calling the parallel matrix mult. subroutine
  call par_mat_mul(K_core,M,E)

  !Put each portion of K into 1-d vector form so they can be properly gathered by master
  do i = 1,num_rows(my_rank+1)
     do j = 1,N
        K_core_vec((i-1)*N+j) = K_core(i,j)
     end do
  end do

  !Master gathers all of the core's portions of K and puts them into a 1-d form of the full K vector
  call mpi_gatherv(K_core_vec,num_rows(my_rank+1)*N,mpi_double,K_vec,N*num_rows,N*(start_vec-1),&
                   mpi_double,master,mpi_comm_world,ierror)
  
  !Master resizes K_vec in 1-d form to K, which is 2-d form
  if (my_rank == master) then
     do i = 1,N
        do j = 1,N
           K(i,j) = K_vec((i-1)*N+j)
        end do
     end do     
  end if

  !All cores should allocate their portion of K inverse in 2-d and 1-d form
  allocate(K_core_inv(num_rows(my_rank+1),N))
  allocate(K_core_inv_vec(num_rows(my_rank+1)*N))

  !Master core should do the following:
  if (my_rank == master) then

     !Allocate the master's K inverse as a 1-d vector so it can be sent and recieved properly
     allocate(K_inv_vec(N*N))

     !Use the OpenMP row reduction subroutine to find K inverse
     call row_red_omp(K,K_inv,N)

     !ADD CALL TO SUBROUTINE FOR MPI IN ADDITION TO OPENMP VERSION!

     !Populate the K inverse 1-d vector
     do i = 1,N
        do j = 1,N
           K_inv_vec((i-1)*N+j) = K_inv(i,j)
        end do
     end do
     
  end if
  
  !Must make sure all cores are caught up before making the scatterv call
  call mpi_barrier(mpi_comm_world,ierror)

  !Scatter the full K inv vector from master to all cores and call each core's portion K_core_inv_vec
  call mpi_scatterv(K_inv_vec,num_rows*N,(start_vec-1)*N,mpi_double,&
       K_core_inv_vec,num_rows(my_rank+1)*N,mpi_double,master,mpi_comm_world,ierror)

  !Allocate the full encoded message matrix, that in 1-d vector form, as well as the portions
  !of the encoded message matrix E
  allocate(E_whole(msg_rows,msg_cols))
  allocate(E_whole_vec(msg_rows*msg_cols))
  allocate(E_vec(num_rows(my_rank+1)*msg_cols))

  !Now that the scatter call has been made, each core's portion of E can be put back into 2-d form
  do i = 1,num_rows(my_rank+1)
     do j = 1,N
        K_core_inv(i,j) = K_core_inv_vec((i-1)*N+j)
     end do
  end do

  !Put the encoded message portion (for each core) into 1-d vector form so they can be gathered
  do i = 1,num_rows(my_rank+1)
     do j = 1,msg_cols
        E_vec((i-1)*msg_cols+j) = E(i,j)
     end do
  end do  
  
  !Gather all portions of E in vector form and put them in E_whole_vec (1-d), which each core will have
  call mpi_allgatherv(E_vec,num_rows(my_rank+1)*msg_cols,mpi_double,E_whole_vec,num_rows*msg_cols,&
       (start_vec-1)*msg_cols,mpi_double,mpi_comm_world,ierror)

  !Now we put the full encoded message matrix into 2-d matrix form
  do i = 1,msg_rows
     do j = 1,msg_cols
        E_whole(i,j) = E_whole_vec((i-1)*msg_cols+j)
     end do
  end do

  !Not sure we need the following barrier anymore:
  call mpi_barrier(mpi_comm_world,ierror)

  !Call the MPI parallel matrix mult. subroutine to get each core's portion of message matrix M
  call par_mat_mul(K_core_inv,E_whole,M)

  allocate(M_vec(num_rows(my_rank+1)*msg_cols))

  !Put each core's portion of M into vector, 1-d, form so they can be gathered by master
  do i = 1,num_rows(my_rank+1)
     do j = 1,msg_cols
        M_vec((i-1)*msg_cols+j) = M(i,j)
     end do
  end do

  !Deallocating all arrays that are no longer necessary to finish the program
  !This will free up space for the allocations coming up
  deallocate(E)
  deallocate(E_vec)
  deallocate(E_whole)
  deallocate(E_whole_vec)
  deallocate(K_core)
  deallocate(K_core_inv)

  !Master core now allocates the full message matrix, that in 1-d vector form, as well as C
  !where C is the message matrix converted to ASCII code equivalents (readable letters)
  if (my_rank == master) then
     allocate(M_whole(msg_rows,msg_cols))
     allocate(M_whole_vec(msg_rows*msg_cols))
     allocate(C(msg_rows,msg_cols))    
  end if
  
  !Gather all cores' portions of M and assemble into M_whole_vec
  call mpi_gatherv(M_vec,num_rows(my_rank+1)*msg_cols,mpi_double,M_whole_vec,&
       num_rows*msg_cols,(start_vec-1)*msg_cols,mpi_double,master,mpi_comm_world,ierror)
  
  !Convert the integer decoded matrix of integers into their ASCII equivalent characters
  if (my_rank == master) then
     do i = 1,msg_rows
        do j = 1,msg_cols
           M_whole(i,j) = M_whole_vec((i-1)*msg_cols + j)
        end do
     end do

     !Put message matrix, M, into readable character form, C
     C = char(nint(M_whole))

     !Write the decoded message to the output text file
     open(unit=3,file='output.txt')
     do i=1,msg_rows    
        write(3,*) C(i,:)
     end do
     close(3)
  end if
   
  !Deallocate all arrays used in master core
  if (my_rank == master) then
     deallocate(K_inv)
     deallocate(K)
     deallocate(C)
     deallocate(M_whole)
  end if

  !Deallocate the rest of the arrays that all cores used
  deallocate(start_vec)
  deallocate(num_rows)
  deallocate(M)
  deallocate(M_vec)

  !Finalize MPI
  call mpi_finalize(ierror)
  
end program main
