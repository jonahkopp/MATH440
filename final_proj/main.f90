program main
  
  use mpi
  use omp_lib
  use mod_file

  implicit none

  !Declare variables and parameters
  integer, parameter :: my_kind = kind(0.0d0)
  real(kind=my_kind), allocatable, dimension(:,:) :: K,K_core,K_inv,K_orig,K_core_inv,M,E
  real(kind=my_kind),allocatable,dimension(:,:) :: E1,M_core,M_whole,E_whole,M_final,M_function
  real(kind=my_kind),allocatable,dimension(:) :: K_core_vec,K_vec,K_inv_vec,K_core_inv_vec,E_vec
  real(kind=my_kind),allocatable,dimension(:) :: E_whole_vec,M_vec,M_whole_vec,K_row_vec
  character,allocatable,dimension(:,:) :: C
  integer,dimension(mpi_status_size) :: mpi_status
  integer :: i,j,N,ierror,my_rank,num_cores,master,div,rem,tag,msg_rows,msg_cols,the_core
  real(kind=my_kind) :: timeseed,time_start,time_end
  character(len=20) :: file_name
  integer,allocatable,dimension(:) :: num_rows,start_vec
  character :: omp_or_mpi = 'a'

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,my_rank,ierror)
  call mpi_comm_size(mpi_comm_world,num_cores,ierror)

  master = 0
  tag = 0

  !Set N to be the number of rows/cols of the key matrix, the number of rows of rows of the message matrix/encoded matrix
  N = 3000

  if (my_rank == master) then
  
     !Prompt user for sequential or parallel version
     print *, "Would you like to run OpenMP inversion or MPI inversion? 'o' for OpenMP or 'm' for MPI:"

     !Keep repeating user input if they didn't enter 'o' or 'm'
     do while (omp_or_mpi .ne. 'o' .and. omp_or_mpi .ne. 'm')
        read(*,*) omp_or_mpi
     end do

  end if
  
  call mpi_bcast(omp_or_mpi,1,mpi_character,master,mpi_comm_world,ierror)
 
  !Used to seed the RNG
  call cpu_time(timeseed)

  !Seed rng with exp(my_rank+1) so each core has different seed
  call srand(int(exp(real(my_rank+1)))*nint(timeseed))

  !Master core does the following:
  if (my_rank == master) then
  
     !Prompt user for the file name containing the message to be encoded
     print *, "Enter the message file name:"
     read(*,*) file_name



     time_start = omp_get_wtime()
     
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
  
  if (omp_or_mpi == 'o') then
  
     !Master core should do the following:
     if (my_rank == master) then

        !Allocate the master's K inverse as a 1-d vector so it can be sent and recieved properly
        allocate(K_inv_vec(N*N))
     
        !Use the OpenMP row reduction subroutine to find K inverse
        call row_red_omp(K,K_inv,N)

        !Populate the K inverse 1-d vector
        do i = 1,N
           do j = 1,N
              K_inv_vec((i-1)*N+j) = K_inv(i,j)
           end do
        end do

     end if

  else

    allocate(K_row_vec(N))

    K_core_inv(:,:) = real(0,my_kind)

    do i=1,num_rows(my_rank+1)
       K_core_inv(i,start_vec(my_rank+1)+i-1) = real(1,my_kind)
    end do


    !clearing the lower triangle
    
    i = 1
    
    do while (i <= N)

       
       if (i >= start_vec(my_rank+1) .and. i < start_vec(my_rank+1)+num_rows(my_rank+1)) then

          K_core_inv(i-start_vec(my_rank+1)+1,:) = 1/K_core(i-start_vec(my_rank+1)+1,i-start_vec(my_rank+1)+1)*&
               K_core_inv(i-start_vec(my_rank+1)+1,:)
          K_core(i-start_vec(my_rank+1)+1,:) = 1/K_core(i-start_vec(my_rank+1)+1,i-start_vec(my_rank+1)+1)*&
               K_core(i-start_vec(my_rank+1)+1,:)
          
          K_row_vec = K_core(i-start_vec(my_rank+1)+1,:)

       end if

       do j = 0,num_cores-1

          if (i >= start_vec(j+1) .and. i < start_vec(j+1)+num_rows(my_rank+1)) then
             the_core = j
          end if
             
       end do
       
       call mpi_bcast(K_row_vec,N,mpi_double,the_core,mpi_comm_world,ierror)

       if (i >= start_vec(my_rank+1) .and. i < start_vec(my_rank+1)+num_rows(my_rank+1)) then

          do j = i-start_vec(my_rank+1)+1,num_rows(my_rank+1)

             K_core_inv(j,:) = K_core_inv(j,:) - K_row_vec*K_core(j,i)
             K_core(j,:) = K_core(j,:) - K_row_vec*K_core(j,i)

          end do

       else
          
          do j = 1,num_rows(my_rank+1)

             K_core_inv(j,:) = K_core_inv(j,:) - K_row_vec*K_core(j,i)
             K_core(j,:) = K_core(j,:) - K_row_vec*K_core(j,i)

          end do

       end if

          i = i + 1
       
    end do

    
    !clearing the upper triangle
    i = 1
    
    do while (i <= N)

       
       
       if (i >= start_vec(my_rank+1) .and. i < start_vec(my_rank+1)+num_rows(my_rank+1)) then

          K_core_inv(i-start_vec(my_rank+1)+1,:) = 1/K_core(i-start_vec(my_rank+1)+1,i-start_vec(my_rank+1)+1)*&
               K_core_inv(i-start_vec(my_rank+1)+1,:)
          K_core(i-start_vec(my_rank+1)+1,:) = 1/K_core(i-start_vec(my_rank+1)+1,i-start_vec(my_rank+1)+1)*&
               K_core(i-start_vec(my_rank+1)+1,:)
          
          K_row_vec = K_core(i-start_vec(my_rank+1)+1,:)

       end if

       do j = 0,num_cores-2

          if (i >= start_vec(j+1) .and. i < start_vec(j+1)+num_rows(my_rank+1)) then
             the_core = j
          end if
             
       end do
    
       call mpi_bcast(K_row_vec,N,mpi_double,the_core,mpi_comm_world,ierror)
       
       if (i >= start_vec(my_rank+1) .and. i < start_vec(my_rank+1)+num_rows(my_rank+1)) then

          do j = 1,start_vec(my_rank+1)+num_rows(my_rank+1)-(N+i)

             K_core_inv(j,:) = K_core_inv(j,:) - K_row_vec*K_core(j,i)
             K_core(j,:) = K_core(j,:) - K_row_vec*K_core(j,i)

          end do

       else
          
          do j = 1,num_rows(my_rank+1)

             K_core_inv(j,:) = K_core_inv(j,:) - K_row_vec*K_core(j,i)
             K_core(j,:) = K_core(j,:) - K_row_vec*K_core(j,i)

          end do

       end if

       i = i + 1
       
    end do

    deallocate(K_row_vec)
    
 end if


       
  !Must make sure all cores are caught up before making the scatterv call
  call mpi_barrier(mpi_comm_world,ierror)

  !Allocate the full encoded message matrix, that in 1-d vector form, as well as the portions
  !of the encoded message matrix E
  allocate(E_whole(msg_rows,msg_cols))
  allocate(E_whole_vec(msg_rows*msg_cols))
  allocate(E_vec(num_rows(my_rank+1)*msg_cols))

  
  if (omp_or_mpi == 'o') then
     !Scatter the full K inv vector from master to all cores and call each core's portion K_core_inv_vec
     call mpi_scatterv(K_inv_vec,num_rows*N,(start_vec-1)*N,mpi_double,&
          K_core_inv_vec,num_rows(my_rank+1)*N,mpi_double,master,mpi_comm_world,ierror)
 

     !Now that the scatter call has been made, each core's portion of K_inv can be put back into 2-d form
     do i = 1,num_rows(my_rank+1)
        do j = 1,N
           K_core_inv(i,j) = K_core_inv_vec((i-1)*N+j)
        end do
     end do
  end if


  !Put the encoded message portion (for each core) into 1-d vector form so they can be gathered
  do i = 1,num_rows(my_rank+1)
     do j = 1,msg_cols
        E_vec((i-1)*msg_cols+j) = E(i,j)
     end do
  end do
  
  !Gather all portions of E in vector form and put them in E_whole_vec (1-d), which each core will have
  call mpi_allgatherv(E_vec,num_rows(my_rank+1)*msg_cols,mpi_double,E_whole_vec,num_rows*msg_cols,&
       (start_vec-1)*msg_cols,mpi_double,mpi_comm_world,ierror)
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  !Now we put the full encoded message matrix into 2-d matrix form
  do i = 1,msg_rows
     do j = 1,msg_cols
        E_whole(i,j) = E_whole_vec((i-1)*msg_cols+j)
     end do
  end do  
  
 
  !Call the MPI parallel matrix mult. subroutine to get each core's portion of message matrix M
  call par_mat_mul(K_core_inv,E_whole,M_function)


 
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
  deallocate(K_core_vec)
  deallocate(K_core_inv_vec)

  
  !Master core now allocates the full message matrix, that in 1-d vector form, as well as C
  !where C is the message matrix converted to ASCII code equivalents (readable letters)
  if (my_rank == master) then
     deallocate(K)
     deallocate(K_vec)
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

     time_end = omp_get_wtime()

     print *,'time = ',time_end-time_start
     
     !Write the decoded message to the output text file
     open(unit=3,file='output.txt')
     do i=1,msg_rows    
        write(3,*) C(i,:)
     end do
     close(3)
  end if
   
  !Deallocate all arrays used in master core
  if (my_rank == master) then
     
     if (omp_or_mpi == 'o') then
        deallocate(K_inv)
     end if
     
     deallocate(C)
     deallocate(M_whole)
     deallocate(M_whole_vec)
  end if

  !Deallocate the rest of the arrays that all cores used
  deallocate(start_vec)
  deallocate(num_rows)
  deallocate(M)
  deallocate(M_vec)

  !Finalize MPI
  call mpi_finalize(ierror)

end program main
