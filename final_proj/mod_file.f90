module mod_file
use mpi
contains

  !parallel algorithm for row reduction below
  subroutine row_red_omp(K,K_inv,N)

    implicit none
    
    !Declare variables and parameters
    integer, parameter :: my_kind = kind(0.0d0)
    integer, intent(in) :: N
    real(kind=my_kind), allocatable, dimension(:,:), intent(inout) :: K
    real(kind=my_kind), allocatable, dimension(:,:), intent(out) :: K_inv
    integer :: i,j,iam,tot,numthreads
    integer,external :: omp_get_thread_num, omp_get_num_threads
    
    allocate(K_inv(N,N))
    
    !Populate identity matrix of size N by N
    K_inv(:,:) = real(0,my_kind)

    do i=1,N
       K_inv(i,i) = real(1,my_kind)
    end do
    
    i=1

    !Set number of threads to 8 (optimal) unless N is less than 8 (then it is inefficient to use omp)
    if (N<8) then
       numthreads=1
    else
       numthreads=8
    end if

    !Set number of threads
    call omp_set_num_threads(numthreads)
    
    !$omp parallel &
    !$omp shared(K,K_inv,N,i) &
    !$omp private(iam,j)

    !iam is the rank of the current thread
    iam = omp_get_thread_num()

    !$omp barrier
    
    do while (i <= N)

       !$omp barrier
       !Scale the leading row
       if (iam == 0) then 

          K_inv(i,:) = (1/K(i,i))*K_inv(i,:)
          K(i,:) = (1/K(i,i))*K(i,:)

       end if
       
       !$omp barrier
       
       !Clear out all entries below leading row's first entry (1 entry)
       !$omp do
       do j=i+1,N
          
          K_inv(j,:) = K_inv(j,:) - K(j,i)*K_inv(i,:)
          K(j,:) = K(j,:) - K(j,i)*K(i,:)
          
       end do
       !$omp end do

       !$omp barrier
       
       !Increment i
       if (iam == 0) then
          i = i + 1
       end if

       !$omp barrier
       
    end do

    !$omp barrier
    
    if (iam == 0) then
       i = 1
    end if

    !$omp barrier
    
    !Repeat but to clear the upper triangle this time; result is K is now the identity
    do while (i <= N)
       
       !$omp do
       do j=i,N-1
          K_inv(N-j,:) = K_inv(N-j,:) - (K(N-j,N-i+1))*(K_inv(N-i+1,:))
          K(N-j,:) = K(N-j,:) - (K(N-j,N-i+1))*(K(N-i+1,:))
       end do
       !$omp end do

       !$omp barrier
       
       if (iam == 0) then
          i=i+1
       end if

       !$omp barrier
       
    end do
    
    !$omp barrier
    !$omp end parallel
    
  end subroutine row_red_omp

  !Generates message matrix
  subroutine gen_msg_mat(N,file_name,M)

    implicit none

    !Declare variables
    integer, parameter :: my_kind = kind(0.0d0)
    integer, intent(in) :: N
    real(kind=my_kind), allocatable, dimension(:,:), intent(out) :: M
    character(len=20), intent(in) :: file_name
    character(len=2000) :: test_str
    integer :: i, j, num_cols
    integer, allocatable, dimension(:) :: int_vec
    real(kind=my_kind) :: rand_num

    !Open input file and read in message to test_str
    open(unit=1,file=file_name,action='read')

    read(1,'(A)') test_str
       
    close(1)

    !Eliminate trailing whitespace
    test_str = trim(test_str)

    !int_vec will hold the ascii equivalent integers for each char in msg
    allocate(int_vec(len(test_str)))

    do i=1,len(test_str)

       int_vec(i) = ichar(test_str(i:i))

    end do

    !Number of columns for msg matrix will be the msg size divided by N
    num_cols = size(int_vec)/N

    !Add one to number of columns if mod of msg size and N is 0
    if (mod(size(int_vec),N) .ne. 0) then
       
       num_cols = num_cols + 1

    end if
    
    allocate(M(N,num_cols))

    !Populate msg matrix M row by row using int_vec and then random nonsense symbols for the remaining space
    do i = 1,N
       do j = 1,num_cols

          if ((i-1)*num_cols+j <= size(int_vec)) then

             M(i,j) = real(int_vec((i-1)*num_cols+j),my_kind)

          else

             call random_number(rand_num)

             M(i,j) = real(rand_num*(47-34)+34,my_kind)

          end if
                    
       end do
    end do

    deallocate(int_vec)
    
  end subroutine gen_msg_mat

  subroutine par_mat_mul(K,M,A)

    implicit none

    !Declare variables
    integer, parameter :: my_kind = kind(0.0d0)
    real(kind=my_kind), allocatable, dimension(:,:), intent(in) :: K,M
    real(kind=my_kind), allocatable, dimension(:,:), intent(out) :: A
    integer :: i,j,p,iam
    integer, external :: omp_get_thread_num
    
    allocate(A(size(K(:,1)),size(M(1,:))))

    !Initialize solution matrix to zeros
    A = 0
    
    !$omp parallel &
    !$omp shared(K,M,A) &
    !$omp private(i,j,p)

    !$omp do
    
    !Each thread will handle certain rows of the matrix for dot products in matrix mult
    do i=1,size(K(:,1))
       do j=1,size(M(1,:))
          do p=1,size(K(1,:))
             A(i,j) = A(i,j) + K(i,p)*M(p,j)
          end do
       end do
    end do

    !$omp end do

    !$omp end parallel
    
  end subroutine par_mat_mul

end module mod_file
