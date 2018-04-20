module mod_file
contains

  subroutine row_red(K,K_inv)

    implicit none
    
    !Declare variables and parameters
    integer, parameter :: my_kind = selected_real_kind(16,300)
    real(kind=my_kind), allocatable, dimension(:,:), intent(inout) :: K
    real(kind=my_kind), allocatable, dimension(:,:), intent(out) :: K_inv
    integer :: i,j,N

    N = size(K(1,:))
    
    allocate(K_inv(N,N))

    K_inv(:,:) = real(0,my_kind)

    do i=1,N
       K_inv(i,i) = real(1,my_kind)
    end do

    do i=1,N
       K_inv(i,:) = (1/K(i,i))*K_inv(i,:)
       K(i,:) = (1/K(i,i))*K(i,:)
       do j=i+1,N
          K_inv(j,:) = K_inv(j,:) - K(j,i)*K_inv(i,:)
          K(j,:) = K(j,:) - K(j,i)*K(i,:)
       end do
    end do

    do i=0,N-2
       do j=i+1,N-1
          K_inv(N-j,:) = K_inv(N-j,:) - (K(N-j,N-i))*(K_inv(N-i,:))
          K(N-j,:) = K(N-j,:) - (K(N-j,N-i))*(K(N-i,:))
       end do
    end do

  end subroutine row_red




  !parallel algorithm below

  subroutine row_red_omp(K,K_inv)

    implicit none
    
    !Declare variables and parameters
    integer, parameter :: my_kind = selected_real_kind(16,300)
    real(kind=my_kind), allocatable, dimension(:,:), intent(inout) :: K
    real(kind=my_kind), allocatable, dimension(:,:), intent(out) :: K_inv
    integer :: i,j,N,iam,tot,numthreads
    integer,external :: omp_get_thread_num, omp_get_num_threads

    N = size(K(1,:))
    
    allocate(K_inv(N,N))

    K_inv(:,:) = real(0,my_kind)

    do i=1,N
       K_inv(i,i) = real(1,my_kind)
    end do

    i=1

    if (N<8) then
       numthreads=1
    else
       numthreads=8
    end if
    
    
    call omp_set_num_threads(numthreads)
    
    !$omp parallel &
    !$omp shared(K,K_inv,N,i) &
    !$omp private(iam,j)

    iam = omp_get_thread_num()
    
    do while (i <= N)
       
       if (iam == 0) then 

          K_inv(i,:) = (1/K(i,i))*K_inv(i,:)
          K(i,:) = (1/K(i,i))*K(i,:)

       end if

       !$omp barrier
    
       !$omp do
       do j=i+1,N
          
          K_inv(j,:) = K_inv(j,:) - K(j,i)*K_inv(i,:)
          K(j,:) = K(j,:) - K(j,i)*K(i,:)
          
       end do
       !$omp end do

       !$omp barrier
       
       if (iam == 0) then
          i = i + 1
       end if

       !$omp barrier
       
    end do    
    
    if (iam == 0) then
       i = 1
    end if

    !$omp barrier
    
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

  subroutine gen_msg_mat(N,file_name,M)

    implicit none

    !Declare variables
    integer, intent(in) :: N
    integer, allocatable, dimension(:,:), intent(out) :: M
    character(len=20), intent(in) :: file_name
    character(len=2000) :: test_str
    integer :: i, j, num_cols
    integer, allocatable, dimension(:) :: int_vec
    
    open(unit=1,file=file_name)

    read(1,'(A)') test_str
       
    close(1)

    test_str = trim(test_str)

    !print *, len(test_str)

    allocate(int_vec(len(test_str)))

    do i=1,len(test_str)

       int_vec(i) = ichar(test_str(i:i))

    end do

    num_cols = size(int_vec)/N

    if (mod(size(int_vec),N) .ne. 0) then
       
       num_cols = num_cols + 1

    end if
    
    allocate(M(N,num_cols))

    do i = 1,N
       do j = 1,num_cols

          if ((i-1)*num_cols+j <= size(int_vec)) then

             M(i,j) = int_vec((i-1)*num_cols+j)

          else

             M(i,j) = rand()*(47-34)+34

          end if
                    
       end do
    end do

    deallocate(int_vec)
    
  end subroutine gen_msg_mat

  subroutine par_mat_mul(A,B,N)

    implicit none

    !Declare variables
     

  end subroutine par_mat_mul
  
    
end module mod_file
