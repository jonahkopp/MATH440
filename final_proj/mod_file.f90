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

       !$omp parallel &
       !$omp shared(K,K_inv,N,i) &
       !$omp private(j)

       !$omp do
       do j=i+1,N
          K_inv(j,:) = K_inv(j,:) - K(j,i)*K_inv(i,:)
          K(j,:) = K(j,:) - K(j,i)*K(i,:)
       end do
       !$omp end do

       !$omp end parallel
       
       end do

    do i=0,N-2

       !$omp parallel &
       !$omp shared(K,K_inv,N,i) &
       !$omp private(j)

       !$omp do
       do j=i+1,N-1
          K_inv(N-j,:) = K_inv(N-j,:) - (K(N-j,N-i))*(K_inv(N-i,:))
          K(N-j,:) = K(N-j,:) - (K(N-j,N-i))*(K(N-i,:))
       end do
       !$omp end do

       !$omp end parallel
              
    end do

  end subroutine row_red_omp

  subroutine gen_msg_mat(N,file_name,M)

    implicit none

    !Declare variables
    integer, intent(in) :: N
    integer, allocatable, dimension(:,:), intent(out) :: M
    character(len=20), intent(in) :: file_name
    character(len=2000) :: test_str
    !character(len=20)  :: test_str2
    !character(len=20) :: test_str3
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

    !print *, M
       
    !print *, int_vec

    !do i = 1,len(test_str)

       !test_str2(i:i) = char(int_vec(i))

    !end do
    
    !print *, test_str2

    deallocate(int_vec)
    
    !test_str2 = ichar(test_str)

    !test_str3 = char(test_str2)
    !print *, test_str3
    
  end subroutine gen_msg_mat
  
    
end module mod_file
