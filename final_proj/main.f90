program main

  use mod_file
  use omp_lib
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
  
  call srand(int(timeseed))
  
  print *, "Enter the file name:"
  read(*,*) file_name

  N = 100

  call gen_msg_mat(N,file_name,M)

  !print *, test_str
  
  allocate(K(N,N),K_orig(N,N))
  allocate(E(size(M(:,1)),size(M(1,:))),E1(size(M(:,1)),size(M(1,:))))

  
  do i=1,N
     do j=1,N
        !K(i,j) = sin(real(i,my_kind)) + cos(real(j,my_kind))
        !K(i,j) = 3.0*(i-1) + j + sqrt((real(i)+1.0)*(J+2))*0.5
        K(i,j) = rand()+1
     end do
  end do

  !print *, K(1,:)
  !print *, K(2,:)
  !print *, K(3,:)

  E = matmul(K,M)

  !print *, E
  
  !allocate(K_orig(N,N))

  K_orig = K

  call cpu_time(time_start)
  call row_red(K,K_inv)
  call cpu_time(time_end)

  print *,"seq time = ",time_end-time_start
  print *,''

  deallocate(K_inv)
  
  K = K_orig

  call cpu_time(time_start)
  call row_red_omp(K,K_inv)
  call cpu_time(time_end)
  
  print *,"parallel time = ",time_end-time_start
  print *,''

  E1 = matmul(K_inv,E)

  !print *, E


  M = nint(E1)

  !print *, M

  C = char(M)

  print *, C(1,:)
  print *, C(2,:)
  print *, C(3,:)
  print *, C(4,:)
  print *, C(5,:)
  print *, C(6,:)
  

! do i = 1,N
 !   do j = 1,size(M(1,:))



!    end do
! end do
  
  
  !print *, K(1,:)
  !print *, K(2,:)
  !print *, K(3,:)

  !K = matmul(K_orig,K_inv)
  
  !print *, K(1,:)
  !print *, K(2,:)
  !print *, K(3,:)

  deallocate(K_inv)
  deallocate(K)
  deallocate(K_orig)
  deallocate(E)
  deallocate(M)

end program main
