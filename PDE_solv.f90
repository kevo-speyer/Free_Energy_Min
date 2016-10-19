program main
implicit none
integer :: i, j, N_tot, dbg_count = 1
real(kind=8) :: L_box, eps, dz, conv_eps = 0.00001**2 !square of conv criteria
real(kind=8) , dimension(:), allocatable :: rho, rho_new
logical :: conver = .False.! Converged to solution True or False

!Read input
call read_input(L_box,eps,N_tot, dz)

allocate(rho(N_tot),rho_new(N_tot))

!open file to save data
open (unit = 73, file = "rho_vs_z.dat")

!Make initial guess of rho
call init_guess(rho,N_tot)

print*, " Starting Iteratrive approach to Solution"
do while (.not.conver)
   !Update solution to diffretian equation
    call update(rho,rho_new,N_tot,eps, dz)

    !Evaluate convergence
    call eval_conv(rho,rho_new,N_tot,conver,conv_eps)    

    ! move New to old variables
    rho = rho_new

end do

print*, " Converged to Solution"

!Save result
call write_result(rho, N_tot, L_box)

!Close file with rho vs L
close (73)


end program

subroutine write_result(rho, N_tot, L_box)
real(kind=8) , dimension(N_tot), intent(in) :: rho
integer, intent(in) :: N_tot
real(kind=8) :: z
real(kind=8), intent(in) :: L_box
integer :: i

do i = 1, N_tot
    z = float((i-1))/float((N_tot-1))*2.*L_box - L_box
    write(73,*) z, rho(i)
end do

write(73,*) ""

end subroutine
 
subroutine eval_conv(rho, rho_new, N_tot, conver, conv_eps)
integer :: i 
real(kind=8), intent(in) :: conv_eps
real(kind=8) , dimension(N_tot), intent(in) :: rho, rho_new
real(kind=8) , dimension(N_tot) ::  diff
integer, intent(in) :: N_tot
logical, intent(inout) :: conver
real(kind=8) :: max_diff_2

diff(:) = (rho_new(:) - rho(:) ) ** 2

max_diff_2 = maxval(diff)

if(max_diff_2.le.conv_eps) then !Apply convergence criteria
    conver = .True.
    !print*,max_diff
end if

end subroutine

subroutine init_guess(rho,N_tot)
integer :: i
real(kind=8) :: inv_Ntot_1
integer, intent(in) :: N_tot
real(kind=8) , dimension(N_tot), intent(inout) :: rho

inv_Ntot_1 = 1. / float((N_tot-1))

do i=2, N_tot-1
    rho(i) = 2. * inv_Ntot_1 *(float(i)-1.) - 1.
end do

!Boundary conditions
rho(1) = -1.
rho(N_tot) = 1.

print*, " Initial Guess done"

end subroutine

subroutine update(rho, rho_new, N_tot, eps, dz ) 
integer :: i
integer, intent(in) :: N_tot
real(kind=8), intent(in) :: eps, dz
real(kind=8) , dimension(N_tot), intent(inout) :: rho, rho_new
real(kind=8) :: dz2, mu

dz2 = dz * dz

!FixedBoundary conditions
rho_new(1) = rho(1)
rho_new(N_tot) = rho(N_tot)

!Update rho
do i = 2, N_tot-1 ! Dont touch boundaries
    mu = -rho(i) + rho(i)**3 -( rho(i+1) + rho(i-1) - 2 * rho(i)) / dz2
    rho_new(i) = rho(i) - eps * mu 
end do

end subroutine

subroutine read_input(L,eps,N_steps,dz)
real(kind=8), intent(out) :: L, eps, dz
integer, intent(out) :: N_steps
open (unit = 53, file = "input.dat", status= "old")

read(53,*) L
read(53,*) N_steps
read(53,*) eps
close(53)

print*, " Reading done"
dz = L / float(N_steps) ! Bin size

N_steps = N_steps * 2 ! Total # of steps is double of input

end subroutine
