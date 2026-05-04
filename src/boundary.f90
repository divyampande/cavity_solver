module boundary
    use input, only: wp
    implicit none
    private
    public :: apply_psi_bc, apply_omega_bc
    
contains
    
    ! Apply boundary conditions for stream function (ψ = 0 on all walls)
    subroutine apply_psi_bc(psi, nx, ny)
        integer, intent(in) :: nx, ny
        real(wp), dimension(0:nx, 0:ny), intent(inout) :: psi
        
        psi(0,:) = 0.0_wp      ! Left wall
        psi(nx,:) = 0.0_wp     ! Right wall
        psi(:,0) = 0.0_wp      ! Bottom wall
        psi(:,ny) = 0.0_wp     ! Top wall (moving lid)
        
    end subroutine apply_psi_bc
    
    
    ! Apply boundary conditions for vorticity at walls
    subroutine apply_omega_bc(omega, psi, nx, ny, dx, dy, u_lid)
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: dx, dy, u_lid
        real(wp), dimension(0:nx, 0:ny), intent(inout) :: omega
        real(wp), dimension(0:nx, 0:ny), intent(in) :: psi
        
        real(wp) :: dx2, dy2
        integer :: i, j
        
        dx2 = dx * dx
        dy2 = dy * dy
        
        ! Bottom wall (y=0): u=0, v=0
        do i = 0, nx
            omega(i,0) = -2.0_wp * psi(i,1) / dy2
        end do
        
        ! Top wall (y=1): u=u_lid, v=0 (moving lid)
        do i = 0, nx
            omega(i,ny) = -2.0_wp * psi(i,ny-1) / dy2 - 2.0_wp * u_lid / dy
        end do
        
        ! Left wall (x=0): u=0, v=0
        do j = 0, ny
            omega(0,j) = -2.0_wp * psi(1,j) / dx2
        end do
        
        ! Right wall (x=1): u=0, v=0
        do j = 0, ny
            omega(nx,j) = -2.0_wp * psi(nx-1,j) / dx2
        end do
        
    end subroutine apply_omega_bc
    
end module boundary
