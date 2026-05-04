module boundary
    implicit none
    private
    public :: apply_psi_bc, apply_omega_bc
    
contains
    
    ! Apply boundary conditions for stream function (ψ = 0 on all walls)
    subroutine apply_psi_bc(psi, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), dimension(0:nx, 0:ny), intent(inout) :: psi
        
        ! All walls: ψ = 0 (no flow through boundaries)
        psi(0,:) = 0.0d0      ! Left wall
        psi(nx,:) = 0.0d0     ! Right wall
        psi(:,0) = 0.0d0      ! Bottom wall
        psi(:,ny) = 0.0d0     ! Top wall (moving lid)
        
    end subroutine apply_psi_bc
    
    
    ! Apply boundary conditions for vorticity at walls
    ! Using no-slip condition: ω = -∂²ψ/∂n² at walls
    subroutine apply_omega_bc(omega, psi, nx, ny, dx, dy, u_lid)
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, dy, u_lid
        real(8), dimension(0:nx, 0:ny), intent(inout) :: omega
        real(8), dimension(0:nx, 0:ny), intent(in) :: psi
        
        real(8) :: dx2, dy2
        integer :: i, j
        
        dx2 = dx * dx
        dy2 = dy * dy
        
        ! Bottom wall (y=0): u=0, v=0
        do i = 0, nx
            omega(i,0) = -2.0d0 * psi(i,1) / dy2
        end do
        
        ! Top wall (y=1): u=u_lid, v=0 (moving lid)
        do i = 0, nx
            omega(i,ny) = -2.0d0 * psi(i,ny-1) / dy2 - 2.0d0 * u_lid / dy
        end do
        
        ! Left wall (x=0): u=0, v=0
        do j = 0, ny
            omega(0,j) = -2.0d0 * psi(1,j) / dx2
        end do
        
        ! Right wall (x=1): u=0, v=0
        do j = 0, ny
            omega(nx,j) = -2.0d0 * psi(nx-1,j) / dx2
        end do
        
    end subroutine apply_omega_bc
    
end module boundary
