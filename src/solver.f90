module solver
    use input, only: wp
    implicit none
    private
    public :: solve_poisson_sor, advance_vorticity, compute_velocity, compute_pressure
    
contains
    
    ! Solve Poisson equation using Successive Over-Relaxation
    subroutine solve_poisson_sor(psi, omega, nx, ny, dx, dy, omega_sor, tol, max_iter, iter_count)
        integer, intent(in) :: nx, ny, max_iter
        real(wp), intent(in) :: dx, dy, omega_sor, tol
        real(wp), dimension(0:nx, 0:ny), intent(inout) :: psi
        real(wp), dimension(0:nx, 0:ny), intent(in) :: omega
        integer, intent(out) :: iter_count
        
        real(wp) :: dx2, dy2, denom, residual, psi_new
        integer :: i, j, iter
        
        dx2 = dx * dx
        dy2 = dy * dy
        denom = 2.0_wp * (dx2 + dy2)
        
        do iter = 1, max_iter
            residual = 0.0_wp
            
            do j = 1, ny-1
                do i = 1, nx-1
                    psi_new = ((psi(i+1,j) + psi(i-1,j)) * dy2 + &
                              (psi(i,j+1) + psi(i,j-1)) * dx2 + &
                              omega(i,j) * dx2 * dy2) / denom
                    
                    psi_new = psi(i,j) + omega_sor * (psi_new - psi(i,j))
                    residual = max(residual, abs(psi_new - psi(i,j)))
                    psi(i,j) = psi_new
                end do
            end do
            
            if (residual < tol) then
                iter_count = iter
                return
            end if
        end do
        
        iter_count = max_iter
        print *, "Warning: SOR did not converge. Residual = ", residual
    end subroutine solve_poisson_sor
    
    
    ! Advance vorticity equation in time
    subroutine advance_vorticity(omega, u, v, nx, ny, dx, dy, dt, Re)
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: dx, dy, dt, Re
        real(wp), dimension(0:nx, 0:ny), intent(inout) :: omega
        real(wp), dimension(0:nx, 0:ny), intent(in) :: u, v
        
        real(wp), dimension(0:nx, 0:ny) :: omega_new
        real(wp) :: dwdx, dwdy, d2wdx2, d2wdy2
        real(wp) :: dx2, dy2, inv_re
        integer :: i, j
        
        dx2 = dx * dx
        dy2 = dy * dy
        inv_re = 1.0_wp / Re
        
        ! Explicit time stepping for vorticity transport equation
        
        do j = 1, ny-1
            do i = 1, nx-1
                ! Convective terms (upwind for stability)
                if (u(i,j) >= 0.0_wp) then
                    dwdx = (omega(i,j) - omega(i-1,j)) / dx
                else
                    dwdx = (omega(i+1,j) - omega(i,j)) / dx
                end if
                
                if (v(i,j) >= 0.0_wp) then
                    dwdy = (omega(i,j) - omega(i,j-1)) / dy
                else
                    dwdy = (omega(i,j+1) - omega(i,j)) / dy
                end if
                
                ! Diffusive terms (central difference)
                d2wdx2 = (omega(i+1,j) - 2.0_wp*omega(i,j) + omega(i-1,j)) / dx2
                d2wdy2 = (omega(i,j+1) - 2.0_wp*omega(i,j) + omega(i,j-1)) / dy2
                
                ! Update
                omega_new(i,j) = omega(i,j) + dt * ( &
                    -u(i,j) * dwdx - v(i,j) * dwdy + &
                    inv_re * (d2wdx2 + d2wdy2))
            end do
        end do
        
        ! Update interior points
        omega(1:nx-1, 1:ny-1) = omega_new(1:nx-1, 1:ny-1)
        
    end subroutine advance_vorticity
    
    
    ! Compute velocity components from stream function
    subroutine compute_velocity(psi, u, v, nx, ny, dx, dy, u_lid)
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: dx, dy, u_lid
        real(wp), dimension(0:nx, 0:ny), intent(in) :: psi
        real(wp), dimension(0:nx, 0:ny), intent(out) :: u, v
        
        integer :: i, j
        
        do j = 1, ny-1
            do i = 0, nx
                u(i,j) = (psi(i,j+1) - psi(i,j-1)) / (2.0_wp * dy)
            end do
        end do
        
        do j = 0, ny
            do i = 1, nx-1
                v(i,j) = -(psi(i+1,j) - psi(i-1,j)) / (2.0_wp * dx)
            end do
        end do
        
        ! Boundaries (no-slip)
        u(0,:) = 0.0_wp
        u(nx,:) = 0.0_wp
        u(:,0) = 0.0_wp
        u(:,ny) = u_lid  ! Top lid velocity
        
        v(:,0) = 0.0_wp
        v(:,ny) = 0.0_wp
        v(0,:) = 0.0_wp
        v(nx,:) = 0.0_wp
        
    end subroutine compute_velocity
    
    
    ! Compute pressure from stream function and vorticity
    subroutine compute_pressure(p, psi, nx, ny, dx, dy, omega_p, tol, max_iter)
        integer, intent(in) :: nx, ny, max_iter
        real(wp), intent(in) :: dx, dy, omega_p, tol
        real(wp), dimension(0:nx, 0:ny), intent(in) :: psi
        real(wp), dimension(0:nx, 0:ny), intent(out) :: p
        
        real(wp), dimension(0:nx, 0:ny) :: rhs
        real(wp) :: d2psidx2, d2psidy2, d2psidxdy
        real(wp) :: dx2, dy2, denom, p_new, residual, p_mean
        integer :: i, j, iter
        
        dx2 = dx * dx
        dy2 = dy * dy
        denom = 2.0_wp * (dx2 + dy2)
        
        ! Compute RHS of pressure Poisson equation
        do j = 1, ny-1
            do i = 1, nx-1
                d2psidx2 = (psi(i+1,j) - 2.0_wp*psi(i,j) + psi(i-1,j)) / dx2
                d2psidy2 = (psi(i,j+1) - 2.0_wp*psi(i,j) + psi(i,j-1)) / dy2
                
                d2psidxdy = (psi(i+1,j+1) - psi(i+1,j-1) - psi(i-1,j+1) + psi(i-1,j-1)) / &
                           (4.0_wp * dx * dy)
                
                rhs(i,j) = 2.0_wp * (d2psidx2 * d2psidy2 - d2psidxdy**2)
            end do
        end do
        
        ! Initialize pressure field
        p = 0.0_wp
        
        ! Solve ∇²p = rhs using SOR with Neumann boundary conditions
        do iter = 1, max_iter
            residual = 0.0_wp
            
            do j = 1, ny-1
                do i = 1, nx-1
                    p_new = ((p(i+1,j) + p(i-1,j)) * dy2 + &
                            (p(i,j+1) + p(i,j-1)) * dx2 - &
                            rhs(i,j) * dx2 * dy2) / denom
                    
                    p_new = p(i,j) + omega_p * (p_new - p(i,j))
                    residual = max(residual, abs(p_new - p(i,j)))
                    p(i,j) = p_new
                end do
            end do
            
            ! Apply Neumann boundary conditions
            p(0,:) = p(1,:)
            p(nx,:) = p(nx-1,:)
            p(:,0) = p(:,1)
            p(:,ny) = p(:,ny-1)
            
            if (residual < tol) exit
        end do
        
        if (iter >= max_iter) then
            print *, "Warning: Pressure solver did not converge. Residual = ", residual
        end if
        
        ! Remove mean pressure
        p_mean = sum(p(1:nx-1, 1:ny-1)) / real((nx-1)*(ny-1), wp)
        p = p - p_mean
        
    end subroutine compute_pressure
    
end module solver
