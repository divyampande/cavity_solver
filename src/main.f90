program lid_driven_cavity
    use input
    use solver
    use boundary
    implicit none
    
    type(simulation_params) :: params
    
    ! Grid variables
    integer :: nx, ny, i, j
    real(wp) :: dx, dy
    real(wp), allocatable, dimension(:,:) :: psi, omega, u, v, p
    real(wp), allocatable, dimension(:,:) :: omega_old
    
    ! Solver variables
    integer :: time_step, iter_count
    real(wp) :: time, residual, max_residual
    
    ! Read input parameters
    call read_input_file('cavity.in', params)
    
    nx = params%nx
    ny = params%ny
    
    ! Allocate arrays
    allocate(psi(0:nx, 0:ny))
    allocate(omega(0:nx, 0:ny))
    allocate(omega_old(0:nx, 0:ny))
    allocate(u(0:nx, 0:ny))
    allocate(v(0:nx, 0:ny))
    allocate(p(0:nx, 0:ny))
    
    ! Initialize grid
    dx = 1.0_wp / real(nx, wp)
    dy = 1.0_wp / real(ny, wp)
    
    print *, "Grid spacing: dx =", dx, ", dy =", dy
    print *
    
    ! Initialize fields
    psi = 0.0_wp
    omega = 0.0_wp
    omega_old = 0.0_wp
    u = 0.0_wp
    v = 0.0_wp
    p = 0.0_wp
    
    ! Apply initial boundary conditions
    call apply_psi_bc(psi, nx, ny)
    call apply_omega_bc(omega, psi, nx, ny, dx, dy, params%u_lid)
    
    print *, "Starting time integration..."
    print *, "============================================"
    
    ! Main time-stepping loop
    time = 0.0_wp
    do time_step = 1, params%max_time_steps
        omega_old = omega
        
        call solve_poisson_sor(psi, omega, nx, ny, dx, dy, params%omega_sor, &
                              params%tol_psi, params%max_iter_psi, iter_count)
        call apply_psi_bc(psi, nx, ny)
        
        call compute_velocity(psi, u, v, nx, ny, dx, dy, params%u_lid)
        
        call advance_vorticity(omega, u, v, nx, ny, dx, dy, params%dt, params%Re)
        
        call apply_omega_bc(omega, psi, nx, ny, dx, dy, params%u_lid)
        
        ! Update time
        time = time + params%dt
        
        ! Check convergence to steady state
        max_residual = 0.0_wp
        do j = 0, ny
            do i = 0, nx
                residual = abs(omega(i,j) - omega_old(i,j))
                max_residual = max(max_residual, residual)
            end do
        end do
        
        ! Output progress
        if (mod(time_step, params%output_interval) == 0) then
            print '(A,I8,A,F10.4,A,E12.4,A,I6)', &
                "Step:", time_step, "  Time:", time, &
                "  Max Δω:", max_residual, "  Poisson iter:", iter_count
        end if
        
        ! Check for steady state
        if (max_residual < params%tol_steady) then
            print *
            print *, "============================================"
            print *, "Steady state reached!"
            print *, "Time step:          ", time_step
            print *, "Physical time:      ", time
            print *, "Max vorticity change:", max_residual
            print *, "============================================"
            exit
        end if
        
    end do
    
    if (time_step >= params%max_time_steps) then
        print *
        print *, "============================================"
        print *, "Warning: Maximum time steps reached"
        print *, "Solution may not be fully converged"
        print *, "============================================"
    end if
    
    ! Compute final pressure field
    ! Compute final pressure field using parameters from input
    call compute_pressure(p, psi, nx, ny, dx, dy, &
                          params%omega_sor_p, &
                          params%tol_pressure, &
                          params%max_iter_pressure)
    
    ! Write output files
    print *
    print *, "Writing output files to results/ directory..."
    
    call write_field('results/psi.dat', psi, nx, ny, dx, dy)
    call write_field('results/omega.dat', omega, nx, ny, dx, dy)
    call write_field('results/u.dat', u, nx, ny, dx, dy)
    call write_field('results/v.dat', v, nx, ny, dx, dy)
    call write_field('results/p.dat', p, nx, ny, dx, dy)
    
    call write_centerline_profiles(u, v, nx, ny, dx, dy)
    
    print *, "Simulation completed successfully!"
    print *
    
    ! Deallocate arrays
    deallocate(psi, omega, omega_old, u, v, p)
    
contains
    
    ! Write 2D field to file
    subroutine write_field(filename, field, nx, ny, dx, dy)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: dx, dy
        real(wp), dimension(0:nx, 0:ny), intent(in) :: field
        
        integer :: unit_num, i, j
        real(wp) :: x, y
        
        open(newunit=unit_num, file=filename, status='replace', action='write')
        
        ! Write header
        write(unit_num, '(A,I6)') '# nx = ', nx+1
        write(unit_num, '(A,I6)') '# ny = ', ny+1
        write(unit_num, '(A)') '# x y value'
        
        ! Write data
        do j = 0, ny
            do i = 0, nx
                x = i * dx
                y = j * dy
                write(unit_num, '(3E20.10)') x, y, field(i,j)
            end do
            write(unit_num, *)
        end do
        
        close(unit_num)
    end subroutine write_field
    
    
    ! Write centerline velocity profiles for validation
    subroutine write_centerline_profiles(u, v, nx, ny, dx, dy)
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: dx, dy
        real(wp), dimension(0:nx, 0:ny), intent(in) :: u, v
        
        integer :: unit_num, i, j, i_mid, j_mid
        real(wp) :: x, y
        
        i_mid = nx / 2  ! x = 0.5
        j_mid = ny / 2  ! y = 0.5
        
        ! u-velocity along vertical centerline (x = 0.5)
        open(newunit=unit_num, file='results/u_centerline.dat', status='replace')
        write(unit_num, '(A)') '# y u(x=0.5,y)'
        do j = 0, ny
            y = j * dy
            write(unit_num, '(2E20.10)') y, u(i_mid, j)
        end do
        close(unit_num)
        
        ! v-velocity along horizontal centerline (y = 0.5)
        open(newunit=unit_num, file='results/v_centerline.dat', status='replace')
        write(unit_num, '(A)') '# x v(x,y=0.5)'
        do i = 0, nx
            x = i * dx
            write(unit_num, '(2E20.10)') x, v(i, j_mid)
        end do
        close(unit_num)
        
    end subroutine write_centerline_profiles
    
end program lid_driven_cavity
