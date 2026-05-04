module input
    implicit none
    
    type :: simulation_params
        integer :: nx, ny                ! Grid points
        integer :: max_iter_psi         ! Max iterations for Poisson solver
        integer :: max_time_steps       ! Max time steps
        integer :: output_interval      ! Output frequency
        real(8) :: Re                   ! Reynolds number
        real(8) :: dt                   ! Time step
        real(8) :: u_lid                ! Lid velocity
        real(8) :: omega_sor            ! SOR relaxation parameter
        real(8) :: tol_psi              ! Tolerance for stream function
        real(8) :: tol_steady           ! Tolerance for steady state
    end type simulation_params
    
    private
    public :: simulation_params, read_input_file
    
contains
    
    subroutine read_input_file(filename, params)
        character(len=*), intent(in) :: filename
        type(simulation_params), intent(out) :: params
        
        integer :: unit_num, ios
        character(len=100) :: line, key
        real(8) :: value
        
        ! Default values
        params%nx = 50
        params%ny = 50
        params%max_iter_psi = 10000
        params%max_time_steps = 100000
        params%output_interval = 1000
        params%Re = 100.0d0
        params%dt = 0.001d0
        params%u_lid = 1.0d0
        params%omega_sor = 1.7d0
        params%tol_psi = 1.0d-6
        params%tol_steady = 1.0d-6
        
        ! Open input file
        open(newunit=unit_num, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Warning: Could not open input file '", trim(filename), "'"
            print *, "Using default parameters"
            return
        end if
        
        ! Read parameters
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! Skip comments and empty lines
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            
            ! Parse key-value pairs
            read(line, *, iostat=ios) key, value
            if (ios /= 0) cycle
            
            select case(trim(key))
                case('nx')
                    params%nx = nint(value)
                case('ny')
                    params%ny = nint(value)
                case('max_iter_psi')
                    params%max_iter_psi = nint(value)
                case('max_time_steps')
                    params%max_time_steps = nint(value)
                case('output_interval')
                    params%output_interval = nint(value)
                case('Re')
                    params%Re = value
                case('dt')
                    params%dt = value
                case('u_lid')
                    params%u_lid = value
                case('omega_sor')
                    params%omega_sor = value
                case('tol_psi')
                    params%tol_psi = value
                case('tol_steady')
                    params%tol_steady = value
                case default
                    print *, "Warning: Unknown parameter '", trim(key), "'"
            end select
        end do
        
        close(unit_num)
        
        ! Print parameters
        print *, "============================================"
        print *, "Simulation Parameters:"
        print *, "============================================"
        print *, "Grid size (nx x ny):        ", params%nx, " x ", params%ny
        print *, "Reynolds number:            ", params%Re
        print *, "Time step:                  ", params%dt
        print *, "Lid velocity:               ", params%u_lid
        print *, "SOR parameter:              ", params%omega_sor
        print *, "Max Poisson iterations:     ", params%max_iter_psi
        print *, "Max time steps:             ", params%max_time_steps
        print *, "Output interval:            ", params%output_interval
        print *, "Stream function tolerance:  ", params%tol_psi
        print *, "Steady state tolerance:     ", params%tol_steady
        print *, "============================================"
        print *
        
    end subroutine read_input_file
    
end module input
