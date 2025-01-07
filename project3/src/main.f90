program md_act
    implicit none
    integer :: i, j, Natoms, total_steps
    character(len=50) :: input_file
    double precision :: dt, total_V, total_T
    double precision, allocatable :: new_coord(:,:), coord(:,:), mass(:), distance(:,:)
    double precision, allocatable :: new_velocity(:,:), velocity(:,:), acceleration(:,:)
    double precision, parameter :: epsilon=0.0661d0, sigma=0.3345d0

    ! Read the value of the filename variable from the user
    print *, "Please, introduce the name of the file: "
    print *, ''
    read(*,'(a)') input_file
    print *, ''
    open(unit=2, file=input_file, status='old', action='read')
    
    Natoms = read_Natoms(2)
    allocate(coord(Natoms,3), new_coord(Natoms,3), mass(Natoms), distance(Natoms, Natoms))
    allocate(new_velocity(Natoms, 3), velocity(Natoms, 3), acceleration(Natoms, 3))
    call read_molecule(2, Natoms, coord, mass)
    close(2)

    call compute_distances(Natoms, coord, distance)
    velocity = 0.0d0
    call compute_acc(Natoms, coord, mass, distance, acceleration)

    open(3, file="traj.xyz", status="new")

    dt = 0.2d0
    total_steps = 1000
    
    do i = 1, total_steps
        new_coord = coord + velocity*dt + acceleration*dt**2/2
        new_velocity = velocity + 0.5d0*dt*acceleration
        call compute_distances(Natoms, new_coord, distance)
        call compute_acc(Natoms, new_coord, mass, distance, acceleration)
        new_velocity = new_velocity + 0.5d0*dt*acceleration
        total_V = V(epsilon, sigma, Natoms, distance)
        total_T = T(Natoms, new_velocity, mass)

        coord = new_coord
        velocity = new_velocity

        print*, i
        print*, total_V
        print*, total_T
        print*, total_T - total_V
        print*, ''
    enddo

    close(3)
    deallocate(coord, mass, distance, velocity, acceleration)
    deallocate(new_coord, new_velocity)

contains

    integer function read_Natoms(input_file) result(Natoms)
        implicit none
        integer, intent(in) :: input_file      

        read(input_file, *) Natoms

    end function read_Natoms
    
    subroutine read_molecule(input_file, Natoms, coord, mass)
        implicit none
        integer :: i, j, ierr
        integer, intent(in) :: input_file, natoms
        double precision, intent(out) :: coord(Natoms,3), mass(Natoms)
    
        ! Read data of every atom
        do i = 1, Natoms
            read(input_file, *, iostat=ierr) (coord(i,j), j=1,3), mass(i)
            ! Check succesfull reading
            if (ierr /= 0) then
                write(*,*) 'Error: Failed to read atom data at line ', i
                stop
            endif
        enddo

    end subroutine read_molecule

    subroutine compute_distances(Natoms, coord, distance)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: coord(Natoms, 3)
        double precision, intent(out) :: distance(Natoms, Natoms)
        double precision :: dist
        integer :: i, j

        distance = 0.d0
        do i=1, Natoms
            do j=i+1, Natoms
                dist = norm2(coord(i, :) - coord(j, :))
                distance(i, j) = dist
                distance(j, i) = dist
            end do
        end do

    end subroutine
    
    double precision function V(epsilon, sigma, Natoms, distance)
        implicit none
        integer :: i
        double precision, intent(in) :: epsilon, sigma
        integer, intent(in) :: Natoms
        double precision, intent(in) :: distance(Natoms,Natoms)

        V = 0.0d0
        do i = 1, Natoms
            do j = i+1, Natoms
                V = V + 4*epsilon*( (sigma/distance(i,j))**12 - &
                (sigma/distance(i,j))**6 )
            enddo 
        enddo

    end function V
    
    double precision function T(Natoms, velocity, mass)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: velocity(Natoms, 3), mass(Natoms)
        integer :: i

        !T = 0.d0
        !do i=1, Natoms
        !    T = T + mass(i)*norm2(velocity(i, :))**2
        !end do
        !print *, T*0.5d0

        T = 0.d0
        T = 0.5d0 * dot_product(mass, norm2(velocity, dim=1)**2)

    end function T
    
    subroutine compute_acc(Natoms, coord, mass, distance, acceleration)
        implicit none
        integer :: i, j
        integer, intent(in) :: Natoms
        double precision, intent(in) :: coord(Natoms,3)
        double precision, intent(in) :: mass(Natoms)
        double precision, intent(in) :: distance(Natoms,Natoms)
        double precision :: U, sum_x, sum_y, sum_z
        double precision, intent(out) :: acceleration(Natoms,3)

        do i = 1, Natoms
            sum_x = 0.0d0
            sum_y = 0.0d0
            sum_z = 0.0d0
            do j = i+1, Natoms
                U = 24*epsilon/distance(i,j)*( (sigma/distance(i,j))**6 - &
                2*(sigma/distance(i,j))**6 )
                sum_x = sum_x + U * (coord(i, 1) - coord(j, 1))/distance(i,j)
                sum_y = sum_y + U * (coord(i, 2) - coord(j, 2))/distance(i,j)
                sum_z = sum_z + U * (coord(i, 3) - coord(j, 3))/distance(i,j)
            enddo
            acceleration(i, 1) = -1.0d0/mass(i) * sum_x
            acceleration(i, 2) = -1.0d0/mass(i) * sum_y
            acceleration(i, 3) = -1.0d0/mass(i) * sum_z
        enddo

    end subroutine compute_acc
    
    subroutine write_xyz(Natoms, V, T, coord, mass)
        integer, intent(in) :: Natoms
        double precision, intent(in) :: V, T, coord(Natoms, 3), mass(Natoms)
        integer :: i

        write(3, *) Natoms
        write(3, *) V, T, V+T

        do i=1, Natoms
            write(3, "(A2, 3F8.4)") from_mass_to_symbol(mass(i)), coord(i, :)
        end do

    end subroutine write_xyz

    character(len=2) function from_mass_to_symbol(mass)
        double precision, intent(in) :: mass

        select case (nint(mass))
        case (1)
            from_mass_to_symbol = "H"
        case (4)
            from_mass_to_symbol = "He"
        case (7)
            from_mass_to_symbol = "Li"
        case (9)
            from_mass_to_symbol = "Be"
        case (11)
            from_mass_to_symbol = "B"
        case (12)
            from_mass_to_symbol = "C"
        case (14)
            from_mass_to_symbol = "N"
        case (16)
            from_mass_to_symbol = "O"
        case (19)
            from_mass_to_symbol = "F"
        case (20)
            from_mass_to_symbol = "Ne"
        case (23)
            from_mass_to_symbol = "Na"
        case (24)
            from_mass_to_symbol = "Mg"
        case (27)
            from_mass_to_symbol = "Al"
        case (28)
            from_mass_to_symbol = "Si"
        case (31)
            from_mass_to_symbol = "P"
        case (32)
            from_mass_to_symbol = "S"
        case (35)
            from_mass_to_symbol = "Cl"
        case (40)
            from_mass_to_symbol = "Ar"
        case default
            from_mass_to_symbol = "??"
    end select

    end function from_mass_to_symbol

end program md_act