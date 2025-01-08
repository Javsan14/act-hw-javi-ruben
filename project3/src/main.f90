program md_act
    implicit none
    ! Variable declaration
    integer :: i, j, Natoms, total_steps
    character(len=50) :: input_file
    logical :: file_exists
    double precision :: dt, total_V, total_T
    double precision, allocatable :: new_coord(:,:), coord(:,:), mass(:), distance(:,:)
    double precision, allocatable :: new_velocity(:,:), velocity(:,:), acceleration(:,:)
    double precision, parameter :: epsilon=0.0661d0, sigma=0.3345d0

    ! Read the value of the filename variable from the user
    print *, "Please, introduce the name of the file: "
    print *, ''
    read(*,'(a)') input_file
    print *, ''
    
    ! Open the input file and obtain the important info: natoms, masses and coordinates
    open(unit=2, file=input_file, status='old', action='read')
    Natoms = read_Natoms(2)
    ! With the number of atoms we can allocate all the necessary arrays
    allocate(coord(Natoms,3), new_coord(Natoms,3), mass(Natoms), distance(Natoms, Natoms))
    allocate(velocity(Natoms, 3), new_velocity(Natoms, 3), acceleration(Natoms, 3))
    call read_molecule(2, Natoms, coord, mass)
    close(2)

    ! Calculate the initial distances and accelerations
    call compute_distances(Natoms, coord, distance)
    call compute_acc(Natoms, coord, mass, distance, acceleration)
    velocity = 0.0d0 ! The initial velocities are just set to zero

    ! If the output trajectory file already exists, delete it first
    inquire(file="traj.xyz", exist=file_exists)
    if (file_exists) then
        open(99, file="traj.xyz", status="old")
        close(99, status='delete')
    end if

    ! Open the output file to print all the trajectory info
    open(3, file="traj.xyz", status="new")

    ! Assign variables of the MD algorithm
    dt = 0.2d0
    total_steps = 1000
    
    ! Main MD iteration
    do i = 1, total_steps
        
        ! Update coordinates based on the Verlet algorithm
        new_coord = coord + velocity*dt + acceleration*dt**2/2
        new_velocity = velocity + 0.5d0*dt*acceleration

        ! For the new snapshot, recalculate distances and acceleration
        call compute_distances(Natoms, new_coord, distance)
        call compute_acc(Natoms, new_coord, mass, distance, acceleration)
        ! Finish to update the velocities
        new_velocity = new_velocity + 0.5d0*dt*acceleration
        ! Evaluate the potential (V) and kinetic (T) energies of the system
        total_V = V(epsilon, sigma, Natoms, distance)
        total_T = T(Natoms, new_velocity, mass)

        ! Prepare variables for the next step
        coord = new_coord
        velocity = new_velocity

        ! Print info in the output trajectory file every 10 steps
        if (mod(i, 10) == 0) then
            call write_xyz(Natoms, total_V, total_T, coord, mass)
        end if
    enddo

    close(3)
    ! Always deallocate at the end
    deallocate(coord, mass, distance, velocity, acceleration)
    deallocate(new_coord, new_velocity)

contains

    ! Function that reads the Natoms from the input file
    integer function read_Natoms(input_file) result(Natoms)
        implicit none
        integer, intent(in) :: input_file      

        ! Read the number of atoms on the first line
        read(input_file, *) Natoms

    end function read_Natoms
    
    ! Subroutine to read the coordinates and masses of the input molecule
    subroutine read_molecule(input_file, Natoms, coord, mass)
        implicit none
        integer :: i, j, ierr
        integer, intent(in) :: input_file, natoms
        double precision, intent(out) :: coord(Natoms,3), mass(Natoms)
    
        ! Read data of every atom (the first line has been already read by read_Natoms)
        do i = 1, Natoms
            read(input_file, *, iostat=ierr) (coord(i,j), j=1,3), mass(i)
            ! Check succesfull reading
            if (ierr /= 0) then
                write(*,*) 'Error: Failed to read atom data at line ', i
                stop
            endif
        enddo

    end subroutine read_molecule

    ! Subroutine to generate a matrix (Natoms x Natoms) of distances between atoms
    subroutine compute_distances(Natoms, coord, distance)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: coord(Natoms, 3)
        double precision, intent(out) :: distance(Natoms, Natoms)
        double precision :: dist
        integer :: i, j

        ! Initialize the distance
        distance = 0.d0
        do i=1, Natoms
            do j=i+1, Natoms
                dist = norm2(coord(i, :) - coord(j, :)) 
                distance(i, j) = dist
                distance(j, i) = dist
            end do
        end do

    end subroutine
    
    ! Function to compute the potential energy of the system
    double precision function V(epsilon, sigma, Natoms, distance)
        implicit none
        integer :: i
        double precision, intent(in) :: epsilon, sigma
        integer, intent(in) :: Natoms
        double precision, intent(in) :: distance(Natoms,Natoms)

        ! Initialize V
        V = 0.0d0
        do i = 1, Natoms
            do j = i+1, Natoms
                V = V + 4*epsilon*( (sigma/distance(i,j))**12 - &
                (sigma/distance(i,j))**6 ) ! Lennard-Jones potential
            enddo 
        enddo

    end function V
    
    ! Function to compute the kinetic energy of the system
    double precision function T(Natoms, velocity, mass)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: velocity(Natoms, 3), mass(Natoms)

        ! The dot product acts as a fancy sum
        T = 0.5d0 * dot_product(mass, norm2(velocity, dim=1)**2)

    end function T
    
    ! Compute the acceleration for each atom 
    subroutine compute_acc(Natoms, coord, mass, distance, acceleration)
        implicit none
        integer :: i, j
        integer, intent(in) :: Natoms
        double precision, intent(in) :: coord(Natoms,3)
        double precision, intent(in) :: mass(Natoms)
        double precision, intent(in) :: distance(Natoms,Natoms)
        double precision :: U, fx, fy, fz
        double precision, intent(out) :: acceleration(Natoms,3)

        ! Initialize acceleration
        acceleration = 0.0d0

        ! Compute pairwise forces
        do i = 1, Natoms
            do j = i+1, Natoms
                if (distance(i, j) > 0.0d0) then
                    U = 24 * epsilon / distance(i, j) * ((sigma / distance(i, j))**6 - &
                    2*(sigma / distance(i, j))**12)
    
                    ! Force components
                    fx = U*(coord(i, 1) - coord(j, 1))/distance(i, j)
                    fy = U*(coord(i, 2) - coord(j, 2))/distance(i, j)
                    fz = U*(coord(i, 3) - coord(j, 3))/distance(i, j)
    
                    ! Update accelerations for both atoms (i and j)
                    acceleration(i, 1) = acceleration(i, 1) - fx/mass(i)
                    acceleration(i, 2) = acceleration(i, 2) - fy/mass(i)
                    acceleration(i, 3) = acceleration(i, 3) - fz/mass(i)
                    
                    acceleration(j, 1) = acceleration(j, 1) + fx/mass(j)
                    acceleration(j, 2) = acceleration(j, 2) + fy/mass(j)
                    acceleration(j, 3) = acceleration(j, 3) + fz/mass(j)
                end if
            end do
        end do

    end subroutine compute_acc
    
    ! Subroutine to write a snapshot of the MD algorithm in xyz format
    subroutine write_xyz(Natoms, V, T, coord, mass)
        integer, intent(in) :: Natoms
        double precision, intent(in) :: V, T, coord(Natoms, 3), mass(Natoms)
        integer :: i

        write(3, *) Natoms ! First line is the number of atoms
        write(3, *) V, T, V+T ! The comment line is used to put the energies (check conservation)

        ! Next Natoms lines for the symbols and coordinates
        do i=1, Natoms
            write(3, "(A2, 3F10.6)") from_mass_to_symbol(mass(i)), coord(i, :)
        end do

    end subroutine write_xyz

    ! Function to assign a symbol based on the mass of the atom
    character(len=2) function from_mass_to_symbol(mass)
        double precision, intent(in) :: mass

        select case (nint(mass)) ! Turn it into an integer to avoid decimals problems
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