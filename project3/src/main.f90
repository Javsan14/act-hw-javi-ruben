program md_act
    implicit none
    integer :: i, j, Natoms
    character(len=50) :: input_file
    double precision :: total_V, total_T
    double precision, allocatable :: coord(:,:), mass(:), distance(:, :)
    double precision, parameter :: epsilon=0.0661d0, sigma=0.3345d0

    ! Read the value of the filename variable from the user
    print *, "Please, introduce the name of the file: "
    print *, ''
    read(*,'(a)') input_file
    print *, ''
    open(unit=2, file=input_file, status='old', action='read')
    
    Natoms = read_Natoms(2)
    allocate(coord(Natoms,3), mass(Natoms), distance(Natoms, Natoms))
    call read_molecule(2, Natoms, coord, mass)
    close(2)

    call compute_distances(Natoms, coord, distance)
    total_V = V(epsilon, sigma, Natoms, distance)
    total_T = T(Natoms, reshape([1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0], [3, 3]), mass)

    deallocate(coord, mass, distance)

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
    
end program md_act