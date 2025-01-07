program md_act
    implicit none
    integer :: i, j, Natoms
    character(len=50) :: input_file
    double precision :: total_V
    double precision, allocatable :: coord(:,:), mass(:), distance(:,:), acceleration(:,:)
    double precision, parameter :: epsilon=0.0661d0, sigma=0.3345d0

    ! Read the value of the filename variable from the user
    print *, "Please, introduce the name of the file: "
    print *, ''
    read(*,'(a)') input_file
    print *, ''
    open(unit=2, file=input_file, status='old', action='read')
    
    Natoms = read_Natoms(2)
    allocate(coord(Natoms,3), mass(Natoms), distance(Natoms, Natoms), acceleration(Natoms, 3))
    call read_molecule(2, Natoms, coord, mass)
    close(2)

    call compute_distances(Natoms, coord, distance)
    total_V = V(epsilon, sigma, Natoms, distance)
    
    call compute_acc(Natoms, coord, mass, distance, acceleration)

    deallocate(coord, mass, distance, acceleration)

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

        print*, acceleration

    end subroutine compute_acc
    
end program md_act