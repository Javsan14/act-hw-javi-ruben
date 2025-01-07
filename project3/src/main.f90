program md_act
    implicit none
    integer :: i, j
    character(len=50) :: input_file
    double precision, allocatable :: coord(:,:), mass(:)

    ! Read the value of the filename variable from the user
    write(*,*) 'Welcome, please enter the name of the file of the molecule:'
    write(*,*)
    read(*,'(a)') input_file
    write(*,*)
    open(unit=2, file=input_file, status='old', action='read')
    allocate(coord(3,3))
    allocate(mass(3))
    call read_molecule(2, 3, coord, mass)
    close(2)

    write(*,*) coord, mass

    deallocate(coord, mass)

contains

    subroutine read_molecule(input_file, natoms, coord, mass)
        implicit none
        integer :: i, j, ierr
        integer, intent(in) :: input_file, natoms
        double precision, intent(out) :: coord(natoms,3), mass(natoms)
    
        ! Skip first line
        read(input_file, *)
        ! Read data of every atom
        do i = 1, natoms
            read(input_file, *, iostat=ierr) (coord(i,j), j=1,3), mass(i)
            ! Check succesfull reading
            if (ierr /= 0) then
                write(*,*) 'Error: Failed to read atom data at line ', i
                stop
            endif
        enddo

    end subroutine read_molecule

end program md_act