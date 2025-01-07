program md_act
    implicit none
    integer :: i, j, Natoms
    character(len=50) :: input_file
    double precision, allocatable :: coord(:,:), mass(:)

    ! Read the value of the filename variable from the user
    print *, "Please, introduce the name of the file: "
    print *,
    read(*,'(a)') input_file
    print *,
    open(unit=2, file=input_file, status='old', action='read')
    
    Natoms = read_Natoms(2)
    allocate(coord(3,3), mass(3))
    call read_molecule(2, Natoms, coord, mass)
    close(2)

    deallocate(coord, mass)

contains

    integer function read_Natoms(input_file) result(Natoms)
        integer, intent(in) :: input_file      

        read(input_file, *) Natoms
        print *, Natoms

    end function read_Natoms
    
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