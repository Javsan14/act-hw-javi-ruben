program md_act
    implicit none
    integer :: Natoms
    character(len=100) :: filename

    print *, "Please, introduce the name of the file: "
    read(*, *) filename

    open(2, file=filename)
    
    Natoms = read_Natoms(2)

    contains

    integer function read_Natoms(input_file) result(Natoms)
        integer, intent(in) :: input_file      

        read(input_file, *) Natoms
        print *, Natoms

    end function read_Natoms

end program md_act