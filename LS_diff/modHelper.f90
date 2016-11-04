module helper
  use parameters
  !use ode

  implicit none
  character(len=50):: file_popu, file_reac


contains

  !**********************************
  !*                                *
  !* Set Parameters of The System   *
  !*                                *
  !**********************************  
  subroutine set_evn()
    implicit none 
    integer:: i = 0, sim_num
    character(len=30):: buffer, path, evn
   

    ! set record file path and name
    path = '../Data/'

    call getarg(1, buffer)
    read (buffer, *), sim_num

    write(file_popu, "(A4, I0, A4)") "TraH", sim_num, ".txt"
    file_popu = trim(path) // trim(file_popu)

    write(file_reac, "(A5, I0, A4)") "FireH", sim_num, ".txt"
    file_reac = trim(path) // trim(file_reac)

    print *, '*****************************************'
    print *, ''

  end subroutine set_evn

  subroutine to_lower(str)
     character(*), intent(in out) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do
   end subroutine to_lower


  !=======================!
  !=                     =!
  != Recording Functions =!
  !=                     =!
  !=======================!
  subroutine write_file(te)
    implicit none   
    integer:: i,j
    integer:: mid = 20
    double precision, intent(in):: te 

    open(unit = mid, file = file_popu, position = "append", action = "write")
    write(unit = mid, fmt = *) te, p
    
    close(mid)


  end subroutine write_file

  subroutine write_fire()
    implicit none   
    integer:: lat = 30

    open(unit = lat, file = file_reac, position = "append", action = "write")
    write(unit = lat, fmt = *) fire
    close(lat)

  end subroutine write_fire


end module helper


