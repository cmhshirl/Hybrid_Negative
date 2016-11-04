module deterministic
  use parameters
  use rand
  integer, parameter:: sys_size = SNUM_ODE + 1
  integer, parameter:: neq = sys_size
  double precision, dimension(neq) :: y

contains

  !****************************
  !*                          *
  !* Initialize ode           *
  !*                          *
  !****************************
  subroutine init_ode()
    implicit none
    integer :: i, j
    double precision:: random

    do i = 1, SNUM_ODE
	y(i) = 0.0
    end do

    call random_number(random)
    y(neq) = log(random)

  end subroutine init_ode


  !****************************
  !*                          *
  !* ODE System of Cell Cycle *
  !*                          *
  !****************************

  subroutine ode_func( t, y, dy)
    implicit none
    integer :: i, j, ileft, iright
    double precision, intent(in) ::t, y(neq)
    double precision, intent(out)::dy(neq)


    call propensity_update_ode()

    !print *, 't:', t, 'a0:', a0, 'y(neq):', y(neq)

    dy(1) = - K_f*y(1) + K_fr*y(2)
    dy(2) = + K_f*y(1) - K_fr*y(2) - K_f*y(2) + K_fr*y(3)
    dy(3) = + K_f*y(2) - K_fr*y(3) - K_f*y(3) + K_fr*y(4)
    dy(4) = + K_f*y(3) - K_fr*y(4) - K_f*y(4) + K_fr*y(5)
    dy(5) = + K_f*y(4) - K_fr*y(5) 
    dy(neq) = a0



  end subroutine ode_func



  subroutine propensity_update_ode()
    implicit none
    integer :: i, j, ita
    double precision :: propensity

    do i=1, SNUM_ODE
        p(i) = y(i)
    end do

    do i = 1, DPNUM_ODE
        ita = ODEDEPEND(i)
        a0 = a0 - a(ita)
        call cal_propensity_ode(ita, propensity)
        a(ita) = propensity
        a0 = a0 + a(ita)
    end do


  end subroutine propensity_update_ode



  subroutine cal_propensity_ode(ita, propensity)

    implicit none
    integer, intent(in) :: ita
    double precision, intent(out) :: propensity
    integer :: j
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0

    select case (react_type)

      case (deg)
        propensity = rt*p(react1_ita)

      case (hill)
	auxt = p(aux_ita)*p(aux_ita)*p(aux_ita)*p(aux_ita)
        temp = auxt/(auxt+Kmdl)
        propensity = rt*temp + K_o 

      case default
        print *, "Reaction_type in ita: ", ita, " not found !"
        print *, "Reaction_type in cal-propensity: ", react_type, " not found !"
        stop -1
    
    end select


!if (propensity .GT. 500.0)  then
!    print *, 'react_type: ', react_type, ', react1_ita:', react1_ita, ', react2_ita: ', react2_ita, ', aux_ita: ', aux_ita, ', rate: ', rt, ', h: ', h, ', temp: ', temp, ', react1_pop: ', p(react1_ita), ', aux_pop: ', p(aux_ita)
!do j=1, MM
!print *, 'react1_pp: ', pp(j,react1_ita)
!print *, 'aux_pp:    ', pp(j,aux_ita)
!print *, 'prePro:    ', prePro(j,ita)
!end do
!stop -1
!end if

  end subroutine cal_propensity_ode



end module deterministic
