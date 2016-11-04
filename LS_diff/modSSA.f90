module stochastic
  use deterministic
  integer::bin_diff

contains

  !****************************
  !*                          *
  !* ssa System of Cell Cycle *
  !*                          *
  !****************************
  subroutine ssa_alg()
    implicit none

    integer :: i, ita, bin_ita, inf_ita, j,k
    double precision :: sum, sum_a
    double precision:: random, propensity
    
    !***************************************
    !choose reaction channel
    !***************************************
    call random_number(random)
    sum = a0 * random
    !print *, 'sum:',sum
    ita=1
    sum_a=a(ita)
    do while(sum_a<sum)
      ita = ita + 1
      sum_a = sum_a + a(ita)
    end do

    !***************************************
    !update population
    !***************************************
    call population_update(ita)


    do i=2, SSADEPEND(1,ita)+1
      !cout<<", i:"<<i<<", ssadepend:"<<SSADEPEND(ita,i)<<", type:"<<NETWORK(ita,TYPE)<<endl
      inf_ita = SSADEPEND(i,ita)
      a0 = a0 - a(inf_ita)
      call propensity_update_ssa(inf_ita, propensity)
      a(inf_ita) = propensity
      a0 = a0 + a(inf_ita)
    end do
  end subroutine ssa_alg


  !****************************
  !*                          *
  !* Initialize ssa           *
  !*                          *
  !****************************
  subroutine init_ssa()
    implicit none
    integer :: i, j

    do i=1, SNUM_ODE
        p(i) = y(i)
    end do

    do i=SNUM_ODE+1, SNUM
        call initZero(i)
    end do

  end subroutine init_ssa

 
  subroutine initZero(ita)
    implicit none
    integer :: i
    integer, intent(in) :: ita

    p(ita) = 0.0
  end subroutine initZero
  
  subroutine clearFire()
    implicit none
    integer :: i

    do i=1, 5
      fire(i) = 0
    end do 
  end subroutine clearFire
  

  !/* update propensity influneced by SSA */
  subroutine propensity_update_ssa(ita, propensity)
    implicit none
    integer :: react_type, react1_ita, react2_ita, aux_ita
    double precision :: rt, temp, avg, tempH, propH
    integer, intent(in) :: ita
    double precision, intent(out) :: propensity

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    aux_ita = NETWORK(AUX, ita)
    rt = RATE(ita)
    temp = 0.0

    select case (react_type)

        case (syn)
            propensity= rt 
            
        case (deg)
            propensity= rt * p(react1_ita)

        case (hill)
            avg = p(aux_ita)
            auxt = avg*avg*avg*avg
            propensity = rt*(auxt/(auxt+Kmdl)) + K_o

        case default
            print *, "Reaction_type in propensity_update_ssa: ", react_type, " not found "
            stop -1

    end select

  end subroutine propensity_update_ssa


  subroutine population_update(ita)
    implicit none
    integer :: react_type, react1_ita, react2_ita, prod1_ita, prod2_ita, aux_ita
    double precision :: ADD = 1.0, REMOVE = -1.0
    integer, intent(in) :: ita
    integer:: i,j

    react_type = NETWORK(RTYPE, ita)
    react1_ita = NETWORK(REACT1, ita)
    react2_ita = NETWORK(REACT2, ita)
    prod1_ita = NETWORK(PROD1, ita)
    prod2_ita = NETWORK(PROD2, ita)
    aux_ita = NETWORK(AUX, ita)
    bin_diff = 0
    fire(react_type) = fire(react_type) + 1;

    select case (react_type)

      case (syn) ! null -> A 
        p(prod1_ita) = p(prod1_ita) + ADD
        y(prod1_ita) = y(prod1_ita) + ADD

      case (deg) ! A -> null
        p(react1_ita) = p(react1_ita) + REMOVE
        if (react1_ita .LE. SNUM_ODE) then
          y(react1_ita) = y(react1_ita) + REMOVE
        end if

      case (hill)
        p(prod1_ita) = p(prod1_ita) + ADD
       

      case default
        print *, "Reaction type in population_update: ", react_type, " not found !"
        stop -1
    
    end select


  end subroutine population_update




end module stochastic

