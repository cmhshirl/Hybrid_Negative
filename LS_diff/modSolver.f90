module solver
  use deterministic

  !integer, parameter:: neq = sys_size
  integer, parameter:: ng = 1
  integer:: istate, itask
  integer:: liw, lrw, jt, itol, iopt
  integer:: lrn, lrs
  integer, allocatable, dimension(:):: iwork, jroot
  double precision:: rtol, atol
  double precision, allocatable, dimension(:):: rwork
  !double precision, parameter:: delta = 1e-11
  external dlsodar

contains
  subroutine fex(neq, t, y, ydot)
    implicit none 
    integer::neq
    double precision :: t, y(neq), ydot(neq)

    call ode_func(t, y, ydot)

  end subroutine fex

  subroutine gex(neq, t, y, ng, gout)
    implicit none 
    integer, intent(in):: ng
    integer, intent(in):: neq
    double precision, intent(in):: t, y(neq)
    double precision, intent(inout)::gout(ng)

    gout(ng) = y(neq)
  end subroutine gex

  subroutine jac(neq, t, y, ydot)    
  end subroutine jac

  subroutine init_solver_parameters()
    implicit none

    lrn = 20 + 16*neq + 3*ng
    lrs = 22 + 9*neq + neq*neq + 3*ng
    lrw = max(lrn, lrs)
    liw = 20 + neq
    allocate(jroot(ng), iwork(liw), rwork(lrw))

    itol  = 1    ! atol: 1--scalar, 2--array
    rtol  = 1d-3
    atol  = 1d-6
    itask = 1    ! input: 1--normal, 2--one step, 3,4,5...
    iopt  = 0    ! input: 0--no optional inputs, 1--optional inputs
    jt    = 2
    !IWORK(6) = 1000

    !----istate
    ! input:  1--first call for the problem
    !         2--not first call, continue with tout and itask changes
    !         3--not first call, continue with parameters change except above
    ! output: 1--nothing was done
    !         2--integration was performed successfully, no root
    !         3--integration successful, roots were found

    !----jroot
    ! output: jroot(i) means g(i) = 1

  end subroutine init_solver_parameters

  subroutine integrate(ts, te, istate, jroot)
    implicit none 
    integer, intent(inout):: istate 
    integer, dimension(:), intent(out):: jroot
    double precision, intent(in):: ts, te
    call dlsodar(fex, neq, y, ts, te, itol, rtol, &
                 atol, itask, istate, iopt, rwork, lrw, &
                 iwork, liw, jac, jt, gex, ng, jroot)

    if (istate == 1) then
       print *, "Nothing to integrate!"
    elseif (istate < 0) then
       print *, "Error happens with id", istate
    end if

  end subroutine integrate

  subroutine deallocate_parameters()
    implicit none 
    deallocate(jroot, iwork, rwork)
  end subroutine deallocate_parameters

end module solver
