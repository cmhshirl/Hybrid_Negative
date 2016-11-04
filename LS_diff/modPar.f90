module parameters
  implicit none

  !************************
  !*                      *
  !* Parameters for system*
  !*                      *
  !************************
  integer, parameter:: SNUM = 6
  integer, parameter:: SNUM_ODE = 5
  integer, parameter:: RNUM = 8
  integer, parameter:: DPNUM_ODE = 6
  integer, parameter:: DPNUM_SSA = 3

  integer, parameter:: syn = 1
  integer, parameter:: deg = 2
  integer, parameter:: first_order= 3
  integer, parameter:: decomposition = 4
  integer, parameter:: hill = 5

  !************************
  !*                      *
  !* Parameters for SPECIE*
  !*                      *
  !************************
  integer, parameter:: RSPACE = 6
  integer, parameter:: RTYPE  = 1
  integer, parameter:: REACT1 = 2
  integer, parameter:: REACT2 = 3
  integer, parameter:: PROD1  = 4
  integer, parameter:: PROD2  = 5
  integer, parameter:: AUX    = 6
  
  integer, parameter :: E1   = 1
  integer, parameter :: E2   = 2
  integer, parameter :: E3   = 3
  integer, parameter :: E4   = 4
  integer, parameter :: E5   = 5
  integer, parameter :: P1   = 6

  !************************
  !*                      *
  !*Para for ODE RATE     *
  !*                      *
  !************************
  double precision, parameter :: K_f   = 1.0e+0
  double precision, parameter :: K_fr  = 1.1e+0
  double precision, parameter :: K_M   = 2.6e-1
  double precision, parameter :: K_n   = 4.0e+0
  double precision, parameter :: K_o   = 1.0e+0
  double precision, parameter :: K_se  = 2.5e-1
  double precision, parameter :: K_sp  = 5.0e+1
  double precision, parameter :: K_de  = 5.0e-2
  double precision, parameter :: K_dp  = 1.0e+1
  double precision, parameter :: Kmdl  = K_M*K_M*K_M*K_M
  !**********************
  !*                    *
  !* changing variables *
  !*                    *
  !**********************
  double precision :: a0
  double precision, dimension(RNUM) :: a
  double precision :: auxt
  double precision, dimension(SNUM) :: p
  double precision:: mu
  integer(KIND=8), dimension(5) :: fire

  !**********************
  !*                    *
  !* Param for SSA rate *
  !*                    *
  !**********************
  double precision, dimension(RNUM) :: RATE
  
  integer, dimension(DPNUM_ODE), parameter :: ODEDEPEND = (/1, 2, 3, 4, 5, 8/)
    
  integer, dimension(RSPACE, RNUM), parameter :: NETWORK = reshape( (/ &
2, 1, 0,  0,  0,  0, &
2, 2, 0,  0,  0,  0, &
2, 3, 0,  0,  0,  0, &
2, 4, 0,  0,  0,  0, &
2, 5, 0,  0,  0,  0, &
1, 0, 0,  1,  0,  0, &
2, 6, 0,  0,  0,  0, &
5, 0, 0,  6,  0,  5/), (/ RSPACE, RNUM /) )
  
  integer, dimension(DPNUM_SSA,RNUM) :: SSADEPEND


contains



  subroutine init_para()

    RATE = (/ K_de, K_de, K_de, K_de, K_de, K_se, K_dp, K_sp/)


    SSADEPEND(:,1) = (/ 1, 1, 0/)
    SSADEPEND(:,2) = (/ 1, 2, 0/)
    SSADEPEND(:,3) = (/ 1, 3, 0/)
    SSADEPEND(:,4) = (/ 1, 4, 0/)
    SSADEPEND(:,5) = (/ 2, 5, 8/)
    SSADEPEND(:,6) = (/ 1, 1, 0/)
    SSADEPEND(:,7) = (/ 1, 7, 0/)
    SSADEPEND(:,8) = (/ 1, 7, 0/)

  end subroutine init_para

end module parameters
