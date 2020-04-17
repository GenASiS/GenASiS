#include "Preprocessor"

submodule ( LaplacianMultipole_Template ) LaplacianMultipole_Kernel

  use Basics

  implicit none

contains

  
  module procedure ComputeMomentContributionsKernel

    integer ( KDI ) :: &
      iE  !-- iEquation
    real ( KDR ), dimension ( nE ) :: &
      Source_dV

    Source_dV  =  [ ( Source ( iaSource ( iE ) )  *  Volume, &
                      iE = 1, nE ) ] 
!call Show ( Source_dV, 'Source_dV' )

    call ComputeMomentContributions_MR_MI_Kernel &
           ( MyM_RC, MyM_IC, SH_RC, SH_IC, Source_dV, nE, nA, iR )
    if ( L > 0 ) &
      call ComputeMomentContributions_MR_MI_Kernel &
             ( MyM_RS, MyM_IS, SH_RS, SH_IS, Source_dV, nE, nA, iR )

  end procedure ComputeMomentContributionsKernel


  module procedure ComputeSolidHarmonicsKernel

    real ( KDR ) :: &
      X, Y, Z

!    call Show ( '>>> Hello world from ComputeSolidHarmonicsKernel' )

    select case ( trim ( CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      X  =  Position ( 1 )  -  Origin ( 1 )
      Y  =  Position ( 2 )  -  Origin ( 2 )
      Z  =  Position ( 3 )  -  Origin ( 3 )
    case ( 'CYLINDRICAL' )
      if ( nDimensions < 3 ) then
        X  =  Position ( 1 )
        Y  =  0.0_KDR
      else
        X  =  Position ( 1 )  *  cos ( Position ( 3 ) )
        Y  =  Position ( 1 )  *  sin ( Position ( 3 ) )
      end if
      Z  =  Position ( 2 )  -  Origin ( 2 )
    case ( 'SPHERICAL' )
      if ( nDimensions < 3 ) then
        X  =  Position ( 1 )  *  sin ( Position ( 2 ) )
        Y  =  0.0_KDR
      else
        X  =  Position ( 1 )  *  sin ( Position ( 2 ) )  &
                              *  cos ( Position ( 3 ) )
        Y  =  Position ( 1 )  *  sin ( Position ( 2 ) )  &
                              *  sin ( Position ( 3 ) )
      end if
      Z  =  Position ( 1 )  *  cos ( Position ( 2 ) )
    end select !-- CoordinateSystem

    if ( nDimensions < 3 ) then
      call ComputeSolidHarmonics_C_M_0_Kernel ( X, Z, L, R_C, I_C )
    else
      call ComputeSolidHarmonics_C_S_Kernel ( X, Y, Z, L, R_C, I_C, R_S, I_S )
    end if

    R  =  sqrt ( X * X  +  Y * Y  +  Z * Z )

    GridError = .false.
    if ( R > RadialEdge ( size ( RadialEdge ) ) ) then
      call Show ( 'Radial grid not large enough', CONSOLE % ERROR )
      call Show ( R, 'R', CONSOLE % ERROR )
      call Show ( RadialEdge ( size ( RadialEdge ) ), 'R_Max', &
                  CONSOLE % ERROR )
      call Show ( 'ComputeSolidHarmonicsKernel', 'subroutine', &
                  CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_Kernel', 'submodule', CONSOLE % ERROR )
      GridError = .true.
    end if

    call Search ( RadialEdge, R, iR ) 
!call Show ( R, 'R', nLeadingLinesOption = 1 )
!call Show ( iR, 'iR' )
!call Show ( RadialEdge ( iR ), 'R_in' )
!call Show ( RadialEdge ( iR + 1 ), 'R_out' )

  end procedure ComputeSolidHarmonicsKernel


  module procedure ComputeMomentContributions_MR_MI_Kernel

    integer ( KDI ) :: &
      iA, &  !-- iAngular
      iE     !-- iEquation   

    do iE = 1, nE
      do iA = 1, nA
        MyMR ( iA, iRS, iE )  &
          =  MyMR ( iA, iRS, iE )  +  SH_R ( iA )  *  Source_dV ( iE ) 
        MyMI ( iA, iRS, iE )  &
          =  MyMI ( iA, iRS, iE )  +  SH_I ( iA )  *  Source_dV ( iE )
      end do
    end do

  end procedure ComputeMomentContributions_MR_MI_Kernel


  module procedure ComputeSolidHarmonics_C_M_0_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM, &
      iPD   !-- iPreviousDiagonal
    real ( KDR ) :: &
      D_2

    D_2  =  X * X  +  Z * Z 

    iV = 0
    iM = 0

    !-- ( L, M ) = ( iM, iM )
    !-- Note iL = iM
    iV = iV + 1
    if ( iM == 0 ) then
      iV = 1
      R_C ( iV ) = 1.0_KDR
      I_C ( iV ) = 1.0_KDR / sqrt ( D_2 )
    else
      R_C ( iV ) = - ( X * R_C ( iPD ) ) / ( 2 * iM )
      I_C ( iV ) = - ( 2 * iM - 1 ) &
                     * ( X * I_C ( iPD ) ) / D_2
    end if
    iPD = iV

    if ( iM == L ) return

    !-- ( L, M ) = ( iM + 1, iM )
    !-- Note iL = iM + 1
    iV = iV + 1
    R_C ( iV ) = Z * R_C ( iV - 1 )  
    I_C ( iV ) = ( 2 * ( iM + 1 ) - 1 ) * Z * I_C ( iV - 1 ) / D_2  

    do iL = iM + 2, L
      !-- ( L, M ) = ( iL, iM )
      iV = iV + 1
      R_C ( iV ) &
        = ( ( 2 * iL - 1 ) * Z * R_C ( iV - 1 )  -  D_2 * R_C ( iV - 2 ) ) &
          / ( ( iL + iM ) * ( iL - iM ) )
      I_C ( iV ) &
        = ( ( 2 * iL - 1 ) * Z * I_C ( iV - 1 )  &
            -  ( ( iL - 1 )**2 - iM**2 ) * I_C ( iV - 2 ) ) &
          / D_2
    end do !-- iL

  end procedure ComputeSolidHarmonics_C_M_0_Kernel 


  module procedure ComputeSolidHarmonics_C_S_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iL, &
      iM, &
      iPD   !-- iPreviousDiagonal
    real ( KDR ) :: &
      D_2

    D_2  =  X * X  +  Y * Y  +  Z * Z 

    iV = 0
    do iM = 0, L

      !-- ( L, M ) = ( iM, iM )
      !-- Note iL = iM
      iV = iV + 1
      if ( iM == 0 ) then
        iV = 1
        R_C ( iV ) = 1.0_KDR
        R_S ( iV ) = 0.0_KDR
        I_C ( iV ) = 1.0_KDR / sqrt ( D_2 )
        I_S ( iV ) = 0.0_KDR
      else
        R_C ( iV ) = - ( X * R_C ( iPD ) - Y * R_S ( iPD ) ) / ( 2 * iM )
        R_S ( iV ) = - ( Y * R_C ( iPD ) + X * R_S ( iPD ) ) / ( 2 * iM )
        I_C ( iV ) = - ( 2 * iM - 1 ) &
                       * ( X * I_C ( iPD ) - Y * I_S ( iPD ) ) / D_2
        I_S ( iV ) = - ( 2 * iM - 1 ) &
                       * ( Y * I_C ( iPD ) + X * I_S ( iPD ) ) / D_2
      end if
      iPD = iV

      if ( iM == L ) exit

      !-- ( L, M ) = ( iM + 1, iM )
      !-- Note iL = iM + 1
      iV = iV + 1
      R_C ( iV ) = Z * R_C ( iV - 1 )  
      R_S ( iV ) = Z * R_S ( iV - 1 )  
      I_C ( iV ) = ( 2 * ( iM + 1 ) - 1 ) * Z * I_C ( iV - 1 ) / D_2  
      I_S ( iV ) = ( 2 * ( iM + 1 ) - 1 ) * Z * I_S ( iV - 1 ) / D_2 

      do iL = iM + 2, L
        !-- (L,M) = (iL,iM)
        iV = iV + 1
        R_C ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * R_C ( iV - 1 )  -  D_2 * R_C ( iV - 2 ) )&
            / ( ( iL + iM ) * ( iL - iM ) )
        R_S ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * R_S ( iV - 1 )  -  D_2 * R_S ( iV - 2 ) )&
            / ( ( iL + iM ) * ( iL - iM ) )
        I_C ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * I_C ( iV - 1 )  &
              -  ( ( iL - 1 )**2 - iM**2 ) * I_C ( iV - 2 ) ) &
            / D_2
        I_S ( iV ) &
          = ( ( 2 * iL - 1 ) * Z * I_S ( iV - 1 )  &
              -  ( ( iL - 1 )**2 - iM**2 ) * I_S ( iV - 2 ) ) &
            / D_2
      end do !-- iL

    end do !-- iM

  end procedure ComputeSolidHarmonics_C_S_Kernel


end submodule LaplacianMultipole_Kernel
