#include "Preprocessor"

submodule ( Poisson_ASC__Form ) Poisson_ASC__Kernel

  use Basics
  use LaplacianMultipoleOld_Template
  
  implicit none

contains


  module procedure CombineMoment_CSL_S_Kernel

    !-- CombineMoment_ChartSingleLevel_Spherical_Kernel

    integer ( KDI ) :: &
      iR, iT, iP, &  !-- iRadius, iTheta, iPhi
      iE, &          !-- iEquation
      iR_G, &        !-- iR_Global
      iFirst, iLast
    real ( KDR ) :: &
      F, &  !-- Fraction
      M_RC_C, M_IC_C, &  !-- M_RC_Center, etc.
      M_RS_C, M_IS_C
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    iFirst  =  1
    iLast   =  nC ( 1 )

    if ( UseDevice ) then

      if ( IsFirstShell ) then

        iFirst  =  iFirst + 1

        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iT, iP, iE, F, M_RC_C, M_IC_C, M_RS_C, M_IS_C )
        do iE  =  1, nE
          do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
            do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )

              !-- first radial shell

                iR  =  oC ( 1 )  +  1
              iR_G  =  oR  -  oC ( 1 )  +  iR
  
                 F  =    ( R_C ( iR, iT, iP )  -  R_I ( iR_G ) )  &
                       / ( R_I ( iR_G + 1 )  -  R_I ( iR_G ) )
  
              M_RC_C  =                  F    *  0.0_KDR  &
                         +   ( 1.0_KDR - F )  *  M_RC ( iR_G,     iE )
  
              M_IC_C  =                  F    *  M_IC ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IC ( iR_G + 1, iE )
  
              M_RS_C  =                  F    *  0.0_KDR  &
                         +   ( 1.0_KDR - F )  *  M_RS ( iR_G    , iE )
  
              M_IS_C  =                  F    *  M_IS ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IS ( iR_G + 1, iE )

              S ( iR, iT, iP, iE )  &
                =  S ( iR, iT, iP, iE )  &
                   +  Delta_M_FourPi  *  (    M_RC_C  *  SH_IC ( iR, iT, iP )  &
                                           +  M_IC_C  *  SH_RC ( iR, iT, iP )  &
                                           +  M_RS_C  *  SH_IS ( iR, iT, iP )  &
                                           +  M_IS_C  *  SH_RS ( iR, iT, iP ) )

            end do !-- iT
          end do !-- iP
        end do !-- iE
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      end if !-- IsFirstShell


      if ( IsLastShell ) then

        iLast  =  iLast - 1

        !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& private ( iT, iP, iE, F, M_RC_C, M_IC_C, M_RS_C, M_IS_C )
        do iE  =  1, nE
          do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
            do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )

              !-- last radial shell
  
                iR  =  oC ( 1 )  +  nC ( 1 )
              iR_G  =  oR  -  oC ( 1 )  +  iR

                 F  =    ( R_C ( iR, iT, iP )  -  R_I ( iR_G ) )  &
                       / ( R_I ( iR_G + 1 )  -  R_I ( iR_G ) )

              M_RC_C  =                  F    *  M_RC ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RC ( iR_G    , iE )

              M_IC_C  =                  F    *  M_IC ( iR_G   ,  iE )  &
                         +   ( 1.0_KDR - F )  *  0.0_KDR

              M_RS_C  =                  F    *  M_RS ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RS ( iR_G    , iE )

              M_IS_C  =                  F    *  M_IS ( iR_G    , iE )  &
                         +   ( 1.0_KDR - F )  *  0.0_KDR

              S ( iR, iT, iP, iE )  &
                =  S ( iR, iT, iP, iE )  &
                   +  Delta_M_FourPi  *  (    M_RC_C  *  SH_IC ( iR, iT, iP )  &
                                           +  M_IC_C  *  SH_RC ( iR, iT, iP )  &
                                           +  M_RS_C  *  SH_IS ( iR, iT, iP )  &
                                           +  M_IS_C  *  SH_RS ( iR, iT, iP ) )

            end do !-- iT
          end do !-- iP
        end do !-- iE
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      end if !-- IsLastShell

      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 4 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iR, iT, iP, iE, F, M_RC_C, M_IC_C, M_RS_C, M_IS_C )
      do iE  =  1, nE
        do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
          do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
            do iR  =  oC ( 1 )  +  iFirst,  oC ( 1 )  +  iLast

              iR_G  =  oR  -  oC ( 1 )  +  iR

                 F  =    ( R_C ( iR, iT, iP )  -  R_I ( iR_G ) )  &
                       / ( R_I ( iR_G + 1 )  -  R_I ( iR_G ) )

              M_RC_C  =                  F    *  M_RC ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RC ( iR_G    , iE )

              M_IC_C  =                  F    *  M_IC ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IC ( iR_G + 1, iE )

              M_RS_C  =                  F    *  M_RS ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RS ( iR_G    , iE )

              M_IS_C  =                  F    *  M_IS ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IS ( iR_G + 1, iE )

              S ( iR, iT, iP, iE )  &
                =  S ( iR, iT, iP, iE )  &
                   +  Delta_M_FourPi  *  (    M_RC_C  *  SH_IC ( iR, iT, iP )  &
                                           +  M_IC_C  *  SH_RC ( iR, iT, iP )  &
                                           +  M_RS_C  *  SH_IS ( iR, iT, iP )  &
                                           +  M_IS_C  *  SH_RS ( iR, iT, iP ) )

            end do !-- iR
          end do !-- iT
        end do !-- iP
      end do !-- iE
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else !-- use host

      if ( IsFirstShell ) then

        iFirst  =  iFirst + 1

        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iT, iP, iE, F, M_RC_C, M_IC_C, M_RS_C, M_IS_C )
        do iE  =  1, nE
          do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
            do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )

              !-- first radial shell

                iR  =  oC ( 1 )  +  1
              iR_G  =  oR  -  oC ( 1 )  +  iR
  
                 F  =    ( R_C ( iR, iT, iP )  -  R_I ( iR_G ) )  &
                       / ( R_I ( iR_G + 1 )  -  R_I ( iR_G ) )
  
              M_RC_C  =                  F    *  0.0_KDR  &
                         +   ( 1.0_KDR - F )  *  M_RC ( iR_G,     iE )
  
              M_IC_C  =                  F    *  M_IC ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IC ( iR_G + 1, iE )
  
              M_RS_C  =                  F    *  0.0_KDR  &
                         +   ( 1.0_KDR - F )  *  M_RS ( iR_G    , iE )
  
              M_IS_C  =                  F    *  M_IS ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IS ( iR_G + 1, iE )

              S ( iR, iT, iP, iE )  &
                =  S ( iR, iT, iP, iE )  &
                   +  Delta_M_FourPi  *  (    M_RC_C  *  SH_IC ( iR, iT, iP )  &
                                           +  M_IC_C  *  SH_RC ( iR, iT, iP )  &
                                           +  M_RS_C  *  SH_IS ( iR, iT, iP )  &
                                           +  M_IS_C  *  SH_RS ( iR, iT, iP ) )

            end do !-- iT
          end do !-- iP
        end do !-- iE
        !$OMP  end parallel do

      end if !-- IsFirstShell


      if ( IsLastShell ) then

        iLast  =  iLast - 1

        !$OMP  parallel do collapse ( 3 ) &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& private ( iT, iP, iE, F, M_RC_C, M_IC_C, M_RS_C, M_IS_C )
        do iE  =  1, nE
          do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
            do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )

              !-- last radial shell
  
                iR  =  oC ( 1 )  +  nC ( 1 )
              iR_G  =  oR  -  oC ( 1 )  +  iR

                 F  =    ( R_C ( iR, iT, iP )  -  R_I ( iR_G ) )  &
                       / ( R_I ( iR_G + 1 )  -  R_I ( iR_G ) )

              M_RC_C  =                  F    *  M_RC ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RC ( iR_G    , iE )

              M_IC_C  =                  F    *  M_IC ( iR_G   ,  iE )  &
                         +   ( 1.0_KDR - F )  *  0.0_KDR

              M_RS_C  =                  F    *  M_RS ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RS ( iR_G    , iE )

              M_IS_C  =                  F    *  M_IS ( iR_G    , iE )  &
                         +   ( 1.0_KDR - F )  *  0.0_KDR

              S ( iR, iT, iP, iE )  &
                =  S ( iR, iT, iP, iE )  &
                   +  Delta_M_FourPi  *  (    M_RC_C  *  SH_IC ( iR, iT, iP )  &
                                           +  M_IC_C  *  SH_RC ( iR, iT, iP )  &
                                           +  M_RS_C  *  SH_IS ( iR, iT, iP )  &
                                           +  M_IS_C  *  SH_RS ( iR, iT, iP ) )

            end do !-- iT
          end do !-- iP
        end do !-- iE
        !$OMP  end parallel do
      end if !-- IsLastShell

      !$OMP  parallel do collapse ( 4 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iR, iT, iP, iE, F, M_RC_C, M_IC_C, M_RS_C, M_IS_C )
      do iE  =  1, nE
        do iP  =  oC ( 3 )  +  1,  oC ( 3 )  +  nC ( 3 )
          do iT  =  oC ( 2 )  +  1,  oC ( 2 )  +  nC ( 2 )
            do iR  =  oC ( 1 )  +  iFirst,  oC ( 1 )  +  iLast

              iR_G  =  oR  -  oC ( 1 )  +  iR

                 F  =    ( R_C ( iR, iT, iP )  -  R_I ( iR_G ) )  &
                       / ( R_I ( iR_G + 1 )  -  R_I ( iR_G ) )

              M_RC_C  =                  F    *  M_RC ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RC ( iR_G    , iE )

              M_IC_C  =                  F    *  M_IC ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IC ( iR_G + 1, iE )

              M_RS_C  =                  F    *  M_RS ( iR_G - 1, iE )  &
                         +   ( 1.0_KDR - F )  *  M_RS ( iR_G    , iE )

              M_IS_C  =                  F    *  M_IS ( iR_G,     iE )  &
                         +   ( 1.0_KDR - F )  *  M_IS ( iR_G + 1, iE )

              S ( iR, iT, iP, iE )  &
                =  S ( iR, iT, iP, iE )  &
                   +  Delta_M_FourPi  *  (    M_RC_C  *  SH_IC ( iR, iT, iP )  &
                                           +  M_IC_C  *  SH_RC ( iR, iT, iP )  &
                                           +  M_RS_C  *  SH_IS ( iR, iT, iP )  &
                                           +  M_IS_C  *  SH_RS ( iR, iT, iP ) )

            end do !-- iR
          end do !-- iT
        end do !-- iP
      end do !-- iE
      !$OMP  end parallel do

    end if !-- UseDevice

  end procedure CombineMoment_CSL_S_Kernel


  module procedure SolveCells_CSL_Kernel

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iR, &  !-- iRadius
      iE     !-- iEquation
    real ( KDR ) :: &
      R, &  !-- Radius
      S     !-- Solution
    real ( KDR ), dimension ( nAngularMomentCells ) :: &
      Zero

    Zero = 0.0_KDR

    !$OMP  parallel do &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
    !$OMP& private ( iC, iR, iE, R, S, SH_RC, SH_IC, SH_RS, SH_IS )
    do iC = 1, nCells

      if ( .not. IsProperCell ( iC ) ) &
        cycle

      call ComputeSolidHarmonicsKernel &
             ( CoordinateSystem, Center ( iC, : ), Origin, RadialEdge, &
               MaxDegree, nDimensions, GridError, SH_RC, SH_IC, SH_RS, SH_IS, &
               R, iR )

      do iE = 1, nEquations

        S = 0.0_KDR

        if ( iR > 1 .and. iR < nRadialCells ) then
          call AssembleSolutionKernel &
                 ( S, M_RC ( :, iR - 1, iE ), M_IC ( :, iR, iE ), &
                   M_RC ( :, iR, iE ), M_IC ( :, iR + 1, iE ), &
                   SH_RC, SH_IC, Delta, RadialEdge ( iR ), &
                   RadialEdge ( iR + 1 ), R, nAngularMomentCells )
        else if ( iR == 1 ) then
          call AssembleSolutionKernel &
                 ( S, Zero, M_IC ( :, iR, iE ), &
                   M_RC ( :, iR, iE ), M_IC ( :, iR + 1, iE ), &
                   SH_RC, SH_IC, Delta, RadialEdge ( iR ), &
                   RadialEdge ( iR + 1 ), R, nAngularMomentCells )
        else if ( iR == nRadialCells ) then
          call AssembleSolutionKernel &
                 ( S, M_RC ( :, iR - 1, iE ), M_IC ( :, iR, iE ), &
                   M_RC ( :, iR, iE ), Zero, &
                   SH_RC, SH_IC, Delta, RadialEdge ( iR ), &
                   RadialEdge ( iR + 1 ), R, nAngularMomentCells )
        end if !-- iR

        if ( MaxOrder > 0 ) then
          if ( iR > 1 .and. iR < nRadialCells ) then
            call AssembleSolutionKernel &
                   ( S, M_RS ( :, iR - 1, iE ), M_IS ( :, iR, iE ), &
                     M_RS ( :, iR, iE ), M_IS ( :, iR + 1, iE ), &
                     SH_RS, SH_IS, Delta, RadialEdge ( iR ), &
                     RadialEdge ( iR + 1 ), R, nAngularMomentCells )
          else if ( iR == 1 ) then
            call AssembleSolutionKernel &
                   ( S, Zero, M_IS ( :, iR, iE ), &
                     M_RS ( :, iR, iE ), M_IS ( :, iR + 1, iE ), &
                     SH_RS ( : ), SH_IS ( : ), Delta, RadialEdge ( iR ), &
                     RadialEdge ( iR + 1 ), R, nAngularMomentCells )
          else if ( iR == nRadialCells ) then
            call AssembleSolutionKernel &
                   ( S, M_RS ( :, iR - 1, iE ), M_IS ( :, iR, iE ), &
                     M_RS ( :, iR, iE ), Zero, &
                     SH_RS ( : ), SH_IS ( : ), Delta, RadialEdge ( iR ), &
                     RadialEdge ( iR + 1 ), R, nAngularMomentCells )
          end if !-- iR
        end if !-- MaxOrder

        Solution ( iC, iaSolution ( iE ) )  =  - S / FourPi

      end do !-- iE

    end do !-- iC
    !$OMP  end parallel do

  end procedure SolveCells_CSL_Kernel 


  module procedure AssembleSolutionKernel

    integer ( KDI ) :: &
      iA  !-- iAngular
    real ( KDR ) :: &
      F, &  !-- Fraction
      M_R, M_I

    F  =  ( R - R_In ) / ( R_Out - R_in )

    do iA = 1, nA

      M_R  =  F  *  M_R_In ( iA )   +   ( 1.0_KDR - F )  *  M_R_Out ( iA )
      M_I  =  F  *  M_I_In ( iA )   +   ( 1.0_KDR - F )  *  M_I_Out ( iA )

      S  =  S   +   Delta ( iA )  &
                    *  ( M_R  *  SH_I ( iA )  +  M_I  *  SH_R ( iA ) )

    end do !-- iA

  end procedure AssembleSolutionKernel


end submodule Poisson_ASC__Kernel
