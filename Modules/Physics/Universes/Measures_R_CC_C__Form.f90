module Measures_R_CC_C__Form

  !-- Measures_Radiation_CentralCore_Collected__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Measures_F_CC__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_MEASURES_NM = 2

  type, public, extends ( Measures_F_CC_Form ) :: Measures_R_CC_C_Form
    integer ( KDI ) :: &
      N_MEASURES_NM = N_MEASURES_NM
    integer ( KDI ) :: &
      nRadiations, &
      nMeasures_R
    real ( KDR ), dimension ( : ), allocatable :: &
      Luminosity, &
      EnergyAverage
    type ( FieldSet_BM_Pointer ), dimension ( : ), allocatable :: &
      Radiation_SA_1D
    class ( Units_R_Form ), pointer :: &
      Units_R => null ( )
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Measures_R_CC_C_Form

    type ( Real_1D_Form ), dimension ( : ), allocatable, private :: &
      H, &
      EA
    type ( CollectiveOperation_R_Form ), allocatable, private :: &
      CO


contains


  subroutine Initialize_R_CC_C &
               ( M, R_SA_1D, F, F_SA, G_SA, A_SA, Units_R, Units_F )

    class ( Measures_R_CC_C_Form ), intent ( inout ) :: &
      M
    class ( FieldSet_BM_Pointer ), dimension ( : ), intent ( in ) :: &
      R_SA_1D
    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      F, &
      F_SA, &
      G_SA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A_SA
    class ( Units_R_Form ), intent ( in ), target :: &
      Units_R
    class ( Units_F_Form ), intent ( in ), target :: &
      Units_F

    integer ( KDI ) :: &
      oV, &
      iR

    M % nRadiations  =  size ( R_SA_1D )
    M % nMeasures_R  =  M % nRadiations  *  M % N_MEASURES_NM

    call M % Measures_F_CC_Form % Initialize &
           ( F, F_SA, G_SA, A_SA, Units_F, &
             nMeasuresAddOption = M % nMeasures_R )

    M % Units_R  =>  Units_R

    associate &
      ( UR   =>  M % Units_R, &
        nMF  =>  M % nMeasures_F, &
        nMN  =>  M % N_MEASURES_NM, &
        nR   =>  M % nRadiations )

    allocate ( M % Luminosity ( nR ) )
    allocate ( M % EnergyAverage ( nR ) )

    allocate ( M % Radiation_SA_1D ( nR ) )
    do iR  =  1, nR

      M % Radiation_SA_1D ( iR ) % Pointer  =>  R_SA_1D ( iR ) % Pointer
      associate ( R  =>  M % Radiation_SA_1D ( iR ) % Pointer )

      oV  =  nMF  +  ( iR - 1 ) * nMN

      M % Name ( oV + 1 )  =  'Luminosity_' // trim ( R % Name )
      M % Name ( oV + 2 )  =  'EnergyAve_'  // trim ( R % Name )

      M % Unit ( oV + 1 )  =  UR % Luminosity
      M % Unit ( oV + 2 )  =  UR % EnergyAverage
  
      end associate !-- R

    end do !-- iR

    end associate !-- UR, etc.

  end subroutine Initialize_R_CC_C


  subroutine Compute ( M )

    class ( Measures_R_CC_C_Form ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      iB, &  !-- iBrick
      iC, &  !-- iCell
      iM, &  !-- iMeasure
      oC, &  !-- oCell
      oI, &  !-- oIncoming
      oV, &
      iR, &   !-- iRadiation
      oF, &   !-- oField
      nFR, &  !-- nFieldsRadiation
      nF      !-- nFields
    real ( KDR ) :: &
      FourPi_C
    real ( KDR ), dimension ( :, : ), pointer :: &
      Outgoing_2D, &
      Incoming_2D

    call M % Measures_F_CC_Form % Compute ( )

    FourPi_C  =  4.0_KDR  *  CONSTANT % PI  *  CONSTANT % SPEED_OF_LIGHT

    select type ( A_SA  =>  M % Atlas_SA )
      class is ( Atlas_SCG_Form )
    associate &
      ( C_SA    =>  A_SA % Chart_GS )
    associate &
      ( nGL   =>  C_SA % nGhostLayers ( 1 ), &
        nC    =>  C_SA % nCells ( 1 ), &
        nCB   =>  C_SA % nCellsBrick ( 1 ), &
        nCBG  =>  C_SA % nCellsBrickGlobal ( 1 ), &
        nB    =>  C_SA % nBricks ( 1 ) )
    associate &
      ( nMF  =>  M % nMeasures_F, &
        nMN  =>  M % N_MEASURES_NM, &
        nR   =>  M % nRadiations )
        
    !-- Gather spherically averaged radiation fields
    !   (assume decomposition in spherical shells)

    nFR  =  2
    nF   =  nFR * nR

    if ( .not. allocated ( CO ) ) then
      allocate ( CO )
      call CO % Initialize &
             ( C_SA % Communicator, &
               nOutgoing = [ nF * nCB ], nIncoming = nF * nCBG % Value )
    end if

    Outgoing_2D ( 1 : nCB, 1 : nF )  =>  CO % Outgoing % Value
    do iR  =  1, nR
      select type ( R_SA  =>  M % Radiation_SA_1D ( iR ) % Pointer )
        class is ( NeutrinoMoments_G_Form )
      associate &
        ( R_SA_V  =>  R_SA % Storage_GS % Value )
      associate &
        (  H_P  =>  R_SA_V ( :, R_SA % MOMENTUM_DENSITY_C_U_1 ), &
          EA_P  =>  R_SA_V ( :, R_SA % ENERGY_AVERAGE ) )
      oF  =  ( iR - 1 ) * nFR
      Outgoing_2D ( 1 : nCB, oF + 1 )  =   H_P ( nGL + 1 : nGL + nCB )
      Outgoing_2D ( 1 : nCB, oF + 2 )  =  EA_P ( nGL + 1 : nGL + nCB )
      end associate !-- H_P, etc.
      end associate !-- R_SA_V
      end select !-- R_SA
    end do !-- iR

    call CO % Gather_V ( )

    if ( .not. allocated ( H ) ) then
      allocate ( H ( nR ), EA ( nR ) )
      do iR  =  1, nR
        call  H ( iR ) % Initialize ( nC )
        call EA ( iR ) % Initialize ( nC )
      end do !-- iR
    end if

    do iB  =  1,  nB
      if ( iB == 1 ) then
        oC = 0
      else
        oC  = oC + nCBG % Value ( iB - 1 )
      end if
      oI  =  oC * nF
      associate ( nCBG_V => nCBG % Value ( iB ) )
      Incoming_2D ( 1 : nCBG_V, 1 : nF )  &
        =>  CO % Incoming % Value ( oI + 1 : oI + nCBG_V * nF )
      do iR  =  1, nR
        oF  =  ( iR - 1 ) * nFR
         H ( iR ) % Value ( oC + 1 : oC + nCBG_V )  &
           =  Incoming_2D ( : , oF + 1 )
        EA ( iR ) % Value ( oC + 1 : oC + nCBG_V )  &
           =  Incoming_2D ( : , oF + 2 )
      end do !--iR
      end associate   !-- nCBG_V
    end do !-- iB

    !-- Values at 500 km

    do iC  =  1, nC
      if ( M % Radius ( iC )  >  500.0_KDR * UNIT % KILOMETER ) then
        do iR  =  1, nR
          M % Luminosity ( iR )  &
            =  H ( iR ) % Value ( iC )  *  FourPi_C  *  M % Radius ( iC ) ** 2
          M % EnergyAverage ( iR )  &
            =  EA ( iR ) % Value ( iC )
        end do !-- iR
        exit
      end if
    end do

    !-- Record

    do iR  =  1, nR
      oV  =  nMF  +  ( iR - 1 ) * nMN
      M % Value ( oV + 1 )  =  M % Luminosity ( iR )
      M % Value ( oV + 2 )  =  M % EnergyAverage ( iR ) 
    end do !-- iR

    !-- Display

    call Show ( 'Radiation_CentralCore Measures' )
    do iM  =  nMF  +  1,  nMF  +  M % nMeasures_R
      call Show ( M % Value ( iM ), M % Unit ( iM ), M % Name ( iM ) )
    end do !-- iM

    !-- Cleanup

    end associate !-- nMF, etc.
    end associate !-- nGL, etc.
    end associate !-- C_SA
    end select !-- A_SA

  end subroutine Compute


  impure elemental subroutine Finalize ( M )

    type ( Measures_R_CC_C_Form ), intent ( inout ) :: &
      M

    if ( allocated ( M % Radiation_SA_1D ) ) &
      deallocate ( M % Radiation_SA_1D )
    if ( allocated ( M % EnergyAverage ) ) &
      deallocate ( M % EnergyAverage )
    if ( allocated ( M % Luminosity ) ) &
      deallocate ( M % Luminosity )

    if ( allocated ( CO ) ) &
      deallocate ( CO )
    if ( allocated ( EA ) ) &
      deallocate ( EA )
    if ( allocated ( H ) ) &
      deallocate ( H )

    nullify ( M % Units_R )

  end subroutine Finalize


end module Measures_R_CC_C__Form
