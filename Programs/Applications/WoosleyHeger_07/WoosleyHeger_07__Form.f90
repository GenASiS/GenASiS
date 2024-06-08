module WoosleyHeger_07__Form

  use GenASiS

  implicit none
  private

!  type, public, extends ( Universe_R_CC_Form ) :: WoosleyHeger_07_Form
  type, public, extends ( Universe_R_CC_C_Form ) :: WoosleyHeger_07_Form
  contains
    procedure, public, pass :: &
      SetFluid
    final :: &
      Finalize
  end type WoosleyHeger_07_Form

    private :: &
      PrepareInterpolation

    integer ( KDI ), private, parameter :: &
      iRADIUS_TS            = 2, &  !-- must match the profile file columns
      iRADIAL_VELOCITY_TS   = 3, &
      iDENSITY_TS           = 4, &
      iTEMPERATURE_TS       = 5, &
      iSPECIFIC_ENERGY_TS   = 10, &
      iELECTRON_FRACTION_TS = 11
    integer ( KDI ), private, parameter :: &
      iRADIAL_VELOCITY_I   = 1, &  !-- interpolation
      iDENSITY_I           = 2, &
      iTEMPERATURE_I       = 3, &
      iSPECIFIC_ENERGY_I   = 4, &
      iELECTRON_FRACTION_I = 5

contains


  impure elemental subroutine Finalize ( WH )
    
    type ( WoosleyHeger_07_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine SetFluid ( WH )

    class ( WoosleyHeger_07_Form ), intent ( inout ) :: &
      WH

    integer ( KDI ) :: &
      iV!, &  !-- iValue
!       iD     !-- iDimension
    real ( KDR ) :: &
      MD, &  !-- MassDensity
      SE     !-- SpecificEnergy
    type ( InterpolationForm ), dimension ( 5 ) :: &
      I_WH
    
    call Show ( 'Setting WoosleyHeger_07 fluid' )

    call PrepareInterpolation ( I_WH )

    select type ( I  =>  WH % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_HN_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( G  =>  F % Geometry, &
        C  =>  A % Chart_GS )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    associate &
      (     N => FV ( :, F % BARYON_DENSITY_C ), &
            T => FV ( :, F % TEMPERATURE ), &
          V_1 => FV ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => FV ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => FV ( :, F % VELOCITY_U ( 3 ) ), &
            E => FV ( :, F % ENERGY_DENSITY_C ), &
          Y_E => FV ( :, F % ELECTRON_FRACTION ), &
            R => GV ( :, G % CENTER_U ( 1 ) ) )

    do iV = 1, size ( N )

      if ( .not. C % ProperCell ( iV ) ) &
        cycle

      call I_WH ( iDENSITY_I ) % Evaluate ( R ( iV ), MD ) 
      call I_WH ( iTEMPERATURE_I ) % Evaluate ( R ( iV ) , T ( iV ) ) 
      call I_WH ( iRADIAL_VELOCITY_I ) % Evaluate ( R ( iV ), V_1 ( iV ) ) 
      call I_WH ( iSPECIFIC_ENERGY_I ) % Evaluate ( R ( iV ), SE ) 
      call I_WH ( iELECTRON_FRACTION_I ) % Evaluate ( R ( iV ), Y_E ( iV ) ) 

      N ( iV )  =  MD  /  CONSTANT % ATOMIC_MASS_UNIT
      E ( iV )  =  SE  *  MD

    end do !-- iV

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

!     do iD = 1, PS % nDimensions
!       associate &
!         ( iaI => PS % Connectivity % iaInner ( iD ), &
!           iaO => PS % Connectivity % iaOuter ( iD ) )
!       call PS % ApplyBoundaryConditions ( F, iD, iaI )
!       call PS % ApplyBoundaryConditions ( F, iD, iaO )
!       end associate !-- iaI, etc.
!     end do !-- iD
    
    end associate !-- N, etc.
    end associate !-- FV, etc.
    end associate !-- G, etc.
    end select !-- A
    end select !-- F
    end select !-- I
     
  end subroutine SetFluid


  subroutine PrepareInterpolation ( I_WH )

    type ( InterpolationForm ), dimension ( 5 ), intent ( inout ) :: &
      I_WH

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Slope_N, &
      Slope_T, &
      Slope_SE, &
      Slope_Y_E
    real ( KDR ), dimension ( : ), allocatable :: &
      RC, &   !-- RadiusCenter
      dRC, &  !-- WidthCenter
      Radius, &
      RadialVelocity, &
      Density, &
      Temperature, &
      SpecificEnergy, &
      ElectronFraction
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Profile
    character ( LDF ) :: &
      Path, &
      Filename
    type ( TableStreamForm ) :: &
      TS

    call Show ( 'Preparing Interpolation' )

    Path = '../Parameters/'
    Filename = 'WH07_S12_08.d.stripped'

    call TS % Initialize &
           ( Filename, PROGRAM_HEADER % Communicator % Rank, &
             PathOption = Path )
    call TS % Read ( Profile, oRowOption = 2 )

    !-- Set "edge" values

    associate &
      (   R => Profile ( :, iRADIUS_TS ), &             !-- cell outer edge
        V_R => Profile ( :, iRADIAL_VELOCITY_TS ), &    !-- cell outer edge
          N => Profile ( :, iDENSITY_TS ), &            !-- cell center
          T => Profile ( :, iTEMPERATURE_TS ), &        !-- cell center
        Y_E => Profile ( :, iELECTRON_FRACTION_TS ), &  !-- cell center
         SE => Profile ( :, iSPECIFIC_ENERGY_TS ), &    !-- cell center
        nProfile => size ( Profile, dim = 1 ) )

    allocate ( Radius ( nProfile + 1 ) )
    allocate ( RadialVelocity ( nProfile + 1 ) )
    allocate ( dRC ( nProfile ) )
    allocate ( RC ( nProfile ) )
    Radius ( 1 )          =  0.0_KDR
    RadialVelocity ( 1 )  =  0.0_KDR
    do iV = 2, nProfile + 1
      Radius         ( iV )  =  R ( iV - 1 )
      RadialVelocity ( iV )  =  V_R ( iV - 1 )
      dRC ( iV - 1 )  =  Radius ( iV )  -  Radius ( iV - 1 )
      RC  ( iV - 1 )  =  Radius ( iV - 1 )  +  0.5_KDR * dRC ( iV - 1 )
    end do

    allocate ( Density ( nProfile + 1 ) )
    allocate ( Temperature ( nProfile + 1 ) )
    allocate ( SpecificEnergy ( nProfile + 1 ) )
    allocate ( ElectronFraction ( nProfile + 1 ) )

    !-- First edge extrapolated
    Slope_N    =  ( N ( 2 )  -  N ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_T    =  ( T ( 2 )  -  T ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_SE   =  ( SE ( 2 )  -  SE ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_Y_E  =  ( Y_E ( 2 )  -  Y_E ( 1 ) )  &
                  /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )

    Density          ( 1 )  =  N ( 1 )  &
                               +  Slope_N   * ( Radius ( 1 )  -  RC ( 1 ) )
    Temperature      ( 1 )  =  T ( 1 )  &
                               +  Slope_T   * ( Radius ( 1 )  -  RC ( 1 ) )
    SpecificEnergy   ( 1 )  =  SE ( 1 )  &
                               +  Slope_SE  * ( Radius ( 1 )  -  RC ( 1 ) )
    ElectronFraction ( 1 )  =  Y_E ( 1 )  &
                               +  Slope_Y_E * ( Radius ( 1 )  -  RC ( 1 ) )

    do iV = 2, nProfile + 1

      if ( iV <= nProfile ) then
        Slope_N    =  ( N ( iV )  -  N ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_T    =  ( T ( iV )  -  T ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_SE   =  ( SE ( iV )  -  SE ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_Y_E  =  ( Y_E ( iV )  -  Y_E ( iV - 1 ) )  &
                      /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
      else
        !-- Last edge extrapolated with same slope
      end if

      Density ( iV )  &
        =  N   ( iV - 1 )  +  Slope_N   * ( Radius ( iV )  -  RC ( iV - 1 ) )
      Temperature ( iV )  &
        =  T   ( iV - 1 )  +  Slope_T   * ( Radius ( iV )  -  RC ( iV - 1 ) )
      SpecificEnergy ( iV )  &
        =  SE  ( iV - 1 )  +  Slope_SE  * ( Radius ( iV )  -  RC ( iV - 1 ) )
      ElectronFraction ( iV )  &
        =  Y_E ( iV - 1 )  +  Slope_Y_E * ( Radius ( iV )  -  RC ( iV - 1 ) )

    end do !-- iV

    end associate !-- R, etc.

    call Show ( 'First few values' )
    call Show ( Profile ( 1 : 5, iRADIUS_TS ), 'RadiusTable' )
    call Show ( Radius ( 1 : 5 ), 'RadiusEdge' )
    call Show ( Profile ( 1 : 5, iRADIAL_VELOCITY_TS ), 'RadialVelocityTable' )
    call Show ( RadialVelocity ( 1 : 5 ), 'RadialVelocityEdge' )
    call Show ( Profile ( 1 : 5, iDENSITY_TS ), 'DensityTable' )
    call Show ( Density ( 1 : 5 ), 'DensityEdge' )
    call Show ( Profile ( 1 : 5, iTEMPERATURE_TS ), 'TemperatureTable' )
    call Show ( Temperature ( 1 : 5 ), 'TemperatureEdge' )
    call Show ( Profile ( 1 : 5, iSPECIFIC_ENERGY_TS ), 'SpecificEnergyTable' )
    call Show ( SpecificEnergy ( 1 : 5 ), 'SpecificEnergyEdge' )
    call Show ( Profile ( 1 : 5, iELECTRON_FRACTION_TS ), &
                'ElectronFractionTable' )
    call Show ( ElectronFraction ( 1 : 5 ), 'ElectronFractionEdge' )
    
    call Show ( UNIT % CENTIMETER, 'UNIT % CENTIMETER' )
    call Show ( UNIT % SECOND    , 'UNIT % SECOND ' )
    call Show ( UNIT % CENTIMETER / UNIT % SECOND, 'Centimeter per second' )
    
    Radius         =  Radius          *  UNIT % CENTIMETER
    Density        =  Density         *  UNIT % MASS_DENSITY_CGS
    Temperature    =  Temperature     *  UNIT % KELVIN
    !-- FIXME: Extraneous parenthesis needed to avoid segfault with XL
    SpecificEnergy =  SpecificEnergy  *  ( UNIT % ERG / UNIT % GRAM )
    RadialVelocity =  RadialVelocity  *  ( UNIT % CENTIMETER / UNIT % SECOND )

    !-- Interpolation initialization
    
    call Show ( 'Initializing Interpolation', CONSOLE % INFO_5 )
    
    call I_WH ( iRADIAL_VELOCITY_I ) % Initialize &
           ( Radius, RadialVelocity, VerbosityOption = CONSOLE % INFO_3 )
    call I_WH ( iDENSITY_I ) % Initialize &
           ( Radius, Density, VerbosityOption = CONSOLE % INFO_3 )
    call I_WH ( iTEMPERATURE_I ) % Initialize &
           ( Radius, Temperature, VerbosityOption = CONSOLE % INFO_3 )
    call I_WH ( iSPECIFIC_ENERGY_I ) % Initialize &
           ( Radius, SpecificEnergy, VerbosityOption = CONSOLE % INFO_3 )
    call I_WH ( iELECTRON_FRACTION_I ) % Initialize &
           ( Radius, ElectronFraction, VerbosityOption = CONSOLE % INFO_3 )

  end subroutine PrepareInterpolation


end module WoosleyHeger_07__Form
