module WoosleyHeger_07__Template

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ), abstract :: &
    WoosleyHeger_07_Template
  contains
    procedure, public, pass :: &
      SetFluid
  end type WoosleyHeger_07_Template

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
      iRADIAL_VELOCITY_SI   = 1, &  !-- spline interpolation
      iDENSITY_SI           = 2, &
      iTEMPERATURE_SI       = 3, &
      iSPECIFIC_ENERGY_SI   = 4, &
      iELECTRON_FRACTION_SI = 5

contains


  subroutine SetFluid ( WH )

    class ( WoosleyHeger_07_Template ), intent ( inout ) :: &
      WH

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      MD, &  !-- MassDensity
      SE     !-- SpecificEnergy
    type ( SplineInterpolationForm ), dimension ( 5 ) :: &
      SI
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_HN_Form ), pointer :: &
      F

    call Show ( 'Setting initial conditions' )

    call PrepareInterpolation ( SI )

    select type ( I => WH % Integrator )
    class is ( Integrator_C_PS_Template )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_HN ( )

    select type ( PS => WH % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (     N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
            T => F % Value ( :, F % TEMPERATURE ), &
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
            E => F % Value ( :, F % INTERNAL_ENERGY ), &
          Y_E => F % Value ( :, F % PROTON_FRACTION ), &
            R => G % Value ( :, G % CENTER_U ( 1 ) ) )

    do iV = 1, size ( N )

      call SI ( iDENSITY_SI ) % Evaluate ( R ( iV ), MD ) 
      call SI ( iTEMPERATURE_SI ) % Evaluate ( R ( iV ) , T ( iV ) ) 
      call SI ( iRADIAL_VELOCITY_SI ) % Evaluate ( R ( iV ), V_1 ( iV ) ) 
      call SI ( iSPECIFIC_ENERGY_SI ) % Evaluate ( R ( iV ), SE ) 
      call SI ( iELECTRON_FRACTION_SI ) % Evaluate ( R ( iV ), Y_E ( iV ) ) 

      N ( iV )  =  MD  /  CONSTANT % ATOMIC_MASS_UNIT
      E ( iV )  =  SE  *  MD

    end do !-- iV

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    call F % ComputeFromPrimitive ( G )

    end associate !-- N, etc.
    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( F, G )
     
  end subroutine SetFluid


  subroutine PrepareInterpolation ( SI )

    type ( SplineInterpolationForm ), dimension ( 5 ), intent ( inout ) :: &
      SI

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

    Radius          =  Radius          *  UNIT % CENTIMETER
    RadialVelocity  =  RadialVelocity  *  UNIT % CENTIMETER / UNIT % SECOND
    Density         =  Density         *  UNIT % MASS_DENSITY_CGS
    Temperature     =  Temperature     *  UNIT % KELVIN
    SpecificEnergy  =  SpecificEnergy  *  UNIT % ERG / UNIT % GRAM

    !-- SplineInterpolation initialization

    call SI ( iRADIAL_VELOCITY_SI ) % Initialize &
           ( Radius, RadialVelocity, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iDENSITY_SI ) % Initialize &
           ( Radius, Density, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iTEMPERATURE_SI ) % Initialize &
           ( Radius, Temperature, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iSPECIFIC_ENERGY_SI ) % Initialize &
           ( Radius, SpecificEnergy, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iELECTRON_FRACTION_SI ) % Initialize &
           ( Radius, ElectronFraction, VerbosityOption = CONSOLE % INFO_3 )

  end subroutine PrepareInterpolation


end module WoosleyHeger_07__Template
