module Thermalization_G_S__Form

  !-- Thermalization_Generic_Grey_Form

  use GenASiS
  use Interactions_G_S_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: Thermalization_G_S_Form
  contains
    procedure, public, pass :: &
      Initialize
  end type Thermalization_G_S_Form

    private :: &
      SetFluid, &
      SetRadiation

    real ( KDR ), private :: &
      TemperatureMin, &
      TemperatureMax, &
      OpacityAbsorption

contains


  subroutine Initialize ( T, Name )

    class ( Thermalization_G_S_Form ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      TimeScale


    if ( T % Type == '' ) &
      T % Type = 'a Thermalization_G_S'

!    Thermalization => T

    call T % InitializeTemplate ( Name )

    TemperatureMin  =   0.1_KDR
    TemperatureMax  =  10.0_KDR
    call PROGRAM_HEADER % GetParameter ( TemperatureMin, 'TemperatureMin' )
    call PROGRAM_HEADER % GetParameter ( TemperatureMax, 'TemperatureMax' )

    OpacityAbsorption = 1.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( OpacityAbsorption, 'OpacityAbsorption' )

    associate &
      ( c       => CONSTANT % SPEED_OF_LIGHT, &
        Kappa_A => OpacityAbsorption )
    TimeScale   =  1.0 / ( c * Kappa_A )
    end associate !-- c, etc.


    !-- Integrator

    allocate ( SpectralRadiationBoxForm :: T % Integrator )
    select type ( SRB => T % Integrator )
    type is ( SpectralRadiationBoxForm )
    call SRB % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'GENERIC' ], &
             Name = Name, &
             ApplyStreamingOption = .false., &
             EvolveFluidOption = .false., &
             FinishTimeOption = 10.0_KDR * TimeScale )
!     SRB % SetReference => SetReference

    select type ( PS => SRB % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( MS => SRB % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    select type ( RMB => SRB % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    class is ( RadiationMoments_BSLL_ASC_CSLD_Form )


    !-- Interactions

    allocate ( Interactions_G_S_BSLL_ASC_CSLD_Form :: &
                 SRB % Interactions_BSLL_ASC_CSLD )
    select type ( IB => SRB % Interactions_BSLL_ASC_CSLD )
    class is ( Interactions_G_S_BSLL_ASC_CSLD_Form )
    call IB % Initialize ( MS, 'GENERIC_SPECTRAL' )
    call IB % Set ( OpacityAbsorption = OpacityAbsorption )
    call RMB % SetInteractions ( IB )
    end select !-- IB


    ! !-- Diagnostics

    ! EnergyDensityUnit  =  UNIT % MEGA_ELECTRON_VOLT ** 4 / UNIT % HBAR_C ** 3

    ! allocate ( T % Reference_ASC )
    ! allocate ( T % FractionalDifference_ASC )
    ! call T % Reference_ASC % Initialize &
    !        ( PS, 'GENERIC', NameShortOption = 'Reference', &
    !          EnergyDensityUnitOption = EnergyDensityUnit, &
    !          AllocateSourcesOption = .false., &
    !          IgnorabilityOption = CONSOLE % INFO_5 )
    ! call T % FractionalDifference_ASC % Initialize &
    !        ( PS, 'GENERIC', NameShortOption = 'FractionalDifference', &
    !          AllocateSourcesOption = .false., &
    !          IgnorabilityOption = CONSOLE % INFO_5 )
    ! T % SetReference => SetReference


    !-- Initial conditions

    call SetFluid ( SRB )
    call SetRadiation ( SRB )


    !-- Cleanup

    end select !-- RMB
    end select !-- MS
    end select !-- PS
    end select !-- SRB

  end subroutine Initialize


  subroutine SetFluid ( SRB )

    class ( SpectralRadiationBoxForm ), intent ( inout ) :: &
      SRB

    real ( KDR ) :: &
      R_Min, R_Max, &
      T_Min, T_Max
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( FA => SRB % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    select type ( PS => SRB % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    associate &
      ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
        L  => ( C % MaxCoordinate - C % MinCoordinate ), &
        X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER_U ( 3 ) ), &
        Temperature => F % Value ( :, F % TEMPERATURE ), &
        Density     => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        Velocity_1  => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
        Velocity_2  => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        Velocity_3  => F % Value ( :, F % VELOCITY_U ( 3 ) ) )
    associate &
      ( R => sqrt ( ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
                    +  ( Z - R0 ( 3 ) ) ** 2 ) )

    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = minval ( R )
    CO % Outgoing % Value ( 2 ) = 1.0_KDR / maxval ( R )
    call CO % Reduce ( REDUCTION % MIN )

    R_Min = CO % Incoming % Value ( 1 )
    R_Max = 1.0_KDR / CO % Incoming % Value ( 2 )

    T_Min = TemperatureMin
    T_Max = TemperatureMax

!    Temperature &
!      = T_Max  -  ( T_Max - T_Min ) / ( R_Max - R_Min ) * ( R - R_Min )
    Temperature = T_Max * ( ( R_Min / R ) ** ( log10 ( T_Max / T_Min ) &
                                               / log10 ( R_Max / R_Min ) ) )

    !-- Only temperature is relevant, but we need values for the following
    !   to avoid temperature getting blitzed by ComputeFromConserved
    Density    = 1.0_KDR
    Velocity_1 = 0.0_KDR
    Velocity_2 = 0.0_KDR
    Velocity_3 = 0.0_KDR

    call F % ComputeFromTemperature ( F % Value, G, G % Value )

    end associate !-- R
    end associate !-- R0, etc.
    end select !-- C
    end select !-- PS
    end select !-- FA

    nullify ( G, F )

  end subroutine SetFluid


  subroutine SetRadiation ( SRB )

    class ( SpectralRadiationBoxForm ), intent ( inout ) :: &
      SRB

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iE     !-- iEnergy  
    real ( KDR ) :: &
      Amplitude, &
      Perturbation
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( RadiationMomentsForm ), pointer :: &
      RM

    select type ( RMB => SRB % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    select type ( FA => SRB % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    select type ( MS => SRB % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    G => MS % Base_CSLD % Geometry ( )

    call InitializeRandomSeed ( PROGRAM_HEADER % Communicator )

    !$OMP parallel do private ( iF )
    do iF = 1, MS % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )
      RM => RMB % RadiationMoments ( iF )
      associate &
        ( J     =>  RM % Value ( :, RM % COMOVING_ENERGY ), &
          H_1   =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
          H_2   =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
          H_3   =>  RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
          T     =>  F % Value ( iBC, F % TEMPERATURE ), &
          E     =>  RMB % Energy )

      call SetPlanckSpectrum ( E, T, J )

      Amplitude = 0.9_KDR
      do iE = 1, RMB % nEnergyValues

        call random_number ( Perturbation )
        Perturbation = Amplitude * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        J ( iE ) = ( 1.0_KDR + Perturbation )  *  J ( iE )

        call random_number ( Perturbation )
        Perturbation &
          = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_1 ( iE ) = Perturbation

        call random_number ( Perturbation )
        Perturbation &
          = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_2 ( iE ) = Perturbation

        call random_number ( Perturbation )
        Perturbation &
          = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_3 ( iE ) = Perturbation
        
      end do !-- iE

      call RM % ComputeFromPrimitive ( iBC, G )

      end associate !-- J, etc.
      end associate !-- iBC
    end do !-- iF
    !$OMP end parallel do

    call RMB % LoadSections ( )

    end select !-- MS
    end select !-- FA
    end select !-- RMB
    nullify ( G, F, RM )

  end subroutine SetRadiation


end module Thermalization_G_S__Form
