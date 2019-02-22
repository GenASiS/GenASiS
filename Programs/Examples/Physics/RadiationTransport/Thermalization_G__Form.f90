module Thermalization_G__Form

  !-- Thermalization_Generic_Form

  use GenASiS
  use InteractionsExamples_ASC__Form
  use InteractionsExamples_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( RadiationBoxForm ) :: Thermalization_G_Form
  contains
    procedure, private, pass :: &
      Initialize_T
    generic, public :: &
      Initialize => Initialize_T
  end type Thermalization_G_Form

    private :: &
      InitializeRadiationBox, &
      InitializeInteractions, &
!      InitializeDiagnostics, &
      SetFluid, &
      SetRadiation

    real ( KDR ), private :: &
      TemperatureMin, &
      TemperatureMax, &
      OpacityAbsorption, &
      TimeScale

contains


  subroutine Initialize_T ( T, MomentsType, Name )

    class ( Thermalization_G_Form ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name


    if ( T % Type == '' ) &
      T % Type = 'a Thermalization_G'

!    Thermalization => T

    !-- Parameters

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

    !-- Initialization

    call InitializeRadiationBox ( T, MomentsType, Name )
    call InitializeInteractions ( T )
!    call InitializeDiagnostics ( T )
    call SetFluid ( T )
    call SetRadiation ( T )

  end subroutine Initialize_T


  subroutine InitializeRadiationBox ( T, MomentsType, Name )

    class ( Thermalization_G_Form ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    call T % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'GENERIC' ], &
             MomentsType = MomentsType, &
             Name = Name, &
             ApplyStreamingOption = .false., &
             EvolveFluidOption = .false., &
             FinishTimeOption = 10.0_KDR * TimeScale )
!    T % Integrator % SetReference => SetReference

  end subroutine InitializeRadiationBox


  subroutine InitializeInteractions ( T )

    class ( Thermalization_G_Form ), intent ( inout ) :: &
      T

    select type ( I => T % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )

      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )

      select type ( PS => I % PositionSpace )
      class is ( Atlas_SC_Form )

      allocate ( InteractionsExamples_ASC_Form :: T % Interactions_ASC )
      select type ( IA => T % Interactions_ASC )
      class is ( InteractionsExamples_ASC_Form )
      call IA % Initialize &
             ( PS, InteractionsType = 'CONSTANT', MomentsType = 'GREY' )
      call IA % Set ( FA, OpacityAbsorption = OpacityAbsorption )
      call RMA % SetInteractions ( IA )
      end select !-- IA

      end select !-- PS
      end select !-- FA
      end select !-- RMA

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      class is ( RadiationMoments_BSLL_ASC_CSLD_Form )

      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )

      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )

      allocate ( InteractionsExamples_BSLL_ASC_CSLD_Form :: &
                   T % Interactions_BSLL_ASC_CSLD )
      select type ( IB => T % Interactions_BSLL_ASC_CSLD )
      class is ( InteractionsExamples_BSLL_ASC_CSLD_Form )
      call IB % Initialize ( MS, InteractionsType = 'CONSTANT' )
      call IB % Set ( FA, OpacityAbsorption = OpacityAbsorption )
      call RMB % SetInteractions ( IB )
      end select !-- IB

      end select !-- PS
      end select !-- FA
      end select !-- RMA

    end select !-- Integrator

  end subroutine InitializeInteractions


  subroutine SetFluid ( T )

    class ( Thermalization_G_Form ), intent ( inout ) :: &
      T

    real ( KDR ) :: &
      R_Min, R_Max, &
      T_Min, T_Max
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( I => T % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    select type ( PS => I % PositionSpace )
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
    end select !-- I

    nullify ( G, F )

  end subroutine SetFluid


  subroutine SetRadiation ( T )

    class ( Thermalization_G_Form ), intent ( inout ) :: &
      T

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iF, &  !-- iFiber
      iE     !-- iEnergy  
    real ( KDR ) :: &
      a, &
      Amplitude, &
      Perturbation
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F
    class ( RadiationMomentsForm ), pointer :: &
      RM

    call InitializeRandomSeed ( PROGRAM_HEADER % Communicator )

    !-- Fluid and Geometry

    select type ( I => T % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )
      select type ( FA => I % Current_ASC )
      class is ( Fluid_ASC_Form )
        F => FA % Fluid_P_I ( )
      end select !-- FA
      select type ( PS => I % PositionSpace )
      class is ( Atlas_SC_Form )
        G => PS % Geometry ( )
      end select !-- PS
    end select !-- I

    !-- Radiation

    select type ( I => T % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      type is ( RadiationMoments_ASC_Form )
        RM => RMA % RadiationMoments ( )
      select type ( PS => I % PositionSpace )
      class is ( Atlas_SC_Form )
      select type ( C => PS % Chart )
      class is ( Chart_SL_Template )

              a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN
      Amplitude  =  0.5_KDR

      !$OMP parallel do private ( iV )
      do iV = 1, RM % nValues
        if ( .not. C % IsProperCell ( iV ) ) &
          cycle
        associate &
          ( J     =>  RM % Value ( iV, RM % COMOVING_ENERGY ), &
            H_1   =>  RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
            H_2   =>  RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
            H_3   =>  RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
            T     =>  F % Value ( iV, F % TEMPERATURE ) )

        J  =  a  *  T ** 4
 
        call random_number ( Perturbation )
        Perturbation = Amplitude * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        J = ( 1.0_KDR + Perturbation ) * J

        call random_number ( Perturbation )
        Perturbation = 0.01_KDR * J * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_1 = Perturbation

        call random_number ( Perturbation )
        Perturbation = 0.01_KDR * J * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_2 = Perturbation

        call random_number ( Perturbation )
        Perturbation = 0.01_KDR * J * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_3 = Perturbation
        
        end associate !-- J, etc.
      end do !-- iV
      !$OMP end parallel do

      call RM % ComputeFromPrimitive ( G )

      end select !-- C
      end select !-- PS
      end select !-- RMA

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      type is ( RadiationMoments_BSLL_ASC_CSLD_Form )
      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )

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
      end select !-- RMB

    end select !-- I
    nullify ( G, F, RM )

  end subroutine SetRadiation


end module Thermalization_G__Form
