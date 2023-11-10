module Thermalization_Form

  use GenASiS
  use Interactions_C__Form

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: ThermalizationForm
    real ( KDR ) :: &
      TemperatureMin, &
      TemperatureMax, &
      TimeScale
    real ( KDR ), dimension ( 2 ) :: &
      OpacityAbsorption
    type ( PhotonMoments_G_Form ), allocatable :: &
      FractionalDifference
  contains
    procedure, private, pass :: &
      Initialize_T
    generic, public :: &
      Initialize => Initialize_T
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type ThermalizationForm

    private :: &
      InitializeUniverse, &
      InitializeDiagnostics, &
      SetInitial, &
      SetReference, &
      Compute_dT_Local

      private :: &
        SetFluid, &
        SetRadiation

        private :: &
          SetFluidKernel, &
          SetPerturbationKernel

contains


  subroutine Initialize_T ( T, FormalismType, Name )

    class ( ThermalizationForm ), intent ( inout ), target :: &
      T
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( T % Type  ==  '' ) &
      T % Type  =  'a Thermalization'

    call InitializeUniverse ( T, FormalismType, Name )
    call InitializeDiagnostics ( T )
 
  end subroutine Initialize_T


  impure elemental subroutine Finalize ( T )

    type ( ThermalizationForm ), intent ( inout ), target :: &
      T

    if ( allocated ( T % FractionalDifference ) ) &
      deallocate ( T % FractionalDifference )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( ThermalizationForm ), intent ( in ) :: &
      U

    call U % Universe_R_B_Form % ShowParameters ( )

    call Show ( U % TemperatureMin,    'TemperatureMin' )
    call Show ( U % TemperatureMax,    'TemperatureMax' )
    call Show ( U % OpacityAbsorption, 'OpacityAbsorption' )
    call Show ( U % TimeScale,         'TimeScale' )

  end subroutine ShowParameters


  subroutine InitializeUniverse ( T, FormalismType, Name )

    class ( ThermalizationForm ), intent ( inout ), target :: &
      T
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    integer ( KDI ) :: &
      iD

    allocate ( Interactions_C_Form :: T % Interactions_BM )

    call T % Initialize &
           ( RadiationName = [ 'Radiation_1', 'Radiation_2' ], &
             RadiationType = [ 'PHOTONS', 'PHOTONS' ], &
             FormalismType = FormalismType, &
             Name = Name, &
             ApplyStreamingOption = .false., &
             EvolveFluidOption = .false., &
             nCellsPositionOption = [ 128, 128, 128 ] )

    select type ( I  =>  T % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    select type ( I  =>  T % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    associate &
      ( R  =>  I % CurrentSet_X_1D )
    do iD  =  1, 3
      call R % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- R
    end select !-- I
             
    T % Integrator % SetInitial    =>  SetInitial
    T % Integrator % SetReference  =>  SetReference
    T % Integrator % System        =>  T
    
    T % Integrator % Compute_dT_Local  =>  Compute_dT_Local
    
  end subroutine InitializeUniverse


  subroutine InitializeDiagnostics ( T )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T

    character ( LDL ) :: &
      DifferenceName

    allocate &
      ( T % FractionalDifference )
    associate &
      ( R_FD  =>  T % FractionalDifference, &
        G     =>  T % Integrator % Geometry_X, &
        S     =>  T % Integrator % Checkpoint_X )

    select case ( T % iRadiation )
    case ( 1 )
      DifferenceName  =  'FractionalDifference_1'
    case ( 2 )
      DifferenceName  =  'FractionalDifference_2'
    end select
    
    call R_FD % Initialize ( G, T % Units_R, NameOption = DifferenceName )
    call R_FD % SetStream ( S )

    end associate !-- R_R, etc.

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( T  =>  I % System )
      class is ( ThermalizationForm )
    select type ( I )
      class is ( Integrator_CS_1D_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )

    T % TemperatureMin  =   1.0_KDR
    T % TemperatureMax  =  10.0_KDR
    call PROGRAM_HEADER % GetParameter ( T % TemperatureMin, 'TemperatureMin' )
    call PROGRAM_HEADER % GetParameter ( T % TemperatureMax, 'TemperatureMax' )

    T % OpacityAbsorption  =  [ 1.0_KDR, 2.0_KDR ]
    call PROGRAM_HEADER % GetParameter &
           ( T % OpacityAbsorption, 'OpacityAbsorption' )

    associate &
      ( c        =>  CONSTANT % SPEED_OF_LIGHT, &
        Kappa_A  =>  minval ( T % OpacityAbsorption ) )
    T % TimeScale   =  1.0 / ( c * Kappa_A )
    end associate !-- c, etc.

    I % T_Finish  =  10.0_KDR  *  T % TimeScale

    call InitializeRandomSeed ( I % Communicator )

    !-- Fluid

    call SetFluid ( T, F )
    call F % SetUseInitialTemperature ( .true. )

    !-- Radiation

    select type ( I )
    class is ( Integrator_CS_1D_BM_CS_Form )

      select type ( R  =>  I % CurrentSet_X_1D )
        class is ( PhotonMoments_G_Form )

      call SetRadiation ( T, R, F )

      end select !-- R

    class default
      call Show ( 'Integrator type not recognized', CONSOLE % ERROR )
      call Show ( 'ThermalizationForm', 'module', CONSOLE % ERROR )
      call Show ( 'SetInitial', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I    

    !-- Interactions

    select type ( Intrctns  =>  T % Interactions_BM )
    class is ( Interactions_C_Form )
       call Intrctns % SetOpacityAbsorption &
              ( T % OpacityAbsorption ( T % iRadiation ) )
    end select !-- Intrctns
    
    end select !-- F
    end select !-- I
    end select !-- T

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( T  =>  I % System )
      class is ( ThermalizationForm )
    select type ( I  =>  T % Integrator )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( PhotonMoments_G_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )
    select type ( A  =>  R % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( R_FD  =>  T % FractionalDifference )
    associate &
      ( C      =>  A % Chart_GS, &
        FV     =>  F    % Storage_GS % Value, &
        RV     =>  R    % Storage_GS % Value, &
        R_FDV  =>  R_FD % Storage_GS % Value )

    call R % ComputeSpectralParameters ( )
    call R % ComputeEquilibrium ( )

    call ComputeFractionalDifferenceKernel &
           ( J_FD = R_FDV ( :, R % ENERGY_DENSITY_C ), &
             J    = RV    ( :, R % ENERGY_DENSITY_C ), &
             J_R  = RV    ( :, R % ENERGY_DENSITY_C_EQ ), & 
             ProperCell = C % ProperCell )

    end associate !-- C, etc.
    end associate !-- R_R, etc.
    end select !-- A
    end select !-- F
    end select !-- R
    end select !-- I
    end select !-- T

  end subroutine SetReference


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option  )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( I )
      class is ( Integrator_CS_1D_BM_CS_Form )
    select type ( T => I % System )
      class is ( ThermalizationForm )

    dT_Candidate  =  huge ( 1.0_KDR )

    associate &
      ( c        =>  CONSTANT % SPEED_OF_LIGHT, &
        Kappa_A  =>  T % OpacityAbsorption ( T % iRadiation ) )

    dT_Candidate ( 3 )  &
      =  T % InteractionFactor  *  1.0 / ( c * Kappa_A )

    end associate !-- c, etc.

    !-- Reduce across CS_1D

    call CO % Initialize &
           ( I % Communicator_X_1D, nOutgoing = [ 1 ], &
             nIncoming = [ 1 ] )

    CO % Outgoing % Value ( 1 )  =  I % dT_Candidate ( 3 )
    call CO % Reduce ( REDUCTION % MIN )
    I % dT_Candidate ( 3 )  =  CO % Incoming % Value ( 1 )

    end select !-- I
    end select !-- T

  end subroutine Compute_dT_Local


  subroutine SetFluid ( T, F )

    class ( ThermalizationForm ), intent ( in ) :: &
      T
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    real ( KDR ), dimension ( : ), allocatable :: &
      R
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    associate &
      ( R0  =>  ( C % MaxCoordinate  +  C % MinCoordinate ) / 2.0_KDR, &
        X   =>  GV ( :, G % CENTER_U_1 ), & 
        Y   =>  GV ( :, G % CENTER_U_2 ), & 
        Z   =>  GV ( :, G % CENTER_U_3 ) )

    allocate ( R, source = sqrt (      ( X  -  R0 ( 1 ) ) ** 2  &
                                    +  ( Y  -  R0 ( 2 ) ) ** 2  &
                                    +  ( Z  -  R0 ( 3 ) ) ** 2 ) )

    call CO % Initialize ( T % Integrator % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = minval ( R )
    CO % Outgoing % Value ( 2 ) = 1.0_KDR / maxval ( R )
    call CO % Reduce ( REDUCTION % MIN )

    call SetFluidKernel &
           ( T     = FV ( :, F % TEMPERATURE ), &
             N     = FV ( :, F % BARYON_DENSITY_C ), &
             V_1   = FV ( :, F % VELOCITY_U_1 ), &
             V_2   = FV ( :, F % VELOCITY_U_2 ), &
             V_3   = FV ( :, F % VELOCITY_U_3 ), &
             R     = R, &
             T_Min = T % TemperatureMin, &
             T_Max = T % TemperatureMax, &
             R_Min = CO % Incoming % Value ( 1 ), &
             R_Max = 1.0_KDR  /  CO % Incoming % Value ( 2 ) )

    end associate !-- X, etc.
    end associate !-- FV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


  subroutine SetRadiation ( T, R, F )

    class ( ThermalizationForm ), intent ( in ) :: &
      T
    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      R
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      F

    real ( KDR ) :: &
      Amplitude

    select type ( A  =>  R % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        RV  =>  R % Storage_GS % Value )

    !-- Grey, equilibrium

    call R % ComputeEquilibrium ( )
    call Copy ( RV ( :, R % ENERGY_DENSITY_C_EQ ), &
                RV ( :, R % ENERGY_DENSITY_C ) )

    !-- Perturbations on equilibrium

    Amplitude  =  0.5_KDR

    call SetPerturbationKernel &
           ( J   = RV ( :, R % ENERGY_DENSITY_C ), &
             H_1 = RV ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
             H_2 = RV ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
             H_3 = RV ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
             ProperCell = C % ProperCell, &
             A   = Amplitude )

    end associate !-- C, etc.
    end select !-- A

  end subroutine SetRadiation


  subroutine SetFluidKernel &
               ( T, N, V_1, V_2, V_3, R, T_Min, T_Max, R_Min, R_Max )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      N, &
      V_1, V_2, V_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R
    real ( KDR ), intent ( in ) :: &
      T_Min, T_Max, &
      R_Min, R_Max

    T  =  T_Max * ( ( R_Min / R ) ** ( log10 ( T_Max / T_Min ) &
                                       / log10 ( R_Max / R_Min ) ) )

    !-- Only temperature is relevant, but we need a nonzero value for density
    !   to avoid temperature getting blitzed by ComputeFromConserved
    N    =  1.0_KDR
    V_1  =  0.0_KDR
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR

  end subroutine SetFluidKernel


  subroutine SetPerturbationKernel ( J, H_1, H_2, H_3, ProperCell, A )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), intent ( in ) :: &
      A  !-- Amplitude

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
    real ( KDR ) :: &      
      P  !-- Perturbation

    nV  =  size ( J )

    do iV  =  1, nV

      if ( .not. ProperCell ( iV ) ) &
        cycle

      call random_number ( P )
      P         =  A * 2.0_KDR * ( P - 0.5_KDR ) 
      J ( iV )  =  ( 1.0_KDR + P )  *  J ( iV )

      call random_number ( P )
      P           =  0.01_KDR  *  J ( iV )  *  2.0_KDR * ( P - 0.5_KDR ) 
      H_1 ( iV )  =  P

      call random_number ( P )
      P           =  0.01_KDR  *  J ( iV )  *  2.0_KDR * ( P - 0.5_KDR ) 
      H_2 ( iV )  =  P

      call random_number ( P )
      P           =  0.01_KDR  *  J ( iV )  *  2.0_KDR * ( P - 0.5_KDR ) 
      H_3 ( iV )  =  P

    end do

  end subroutine SetPerturbationKernel


  subroutine ComputeFractionalDifferenceKernel ( J_FD, J, J_R, ProperCell )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J_FD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, J_R
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
 
    nV  =  size ( J )

    !$OMP parallel do private ( iV )
    do iV  =  1,  nV

      if ( .not. ProperCell ( iV ) ) &
        cycle

      J_FD ( iV )  =  ( J ( iV )  -  J_R ( iV ) )  /  J_R ( iV )

    end do !-- iV
    !$OMP end parallel do
    
  end subroutine ComputeFractionalDifferenceKernel


end module Thermalization_Form
