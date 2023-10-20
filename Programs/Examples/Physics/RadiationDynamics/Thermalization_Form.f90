module Thermalization_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: ThermalizationForm
    real ( KDR ) :: &
      TemperatureMin, &
      TemperatureMax, &
      OpacityAbsorption, &
      TimeScale
    type ( RadiationMoments_BM_Form ), allocatable :: &
      Reference, &
      FractionalDifference
  contains
    procedure, private, pass :: &
      Initialize_T
    generic, public :: &
      Initialize => Initialize_T
    final :: &
      Finalize
  end type ThermalizationForm

    private :: &
      InitializeUniverse!, &
      ! InitializeDiagnostics, &
      ! SetInitial, &
      ! SetReference

      private :: &
        SetFluid
      !   SetRadiation

        private :: &
          SetFluidKernel!, &
!          SetRadiationKernel

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
!    call InitializeDiagnostics ( T )
 
  end subroutine Initialize_T


  impure elemental subroutine Finalize ( T )

    type ( ThermalizationForm ), intent ( inout ), target :: &
      T

    if ( allocated ( T % FractionalDifference ) ) &
      deallocate ( T % FractionalDifference )
    if ( allocated ( T % Reference ) ) &
      deallocate ( T % Reference )

  end subroutine Finalize


  subroutine InitializeUniverse ( T, FormalismType, Name )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    ! integer ( KDI ) :: &
    !   iD

    call T % Initialize &
           ( RadiationName = [ 'Radiation_1', 'Radiation_2' ], &
             RadiationType = [ 'GENERIC', 'GENERIC' ], &
             FormalismType = FormalismType, &
             Name = Name, &
             ApplyStreamingOption = .false., &
             EvolveFluidOption = .false., &
             nCellsPositionOption = [ 128, 128, 128 ] )

    ! select type ( I  =>  T % Integrator )
    !   class is ( Integrator_CS_Form )
    ! associate &
    !   ( F  =>  I % CurrentSet_X )
    ! do iD  =  1, 3
    !   call F % SetBoundaryConditionsFace &
    !          ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    ! end do !-- iD
    ! end associate !-- F
    ! end select !-- I
             
    ! select type ( I  =>  T % Integrator )
    !   class is ( Integrator_CS_1D_BM_CS_Form )
    ! associate &
    !   ( R  =>  I % CurrentSet_X_1D )
    ! do iD  =  1, 3
    !   call R % SetBoundaryConditionsFace &
    !          ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    ! end do !-- iD
    ! end associate !-- R
    ! end select !-- I
             
    ! T % Integrator % SetInitial    =>  SetInitial
    ! T % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


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
    end associate !-- RV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


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


end module Thermalization_Form
