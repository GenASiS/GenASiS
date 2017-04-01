module Thermalization_Form

  use Basics
  use Mathematics
  use Fluid_P_NR__Form
  use Fluid_ASC__Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  ! use Interactions_F__Form
  ! use SetFermiDiracSpectrum_Command
  use RadiationMoments_BSLL_ASC_CSLD__Form
  use Interactions_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( Integrator_C_MS_C_PS_Template ) :: &
    ThermalizationForm
      type ( RadiationMoments_ASC_Form ), allocatable :: &
        Reference_ASC, &
        FractionalDifference_ASC
      type ( Interactions_BSLL_ASC_CSLD_Form ), allocatable :: &
        Interactions_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
  !   procedure, private, pass :: &
  !     ComputeCycle
  !   procedure, private, pass :: &
  !     ComputeTimeStepLocal
  end type ThermalizationForm

    private :: &
!       SetMatter, &
!       SetRadiation, &
!       SetInteractions, &
      SetReference!, &
!       ComputeCycle_BSLL_ASC_CSLD

    real ( KDR ), private :: &
      TemperatureMin, &
      TemperatureMax, &
      EnergyScale, &
      EffectiveOpacity, &
      TransportOpacity, &
      TimeScale

contains


  subroutine Initialize ( T, Name )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nEnergyCells
    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ) :: &
      EnergyDensityUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    TemperatureMin  =   0.1_KDR  *  UNIT % MEV
    TemperatureMax  =  10.0_KDR  *  UNIT % MEV
    EnergyScale     =   3.0_KDR  *  UNIT % MEV
    call PROGRAM_HEADER % GetParameter ( TemperatureMin, 'TemperatureMin' )
    call PROGRAM_HEADER % GetParameter ( TemperatureMax, 'TemperatureMax' )
    call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

    EffectiveOpacity = 1.0_KDR
    TransportOpacity = 1.1_KDR
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )
    call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    associate &
      ( c     => CONSTANT % SPEED_OF_LIGHT, &
        Kappa => TransportOpacity )
    TimeScale  =  1.0 / ( c * Kappa )
    end associate !-- c, etc.
    
    !-- PositionSpace

    allocate ( Atlas_SC_Form :: T % PositionSpace )
    select type ( PS => T % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( T % Geometry_ASC )
    associate ( GA => T % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: T % MomentumSpace )
    select type ( MS => T % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, Name )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = EnergyScale

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEV

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    nEnergyCells = 20
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Radiation

    EnergyDensityUnit  =  UNIT % MEV / UNIT % HBAR_C ** 3

    allocate ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
               T % Current_BSLL_ASC_CSLD )
    select type ( RMB => T % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )
    call RMB % Initialize &
           ( MS, 'GENERIC', EnergyDensityUnitOption = EnergyDensityUnit )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: T % Current_ASC )
    select type ( FA => T % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize &
           ( PS, 'NON_RELATIVISTIC', TemperatureUnitOption = UNIT % MEV )

    !-- Interactions

    allocate ( T % Interactions_BSLL_ASC_CSLD )
    associate ( IB => T % Interactions_BSLL_ASC_CSLD )
    call IB % Initialize &
           ( MS, 'FIXED', EnergyDensityUnitOption = EnergyDensityUnit )
    call RMB % SetInteractions ( IB )
    end associate !-- IB

    !-- Step

    allocate ( Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form :: T % Step )
    select type ( S => T % Step )
    class is ( Step_RK2_C_BSLL_ASC_CSLD_C_ASC_Form )
    call S % Initialize ( RMB, FA, Name )
    S % ApplyDivergence   % Pointer => null ( )  !-- Disable fluid evolution
    S % ApplyDivergence_S % Pointer => null ( )  !-- Disable spatial 
                                                 !   section evolution
    S % ApplyRelaxation_F % Pointer => ApplyRelaxation_Interactions
    end select !-- S

    !-- Diagnostics

    EnergyDensityUnit  =  UNIT % MEV ** 4 / UNIT % HBAR_C ** 3

    allocate ( T % Reference_ASC )
    allocate ( T % FractionalDifference_ASC )
    call T % Reference_ASC % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Reference', &
             EnergyDensityUnitOption = EnergyDensityUnit )
    call T % FractionalDifference_ASC % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'FractionalDifference' )
    T % SetReference => SetReference

!     !-- Initial conditions

!     call SetMatter ( T )
!     call SetRadiation ( T )
!     call SetInteractions ( T )

!     !-- Initialize template

!     call T % InitializeTemplate &
!            ( Name, FinishTimeOption = 10.0_KDR * TimeScale )

    !-- Cleanup

    end select !-- FA
    end select !-- RMB
    end select !-- MS
    end select !-- PS

  end subroutine Initialize


  subroutine Finalize ( T )

    type ( ThermalizationForm ), intent ( inout ) :: &
      T

    if ( allocated ( T % Interactions_BSLL_ASC_CSLD ) ) &
      deallocate ( T % Interactions_BSLL_ASC_CSLD )
    if ( allocated ( T % FractionalDifference_ASC ) ) &
      deallocate ( T % FractionalDifference_ASC )
    if ( allocated ( T % Reference_ASC ) ) &
      deallocate ( T % Reference_ASC )

    call T % FinalizeTemplate_C_MS_C_PS ( )

  end subroutine Finalize


!   subroutine ComputeCycle ( I )

!     class ( ThermalizationForm ), intent ( inout ) :: &
!       I

!     associate ( Timer => PROGRAM_HEADER % Timer ( I % iTimerComputeCycle ) )
!     call Timer % Start ( )

!     select type ( MS => I % MomentumSpace )
!     class is ( Bundle_SLL_ASC_CSLD_Form )
!       call ComputeCycle_BSLL_ASC_CSLD ( I, MS )
!     end select !-- MS

!     call Timer % Stop ( )
!     end associate !-- Timer

!   end subroutine ComputeCycle


!   subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

!     class ( ThermalizationForm ), intent ( in ) :: &
!       I
!     real ( KDR ), dimension ( : ), intent ( inout ) :: &
!       TimeStepCandidate

!     TimeStepCandidate ( 1 ) = 0.01_KDR * TimeScale

!   end subroutine ComputeTimeStepLocal


!   subroutine SetMatter ( T )

!     class ( ThermalizationForm ), intent ( inout ) :: &
!       T

!     real ( KDR ) :: &
!       R_Min, R_Max, &
!       T_Min, T_Max
!     type ( CollectiveOperation_R_Form ) :: &
!       CO
!     class ( GeometryFlatForm ), pointer :: &
!       G
!     class ( MatterForm ), pointer :: &
!       M

!     M => T % Matter_ASC % Matter ( )

!     select type ( PS => T % PositionSpace )
!     class is ( Atlas_SC_Form )
!     G => PS % Geometry ( )

!     select type ( C => PS % Chart )
!     class is ( Chart_SL_Template )

!     associate &
!       ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
!         L  => ( C % MaxCoordinate - C % MinCoordinate ), &
!         X  => G % Value ( :, G % CENTER ( 1 ) ), &
!         Y  => G % Value ( :, G % CENTER ( 2 ) ), &
!         Z  => G % Value ( :, G % CENTER ( 3 ) ), &
!         Temperature => M % Value ( :, M % TEMPERATURE ), &
!         ChemicalPotential => M % Value ( :, M % CHEMICAL_POTENTIAL ) )
!     associate &
!       ( R => sqrt ( ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
!                     +  ( Z - R0 ( 3 ) ) ** 2 ) )

!     call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
!     CO % Outgoing % Value ( 1 ) = minval ( R )
!     CO % Outgoing % Value ( 2 ) = 1.0_KDR / maxval ( R )
!     call CO % Reduce ( REDUCTION % MIN )

!     R_Min = CO % Incoming % Value ( 1 )
!     R_Max = 1.0_KDR / CO % Incoming % Value ( 2 )

!     T_Min = TemperatureMin
!     T_Max = TemperatureMax

! !    Temperature &
! !      = T_Max  -  ( T_Max - T_Min ) / ( R_Max - R_Min ) * ( R - R_Min )
!     Temperature = T_Max * ( ( R_Min / R ) ** ( log10 ( T_Max / T_Min ) &
!                                                / log10 ( R_Max / R_Min ) ) )
!     ChemicalPotential = 0.0_KDR

!     end associate !-- R
!     end associate !-- R0, etc.
!     end select !-- C
!     end select !-- PS

!     nullify ( G, M )

!   end subroutine SetMatter


!   subroutine SetRadiation ( T )

!     class ( ThermalizationForm ), intent ( inout ) :: &
!       T

!     integer ( KDI ) :: &
!       iF, &  !-- iFiber
!       iE     !-- iEnergy  
!     real ( KDR ) :: &
!       Amplitude, &
!       Perturbation
!     class ( RadiationMomentsForm ), pointer :: &
!       RM
!     class ( GeometryFlatForm ), pointer :: &
!       G
!     class ( MatterForm ), pointer :: &
!       M

!     associate ( RMB => T % RadiationMoments_BSLL_ASC_CSLD )

!     select type ( MS => T % MomentumSpace )
!     class is ( Bundle_SLL_ASC_CSLD_Form )

!     G => MS % Base_CSLD % Geometry ( )
!     M => T % Matter_ASC % Matter ( )

!     call InitializeRandomSeed ( PROGRAM_HEADER % Communicator )

!     do iF = 1, MS % nFibers
!       associate ( iBC => MS % iaBaseCell ( iF ) )
!       RM => RMB % RadiationMomentsFiber ( iF )
!       associate &
!         ( J   => RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
!           H_1 => RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
!           H_2 => RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
!           H_3 => RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
!           T   => M % Value ( iBC, M % TEMPERATURE ), &
!           Mu  => M % Value ( iBC, M % CHEMICAL_POTENTIAL ), &
!           E   => RMB % Energy )

!       call SetFermiDiracSpectrum ( E, T, Mu, J )

!       Amplitude = 0.9_KDR
!       do iE = 1, RMB % nEnergyValues

!         call random_number ( Perturbation )
!         Perturbation = Amplitude * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
!         J ( iE ) = ( 1.0_KDR + Perturbation )  *  J ( iE )

!         call random_number ( Perturbation )
!         Perturbation &
!           = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
!         H_1 ( iE ) = Perturbation

!         call random_number ( Perturbation )
!         Perturbation &
!           = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
!         H_2 ( iE ) = Perturbation

!         call random_number ( Perturbation )
!         Perturbation &
!           = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
!         H_3 ( iE ) = Perturbation
        
!       end do !-- iE

!       call RM % ComputeFromPrimitive ( iBC, G )

!       end associate !-- J, etc.
!       end associate !-- iBC
!     end do !-- iF

!     end select !-- MS
!     end associate !-- RMB
!     nullify ( G, M, RM )

!   end subroutine SetRadiation


!   subroutine SetInteractions ( T )

!     class ( ThermalizationForm ), intent ( inout ) :: &
!       T

!     class ( Interactions_F_Form ), pointer :: &
!       I
!     class ( MatterForm ), pointer :: &
!       M

!     integer ( KDI ) :: &
!       iF  !-- iFiber

!     associate &
!       ( IB => T % Interactions_BSLL_ASC_CSLD, &
!         RMB => T % RadiationMoments_BSLL_ASC_CSLD )

!     select type ( MS => T % MomentumSpace )
!     class is ( Bundle_SLL_ASC_CSLD_Form )

!     M => T % Matter_ASC % Matter ( )

!     do iF = 1, MS % nFibers
!       associate ( iBC => MS % iaBaseCell ( iF ) )
!       I => IB % InteractionsFiber_F ( iF )
!       associate &
!         ( J_Eq  => I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
!           Chi   => I % Value ( :, I % EFFECTIVE_OPACITY ), &
!           Kappa => I % Value ( :, I % TRANSPORT_OPACITY ), &
!           T     => M % Value ( iBC, M % TEMPERATURE ), &
!           Mu    => M % Value ( iBC, M % CHEMICAL_POTENTIAL ), &
!           E     => RMB % Energy )

!       call SetFermiDiracSpectrum ( E, T, Mu, J_Eq )

!       Chi   = EffectiveOpacity
!       Kappa = TransportOpacity

!       end associate !-- J_Eq, etc.
!       end associate !-- iBC
!     end do !-- iF

!     end select !-- MS
!     end associate !-- IB, etc.
!     nullify ( M, I )

!   end subroutine SetInteractions


  subroutine SetReference ( T )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      T

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iF     !-- iFiber
    class ( RadiationMomentsForm ), pointer :: &
      RM_R, &  !-- RM_Reference
      RM_C, &  !-- RM_Computed
      RM_FD    !-- RM_FractionalDifference
    class ( Fluid_P_NR_Form ), pointer :: &
      F

    select type ( T )
    class is ( ThermalizationForm )

    select type ( MS => T % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    select type ( RMB => T % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    RM_R  => T % Reference_ASC % RadiationMoments ( )
    RM_C  => RMB % EnergyIntegral % RadiationMoments ( )
    RM_FD => T % FractionalDifference_ASC % RadiationMoments ( )

    select type ( FA => T % Current_ASC )
    class is ( Fluid_ASC_Form )

    F => FA % Fluid_P_NR ( )

    !-- Reference blackbody

    associate &
      ( J_R  =>  RM_R % Value ( :, RM_R % COMOVING_ENERGY_DENSITY ), &
        T    =>  F % Value ( :, F % TEMPERATURE ), &
        a    =>  CONSTANT % RADIATION )

    do iV = 1, RM_R % nValues
      if ( T ( iV ) > 0.0_KDR ) then
        J_R ( iV )  =  a  *  T ( iV ) ** 4
      else
        J_R ( iV ) = huge ( 0.0_KDR )
      end if
    end do !-- iV

    !-- Integrated spectrum and difference

    call RMB % ComputeEnergyIntegral ( )

    associate &
      ( J_C  => RM_C  % Value ( :, RM_C  % COMOVING_ENERGY_DENSITY ), &
        J_FD => RM_FD % Value ( :, RM_FD % COMOVING_ENERGY_DENSITY ) )

    do iF = 1, MS % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )
      J_FD ( iBC )  =  ( J_C ( iBC )  -  J_R ( iBC ) )  /  J_R ( iBC )
      end associate !-- iBC
    end do !-- iF

    end associate !-- J_C, etc.
    end associate !-- J_R, etc.
    end select !-- FA
    end select !-- RMB
    end select !-- MS
    end select !-- T

    nullify ( F, RM_R, RM_C, RM_FD )

  end subroutine SetReference


!   subroutine ComputeCycle_BSLL_ASC_CSLD ( T, MS )

!     class ( ThermalizationForm ), intent ( inout ) :: &
!       T
!     class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
!       MS

!     integer ( KDI ) :: &
!       iF     !-- iFiber
!     real ( KDR ) :: &
!       TimeNew
!     class ( GeometryFlatForm ), pointer :: &
!       G
!     class ( RadiationMomentsForm ), pointer :: &
!       RM

!     call T % ComputeNewTime ( TimeNew )

!     associate &
!       ( S   => T % Step, &
!         RMB => T % RadiationMoments_BSLL_ASC_CSLD, &
!         CF  => MS % Fiber_CSLL, &
!         TimeStep => TimeNew - T % Time )    

!     G => MS % Base_CSLD % Geometry ( )

!     do iF = 1, MS % nFibers
!       associate ( iBC => MS % iaBaseCell ( iF ) )
!       RM => RMB % RadiationMomentsFiber ( iF )
!       call S % Compute &
!              ( RM, CF, T % Time, TimeStep, GeometryOption = G, &
!                iGeometryValueOption = iBC )
!       end associate !-- iBC
!     end do !-- iF

!     T % iCycle = T % iCycle + 1
!     T % Time = T % Time + TimeStep
!     if ( T % Time == T % WriteTime ) &
!       T % IsCheckpointTime = .true.

!     end associate !-- S, etc.
!     nullify ( G, RM )

!   end subroutine ComputeCycle_BSLL_ASC_CSLD


end module Thermalization_Form
