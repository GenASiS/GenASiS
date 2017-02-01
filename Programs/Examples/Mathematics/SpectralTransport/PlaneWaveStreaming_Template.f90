module PlaneWaveStreaming_Template

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use RadiationMoments_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( Integrator_C_Template ), abstract :: &
    PlaneWaveStreamingTemplate
      type ( RadiationMoments_ASC_Form ), dimension ( : ), allocatable :: &
        Section_ASC, &
        Reference_ASC, &
        Difference_ASC        
      type ( RadiationMoments_BSLL_ASC_CSLD_Form ), allocatable :: &
        RadiationMoments_BSLL_ASC_CSLD
  contains 
    procedure, public, pass :: &
      InitializeTemplate_PWS
    procedure, public, pass :: &
      ComputeError
    procedure, public, pass :: &
      FinalizeTemplate_PWS
    procedure ( Waveform ), private, pass, deferred :: &
      Waveform
    procedure, private, pass :: &
      ComputeCycle
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type PlaneWaveStreamingTemplate

  abstract interface
    elemental function Waveform ( PWS, X ) result ( W )
      !-- Waveform with a full period in the range 0 < X < 1
      use Basics
      import PlaneWaveStreamingTemplate
      class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
        PWS
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ) :: &
        W
    end function Waveform
  end interface

    private :: &
      SetReference, &
      SetRadiation, &
      ComputeCycle_BSLL_ASC_CSLD

      private :: &
        SetWave

        private :: &
          SetWaveKernel

contains


  subroutine InitializeTemplate_PWS ( PWS, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iE, &  !-- iEnergy
      nEnergyCells
    integer ( KDI ), dimension ( 3 ) :: &
      nCellsPosition
    real ( KDR ) :: &
      Period
    character ( 1 + 2 ) :: &
      EnergyNumber

    if ( PWS % Type == '' ) &
      PWS % Type = 'a PlaneWaveStreaming'

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: PWS % PositionSpace )
    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    nCellsPosition = [ 128, 128, 128 ]
    call PROGRAM_HEADER % GetParameter ( nCellsPosition, 'nCellsPosition' )

    call PS % CreateChart ( nCellsOption = nCellsPosition )

    !-- Geometry of PositionSpace

    allocate ( PWS % Geometry_ASC )
    associate ( GA => PWS % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: PWS % MomentumSpace )
    select type ( MS => PWS % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, Name )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    nEnergyCells = 4
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( CoordinateSystemOption = 'SPHERICAL', &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Radiation

    allocate ( PWS % RadiationMoments_BSLL_ASC_CSLD )
    associate ( RMB => PWS % RadiationMoments_BSLL_ASC_CSLD )
    call RMB % Initialize ( MS, 'GENERIC' )

    !-- Step

    allocate ( Step_RK2_C_Form :: PWS % Step )
    select type ( S => PWS % Step )
    class is ( Step_RK2_C_Form )
    call S % Initialize ( Name )
    end select !-- S

    !-- Diagnostics

    allocate ( PWS % Section_ASC ( RMB % nEnergyValues ) )
    allocate ( PWS % Reference_ASC ( RMB % nEnergyValues ) )
    allocate ( PWS % Difference_ASC ( RMB % nEnergyValues ) )
    do iE = 1, RMB % nEnergyValues
      write ( EnergyNumber, fmt = '(a1,i2.2)' ) '_', iE
      call PWS % Section_ASC ( iE ) % Initialize &
             ( PS, 'GENERIC', NameOutputOption = 'Section' // EnergyNumber )
      call PWS % Reference_ASC ( iE ) % Initialize &
             ( PS, 'GENERIC', NameOutputOption = 'Reference' // EnergyNumber )
      call PWS % Difference_ASC ( iE ) % Initialize &
             ( PS, 'GENERIC', NameOutputOption = 'Difference' // EnergyNumber )
    end do !-- iE
    PWS % SetReference => SetReference

    !-- Initial conditions

    call SetRadiation ( PWS, Time = 0.0_KDR, Period = Period )

    !-- Initialize template

    call PWS % InitializeTemplate_C ( Name, FinishTimeOption = Period )

    !-- Cleanup

    end associate !-- RMB
    end select !-- MS
    end select !-- PS

  end subroutine InitializeTemplate_PWS


  subroutine ComputeError ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    
    integer ( KDI ) :: &
      iE  !-- iEnergy
    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO
    class ( RadiationMomentsForm ), pointer :: &
      RS_D     !-- RM_Difference

    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( CB => PS % Chart )
    class is ( Chart_SL_Template )

    associate ( RMB => PWS % RadiationMoments_BSLL_ASC_CSLD )

    do iE = 1, RMB % nEnergyValues
      associate ( RSA_D => PWS % Difference_ASC ( iE ) )
      RS_D => RSA_D % RadiationMoments ( )
      allocate ( CO )

      associate &
        ( Difference => RS_D % Value ( :, RS_D % COMOVING_ENERGY_DENSITY ) )
      call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
      CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference ), &
                                          mask = CB % IsProperCell )
      CO % Outgoing % Value ( 2 ) = CB % nProperCells
      call CO % Reduce ( REDUCTION % SUM )
      end associate !-- Difference

      associate &
        ( DifferenceSum => CO % Incoming % Value ( 1 ), &
          nValues => CO % Incoming % Value ( 2 ) )
      L1 = DifferenceSum / nValues
      end associate 

      call Show ( iE, '*** iEnergyBin', nLeadingLinesOption  = 2 ) 
      call Show ( L1, '*** L1 error',   nTrailingLinesOption = 2 )

      deallocate ( CO )
      end associate !-- RSA_D
    end do !-- iE

    end associate !-- RMB
    end select !-- CB
    end select !-- PS
    nullify ( RS_D )

  end subroutine ComputeError


  subroutine FinalizeTemplate_PWS ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS

    if ( allocated ( PWS % Difference_ASC ) ) &
      deallocate ( PWS % Difference_ASC )
    if ( allocated ( PWS % Reference_ASC ) ) &
      deallocate ( PWS % Reference_ASC )
    if ( allocated ( PWS % Section_ASC ) ) &
      deallocate ( PWS % Section_ASC )
    if ( allocated ( PWS % RadiationMoments_BSLL_ASC_CSLD ) ) &
      deallocate ( PWS % RadiationMoments_BSLL_ASC_CSLD )

    call PWS % FinalizeTemplate_C ( )

  end subroutine FinalizeTemplate_PWS


  subroutine ComputeCycle ( I )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      I

    associate ( Timer => PROGRAM_HEADER % Timer ( I % iTimerComputeCycle ) )
    call Timer % Start ( )

    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
      call ComputeCycle_BSLL_ASC_CSLD ( I, MS )
    end select !-- MS

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeCycle


  subroutine ComputeTimeStepLocal ( I, TimeStep )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      I
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    integer ( KDI ) :: &
      iE     !-- iEnergy
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( RadiationMomentsForm ), allocatable :: &
      RadiationSection

    associate ( RMB => I % RadiationMoments_BSLL_ASC_CSLD )

    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    G => MS % Base_CSLD % Geometry ( )

    allocate ( RadiationSection )
    associate ( RS => RadiationSection )
    call RS % Initialize &
           ( RMB % Velocity_U_Unit, RMB % MomentumDensity_U_Unit, &
             RMB % MomentumDensity_D_Unit, RMB % EnergyDensityUnit, &
             G % nValues, ClearOption = .true. )
    
    do iE = 1, RMB % nEnergyValues
      call MS % LoadSection ( RS, RMB, iE )
      call I % ComputeTimeStepKernel_CSL &
             ( MS % Base_CSLD % IsProperCell, &
               RS % Value ( :, RS % FAST_EIGENSPEED_PLUS ( 1 ) ), &
               RS % Value ( :, RS % FAST_EIGENSPEED_PLUS ( 2 ) ), &
               RS % Value ( :, RS % FAST_EIGENSPEED_PLUS ( 3 ) ), &
               RS % Value ( :, RS % FAST_EIGENSPEED_MINUS ( 1 ) ), &
               RS % Value ( :, RS % FAST_EIGENSPEED_MINUS ( 2 ) ), &
               RS % Value ( :, RS % FAST_EIGENSPEED_MINUS ( 3 ) ), &
               G % Value ( :, G % WIDTH ( 1 ) ), &
               G % Value ( :, G % WIDTH ( 2 ) ), & 
               G % Value ( :, G % WIDTH ( 3 ) ), &
               MS % Base_CSLD % nDimensions, TimeStep )
    end do !-- iE

    TimeStep = I % CourantFactor * TimeStep

    end associate !-- RS
    end select !-- MS
    end associate !-- RMB
    nullify ( G )

  end subroutine ComputeTimeStepLocal


  subroutine SetReference ( PWS )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      PWS

    integer ( KDI ) :: &
      iE  !-- iEnergy
    class ( RadiationMomentsForm ), pointer :: &
      RS, &
      RS_R, &  !-- RM_Reference
      RS_D     !-- RM_Difference

    select type ( PWS )
    class is ( PlaneWaveStreamingTemplate )

    select type ( MS => PWS % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    associate ( RMB => PWS % RadiationMoments_BSLL_ASC_CSLD )

    do iE = 1, RMB % nEnergyValues
      associate &
        ( RSA   => PWS % Section_ASC ( iE ), &
          RSA_R => PWS % Reference_ASC ( iE ), &
          RSA_D => PWS % Difference_ASC ( iE ) )
      RS   => RSA   % RadiationMoments ( )
      RS_R => RSA_R % RadiationMoments ( )
      RS_D => RSA_D % RadiationMoments ( )

      !-- Computed
      call MS % LoadSection ( RS, RMB, iE )

      !-- Reference
      call SetWave ( PWS, RS_R, MS % Base_CSLD, PWS % Time, nWavelengths = iE )

      !-- Difference
! !    RS_D % Value  =  RS % Value  -  RS_R % Value
      call MultiplyAdd ( RS % Value, RS_R % Value, -1.0_KDR, RS_D % Value )
      end associate !-- RSA, etc.
    end do !-- iE

    end associate !-- RMB
    end select !-- MS
    end select !-- PWS
    nullify ( RS, RS_R, RS_D )

  end subroutine SetReference


  subroutine SetRadiation ( PWS, Time, Period )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS
    real ( KDR ), intent ( in ) :: &
      Time
    real ( KDR ), intent ( out ) :: &
      Period
    
    integer ( KDI ) :: &
      iE  !-- iEnergy
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( RadiationMomentsForm ), allocatable :: &
      RadiationSection

    associate ( RMB => PWS % RadiationMoments_BSLL_ASC_CSLD )

    select type ( MS => PWS % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    G => MS % Base_CSLD % Geometry ( )

    allocate ( RadiationSection )
    associate ( RS => RadiationSection )
    call RS % Initialize &
           ( RMB % Velocity_U_Unit, RMB % MomentumDensity_U_Unit, &
             RMB % MomentumDensity_D_Unit, RMB % EnergyDensityUnit, &
             G % nValues, ClearOption = .true. )

    iE = 1
    call SetWave &
           ( PWS, RS, MS % Base_CSLD, Time, nWavelengths = iE, &
             PeriodOption = Period )
    call MS % StoreSection ( RMB, RS, iE )
    do iE = 2, RMB % nEnergyValues
      call SetWave ( PWS, RS, MS % Base_CSLD, Time, nWavelengths = iE )
      call MS % StoreSection ( RMB, RS, iE )
    end do !-- iE

    end associate !-- RS
    end select !-- MS
    end associate !-- RMB
    nullify ( G )

  end subroutine SetRadiation


  subroutine ComputeCycle_BSLL_ASC_CSLD ( PWS, MS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      MS

    integer ( KDI ) :: &
      iE     !-- iEnergy
    real ( KDR ) :: &
      TimeNew
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( RadiationMomentsForm ), allocatable :: &
      RadiationSection

    select type ( S => PWS % Step )
    class is ( Step_RK_C_Template )

    call PWS % ComputeNewTime ( TimeNew )

    associate &
      ( RMB => PWS % RadiationMoments_BSLL_ASC_CSLD, &
        CB  => MS % Base_CSLD, &
        TimeStep => TimeNew - PWS % Time )    

    G => CB % Geometry ( )

    allocate ( RadiationSection )
    associate ( RS => RadiationSection )
    call RS % Initialize &
           ( RMB % Velocity_U_Unit, RMB % MomentumDensity_U_Unit, &
             RMB % MomentumDensity_D_Unit, RMB % EnergyDensityUnit, &
             G % nValues, ClearOption = .true. )

    do iE = 1, RMB % nEnergyValues
      call MS % LoadSection ( RS, RMB, iE, GhostExchangeOption = .true. )
      call S % Compute ( RS, CB, PWS % Time, TimeStep )
      call MS % StoreSection ( RMB, RS, iE )
    end do !-- iF

    PWS % iCycle = PWS % iCycle + 1
    PWS % Time = PWS % Time + TimeStep
    if ( PWS % Time == PWS % WriteTime ) &
      PWS % IsCheckpointTime = .true.

    end associate !-- RS
    end associate !-- RMB, etc.
    end select !-- S
    nullify ( G )

  end subroutine ComputeCycle_BSLL_ASC_CSLD


  subroutine SetWave ( PWS, RS, PSC, Time, nWavelengths, PeriodOption )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    type ( RadiationMomentsForm ), intent ( inout ) :: &
      RS  !-- RadiationSection
    class ( Chart_SL_Template ), intent ( in ) :: &
      PSC  !-- PositionSpaceChart
    real ( KDR ), intent ( in ) :: &
      Time
    integer ( KDI ), intent ( in ) :: &
      nWavelengths
    real ( KDR ), intent ( out ), optional :: &
      PeriodOption

    real ( KDR ), dimension ( 3 ) :: &
      BoxSize, &
      Wavenumber
    class ( GeometryFlatForm ), pointer :: &
      G

    BoxSize = PSC % MaxCoordinate - PSC % MinCoordinate
    where ( BoxSize > 0.0_KDR )
      Wavenumber = nWavelengths / BoxSize
    elsewhere
      Wavenumber = 0.0_KDR
    end where

    associate &
      ( K => Wavenumber, &
        c => CONSTANT % SPEED_OF_LIGHT )
    associate ( Abs_K => sqrt ( dot_product ( K, K ) ) )

    if ( present ( PeriodOption ) ) then
      PeriodOption = 1.0_KDR / ( c * Abs_K )
      call Show ( PeriodOption, 'Period' )
    end if

    G => PSC % Geometry ( )

    call SetWaveKernel &
           ( PWS, &
             X  = G % Value ( :, G % CENTER ( 1 ) ), &
             Y  = G % Value ( :, G % CENTER ( 2 ) ), &
             Z  = G % Value ( :, G % CENTER ( 3 ) ), &
             K  = K, &
             V  = c, &
             T  = Time, &
             J  = RS % Value ( :, RS % COMOVING_ENERGY_DENSITY ), &
             HX = RS % Value ( :, RS % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
             HY = RS % Value ( :, RS % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
             HZ = RS % Value ( :, RS % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ) )

    call RS % ComputeFromPrimitive ( G )

    end associate !-- Abs_K
    end associate !-- K, etc.

    nullify ( G )

  end subroutine SetWave


  subroutine SetWaveKernel ( PWS, X, Y, Z, K, V, T, J, HX, HY, HZ )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      K
    real ( KDR ), intent ( in ) :: &
      V, &
      T
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J, &
      HX, HY, HZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      Abs_K, &
      VX, VY, VZ

    nValues = size ( X )
    
    Abs_K  =  sqrt ( dot_product ( K, K ) )
       VX  =  V * K ( 1 ) / Abs_K
       VY  =  V * K ( 2 ) / Abs_K
       VZ  =  V * K ( 3 ) / Abs_K

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      J ( iV )  =  PWS % Waveform &
                     (    K ( 1 ) * ( X ( iV )  -  VX * T ) &
                       +  K ( 2 ) * ( Y ( iV )  -  VY * T ) &
                       +  K ( 3 ) * ( Z ( iV )  -  VZ * T ) )
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      HX ( iV )  =  VX  *  J ( iV )
      HY ( iV )  =  VY  *  J ( iV )
      HZ ( iV )  =  VZ  *  J ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetWaveKernel


end module PlaneWaveStreaming_Template
