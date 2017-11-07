module PlaneWaveStreaming_Template

  use Basics
  use Mathematics
  use Fluid_ASC__Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use RadiationMoments_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( Integrator_C_MS_C_PS_Template ), abstract :: &
    PlaneWaveStreamingTemplate
      type ( RadiationMoments_ASC_Form ), dimension ( : ), allocatable :: &
        Reference_ASC, &
        Difference_ASC        
  contains 
    procedure, public, pass :: &
      InitializeTemplate_PWS
    procedure, public, pass :: &
      ComputeError
    procedure, public, pass :: &
      FinalizeTemplate_PWS
    procedure ( Waveform ), private, pass, deferred :: &
      Waveform
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
      SetRadiation

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

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: PWS % PositionSpace )
    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    nCellsPosition = [ 128, 128, 128 ]
    call PROGRAM_HEADER % GetParameter ( nCellsPosition, 'nCellsPosition' )

    call PS % CreateChart ( nCellsOption = nCellsPosition )
    call PS % SetGeometry ( )

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: PWS % MomentumSpace )
    select type ( MS => PWS % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, 'MomentumSpace' )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    nEnergyCells = 4
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( CoordinateSystemOption = 'SPHERICAL', &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Prepare for Currents

    allocate ( PWS % TimeStepLabel ( 2 ) )
    PWS % TimeStepLabel ( 1 ) = 'Radiation'
    PWS % TimeStepLabel ( 2 ) = 'Fluid'

    !-- Radiation

    allocate ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
               PWS % Current_BSLL_ASC_CSLD )
    select type ( RMB => PWS % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )
    call RMB % Initialize ( MS, 'GENERIC' )

    !-- Fluid ( This doesn't do anything, just a Step interface placeholder )

    allocate ( Fluid_ASC_Form :: PWS % Current_ASC )
    select type ( FA => PWS % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'NON_RELATIVISTIC' )

    !-- Step

    allocate ( Step_RK2_C_BSLL_ASC_CSLD_Form :: PWS % Step_MS )
    select type ( S_MS => PWS % Step_MS )
    class is ( Step_RK2_C_BSLL_ASC_CSLD_Form )
    call S_MS % Initialize ( RMB, Name )
    end select !-- S_MS

    allocate ( Step_RK2_C_ASC_Form :: PWS % Step_PS )
    select type ( S_PS => PWS % Step_PS )
    class is ( Step_RK2_C_ASC_Form )
    call S_PS % Initialize ( FA, Name )
    S_PS % ApplyDivergence % Pointer => null ( )  !-- Disable fluid evolution
    end select !-- S

    !-- Diagnostics

    allocate ( PWS % Reference_ASC ( RMB % nEnergyValues ) )
    allocate ( PWS % Difference_ASC ( RMB % nEnergyValues ) )
    do iE = 1, RMB % nEnergyValues
      write ( EnergyNumber, fmt = '(a1,i2.2)' ) '_', iE
      call PWS % Reference_ASC ( iE ) % Initialize &
             ( PS, 'GENERIC', &
               NameShortOption = 'Reference' // EnergyNumber, &
               AllocateSourcesOption = .false., &
               IgnorabilityOption = CONSOLE % INFO_5 )
      call PWS % Difference_ASC ( iE ) % Initialize &
             ( PS, 'GENERIC', &
               NameShortOption = 'Difference' // EnergyNumber, &
               AllocateSourcesOption = .false., &
               IgnorabilityOption = CONSOLE % INFO_5 )
    end do !-- iE
    PWS % SetReference => SetReference

    !-- Initial conditions

    call SetRadiation ( PWS, Time = 0.0_KDR, Period = Period )

    !-- Initialize template

    call PWS % InitializeTemplate_C_MS_C_PS &
           ( Name, FinishTimeOption = Period )

    !-- Cleanup

    end select !-- FA
    end select !-- RMB
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

    select type ( RMB => PWS % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    do iE = 1, RMB % nEnergyValues
      associate ( RSA_D => PWS % Difference_ASC ( iE ) )
      RS_D => RSA_D % RadiationMoments ( )
      allocate ( CO )

      associate &
        ( Difference => RS_D % Value ( :, RS_D % COMOVING_ENERGY ) )
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

    end select !-- RMB
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

    call PWS % FinalizeTemplate_C_MS_C_PS ( )

  end subroutine FinalizeTemplate_PWS


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

    select type ( RMB => PWS % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    do iE = 1, RMB % nEnergyValues
      select type ( RSA => RMB % Section % Atlas ( iE ) % Element )
      class is ( RadiationMoments_ASC_Form )
      associate &
        ( RSA_R => PWS % Reference_ASC ( iE ), &
          RSA_D => PWS % Difference_ASC ( iE ) )
      RS   => RSA   % RadiationMoments ( )
      RS_R => RSA_R % RadiationMoments ( )
      RS_D => RSA_D % RadiationMoments ( )

      !-- Reference
      call SetWave ( PWS, RS_R, MS % Base_CSLD, PWS % Time, nWavelengths = iE )

      !-- Difference
 !    RS_D % Value  =  RS % Value  -  RS_R % Value
      call MultiplyAdd ( RS % Value, RS_R % Value, -1.0_KDR, RS_D % Value )
      end associate !-- RSA_R, etc.
      end select !-- RSA
    end do !-- iE

    end select !-- RMB
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
    class ( RadiationMomentsForm ), pointer :: &
      RS

    call Show ( 'Setting SineWaveStreaming initial conditions' )

    select type ( RMB => PWS % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    select type ( MS => PWS % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    G => MS % Base_CSLD % Geometry ( )

    iE = 1
    select type ( RSA => RMB % Section % Atlas ( 1 ) % Element )
    class is ( RadiationMoments_ASC_Form )
      RS => RSA % RadiationMoments ( )
      call SetWave ( PWS, RS, MS % Base_CSLD, Time, nWavelengths = iE, &
                     PeriodOption = Period )
    end select !-- RSA

    do iE = 2, RMB % nEnergyValues
      select type ( RSA => RMB % Section % Atlas ( iE ) % Element )
      class is ( RadiationMoments_ASC_Form )
        RS => RSA % RadiationMoments ( )
        call SetWave ( PWS, RS, MS % Base_CSLD, Time, nWavelengths = iE )
      end select !-- RSA
    end do !-- iE

    call RMB % StoreSections ( )

    end select !-- MS
    end select !-- RMB
    nullify ( G, RS )

  end subroutine SetRadiation


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
             J  = RS % Value ( :, RS % COMOVING_ENERGY ), &
             HX = RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 1 ) ), &
             HY = RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 2 ) ), &
             HZ = RS % Value ( :, RS % COMOVING_MOMENTUM_U ( 3 ) ) )

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
