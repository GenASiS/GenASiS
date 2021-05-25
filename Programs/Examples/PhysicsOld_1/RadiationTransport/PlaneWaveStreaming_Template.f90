module PlaneWaveStreaming_Template

  use GenASiS

  implicit none
  private

  type, public, extends ( RadiationBoxForm ), abstract :: &
    PlaneWaveStreamingTemplate
      integer ( KDI ) :: &
        nWaves
      real ( KDR ) :: &
        Speed
      real ( KDR ), dimension ( 3 ) :: &
        Wavenumber
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
      use GenASiS
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
      InitializeRadiationBox, &
      InitializeDiagnostics, &
      SetProblem

      private :: &
        SetRadiation

        private :: &
          SetRadiationKernel

contains


  subroutine InitializeTemplate_PWS ( PWS, MomentsType, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    if ( PWS % Type == '' ) &
      PWS % Type = 'a PlaneWaveStreaming'

    call InitializeRadiationBox ( PWS, MomentsType, Name )
    call InitializeDiagnostics ( PWS, MomentsType )
    call SetProblem ( PWS )
 
  end subroutine InitializeTemplate_PWS


  subroutine ComputeError ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS

    integer ( KDI ) :: &
      iW  !-- iWave
    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO
    class ( RadiationMomentsForm ), pointer :: &
      RM_D     !-- RM_Difference

    select type ( PS => PWS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( PSC => PS % Chart )
    class is ( Chart_SL_Template )

    do iW = 1, PWS % nWaves
      associate ( RMA_D => PWS % Difference_ASC ( iW ) )
      RM_D => RMA_D % RadiationMoments ( )
      allocate ( CO )

      associate &
        ( Difference => RM_D % Value ( :, RM_D % COMOVING_ENERGY ) )
      call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
      CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference ), &
                                          mask = PSC % IsProperCell )
      CO % Outgoing % Value ( 2 ) = PSC % nProperCells
      call CO % Reduce ( REDUCTION % SUM )
      end associate !-- Difference

      associate &
        ( DifferenceSum => CO % Incoming % Value ( 1 ), &
          nValues => CO % Incoming % Value ( 2 ) )
      L1 = DifferenceSum / nValues
      end associate 

      call Show ( iW, '*** iWave', nLeadingLinesOption  = 2 ) 
      call Show ( L1, '*** L1 error',   nTrailingLinesOption = 2 )

      deallocate ( CO )
      end associate !-- RSA_D
    end do !-- iW

    end select !-- PSC
    end select !-- PS

  end subroutine ComputeError


  impure elemental subroutine FinalizeTemplate_PWS ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS

    if ( allocated ( PWS % Difference_ASC ) ) &
      deallocate ( PWS % Difference_ASC )
    if ( allocated ( PWS % Reference_ASC ) ) &
      deallocate ( PWS % Reference_ASC )

    call PWS % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_PWS


  subroutine SetReference ( I )
    
    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iW  !-- iWave
    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( PWS => I % Universe )
    class is ( PlaneWaveStreamingTemplate )
    
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( PSC => PS % Chart )
    class is ( Chart_SL_Template )

    do iW = 1, PWS % nWaves

      select type ( I )
      class is ( Integrator_C_1D_PS_C_PS_Form )
        select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
        class is ( RadiationMoments_ASC_Form )
        RM => RMA % RadiationMoments ( )
        end select !-- RMA
      class is ( Integrator_C_1D_MS_C_PS_Form )
        select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
        class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
          select type ( RMA => RMB % Section % Atlas ( iW ) % Element )
          class is ( RadiationMoments_ASC_Form )
          RM => RMA % RadiationMoments ( )
          end select !-- RMA
        end select !-- RMB
      end select !-- I

      associate &
        ( RMA_R => PWS % Reference_ASC ( iW ), &
          RMA_D => PWS % Difference_ASC ( iW ) )
      RM_R => RMA_R % RadiationMoments ( )
      RM_D => RMA_D % RadiationMoments ( )

      !-- Reference
      call SetRadiation ( PWS, RM_R, PSC, nWavelengths = iW )

      !-- Difference
      call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

      end associate !-- RMA_R, etc.
    end do !-- iW

    end select !-- PSC
    end select !-- PS
    end select !-- PWS
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


  subroutine InitializeRadiationBox ( PWS, MomentsType, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS
    character ( * ), intent ( in )  :: &
      MomentsType, &
      Name

    call PWS % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'GENERIC' ], &
             MomentsType = MomentsType, &
             Name = Name, &
             EnergySpacingOption = 'COMPACTIFIED', &
             ApplyInteractionsOption = .false., &
             EvolveFluidOption = .false., &
             nCellsPositionOption = [ 128, 128, 128 ], &
             nCellsEnergyOption = 4 )
    PWS % Integrator % SetReference => SetReference

    select type ( I => PWS % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )
      PWS % nWaves  =  1
    class is ( Integrator_C_1D_MS_C_PS_Form )
      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )
        PWS % nWaves  =  MS % nSections 
      end select !
    end select !-- I

  end subroutine InitializeRadiationBox


  subroutine InitializeDiagnostics ( PWS, MomentsType )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS
    character ( * ), intent ( in ) :: &
      MomentsType

    integer ( KDI ) :: &
      iW  !-- iWave
    character ( 1 + 2 ) :: &
      WaveIndex

    select type ( PS => PWS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( PWS % Reference_ASC ( PWS % nWaves ) )
    allocate ( PWS % Difference_ASC ( PWS % nWaves ) )
    do iW = 1, PWS % nWaves
      if ( PWS % nWaves > 1 ) then
        write ( WaveIndex, fmt = '(a1,i2.2)' ) '_', iW
      else
        WaveIndex = ''
      end if
      call PWS % Reference_ASC ( iW ) % Initialize &
             ( PS, 'GENERIC', MomentsType = MomentsType, Units = PWS % Units, &
               NameShortOption = 'Reference' // WaveIndex, &
               AllocateSourcesOption = .false., &
               IgnorabilityOption = CONSOLE % INFO_5 )
      call PWS % Difference_ASC ( iW ) % Initialize &
             ( PS, 'GENERIC', MomentsType = MomentsType, Units = PWS % Units, &
               NameShortOption = 'Difference' // WaveIndex, &
               AllocateSourcesOption = .false., &
               IgnorabilityOption = CONSOLE % INFO_5 )
    end do !-- iE

    end select !-- PS

  end subroutine InitializeDiagnostics


  subroutine SetProblem ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS
    
    integer ( KDI ) :: &
      iE, &  !-- iEnergy
      nPeriods
    real ( KDR ) :: &
      Period
    class ( RadiationMomentsForm ), pointer :: &
      RM

    select type ( I => PWS % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      select type ( PS => I % PositionSpace )
      class is ( Atlas_SC_Form )

      select type ( PSC => PS % Chart )
      class is ( Chart_SL_Template )

      select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )

      RM => RMA % RadiationMoments ( )
      call SetRadiation &
             ( PWS, RM, PSC, nWavelengths = 1, PeriodOption = Period )

      end select !-- RMA
      end select !-- PSC
      end select !-- PS

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      select type ( MS => I % MomentumSpace )
      class is ( Bundle_SLL_ASC_CSLD_Form )

      select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

      iE = 1
      select type ( RMA => RMB % Section % Atlas ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )
        RM => RMA % RadiationMoments ( )
        call SetRadiation &
               ( PWS, RM, MS % Base_CSLD, nWavelengths = iE, &
                 PeriodOption = Period )
      end select !-- RMA

      do iE = 2, RMB % nEnergyValues
        select type ( RMA => RMB % Section % Atlas ( iE ) % Element )
        class is ( RadiationMoments_ASC_Form )
          RM => RMA % RadiationMoments ( )
          call SetRadiation &
                 ( PWS, RM, MS % Base_CSLD, nWavelengths = iE )
        end select !-- RMA
      end do !-- iE

      call RMB % StoreSections ( )

      end select !-- MS
      end select !-- RMB

    end select !-- I

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )
    call Show ( nPeriods, 'nPeriods' )

    associate ( I => PWS % Integrator )
    I % FinishTime = nPeriods * Period
    call Show ( I % FinishTime, 'Reset FinishTime' )
    end associate !-- I

    nullify ( RM )

  end subroutine SetProblem


  subroutine SetRadiation ( PWS, RM, PSC, nWavelengths, PeriodOption )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    type ( RadiationMomentsForm ), intent ( inout ) :: &
      RM  !-- RadiationSection
    class ( Chart_SL_Template ), intent ( in ) :: &
      PSC  !-- PositionSpaceChart
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

    call SetRadiationKernel &
           ( PWS, &
             X  = G % Value ( :, G % CENTER_U ( 1 ) ), &
             Y  = G % Value ( :, G % CENTER_U ( 2 ) ), &
             Z  = G % Value ( :, G % CENTER_U ( 3 ) ), &
             K  = K, &
             V  = c, &
             T  = PWS % Integrator % Time, &
             J  = RM % Value ( :, RM % COMOVING_ENERGY ), &
             HX = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
             HY = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
             HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ) )

    call RM % ComputeFromPrimitive ( G )

    end associate !-- Abs_K
    end associate !-- K, etc.
    nullify ( G )

  end subroutine SetRadiation


  subroutine SetRadiationKernel ( PWS, X, Y, Z, K, V, T, J, HX, HY, HZ )

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

      HX ( iV )  =  VX  *  J ( iV )
      HY ( iV )  =  VY  *  J ( iV )
      HZ ( iV )  =  VZ  *  J ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


end module PlaneWaveStreaming_Template
