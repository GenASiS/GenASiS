module PlaneWaveStreaming_Template

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ), abstract :: &
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
      SetRadiation

      private :: &
        SetWave

        private :: &
          SetWaveKernel

   class ( PlaneWaveStreamingTemplate ), private, pointer :: &
     PlaneWaveStreaming => null ( )

contains


  subroutine InitializeTemplate_PWS ( PWS, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iE, &  !-- iEnergy
      nPeriods
    real ( KDR ) :: &
      Period
    character ( 1 + 2 ) :: &
      EnergyNumber


    if ( PWS % Type == '' ) &
      PWS % Type = 'a PlaneWaveStreaming'

    PlaneWaveStreaming => PWS

    call PWS % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( SpectralRadiationBoxForm :: PWS % Integrator )
    select type ( SRB => PWS % Integrator )
    type is ( SpectralRadiationBoxForm )
    call SRB % Initialize &
           ( RadiationName = [ 'Radiation' ], RadiationType = [ 'GENERIC' ], &
             Name = Name )
    SRB % SetReference => SetReference

    select type ( PS => SRB % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( RMB => SRB % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    class is ( RadiationMoments_BSLL_ASC_CSLD_Form )


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


    !-- Initial conditions

    call SetRadiation ( PWS, Time = 0.0_KDR, Period = Period )

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )
    call Show ( nPeriods, 'nPeriods' )

    SRB % FinishTime = nPeriods * Period
    call Show ( SRB % FinishTime, 'Reset FinishTime' )


    !-- Cleanup

    end select !-- RMB
    end select !-- PS
    end select !-- SRB

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

    select type ( PS => PWS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( CB => PS % Chart )
    class is ( Chart_SL_Template )

    select type ( I => PWS % Integrator )
    class is ( Integrator_C_1D_MS_C_PS_Template )

    select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
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
    end select !-- I
    end select !-- CB
    end select !-- PS
    ! nullify ( RS_D )

  end subroutine ComputeError


  subroutine FinalizeTemplate_PWS ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS

    if ( allocated ( PWS % Difference_ASC ) ) &
      deallocate ( PWS % Difference_ASC )
    if ( allocated ( PWS % Reference_ASC ) ) &
      deallocate ( PWS % Reference_ASC )

    call PWS % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_PWS


  subroutine SetReference ( SRB )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      SRB

    integer ( KDI ) :: &
      iE  !-- iEnergy
    class ( RadiationMomentsForm ), pointer :: &
      RS, &
      RS_R, &  !-- RM_Reference
      RS_D     !-- RM_Difference

    select type ( SRB )
    class is ( SpectralRadiationBoxForm )

    select type ( MS => SRB % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    select type ( RMB => SRB % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    do iE = 1, RMB % nEnergyValues
      select type ( RSA => RMB % Section % Atlas ( iE ) % Element )
      class is ( RadiationMoments_ASC_Form )
      associate &
        ( RSA_R => PlaneWaveStreaming % Reference_ASC ( iE ), &
          RSA_D => PlaneWaveStreaming % Difference_ASC ( iE ) )
      RS   => RSA   % RadiationMoments ( )
      RS_R => RSA_R % RadiationMoments ( )
      RS_D => RSA_D % RadiationMoments ( )

      !-- Reference
      call SetWave ( PlaneWaveStreaming, RS_R, MS % Base_CSLD, SRB % Time, &
                     nWavelengths = iE )

      !-- Difference
      call MultiplyAdd ( RS % Value, RS_R % Value, -1.0_KDR, RS_D % Value )
      end associate !-- RSA_R, etc.
      end select !-- RSA
    end do !-- iE

    end select !-- RMB
    end select !-- MS
    end select !-- 
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
    class ( RadiationMomentsForm ), pointer :: &
      RS

    select type ( I => PWS % Integrator )
    class is ( Integrator_C_1D_MS_C_PS_Template )

    select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

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
    end select !-- I
    nullify ( RS )

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
             X  = G % Value ( :, G % CENTER_U ( 1 ) ), &
             Y  = G % Value ( :, G % CENTER_U ( 2 ) ), &
             Z  = G % Value ( :, G % CENTER_U ( 3 ) ), &
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
