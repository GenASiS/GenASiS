module PlaneWaveStreaming_Template

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ), abstract :: &
    PlaneWaveStreamingTemplate
      real ( KDR ) :: &
        Speed
      real ( KDR ), dimension ( 3 ) :: &
        Wavenumber
      type ( RadiationMoments_ASC_Form ), allocatable :: &
        Reference, &
        Difference
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
        SetRadiationKernel

   class ( PlaneWaveStreamingTemplate ), private, pointer :: &
     PlaneWaveStreaming => null ( )

contains


  subroutine InitializeTemplate_PWS ( PWS, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period
    class ( RadiationMomentsForm ), pointer :: &
      RM

    if ( PWS % Type == '' ) &
      PWS % Type = 'a PlaneWaveStreaming'

    PlaneWaveStreaming => PWS

    call PWS % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( GreyRadiationBoxForm :: PWS % Integrator )
    select type ( GRB => PWS % Integrator )
    type is ( GreyRadiationBoxForm )
    call GRB % Initialize &
           ( RadiationName = [ 'Radiation' ], RadiationType = [ 'GENERIC' ], &
             Name = Name )
    GRB % SetReference => SetReference

    select type ( PS => GRB % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( RMA => GRB % Current_ASC_1D ( 1 ) % Element )
    class is ( RadiationMoments_ASC_Form )


    !-- Diagnostics

    allocate ( PWS % Reference )
    allocate ( PWS % Difference )
    call PWS % Reference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call PWS % Difference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )


    !-- Parameters

    nWavelengths = 0
    nWavelengths ( 1 : PS % nDimensions ) = 1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    associate ( C => PS % Chart )
    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      PWS % Wavenumber = nWavelengths / BoxSize
    elsewhere
      PWS % Wavenumber = 0.0_KDR
    end where

    PWS % Speed = CONSTANT % SPEED_OF_LIGHT

    associate &
      ( K     => PWS % Wavenumber, &
        Abs_K => sqrt ( dot_product ( PWS % Wavenumber, PWS % Wavenumber ) ), &
        V     => PWS % Speed )
    Period = 1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )
    call Show ( nPeriods, 'nPeriods' )

    GRB % FinishTime = nPeriods * Period
    call Show ( GRB % FinishTime, 'Reset FinishTime' )


    !-- Initial conditions

    RM => RMA % RadiationMoments ( )
    call SetRadiation ( PWS, RM, Time = 0.0_KDR )


    !-- Cleanup

    end associate !-- K, etc.
    end associate !-- BoxSize
    end associate !-- C
    end select !-- RMA
    end select !-- PS
    end select !-- GRB
    nullify ( RM )

  end subroutine InitializeTemplate_PWS


  subroutine ComputeError ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS

    real ( KDR ) :: &
      L1
    class ( RadiationMomentsForm ), pointer :: &
      RM
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( PS => PWS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    RM => PWS % Difference % RadiationMoments ( )

    associate &
      ( Difference => RM % Value ( :, RM % COMOVING_ENERGY ) )
    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )

    associate &
      ( DifferenceSum => CO % Incoming % Value ( 1 ), &
        nValues => CO % Incoming % Value ( 2 ) )
    L1 = DifferenceSum / nValues

    call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- DifferenceSum, etc.
    end associate !-- Difference
    end select !-- C
    end select !-- PS

  end subroutine ComputeError


  impure elemental subroutine FinalizeTemplate_PWS ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS

    if ( allocated ( PWS % Difference ) ) &
      deallocate ( PWS % Difference )
    if ( allocated ( PWS % Reference ) ) &
      deallocate ( PWS % Reference )

    call PWS % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_PWS


  subroutine SetReference ( GRB )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      GRB

    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( GRB )
    class is ( GreyRadiationBoxForm )

    select type ( RMA => GRB % Current_ASC_1D ( 1 ) % Element )
    class is ( RadiationMoments_ASC_Form )
    RM => RMA % RadiationMoments ( )
    end select !-- RMA

    RM_R => PlaneWaveStreaming % Reference % RadiationMoments ( )
    call SetRadiation ( PlaneWaveStreaming, RM_R, GRB % Time )

    RM_D => PlaneWaveStreaming % Difference % RadiationMoments ( )
    call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- GRB
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


  subroutine SetRadiation ( PWS, RM, Time )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      Time
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => PWS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetRadiationKernel &
           ( PWS, &
             X  = G % Value ( :, G % CENTER_U ( 1 ) ), &
             Y  = G % Value ( :, G % CENTER_U ( 2 ) ), &
             Z  = G % Value ( :, G % CENTER_U ( 3 ) ), &
             K  = PWS % Wavenumber, &
             V  = PWS % Speed, &
             T  = Time, &
             J  = RM % Value ( :, RM % COMOVING_ENERGY ), &
             HX = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
             HY = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
             HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ) )

    call RM % ComputeFromPrimitive ( G )

    end select    !-- PS
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
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      HX ( iV )  =  VX  *  J ( iV )
      HY ( iV )  =  VY  *  J ( iV )
      HZ ( iV )  =  VZ  *  J ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


end module PlaneWaveStreaming_Template
