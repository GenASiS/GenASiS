module PlaneWave_Template

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidBoxForm ), abstract :: PlaneWaveTemplate
    real ( KDR ) :: &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
    type ( Fluid_ASC_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, public, pass :: &
      InitializeTemplate_PW
    procedure, public, pass :: &
      ComputeError
    procedure, public, pass :: &
      FinalizeTemplate_PW
    procedure ( Waveform ), private, pass, deferred :: &
      Waveform
  end type PlaneWaveTemplate

  abstract interface
    elemental function Waveform ( PW, X ) result ( W )
      !-- Waveform with a full period in the range 0 < X < 1
      use Basics
      import PlaneWaveTemplate
      class ( PlaneWaveTemplate ), intent ( in ) :: &
        PW
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ) :: &
        W
    end function Waveform
  end interface

    private :: &
      SetReference, &
      InitializeFluidBox, &
      InitializeDiagnostics, &
      SetProblem

      private :: &
        SetFluid

        private :: &
          SetFluidKernel

contains


  subroutine InitializeTemplate_PW ( PW, Name )

    class ( PlaneWaveTemplate ), intent ( inout ) :: &
      PW
    character ( * ), intent ( in )  :: &
      Name

    if ( PW % Type == '' ) &
      PW % Type = 'a PlaneWave'

    call InitializeFluidBox ( PW, Name )
    call InitializeDiagnostics ( PW )
    call SetProblem ( PW )

  end subroutine InitializeTemplate_PW


  subroutine ComputeError ( PW )

    class ( PlaneWaveTemplate ), intent ( in ) :: &
      PW

    real ( KDR ) :: &
      L1
    class ( Fluid_D_Form ), pointer :: &
      F_D
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( PS => PW % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    F_D => PW % Difference % Fluid_D ( )

    associate &
      ( Difference => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ) )
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
    nullify ( F_D )

  end subroutine ComputeError


  impure elemental subroutine FinalizeTemplate_PW ( PW )

    class ( PlaneWaveTemplate ), intent ( inout ) :: &
      PW

    if ( allocated ( PW % Difference ) ) &
      deallocate ( PW % Difference )
    if ( allocated ( PW % Reference ) ) &
      deallocate ( PW % Reference )

  end subroutine FinalizeTemplate_PW


  subroutine SetReference ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_R, &  !-- F_Reference
      F_D     !-- F_Difference

    select type ( I )
    class is ( Integrator_C_PS_Form )

    select type ( PW => I % Universe )
    class is ( PlaneWaveTemplate )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    F_R => PW % Reference % Fluid_D ( )
    call SetFluid ( PW, F_R, I % Time )

    F_D => PW % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

    end select !-- FA
    end select !-- PW
    end select !-- I
    nullify ( F, F_R, F_D )

  end subroutine SetReference


  subroutine InitializeFluidBox ( PW, Name )

    class ( PlaneWaveTemplate ), intent ( inout ) :: &
      PW
    character ( * ), intent ( in )  :: &
      Name

    call PW % Initialize &
           ( FluidType = 'DUST', GeometryType = 'GALILEAN', Name = Name, &
             nCellsOption = [ 128, 128, 128 ] )
    PW % Integrator % SetReference => SetReference

  end subroutine InitializeFluidBox


  subroutine InitializeDiagnostics ( PW )

    class ( PlaneWaveTemplate ), intent ( inout ) :: &
      PW

    select type ( PS => PW % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( PW % Reference )
    allocate ( PW % Difference )
    call PW % Reference % Initialize &
           ( PS, 'DUST', PW % Units, NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call PW % Difference % Initialize &
           ( PS, 'DUST', PW % Units, NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    end select !-- PS

  end subroutine InitializeDiagnostics


  subroutine SetProblem ( PW )

    class ( PlaneWaveTemplate ), intent ( inout ) :: &
      PW

    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period
    class ( Fluid_D_Form ), pointer :: &
      F

    select type ( I => PW % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    nWavelengths = 0
    nWavelengths ( 1 : PS % nDimensions ) = 1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    associate ( PSC => PS % Chart )
    associate ( BoxSize => PSC % MaxCoordinate - PSC % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      PW % Wavenumber = nWavelengths / BoxSize
    elsewhere
      PW % Wavenumber = 0.0_KDR
    end where

    PW % Speed = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PW % Speed, 'Speed' )

    associate &
      ( K     => PW % Wavenumber, &
        Abs_K => sqrt ( dot_product ( PW % Wavenumber, PW % Wavenumber ) ), &
        V     => PW % Speed )
    Period = 1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )
    call Show ( nPeriods, 'nPeriods' )

    I % FinishTime = nPeriods * Period
    call Show ( I % FinishTime, 'Reset FinishTime' )

    F => FA % Fluid_D ( )
    call SetFluid ( PW, F, Time = 0.0_KDR )

    end associate !-- K, etc.
    end associate !-- BoxSize
    end associate !-- PSC
    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( F )

  end subroutine SetProblem


  subroutine SetFluid ( PW, F, Time )

    class ( PlaneWaveTemplate ), intent ( in ) :: &
      PW
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => PW % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetFluidKernel &
           ( PW, &
             X = G % Value ( :, G % CENTER_U ( 1 ) ), &
             Y = G % Value ( :, G % CENTER_U ( 2 ) ), &
             Z = G % Value ( :, G % CENTER_U ( 3 ) ), &
             K = PW % Wavenumber, &
             V = PW % Speed, &
             T = Time, &
             N  = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetFluid


  subroutine SetFluidKernel ( PW, X, Y, Z, K, V, T, N, VX, VY, VZ )

    class ( PlaneWaveTemplate ), intent ( in ) :: &
      PW
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      K
    real ( KDR ), intent ( in ) :: &
      V, &
      T
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      VX, VY, VZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      Abs_K

    nValues = size ( X )
    
    Abs_K = sqrt ( dot_product ( K, K ) )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      VX ( iV ) = V * K ( 1 ) / Abs_K
      VY ( iV ) = V * K ( 2 ) / Abs_K
      VZ ( iV ) = V * K ( 3 ) / Abs_K

      N ( iV ) = PW % Waveform &
                   (    K ( 1 ) * ( X ( iV )  -  VX ( iV ) * T ) &
                     +  K ( 2 ) * ( Y ( iV )  -  VY ( iV ) * T ) &
                     +  K ( 3 ) * ( Z ( iV )  -  VZ ( iV ) * T ) )

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetFluidKernel


end module PlaneWave_Template
