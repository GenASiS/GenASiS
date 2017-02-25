module PlaneWave_Template

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Fluid_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    PlaneWaveTemplate
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
      SetFluid, &
      SetReference

      private :: &
        SetFluidKernel

contains


  subroutine InitializeTemplate_PW ( PW, Name )

    class ( PlaneWaveTemplate ), intent ( inout ), target :: &
      PW
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period
    class ( Fluid_D_Form ), pointer :: &
      F

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: PW % PositionSpace )
    select type ( PS => PW % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( PW % Geometry_ASC )
    associate ( GA => PW % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Fluid (Dust, i.e. pressureless fluid)

    allocate ( Fluid_ASC_Form :: PW % Current_ASC )
    select type ( FA => PW % Current_ASC )  !-- FluidAtlas
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'DUST' )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: PW % Step )
    select type ( S => PW % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( Name )
    end select !-- S

    !-- Diagnostics

    allocate ( PW % Reference )
    allocate ( PW % Difference )
    call PW % Reference % Initialize &
           ( PS, 'DUST', NameOutputOption = 'Reference' )
    call PW % Difference % Initialize &
           ( PS, 'DUST', NameOutputOption = 'Difference' )
    PW % SetReference => SetReference

    !-- Initial conditions

    nWavelengths = 0
    nWavelengths ( 1 : PS % nDimensions ) = 1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    associate ( C => PS % Chart )
    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      PW % Wavenumber = nWavelengths / BoxSize
    elsewhere
      PW % Wavenumber = 0.0_KDR
    end where
    end associate !-- BoxSize
    end associate !-- C

    PW % Speed = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PW % Speed, 'Speed' )

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )

    associate &
      ( K     => PW % Wavenumber, &
        Abs_K => sqrt ( dot_product ( PW % Wavenumber, PW % Wavenumber ) ), &
        V     => PW % Speed )
    Period = 1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )
    end associate !-- K, etc.

    F => FA % Fluid_D ( )
    call SetFluid ( PW, F, Time = 0.0_KDR )

    !-- Initialize template

    call PW % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption = nPeriods * Period )

    !-- Cleanup

    end select !-- FA
    end select !-- PS
    nullify ( F )

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
    
    select type ( PS => PW % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    F_D => PW % Difference % Fluid_D ( )

    associate ( Difference => F_D % Value ( :, F_D % COMOVING_DENSITY ) )
    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )
    end associate !-- Difference

    associate &
      ( DifferenceSum => CO % Incoming % Value ( 1 ), &
        nValues => CO % Incoming % Value ( 2 ) )
    L1 = DifferenceSum / nValues
    end associate

    call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end select !-- C
    end select !-- PS

  end subroutine ComputeError


  impure elemental subroutine FinalizeTemplate_PW ( PW )

    class ( PlaneWaveTemplate ), intent ( inout ) :: &
      PW

    if ( allocated ( PW % Difference ) ) &
      deallocate ( PW % Difference )
    if ( allocated ( PW % Reference ) ) &
      deallocate ( PW % Reference )

    call PW % FinalizeTemplate_C_PS

  end subroutine FinalizeTemplate_PW


  subroutine SetFluid ( PW, F, Time )

    class ( PlaneWaveTemplate ), intent ( in ) :: &
      PW
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => PW % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetFluidKernel &
           ( PW, &
             X = G % Value ( :, G % CENTER ( 1 ) ), &
             Y = G % Value ( :, G % CENTER ( 2 ) ), &
             Z = G % Value ( :, G % CENTER ( 3 ) ), &
             K = PW % Wavenumber, &
             V = PW % Speed, &
             T = Time, &
             N  = F % Value ( :, F % COMOVING_DENSITY ), &
             VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetFluid


  subroutine SetReference ( PW )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      PW

    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_R, &  !-- F_Reference
      F_D     !-- F_Difference

    select type ( PW )
    class is ( PlaneWaveTemplate )

    select type ( FA => PW % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )
    end select !-- FA

    F_R => PW % Reference % Fluid_D ( )
    call SetFluid ( PW, F_R, PW % Time )

    F_D => PW % Difference % Fluid_D ( )
!    F_D % Value  =  F % Value  -  F_R % Value
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

    end select !-- PW
    nullify ( F, F_R, F_D )

  end subroutine SetReference


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
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      N ( iV ) = PW % Waveform &
                   (    K ( 1 ) * ( X ( iV )  -  VX ( iV ) * T ) &
                     +  K ( 2 ) * ( Y ( iV )  -  VY ( iV ) * T ) &
                     +  K ( 3 ) * ( Z ( iV )  -  VZ ( iV ) * T ) )
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetFluidKernel


end module PlaneWave_Template
