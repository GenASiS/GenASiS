module PlaneWaveAdvection_Template

  use Basics
  use ConservationLawEvolution_Template
  use PressurelessFluid_Form

  implicit none
  private

  type, public, extends ( ConservationLawEvolutionTemplate), abstract :: &
    PlaneWaveAdvectionTemplate
  contains
    procedure, private, pass :: &
      Initialize_PWA
    generic, public :: &
      Initialize => Initialize_PWA
    procedure ( WaveformInterface ), public, pass, deferred :: &
      Waveform
  end type

  abstract interface
    elemental function WaveformInterface ( PWA, X ) result ( W )
      !-- Waveform with a full period in the range 0 < X < 1
      use Basics
      import PlaneWaveAdvectionTemplate
      class ( PlaneWaveAdvectionTemplate ), intent ( in ) :: &
        PWA
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ) :: &
        W
    end function WaveformInterface
  end interface

contains


  subroutine Initialize_PWA ( PWA, DensityUnit )

    class ( PlaneWaveAdvectionTemplate ), intent ( inout ) :: &
      PWA
    type ( MeasuredValueForm ), intent ( in ) :: &
      DensityUnit

    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit

    call PWA % Initialize ( PROGRAM_HEADER % Communicator )

    associate ( DM => PWA % DistributedMesh )

    allocate ( PressurelessFluidForm :: PWA % ConservedFields )
    select type ( PF => PWA % ConservedFields )
    type is ( PressurelessFluidForm )

    call PF % Initialize ( DM, NameOption = 'PressurelessFluid' )

    nWavelengths = 0
    nWavelengths ( 1 : DM % nDimensions ) = 1
    Period = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )
    call PROGRAM_HEADER % GetParameter ( Period, 'Period' )

    associate &
      ( X  => DM % Position % Value ( :, 1 ), &
        Y  => DM % Position % Value ( :, 2 ), &
        Z  => DM % Position % Value ( :, 3 ), &
        VX => PF % Value ( :, PF % VELOCITY ( 1 ) ), &
        VY => PF % Value ( :, PF % VELOCITY ( 2 ) ), &
        VZ => PF % Value ( :, PF % VELOCITY ( 3 ) ), &
        N  => PF % Value ( :, PF % COMOVING_DENSITY ) )

    associate ( L => DM % MaxCoordinate - DM % MinCoordinate )
    where ( L > 0.0_KDR )
      Wavenumber = nWavelengths / L
    elsewhere
      Wavenumber = 0.0_KDR
    end where
    end associate !-- L

    associate ( K => Wavenumber )

    !-- Set density according to the waveform specified in an extension of this
    !   template (abstract class)
    N = PWA % Waveform ( K ( 1 ) * X  +  K ( 2 ) * Y  +  K ( 3 ) * Z )
 
    !-- Set velocity such that the wave travels a full wavelength (i.e. returns
    !   to its initial condition) in the specified Period
    VX = K ( 1 ) / ( dot_product ( K, K ) * Period )
    VY = K ( 2 ) / ( dot_product ( K, K ) * Period )
    VZ = K ( 3 ) / ( dot_product ( K, K ) * Period )

    end associate !-- K

    call PF % ComputeAuxiliary ( PF % Value )
    call PF % ComputeConserved ( PF % Value )

    !-- FIXME: Implicit do-loop array needed to workaround Cray compiler
    !          crashes for elemental function argument of array and scalar 
    !-- FIXME: Implied do-loop for iD = 1, 3 fails to compile with NAG
    ! VelocityUnit &
    !   = [ ( DM % CoordinateUnit ( iD ) / PWA % TimeUnit, iD = 1, 3 ) ]
    do iD = 1, 3
      VelocityUnit = DM % CoordinateUnit ( iD ) / PWA % TimeUnit
    end do
    call PF % SetOutputPressureless &
           ( VelocityUnitOption = VelocityUnit, &
             DensityUnitOption = DensityUnit )

    end associate !-- X, etc.
    end select    !-- PF
    end associate !-- DM

  end subroutine Initialize_PWA


end module PlaneWaveAdvection_Template
