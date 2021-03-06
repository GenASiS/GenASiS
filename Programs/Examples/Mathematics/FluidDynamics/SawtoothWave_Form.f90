module SawtoothWave_Form

  use Basics
  use PlaneWave_Template

  implicit none
  private

  type, public, extends ( PlaneWaveTemplate ) :: SawtoothWaveForm
      real ( KDR ) :: &
        Offset, &
        Amplitude
  contains
    procedure, private, pass :: &
      Initialize_SW
    generic, public :: &
      Initialize => Initialize_SW
    procedure, private, pass :: &
      Waveform
    final :: &
      Finalize
  end type SawtoothWaveForm

contains


  subroutine Initialize_SW ( SW, Name )

    class ( SawtoothWaveForm ), intent ( inout ) :: &
      SW
    character ( * ), intent ( in )  :: &
      Name

    SW % Offset    = 2.0_KDR
    SW % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( SW % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( SW % Amplitude, 'Amplitude' )

    call SW % InitializeTemplate_PW ( Name )

  end subroutine Initialize_SW


  elemental function Waveform ( PW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SawtoothWaveForm ), intent ( in ) :: &
      PW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W

    associate &
      ( O => PW % Offset, &
        A => PW % Amplitude, &
        Pi => CONSTANT % PI )

    W  =  O  -  2.0_KDR * A / Pi * atan ( 1.0_KDR / tan ( Pi * X ) )

    end associate !-- O, etc.

  end function Waveform


  impure elemental subroutine Finalize ( SW )
    
    type ( SawtoothWaveForm ), intent ( inout ) :: &
      SW

    call SW % FinalizeTemplate_PW ( )

  end subroutine Finalize


end module SawtoothWave_Form
