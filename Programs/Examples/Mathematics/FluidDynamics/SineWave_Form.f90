module SineWave_Form

  use Basics
  use PlaneWave_Template

  implicit none
  private

  type, public, extends ( PlaneWaveTemplate ) :: SineWaveForm
    real ( KDR ) :: &
      Offset, &
      Amplitude
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Waveform
    final :: &
      Finalize
  end type SineWaveForm

contains


  subroutine Initialize ( SW, Name )

    class ( SineWaveForm ), intent ( inout ) :: &
      SW
    character ( * ), intent ( in )  :: &
      Name

    SW % Offset    = 2.0_KDR
    SW % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( SW % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( SW % Amplitude, 'Amplitude' )

    call SW % InitializeTemplate_PW ( Name )

  end subroutine Initialize


  elemental function Waveform ( PW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveForm ), intent ( in ) :: &
      PW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    associate &
      ( O => PW % Offset, &
        A => PW % Amplitude, &
        TwoPi => 2.0_KDR * CONSTANT % PI )

    W  =  O  +  A * sin ( TwoPi * X )

    end associate !-- O, etc.

  end function Waveform


  impure elemental subroutine Finalize ( SW )
    
    type ( SineWaveForm ), intent ( inout ) :: &
      SW

    call SW % FinalizeTemplate_PW ( )

  end subroutine Finalize


end module SineWave_Form
