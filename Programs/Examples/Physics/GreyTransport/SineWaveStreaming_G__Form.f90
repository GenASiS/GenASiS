module SineWaveStreaming_G__Form

  use Basics
  use PlaneWaveStreaming_G__Template

  implicit none
  private

  type, public, extends ( PlaneWaveStreaming_G_Template ) :: &
    SineWaveStreaming_G_Form
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
  end type SineWaveStreaming_G_Form

contains


  subroutine Initialize ( SWS, Name )

    class ( SineWaveStreaming_G_Form ), intent ( inout ) :: &
      SWS
    character ( * ), intent ( in )  :: &
      Name

    SWS % Offset    = 2.0_KDR
    SWS % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( SWS % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( SWS % Amplitude, 'Amplitude' )

    call SWS % InitializeTemplate_PWS ( Name )

  end subroutine Initialize


  elemental function Waveform ( PWS, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveStreaming_G_Form ), intent ( in ) :: &
      PWS
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    associate &
      ( O => PWS % Offset, &
        A => PWS % Amplitude, &
        TwoPi => 2.0_KDR * CONSTANT % PI )

    W  =  O  +  A * sin ( TwoPi * X )

    end associate !-- O, etc.

  end function Waveform


  impure elemental subroutine Finalize ( SWS )
    
    type ( SineWaveStreaming_G_Form ), intent ( inout ) :: &
      SWS

    call SWS % FinalizeTemplate_PWS ( )

  end subroutine Finalize


end module SineWaveStreaming_G__Form
