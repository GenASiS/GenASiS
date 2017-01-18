module SineWaveStreaming_Form

  use Basics
  use PlaneWaveStreaming_Template

  implicit none
  private

  type, public, extends ( PlaneWaveStreamingTemplate ) :: &
    SineWaveStreamingForm
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
  end type SineWaveStreamingForm

contains


  subroutine Initialize ( SWS, Name )

    class ( SineWaveStreamingForm ), intent ( inout ) :: &
      SWS
    character ( * ), intent ( in )  :: &
      Name

    if ( SWS % Type == '' ) &
      SWS % Type = 'a SineWaveStreaming' 

    SWS % Offset    = 2.0_KDR
    SWS % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( SWS % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( SWS % Amplitude, 'Amplitude' )

    call SWS % InitializeTemplate_PWS ( Name )

  end subroutine Initialize


  elemental function Waveform ( PWS, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveStreamingForm ), intent ( in ) :: &
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
    
    type ( SineWaveStreamingForm ), intent ( inout ) :: &
      SWS

    call SWS % FinalizeTemplate_PWS ( )

  end subroutine Finalize


end module SineWaveStreaming_Form
