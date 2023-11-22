module SineWaveStreaming_Form

  use GenASiS
  use PlaneWaveStreaming_Form

  implicit none
  private

  type, public, extends ( PlaneWaveStreamingForm ) :: SineWaveStreamingForm
    real ( KDR ) :: &
      Offset, &
      Amplitude
  contains
    procedure, public, pass :: &
      Initialize_PWS
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
    procedure, private, pass :: &
      Waveform
  end type SineWaveStreamingForm
 

contains


  subroutine Initialize_PWS ( PWS, FormalismType, Name )

    class ( SineWaveStreamingForm ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( PWS % Type  ==  '' ) &
      PWS % Type  =  'a SineWaveStreaming'

    PWS % Offset    = 2.0_KDR
    PWS % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PWS % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( PWS % Amplitude, 'Amplitude' )

    call PWS % PlaneWaveStreamingForm % Initialize_PWS ( FormalismType, Name )

    PWS % Integrator % System  =>  PWS

  end subroutine Initialize_PWS


  impure elemental subroutine Finalize ( SWS )
    
    type ( SineWaveStreamingForm ), intent ( inout ) :: &
      SWS

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( SineWaveStreamingForm ), intent ( in ) :: &
      U

    call U % PlaneWaveStreamingForm % ShowParameters ( )

    call Show ( U % Offset,    'Offset' )
    call Show ( U % Amplitude, 'Amplitude' )

  end subroutine ShowParameters


  function Waveform ( PWS, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveStreamingForm ), intent ( in ) :: &
      PWS
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    real ( KDR ) :: &
      TwoPi

    associate &
      ( O => PWS % Offset, &
        A => PWS % Amplitude )

    TwoPi  =  2.0_KDR * CONSTANT % PI

    W  =  O  +  A * sin ( TwoPi * X )

    end associate !-- O, etc.

  end function Waveform


end module SineWaveStreaming_Form
