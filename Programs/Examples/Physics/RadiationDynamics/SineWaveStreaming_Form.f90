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
    procedure, private, pass :: &
      Initialize_PWS
    final :: &
      Finalize
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

    call PWS % PlaneWaveStreamingForm % Initialize ( FormalismType, Name )

  end subroutine Initialize_PWS


  impure elemental subroutine Finalize ( SWS )
    
    type ( SineWaveStreamingForm ), intent ( inout ) :: &
      SWS

  end subroutine Finalize


end module SineWaveStreaming_Form
