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
!    procedure, private, pass :: &
!      Initialize_PWS
    final :: &
      Finalize
  end type SineWaveStreamingForm
 

contains


  impure elemental subroutine Finalize ( SWS )
    
    type ( SineWaveStreamingForm ), intent ( inout ) :: &
      SWS

  end subroutine Finalize


end module SineWaveStreaming_Form
