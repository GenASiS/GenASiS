module SineWaveAdvection_Form

  use GenASiS
  use PlaneWaveAdvection_Form

  implicit none
  private

  type, public, extends ( PlaneWaveAdvectionForm ) :: SineWaveAdvectionForm
    real ( KDR ) :: &
      Offset, &
      Amplitude
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      Show => Show_U
    final :: &
      Finalize
    procedure, private, pass :: &
      Waveform
  end type SineWaveAdvectionForm


contains


  subroutine Initialize_H ( U, Name )

    class ( SineWaveAdvectionForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a SineWaveAdvection'

    U % Offset     =  2.0_KDR
    U % Amplitude  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( U % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( U % Amplitude, 'Amplitude' )

    call U % PlaneWaveAdvectionForm % Initialize ( Name = 'SineWaveAdvection' )

    U % Integrator % System  =>  U

  end subroutine Initialize_H


  subroutine Show_U ( U )

    class ( SineWaveAdvectionForm ), intent ( in ) :: &
      U

    call U % PlaneWaveAdvectionForm % Show ( )

    call Show ( U % Offset,    'Offset' )
    call Show ( U % Amplitude, 'Amplitude' )

  end subroutine Show_U


  impure elemental subroutine Finalize ( SWA )
    
    type ( SineWaveAdvectionForm ), intent ( inout ) :: &
      SWA

  end subroutine Finalize


  function Waveform ( PWA, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveAdvectionForm ), intent ( in ) :: &
      PWA
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    real ( KDR ) :: &
      TwoPi

    associate &
      ( O => PWA % Offset, &
        A => PWA % Amplitude )

    TwoPi  =  2.0_KDR * CONSTANT % PI

    W  =  O  +  A * sin ( TwoPi * X )

    end associate !-- O, etc.

  end function Waveform


end module SineWaveAdvection_Form
