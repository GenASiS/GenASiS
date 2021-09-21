module SineWave_Form

  use GenASiS
  use PlaneWave_Form

  implicit none
  private

  type, public, extends ( PlaneWaveForm ) :: SineWaveForm
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
  end type SineWaveForm


contains


  subroutine Initialize_H ( U, NameOption )

    class ( SineWaveForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a SineWave'

    U % Offset     =  2.0_KDR
    U % Amplitude  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( U % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( U % Amplitude, 'Amplitude' )

    call U % PlaneWaveForm % Initialize ( NameOption = 'SineWave' )

    U % Integrator % System  =>  U

  end subroutine Initialize_H


  subroutine Show_U ( U )

    class ( SineWaveForm ), intent ( in ) :: &
      U

    call U % PlaneWaveForm % Show ( )

    call Show ( U % Offset,    'Offset' )
    call Show ( U % Amplitude, 'Amplitude' )

  end subroutine Show_U


  impure elemental subroutine Finalize ( SW )
    
    type ( SineWaveForm ), intent ( inout ) :: &
      SW

  end subroutine Finalize


  function Waveform ( PW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveForm ), intent ( in ) :: &
      PW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    real ( KDR ) :: &
      TwoPi

    associate &
      ( O => PW % Offset, &
        A => PW % Amplitude )

    TwoPi  =  2.0_KDR * CONSTANT % PI

    W  =  O  +  A * sin ( TwoPi * X )

    end associate !-- O, etc.

  end function Waveform


end module SineWave_Form
