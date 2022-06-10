module SawtoothWave_Form

  use GenASiS
  use PlaneWave_Form

  implicit none
  private

  type, public, extends ( PlaneWaveForm ) :: SawtoothWaveForm
    real ( KDR ) :: &
      Offset, &
      Amplitude
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      Show => Show_U
    procedure, private, pass :: &
      Waveform
    final :: &
      Finalize
  end type SawtoothWaveForm


contains


  subroutine Initialize_H ( U, NameOption )

    class ( SawtoothWaveForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a SawtoothWave'

    U % Offset     =  2.0_KDR
    U % Amplitude  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( U % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( U % Amplitude, 'Amplitude' )

    call U % PlaneWaveForm % Initialize ( NameOption = 'SawtoothWave' )
    
    U % Integrator % System  =>  U

  end subroutine Initialize_H


  subroutine Show_U ( U )

    class ( SawtoothWaveForm ), intent ( in ) :: &
      U

    call U % PlaneWaveForm % Show ( )

    call Show ( U % Offset,    'Offset' )
    call Show ( U % Amplitude, 'Amplitude' )

  end subroutine Show_U


  impure elemental subroutine Finalize ( SW )
    
    type ( SawtoothWaveForm ), intent ( inout ) :: &
      SW

  end subroutine Finalize


  function Waveform ( PW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SawtoothWaveForm ), intent ( in ) :: &
      PW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    real ( KDR ) :: &
      Pi

    associate &
      ( O => PW % Offset, &
        A => PW % Amplitude )

    Pi  =  CONSTANT % PI

    W  =  O  -  2.0_KDR * A / Pi * atan ( 1.0_KDR / tan ( Pi * X ) )

    end associate !-- O, etc.

  end function Waveform


end module SawtoothWave_Form
