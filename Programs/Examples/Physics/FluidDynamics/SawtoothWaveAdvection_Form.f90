module SawtoothWaveAdvection_Form

  use GenASiS
  use PlaneWaveAdvection_Form

  implicit none
  private

  type, public, extends ( PlaneWaveAdvectionForm ) :: SawtoothWaveAdvectionForm
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
  end type SawtoothWaveAdvectionForm


contains


  subroutine Initialize_H ( U, Name )

    class ( SawtoothWaveAdvectionForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a SawtoothWaveAdvection'

    U % Offset     =  2.0_KDR
    U % Amplitude  =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( U % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( U % Amplitude, 'Amplitude' )

    call U % PlaneWaveAdvectionForm % Initialize ( Name = 'SawtoothWaveAdvection' )
    
    U % Integrator % System  =>  U

  end subroutine Initialize_H


  subroutine Show_U ( U )

    class ( SawtoothWaveAdvectionForm ), intent ( in ) :: &
      U

    call U % PlaneWaveAdvectionForm % Show ( )

    call Show ( U % Offset,    'Offset' )
    call Show ( U % Amplitude, 'Amplitude' )

  end subroutine Show_U


  impure elemental subroutine Finalize ( SWA )
    
    type ( SawtoothWaveAdvectionForm ), intent ( inout ) :: &
      SWA

  end subroutine Finalize


  function Waveform ( PWA, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SawtoothWaveAdvectionForm ), intent ( in ) :: &
      PWA
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    real ( KDR ) :: &
      Pi

    associate &
      ( O => PWA % Offset, &
        A => PWA % Amplitude )

    Pi  =  CONSTANT % PI

    W  =  O  -  2.0_KDR * A / Pi * atan ( 1.0_KDR / tan ( Pi * X ) )

    end associate !-- O, etc.

  end function Waveform


end module SawtoothWaveAdvection_Form
