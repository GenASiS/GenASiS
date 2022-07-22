module SineWaveAdvection_Form

  use Basics
  use PlaneWaveAdvection_Template

  implicit none
  private

  type, public, extends ( PlaneWaveAdvectionTemplate ) :: SineWaveAdvectionForm
    type ( QuantityForm ) :: &
      Offset, &
      Amplitude
  contains
    procedure, private, pass :: &
      Initialize_SWA
    generic, public :: &
      Initialize => Initialize_SWA
    procedure, public, pass :: &
      Waveform
    final :: &
      Finalize
  end type SineWaveAdvectionForm

contains


  subroutine Initialize_SWA ( SWA )
    
    class ( SineWaveAdvectionForm ), intent ( inout ) :: &
      SWA

    type ( QuantityForm ) :: &
      OffsetUnit, &
      AmplitudeUnit

    if ( SWA % Type == '' ) SWA % Type = 'a SineWaveAdvection' 
    
    OffsetUnit      = UNIT % IDENTITY
    AmplitudeUnit   = UNIT % IDENTITY
    SWA % Offset    = 2.0_KDR * OffsetUnit
    SWA % Amplitude = 1.0_KDR * AmplitudeUnit
    call PROGRAM_HEADER % GetParameter &
           ( SWA % Offset, 'Offset', InputUnitOption = OffsetUnit )
    call PROGRAM_HEADER % GetParameter &
           ( SWA % Amplitude, 'Amplitude', InputUnitOption = AmplitudeUnit )

    if ( OffsetUnit /= AmplitudeUnit ) then
      call Show ( 'Please use the same units for Offset and Amplitude', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call SWA % Initialize ( AmplitudeUnit )

  end subroutine Initialize_SWA


  elemental function Waveform ( PWA, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( SineWaveAdvectionForm ), intent ( in ) :: &
      PWA
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    associate ( TwoPi => 2.0_KDR * CONSTANT % PI )

    W = PWA % Offset + PWA % Amplitude * sin ( TwoPi * X )

    end associate !-- TwoPi

  end function Waveform


  subroutine Finalize ( SWA )

    type ( SineWaveAdvectionForm ), intent ( inout ) :: &
      SWA

    if ( allocated ( SWA % ConservedFields ) ) &
      deallocate ( SWA % ConservedFields )

  end subroutine Finalize


end module SineWaveAdvection_Form
