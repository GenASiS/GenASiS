module SawtoothWaveAdvection_Form

  use Basics
  use PlaneWaveAdvection_Template

  implicit none
  private

  type, public, extends ( PlaneWaveAdvectionTemplate ) :: SawtoothWaveAdvectionForm
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
  end type SawtoothWaveAdvectionForm

contains


  subroutine Initialize_SWA ( SWA )
    
    class ( SawtoothWaveAdvectionForm ), intent ( inout ) :: &
      SWA

    type ( QuantityForm ) :: &
      OffsetUnit, &
      AmplitudeUnit

    if ( SWA % Type == '' ) SWA % Type = 'a SawtoothWaveAdvection' 
    
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

    class ( SawtoothWaveAdvectionForm ), intent ( in ) :: &
      PWA
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    real ( KDR ) :: &
      X_Safe, &
      Pi

    X_Safe = sign ( max ( tiny ( 0.0_KDR ), abs ( X ) ), X )
    
    associate &
      ( O => PWA % Offset, &
        A => PWA % Amplitude)
        ! Pi => CONSTANT % PI ) !-- FIXME: XL bug prevents correct association

    Pi =  CONSTANT % PI
    W  =  O  -  2.0_KDR * A / Pi * atan ( 1.0_KDR / tan ( Pi * X_Safe ) )

    end associate !-- O, etc.

  end function Waveform


  subroutine Finalize ( SWA )

    type ( SawtoothWaveAdvectionForm ), intent ( inout ) :: &
      SWA

    if ( allocated ( SWA % ConservedFields ) ) &
      deallocate ( SWA % ConservedFields )

  end subroutine Finalize


end module SawtoothWaveAdvection_Form
