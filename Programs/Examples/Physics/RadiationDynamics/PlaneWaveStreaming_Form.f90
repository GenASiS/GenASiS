#include "Preprocessor"

module PlaneWaveStreaming_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: PlaneWaveStreamingForm
    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Speed, &
      Period
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
!    type ( Fluid_D_Form ), allocatable :: &
!      Reference, &
!      Difference
  contains
    procedure, private, pass :: &
      Initialize_PWS
    generic, public :: &
      Initialize => Initialize_PWS
    procedure, public, pass :: &
      Show => Show_U
    final :: &
      Finalize
    procedure, private, pass :: &
      Waveform
  end type PlaneWaveStreamingForm

    private :: &
      InitializeUniverse


contains


  subroutine Initialize_PWS ( PWS, FormalismType, Name )

    class ( PlaneWaveStreamingForm ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    if ( PWS % Type  ==  '' ) &
      PWS % Type  =  'a PlaneWaveStreaming'

    call InitializeUniverse ( PWS, FormalismType, Name )

    ! call InitializeRadiationBox ( PWS, FormalismType, Name )
    ! call InitializeDiagnostics ( PWS, FormalismType )
    ! call SetProblem ( PWS )
 
  end subroutine Initialize_PWS


  subroutine Show_U ( U )

    class ( PlaneWaveStreamingForm ), intent ( in ) :: &
      U

    call U % Universe_H_Form % Show ( )

    call Show ( 'PlaneWaveStreaming Parameters' )
    call Show ( U % nWavelengths, 'nWavelengths' )
    call Show ( U % nPeriods,     'nPeriods' )
    call Show ( U % Period,       'Period' )

  end subroutine Show_U


  impure elemental subroutine Finalize ( PWS )

    type ( PlaneWaveStreamingForm ), intent ( inout ) :: &
      PWS

    ! if ( allocated ( PW % Difference ) ) &
    !   deallocate ( PW % Difference )
    ! if ( allocated ( PW % Reference ) ) &
    !   deallocate ( PW % Reference )

  end subroutine Finalize


  function Waveform ( PWA, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( PlaneWaveStreamingForm ), intent ( in ) :: &
      PWA
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    W = huge ( 1.0_KDR ) 
    call Show ( 'Waveform should be overridden', CONSOLE % WARNING )
    call Show ( 'PlaneWaveAdvection_Form', 'module', CONSOLE % WARNING )
    call Show ( 'Waveform', 'function', CONSOLE % WARNING )

  end function Waveform


  subroutine InitializeUniverse ( PWS, FormalismType, Name )

    class ( PlaneWaveStreamingForm ), intent ( inout ) :: &
      PWS
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

!    integer ( KDI ) :: &
!      iD

    call PWS % Initialize &
           ( RadiationName = [ 'Radiation_1', 'Radiation_2' ], &
             RadiationType = [ 'GENERIC', 'GENERIC' ], &
             FormalismType = FormalismType, &
             Name = Name )
             ! EnergySpacingOption = 'COMPACTIFIED', &
             ! ApplyInteractionsOption = .false., &
             ! EvolveFluidOption = .false., &
             ! nCellsPositionOption = [ 128, 128, 128 ], &
             ! nCellsEnergyOption = 4 )

!    call PWS % Initialize &
!           ( FluidType = 'DUST', &
!             GravitationType = 'GALILEO', &
!             NameOption = Name, &
!             nCellsOption = [ 128, 128, 128 ] )

    ! select type ( I  =>  PW % Integrator )
    !   class is ( Integrator_CS_Form )
    ! associate &
    !   ( F  =>  I % CurrentSet_X )
    ! do iD  =  1, 3
    !   call F % SetBoundaryConditionsFace &
    !          ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    ! end do !-- iD
    ! end associate !-- F
    ! end select !-- I
             
    ! PW % Integrator % SetInitial    =>  SetInitial
    ! PW % Integrator % SetReference  =>  SetReference

  end subroutine InitializeUniverse


end module PlaneWaveStreaming_Form
