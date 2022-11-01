#include "Preprocessor"

module PlaneWaveStreaming_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_R_B_Form ) :: PlaneWaveStreamingForm
  contains
    procedure, private, pass :: &
      Initialize_PWS
    generic, public :: &
      Initialize => Initialize_PWS
    final :: &
      Finalize
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


  impure elemental subroutine Finalize ( PWS )

    type ( PlaneWaveStreamingForm ), intent ( inout ) :: &
      PWS

    ! if ( allocated ( PW % Difference ) ) &
    !   deallocate ( PW % Difference )
    ! if ( allocated ( PW % Reference ) ) &
    !   deallocate ( PW % Reference )

  end subroutine Finalize


  subroutine InitializeUniverse ( PWS, FormalismType, Name )

    class ( PlaneWaveStreamingForm ), intent ( inout ) :: &
      PWS
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

!    integer ( KDI ) :: &
!      iD

    call PWS % Initialize &
           ( RadiationName = [ 'Radiation' ], &
             RadiationType = [ 'GENERIC' ], &
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
