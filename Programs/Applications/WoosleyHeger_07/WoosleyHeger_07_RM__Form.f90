module WoosleyHeger_07_RM__Form

  !-- WoosleyHeger_07_RadiationMoments_Form

  use GenASiS
  use WoosleyHeger_07__Form

  implicit none
  private

  type, public, extends ( WoosleyHeger_07_Form ) :: WoosleyHeger_07_RM_Form
  contains
    procedure, private, pass :: &
      Initialize_RM
    generic, public :: &
      Initialize => Initialize_RM
    final :: &
      Finalize
  end type WoosleyHeger_07_RM_Form

    private :: &
      InitializeUniverse, &
      SetInitial


contains


  subroutine Initialize_RM ( U, FormalismType, Name )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a WoosleyHeger_07_RM'

    call InitializeUniverse ( U, FormalismType, Name )

  end subroutine Initialize_RM


  subroutine Finalize ( WH )

    type ( WoosleyHeger_07_RM_Form ), intent ( inout ) :: &
      WH

  end subroutine Finalize


  subroutine InitializeUniverse ( WH, FormalismType, Name )

    class ( WoosleyHeger_07_RM_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      FormalismType, &
      Name

    real ( KDR ) :: &
      FinishTime

    FinishTime  =  0.7_KDR  *  UNIT % SECOND

    !-- Initialization

    call WH % Initialize &
           ( RadiationName = [ 'Neutrinos_E    ', &
                               'Neutrinos_E_Bar' ], &
             RadiationType = [ 'NEUTRINOS_E    ', &
                               'NEUTRINOS_E_BAR' ], &
             FormalismType = FormalismType, &
             FluidType = 'HEAVY_NUCLEUS', &
             GravitationType = 'NEWTON_SG', &
             Name = Name, &
             FinishTimeOption = FinishTime, &
             nCellsPolarOption = 128, &
             nWriteOption = 30 )

    WH % Integrator % SetInitial  =>  SetInitial
    WH % Integrator % System      =>  WH

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( WH  =>  I % System )
      class is ( WoosleyHeger_07_Form )

    call WH % SetFluid ( )

    end select !-- WH

  end subroutine SetInitial


end module WoosleyHeger_07_RM__Form
