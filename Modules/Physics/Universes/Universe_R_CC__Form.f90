module Universe_R_CC__Form

  !-- Universe_Radiation_CentralCore__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Universe_F_CC__Form

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: Universe_R_CC_Form
    integer ( KDI ) :: &
      iRadiation  = 0, &
      nRadiations = 0
    real ( KDR ) :: &
      InteractionFactor
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
  !   type ( CommunicatorForm ), allocatable :: &
  !     Communicator_PS  !-- PositionSpace
  !   type ( CollectiveOperation_R_Form ), dimension ( : ), allocatable :: &
  !     CO_SplitSource
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
  !   class ( Interactions_BM_Form ), allocatable :: &
  !     Interactions_BM
  contains
    procedure, private, pass :: &
      Initialize_R_CC
    generic, public :: &
      Initialize => Initialize_R_CC
    final :: &
      Finalize
  !   procedure, private, pass :: &
  !     SetCommunicator
  !   procedure, private, pass :: &
  !     AllocateIntegrator
  !   procedure, public, pass :: &
  !     InitializePositionSpace
  !   procedure, public, pass :: &
  !     InitializeInteractions
  !   procedure, public, pass :: &
  !     InitializeRadiation
  !   procedure, public, pass :: &
  !     InitializeSteps
  !   procedure, public, pass :: &
  !     InitializeIntegrator
  !   procedure, public, pass :: &
  !     ShowParameters
  !   procedure, public, pass ( U ) :: &
  !     Compute_dT_ET_CGS
  end type Universe_R_CC_Form

    ! private :: &
    !   ResolveCycle_R, &
    !   PrepareStep_F, &
    !   Compute_dT_Local, &
    !   SetSlope_F_P_DFV_SS, &
    !   SetSlope_F_P_SS, &
    !   SetSlope_RM_I, &
    !   SetSlope_RM_DFV_I, &
    !   SetSlope_RM_DFV_I_I

    !   private :: &
    !     Compute_dT_ET_CGS_Kernel

    ! interface
    
    !   module subroutine Compute_dT_ET_CGS_Kernel &
    !            ( dT, ProperCell, Q, E, UseDeviceOption )
    !     use Basics
    !     implicit none
    !     real ( KDR ), intent ( inout ) :: &
    !       dT
    !     logical ( KDL ), dimension ( : ), intent ( in ) :: &
    !       ProperCell
    !     real ( KDR ), dimension ( : ), intent ( in ) :: &
    !       Q, E
    !     logical ( KDL ), intent ( in ), optional :: &
    !       UseDeviceOption
    !   end subroutine Compute_dT_ET_CGS_Kernel

    ! end interface


contains


  subroutine Initialize_R_CC &
               ( U, RadiationName, RadiationType, FormalismType, Name )!, &
!                MinCoordinateOption, MaxCoordinateOption, &
!                 FinishTimeOption, nCellsPositionOption, nWriteOption )

    class ( Universe_R_CC_Form ), intent ( inout ) :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name
    ! logical ( KDL ), intent ( in ), optional :: &
    !   ApplyStreamingOption, &
    !   ApplyInteractionsOption, &
    !   EvolveFluidOption
    ! real ( KDR ), dimension ( : ), intent ( in ), optional :: &
    !   MinCoordinateOption, &
    !   MaxCoordinateOption
    ! real ( KDR ), intent ( in ), optional :: &
    !   FinishTimeOption
    ! integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
    !   nCellsPositionOption
    ! integer ( KDI ), intent ( in ), optional :: &
    !   nWriteOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_CC'

    call U % Universe_H_Form % Initialize ( Name )

    !-- Radiations

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    !-- Units

    if ( .not. allocated ( U % Units_F ) ) then
      allocate ( U % Units_F ( 1 ) )
      call U % Units_F ( 1 ) % Initialize ( TypeOption = 'ASTROPHYSICS' )
    end if

    if ( .not. allocated ( U % Units_R ) ) then
      allocate ( U % Units_R ( 1 ) )
      call U % Units_R ( 1 ) % Initialize ( TypeOption = 'ASTROPHYSICS' )
    end if

call Show ( '>>> 1' )
    ! !-- Initializations

    ! call U % SetCommunicator &
    !        ( )
    ! call U % AllocateIntegrator &
    !        ( )
    ! call U % InitializePositionSpace &
    !        ( MinCoordinateOption = MinCoordinateOption, &
    !          MaxCoordinateOption = MaxCoordinateOption, &
    !          nCellsOption = nCellsPositionOption )
    ! ! call RB % InitializeMomentumSpace &
    ! !        ( EnergySpacingOption = EnergySpacingOption, &
    ! !          MinEnergyOption = MinEnergyOption, &
    ! !          MaxEnergyOption = MaxEnergyOption, &
    ! !          MinWidthEnergyOption = MinWidthEnergyOption, &
    ! !          EnergyScaleOption = EnergyScaleOption, &
    ! !          nCellsEnergyOption = nCellsEnergyOption )
    ! call U % InitializeGravitation &
    !        ( GravitationType = 'GALILEO' )
    ! call U % InitializeFluid &
    !        ( FluidType = 'IDEAL' )
    ! call U % InitializeInteractions &
    !        ( )
    ! call U % InitializeRadiation &
    !        ( )
    ! call U % InitializeSteps &
    !        ( )
    ! call U % InitializeIntegrator &
    !        ( FinishTimeOption = FinishTimeOption, &
    !          nWriteOption = nWriteOption )

    ! !-- Integrator methods

    ! associate ( I  =>  U % Integrator )
    ! I % ResolveCycle      =>  ResolveCycle_R
    ! I % PrepareStep       =>  PrepareStep_F
    ! I % Compute_dT_Local  =>  Compute_dT_Local
    ! end associate !-- I

  end subroutine Initialize_R_CC


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_CC_Form ), intent ( inout ) :: &
      U

    ! if ( allocated ( U % Interactions_BM ) ) &
    !   deallocate ( U % Interactions_BM )
    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
    ! if ( allocated ( U % CO_SplitSource ) ) &
    !   deallocate ( U % CO_SplitSource )
    ! if ( allocated ( U % Communicator_PS ) ) &
    !   deallocate ( U % Communicator_PS )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( Universe_R_CC_Form ), intent ( in ) :: &
      U

    call U % Universe_F_CC_Form % ShowParameters ( )

    call Show ( U % RadiationName,     'RadiationName',     U % IGNORABILITY )
    call Show ( U % RadiationType,     'RadiationType',     U % IGNORABILITY )
    call Show ( U % iRadiation,        'iRadiation',        U % IGNORABILITY )
    call Show ( U % FormalismType,     'FormalismType',     U % IGNORABILITY )
    call Show ( U % InteractionFactor, 'InteractionFactor', U % IGNORABILITY )

  end subroutine ShowParameters


end module Universe_R_CC__Form

