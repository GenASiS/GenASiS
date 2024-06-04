module Universe_R_CC_C__Form

  !-- Universe_Radiation_CentralCore_Collected__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Radiations
  use Universe_F_CC__Form

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: Universe_R_CC_C_Form
!     integer ( KDI ) :: &
!       iRadiation  = 0, &
!       nRadiations = 0
!     real ( KDR ) :: &
!       InteractionFactor = 0.0_KDR
!     character ( LDL ) :: &
!       FormalismType = ''
!     character ( LDL ), dimension ( : ), allocatable :: &
!       RadiationName, &
!       RadiationType
!     type ( CommunicatorForm ), allocatable :: &
!       Communicator_PS  !-- PositionSpace
! !    type ( CollectiveOperation_R_Form ), dimension ( : ), allocatable :: &
! !      CO_SplitSource
!     type ( Units_R_Form ), dimension ( : ), allocatable :: &
!       Units_R
!     type ( Coarsening_C_RM_Form ), allocatable :: &
!       Coarsening_R
!     class ( Interactions_NM_G_Form ), allocatable :: &
!       Interactions_NM_G
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    final :: &
      Finalize
!     procedure, private, pass :: &
!       SetCommunicator
!     procedure, private, pass :: &
!       AllocateIntegrator
!     procedure, public, pass :: &
!       InitializeRadiation
!     procedure, public, pass :: &
!       InitializeInteractions
!     procedure, public, pass :: &
!       SetBoundaryConditions
! !    procedure, public, pass :: &
! !      InitializeSteps
!     procedure, public, pass :: &
!       InitializeStep
!     procedure, public, pass :: &
!       InitializeIntegrator
!     procedure, public, pass :: &
!       ShowParameters
!     procedure, public, pass :: &
!       ShowDiagnostics
!     procedure, public, pass ( U ) :: &
!       Compute_dT_RI_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_RT_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_RK_F_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_RK_R_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_IS_F_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_IS_R_CGS
  end type Universe_R_CC_C_Form


contains


  subroutine Initialize_R_CC_C &
               ( U, RadiationName, RadiationType, FormalismType, FluidType, &
                 GravitationType, Name, UnitsTypeOption, FinishTimeOption, &
                 RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_R_CC_C_Form ), intent ( inout ), target :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      FluidType, &
      GravitationType, &
      Name
    character ( * ), intent ( in ), optional :: &
      UnitsTypeOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

!     if ( U % Type  ==  '' ) &
!       U % Type  =  'a Universe_R_CC'

!     call U % Universe_H_Form % Initialize &
!            ( Name, UnitsTypeOption = UnitsTypeOption )

!     !-- FIXME: This is for a workaround in SetSlope routines below
!     UNIVERSE  =>  U

!     !-- Radiations

!     U % nRadiations  =  size ( RadiationName )

!     allocate ( U % RadiationName ( U % nRadiations ) )
!     allocate ( U % RadiationType ( U % nRadiations ) )
!     U % RadiationName  =  RadiationName
!     U % RadiationType  =  RadiationType

!     U % FormalismType  =  FormalismType

!     if ( trim ( RadiationType ( 1 ) )  ==  'NONE' ) then
!       call U % Universe_F_CC_Form % Initialize &
!              ( FluidType, GravitationType, Name, &
!                FinishTimeOption = FinishTimeOption, &
!                RadiusMaxOption = RadiusMaxOption, &
!                RadiusCoreOption = RadiusCoreOption, &
!                RadialRatioOption = RadialRatioOption, &
!                nCellsPolarOption = nCellsPolarOption, &
!                nWriteOption = nWriteOption )
!       return
!     end if

!     !-- Initializations

!     call U % SetCommunicator &
!            ( )
!     call U % AllocateIntegrator &
!            ( )
!     call U % InitializePositionSpace &
!            ( CommunicatorOption = U % Communicator_PS, &
!              RadiusMaxOption = RadiusMaxOption, &
!              RadiusCoreOption = RadiusCoreOption, &
!              RadialRatioOption = RadialRatioOption, &
!              nCellsPolarOption = nCellsPolarOption )
!     ! call RB % InitializeMomentumSpace &
!     !        ( EnergySpacingOption = EnergySpacingOption, &
!     !          MinEnergyOption = MinEnergyOption, &
!     !          MaxEnergyOption = MaxEnergyOption, &
!     !          MinWidthEnergyOption = MinWidthEnergyOption, &
!     !          EnergyScaleOption = EnergyScaleOption, &
!     !          nCellsEnergyOption = nCellsEnergyOption )
!     call U % InitializeGravitation &
!            ( GravitationType )
!     call U % InitializeFluid &
!            ( FluidType )
!     call U % InitializeRadiation &
!            ( )
!     call U % InitializeInteractions &
!            ( )
!     call U % SetBoundaryConditions &
!            ( )
! !    call U % InitializeSteps &
! !           ( )
!     call U % InitializeStep &
!            ( )
!     call U % InitializeIntegrator &
!            ( GravitationType, &
!              FinishTimeOption = FinishTimeOption, &
!              nWriteOption = nWriteOption )
!     call U % InitializeDiagnostics &
!            ( )

!     call U % SetMeasures ( )

!     !-- Integrator methods

!     associate ( I  =>  U % Integrator )
! !    I % ResolveCycle              =>  ResolveCycle_R
! !    I % PrepareStep               =>  PrepareStep_F
!     I % Compute_dT_Local          =>  Compute_dT_Local
!     I % InitializeSeries          =>  InitializeSeries
!     I % Analyze                   =>  Analyze
!     I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
!     end associate !-- I

  end subroutine Initialize_R_CC_C


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

!     if ( allocated ( U % Interactions_NM_G ) ) &
!       deallocate ( U % Interactions_NM_G )
!     if ( allocated ( U % Coarsening_R ) ) &
!       deallocate ( U % Coarsening_R )
!     if ( allocated ( U % Units_R ) ) &
!       deallocate ( U % Units_R )
! !    if ( allocated ( U % CO_SplitSource ) ) &
! !      deallocate ( U % CO_SplitSource )
!     if ( allocated ( U % Communicator_PS ) ) &
!       deallocate ( U % Communicator_PS )
!     if ( allocated ( U % RadiationType ) ) &
!       deallocate ( U % RadiationType )
!     if ( allocated ( U % RadiationName ) ) &
!       deallocate ( U % RadiationName )

  end subroutine Finalize


end module Universe_R_CC_C__Form
