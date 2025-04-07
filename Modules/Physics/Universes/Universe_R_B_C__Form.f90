module Universe_R_B_C__Form

  !-- Universe_Radiation_Box_Collected__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Universe_F_B__Form

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: Universe_R_B_C_Form
    integer ( KDI ) :: &
      nRadiations = 0
    logical ( KDL ) :: &
      ApplyStreaming, &
      ApplyInteractions, &
      EvolveFluid
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
! !    type ( CollectiveOperation_R_Form ), dimension ( : ), allocatable :: &
! !      CO_SplitSource
!     type ( Units_R_Form ), dimension ( : ), allocatable :: &
!       Units_R
!     class ( Interactions_BM_Form ), allocatable :: &
!       Interactions_BM
  contains
    procedure, private, pass :: &
      Initialize_R_B_C
    generic, public :: &
      Initialize => Initialize_R_B_C
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator
!     procedure, public, pass :: &
!       InitializeRadiation
!     procedure, public, pass :: &
!       InitializeInteractions
! !    procedure, public, pass :: &
! !      InitializeSteps
!     procedure, public, pass :: &
!       InitializeStep
!     procedure, public, pass :: &
!       InitializeIntegrator
!     procedure, public, pass :: &
!       ShowParameters
!     procedure, public, pass ( U ) :: &
!       Compute_dT_R_E_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_ET_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_RK_F_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_RK_R_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_IS_F_CGS
!     procedure, public, pass ( U ) :: &
!       Compute_dT_IS_R_CGS
  end type Universe_R_B_C_Form


contains


  subroutine Initialize_R_B_C &
               ( U, RadiationName, RadiationType, FormalismType, Name, &
                 UnitsTypeOption, ApplyStreamingOption, &
                 ApplyInteractionsOption, EvolveFluidOption, &
                 MinCoordinateOption, MaxCoordinateOption, &
                 FinishTimeOption, nCellsPositionOption, nWriteOption )

    class ( Universe_R_B_C_Form ), intent ( inout ), target :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name
    character ( * ), intent ( in ), optional :: &
      UnitsTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      ApplyStreamingOption, &
      ApplyInteractionsOption, &
      EvolveFluidOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsPositionOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_B_C'

    call U % Universe_H_Form % Initialize &
           ( Name, UnitsTypeOption = UnitsTypeOption )

    !-- Radiations

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    !-- Operators

    U % ApplyStreaming    = .true.
    U % ApplyInteractions = .true.
    U % EvolveFluid       = .true.
    if ( present ( ApplyStreamingOption ) ) &
      U % ApplyStreaming = ApplyStreamingOption
    if ( present ( ApplyInteractionsOption ) ) &
      U % ApplyInteractions = ApplyInteractionsOption
    if ( present ( EvolveFluidOption ) ) &
      U % EvolveFluid = EvolveFluidOption

    !-- Initializations

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsPositionOption )
    call U % InitializeGravitation &
           ( GravitationType = 'GALILEO' )
    call U % InitializeFluid &
           ( FluidType = 'IDEAL' )
!     call U % InitializeRadiation &
!            ( )
!     call U % InitializeInteractions &
!            ( )
!  !   call U % InitializeSteps &
!  !          ( )
!     call U % InitializeStep &
!            ( )
!     call U % InitializeIntegrator &
!            ( FinishTimeOption = FinishTimeOption, &
!              nWriteOption = nWriteOption )

!     !-- Integrator methods

!     associate ( I  =>  U % Integrator )
! !    I % ResolveCycle      =>  ResolveCycle_R
! !    I % PrepareStep       =>  PrepareStep_F
!     I % Compute_dT_Local  =>  Compute_dT_Local
!     end associate !-- I

  end subroutine Initialize_R_B_C

  
  impure elemental subroutine Finalize ( U )

    type ( Universe_R_B_C_Form ), intent ( inout ) :: &
      U

!     if ( allocated ( U % Interactions_BM ) ) &
!       deallocate ( U % Interactions_BM )
!     if ( allocated ( U % Units_R ) ) &
!       deallocate ( U % Units_R )
!    if ( allocated ( U % CO_SplitSource ) ) &
!      deallocate ( U % CO_SplitSource )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine AllocateIntegrator ( U )

    class ( Universe_R_B_C_Form ), intent ( inout ) :: &
      U

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )
      allocate ( Integrator_CS_1D_C_CS_Form :: U % Integrator )
    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_B_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine AllocateIntegrator


end module Universe_R_B_C__Form
