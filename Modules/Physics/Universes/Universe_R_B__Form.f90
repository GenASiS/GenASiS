module Universe_R_B__Form

  !-- Universe_Radiation_Box__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Universe_F_B__Form

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: Universe_R_B_Form
    integer ( KDI ) :: &
      iRadiation  = 0, &
      nRadiations = 0
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
    type ( CommunicatorForm ), allocatable :: &
      Communicator_PS  !-- PositionSpace
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
  contains
    procedure, private, pass :: &
      Initialize_R_B
    generic, public :: &
      Initialize => Initialize_R_B
    final :: &
      Finalize
    procedure, private, pass :: &
      SetCommunicator
    procedure, private, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      ShowParameters
  end type Universe_R_B_Form


contains


  subroutine Initialize_R_B &
               ( U, RadiationName, RadiationType, FormalismType, Name, &
                 MinCoordinateOption, MaxCoordinateOption, &
                 nCellsPositionOption )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsPositionOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_B'

    call U % Universe_H_Form % Initialize ( Name )

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    allocate ( U % Units_F ( 1 ) )
    call U % Units_F ( 1 ) % Initialize ( )

    allocate ( U % Units_R ( 1 ) )
    call U % Units_R ( 1 ) % Initialize ( )

    call U % SetCommunicator &
           ( )
    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsPositionOption )
    ! call RB % InitializeMomentumSpace &
    !        ( EnergySpacingOption = EnergySpacingOption, &
    !          MinEnergyOption = MinEnergyOption, &
    !          MaxEnergyOption = MaxEnergyOption, &
    !          MinWidthEnergyOption = MinWidthEnergyOption, &
    !          EnergyScaleOption = EnergyScaleOption, &
    !          nCellsEnergyOption = nCellsEnergyOption )
    call U % InitializeGravitation &
           ( GravitationType = 'GALILEO' )
    call U % InitializeFluid &
           ( FluidType = 'IDEAL' )
    call U % InitializeRadiation &
           ( )

  end subroutine Initialize_R_B


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
    if ( allocated ( U % Communicator_PS ) ) &
      deallocate ( U % Communicator_PS )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine SetCommunicator ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iPS, &  !-- iPositionSpace
      nProcesses, &
      nProcesses_PS  !-- PositionSpace
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank

    nProcesses  =  U % Communicator % Size
    if ( mod ( nProcesses, U % nRadiations )  /=  0 ) then
      call Show ( 'nRadiations must evenly divide nProcesses', CONSOLE % ERROR )
      call Show ( nProcesses, 'nProcesses', CONSOLE % ERROR )
      call Show ( U % nRadiations, 'nRadiations', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )

      nProcesses_PS  =  nProcesses / U % nRadiations

      do iPS  =  1,  U % nRadiations
        allocate ( Rank, &
                   source =  [ ( iP, iP = ( iPS - 1 ) * nProcesses_PS, &
                                            iPS * nProcesses_PS  -  1 ) ] )
        if ( any ( U % Communicator % Rank  ==  Rank ) ) then
          U % iRadiation = iPS
          allocate ( U % Communicator_PS )
          call U % Communicator_PS % Initialize &
                 ( U % Communicator, Rank, NameOption = 'Communicator_PS' )
        end if
        deallocate ( Rank )
      end do !-- iPS

    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_B_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCommunicator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine SetCommunicator


  subroutine AllocateIntegrator ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iCS  !-- iCurrentSet

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )
      allocate ( Integrator_CS_1D_BM_CS_Form :: U % Integrator )
    case ( 'SPECTRAL' )
      allocate ( Integrator_CS_1D_CB_CS_Form :: U % Integrator )
    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_B_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator_R_B', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

    select type ( I => U % Integrator )
    class is ( Integrator_CS_1D_CS_Form )

      I % iCurrentSet  =  U % iRadiation

      I % N_CURRENT_SETS_1D  =  size ( U % RadiationName )
      allocate ( I % dT_Label &
                   ( 1  +  I % N_CURRENT_SETS_1D  +  I % N_CURRENT_SETS_1D ) )

      I % dT_Label ( 1 )  &
        =  'FluidAdvection'

      do iCS = 1, I % N_CURRENT_SETS_1D
        I % dT_Label ( 1 + iCS )  &
          =  trim ( U % RadiationName ( iCS ) ) // 'Streaming'
      end do !-- iCS

      do iCS = 1, I % N_CURRENT_SETS_1D
        I % dT_Label ( I % N_CURRENT_SETS_1D  +  1  +  iCS )  &
          =  trim ( U % RadiationName ( iCS ) ) // 'Interactions'
      end do !-- iCS

    end select !-- I

    ! select type ( I => U % Integrator )
    ! class is ( Integrator_C_1D_PS_C_PS_Form )
    !   allocate ( I % Current_ASC_1D ( I % N_CURRENT_SETS_1D ) )
    ! class is ( Integrator_C_1D_MS_C_PS_Form )
    !   allocate ( I % Current_BSLL_ASC_CSLD_1D ( I % N_CURRENT_SETS_1D ) )
    ! end select !-- I

  end subroutine AllocateIntegrator


  subroutine InitializePositionSpace &
               ( U, CommunicatorOption, MinCoordinateOption, &
                 MaxCoordinateOption, nCellsOption )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption

    integer ( KDI ), dimension ( 3 ) :: &
      nCellsPosition

    nCellsPosition = [ 128, 128, 128 ]
    if ( present ( nCellsOption ) ) &
      nCellsPosition = nCellsOption
    call PROGRAM_HEADER % GetParameter ( nCellsPosition, 'nCellsPosition' )

    call U % Universe_F_B_Form % InitializePositionSpace &
           ( CommunicatorOption = U % Communicator_PS, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsPosition )

  end subroutine InitializePositionSpace


  subroutine InitializeRadiation ( U )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_BM_CS_Form )

      associate &
        ( G   =>  I % Geometry_X, &
          iR  =>  U % iRadiation )

      select case ( trim ( U % RadiationType ( iR ) ) )
      case ( 'GENERIC' )

        allocate ( RadiationMoments_BM_Form :: I % CurrentSet_X_1D )
        select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )

        call R % Initialize &
               ( G, U % Units_R, NameOption = U % RadiationName ( iR ) )

        end select !-- R

      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( U % RadiationType ( iR ), 'RadiationType', &
                    CONSOLE % ERROR )
        call Show ( 'Universe_R_B__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- FluidType

      end associate !-- G, etc.

    end select !-- I

  end subroutine InitializeRadiation


  subroutine ShowParameters ( U )

    class ( Universe_R_B_Form ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % RadiationName, 'RadiationName', U % IGNORABILITY )
    call Show ( U % RadiationType, 'RadiationType', U % IGNORABILITY )
    call Show ( U % iRadiation,    'iRadiation   ', U % IGNORABILITY )
    call Show ( U % FormalismType, 'FormalismType', U % IGNORABILITY )

  end subroutine ShowParameters


end module Universe_R_B__Form
