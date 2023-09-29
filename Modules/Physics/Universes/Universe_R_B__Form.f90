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
      nRadiations = 0
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName
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
      AllocateIntegrator_R_B
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_R_B
    procedure, public, pass :: &
      ShowParameters
  end type Universe_R_B_Form


contains


  subroutine Initialize_R_B &
               ( U, RadiationName, RadiationType, FormalismType, Name )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      FormalismType, &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_B'

    call U % Universe_H_Form % Initialize ( Name )

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    U % RadiationName  =  RadiationName

    U % FormalismType  =  FormalismType

    allocate ( U % Units_R ( 1 ) )
    call U % Units_R ( 1 ) % Initialize ( )

    call U % AllocateIntegrator &
           ( RadiationName )

  end subroutine Initialize_R_B


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_B_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine AllocateIntegrator_R_B ( U, RadiationName )

    class ( Universe_R_B_Form ), intent ( inout ) :: &
      U
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName

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

      I % N_CURRENT_SETS_1D  =  size ( RadiationName )
      allocate ( I % dT_Label &
                   ( 1  +  I % N_CURRENT_SETS_1D  +  I % N_CURRENT_SETS_1D ) )

      I % dT_Label ( 1 )  &
        =  'FluidAdvection'

      do iCS = 1, I % N_CURRENT_SETS_1D
        I % dT_Label ( 1 + iCS )  &
          =  trim ( RadiationName ( iCS ) ) // 'Streaming'
      end do !-- iCS

      do iCS = 1, I % N_CURRENT_SETS_1D
        I % dT_Label ( I % N_CURRENT_SETS_1D  +  1  +  iCS )  &
          =  trim ( RadiationName ( iCS ) ) // 'Interactions'
      end do !-- iCS

    end select !-- I

    ! select type ( I => U % Integrator )
    ! class is ( Integrator_C_1D_PS_C_PS_Form )
    !   allocate ( I % Current_ASC_1D ( I % N_CURRENT_SETS_1D ) )
    ! class is ( Integrator_C_1D_MS_C_PS_Form )
    !   allocate ( I % Current_BSLL_ASC_CSLD_1D ( I % N_CURRENT_SETS_1D ) )
    ! end select !-- I

  end subroutine AllocateIntegrator_R_B


  subroutine ShowParameters ( U )

    class ( Universe_R_B_Form ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % RadiationName, 'RadiationName', U % IGNORABILITY )
    call Show ( U % FormalismType, 'FormalismType', U % IGNORABILITY )

  end subroutine ShowParameters


end module Universe_R_B__Form
