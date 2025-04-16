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
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
    class ( Interactions_BM_Form ), dimension ( : ), allocatable :: &
      Interactions_BM
  contains
    procedure, private, pass :: &
      Initialize_R_B_C
    generic, public :: &
      Initialize => Initialize_R_B_C
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      InitializeInteractions
    procedure, public, pass :: &
      InitializeStep
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
    call U % InitializeRadiation &
           ( )
    call U % InitializeInteractions &
           ( )
    call U % InitializeStep &
           ( )
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

    if ( allocated ( U % Interactions_BM ) ) &
      deallocate ( U % Interactions_BM )
    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
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


  subroutine InitializeRadiation ( U )

    class ( Universe_R_B_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

      select type ( A  =>  I % X )
      class is ( Atlas_SCG_Form )
      associate ( C  =>  A % Chart_GS )
        allocate ( U % Units_R ( 1 ) )
        call U % Units_R ( 1 ) % Initialize &
               ( C % CoordinateUnit, TypeOption = U % UnitsType )
      end associate !-- C
      end select !-- A

      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_Form )

      associate ( nR  =>  U % nRadiations )

      select case ( trim ( U % RadiationType ( 1 ) ) )
      case ( 'GENERIC' )

        allocate ( RadiationMoments_BM_Form :: I % CurrentSet_X_1D ( nR ) )

        do iR  =  1, nR

          select type ( R  =>  I % CurrentSet_X_1D ( iR ) )
          class is ( RadiationMoments_BM_Form )

            call R % Initialize &
                   ( F, U % Units_R, U % RadiationType ( iR ), &
                     NameOption = U % RadiationName ( iR ) )

          end select !-- R

        end do !-- iR

      case ( 'PHOTONS' )

        allocate ( PhotonMoments_G_Form :: I % CurrentSet_X_1D ( nR ) )

        do iR  =  1, nR

          select type ( R  =>  I % CurrentSet_X_1D ( iR ) )
          class is ( PhotonMoments_G_Form )

            call R % Initialize &
                   ( F, U % Units_R, U % RadiationType ( iR ), &
                     NameOption = U % RadiationName ( iR ) )

          end select !-- R

        end do !-- iR

      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( U % RadiationType ( iR ), 'RadiationType', &
                    CONSOLE % ERROR )
        call Show ( 'Universe_R_B_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- RadiationType

      end associate !-- nR
      end select !-- F

    end associate !-- G
    end select !-- I

  end subroutine InitializeRadiation


  subroutine InitializeInteractions ( U )

    class ( Universe_R_B_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( RadiationMoments_BM_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )
    associate &
      ( nR  =>  U % nRadiations )

   if ( allocated ( U % Interactions_BM ) ) then
      do iR  =  1,  nR
        call U % Interactions_BM ( iR ) % Initialize &
               ( R ( iR ), U % Units_R, F )
      end do !-- iR
    end if

    end associate !-- nR
    end select !-- F
    end select !-- R
    end select !-- I

  end subroutine InitializeInteractions


  subroutine InitializeStep ( U )

    class ( Universe_R_B_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR
    character ( LDL ) :: &
      RiemannSolverType

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_C_CS_Form )

!      allocate ( Step_RK_NM_G_1D_C_Form :: I % Step_X )
!      select type ( S  =>  I % Step_X )
!        class is ( Step_RK_NM_G_1D_C_Form )
      allocate ( Step_RK_CS_1D_C_CS_Form :: I % Step_X )
      select type ( S  =>  I % Step_X )
        class is ( Step_RK_CS_1D_C_CS_Form )

      associate ( nR  =>  U % nRadiations )

      allocate ( S % Step_CS_1D ( nR ) )
      allocate ( S % Step_CS )
      associate &
        ( S_R  =>  S % Step_CS_1D ( : ), &
          S_F  =>  S % Step_CS )

      !-- Radiation
      select type ( R  =>  I % CurrentSet_X_1D )
        class is ( RadiationMoments_BM_Form )

      do iR  =  1, nR
 
        allocate ( DivergencePart_RM_Form :: S_R ( iR ) % DivergenceTotal )
        associate ( DT  =>  S_R ( iR ) % DivergenceTotal )
        call DT % Initialize ( R ( iR ) )
        end associate !-- DT

        allocate ( DiffusionFactor_RM_Form :: S_R ( iR ) % DiffusionFactor )
        select type ( DF  =>  S_R ( iR ) % DiffusionFactor )
        class is ( DiffusionFactor_RM_Form )
          call DF % Initialize ( U % Interactions_BM ( iR ) )
        end select !-- DF

      end do !-- iR

      !-- Fluid
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_HN_Form )

      allocate ( DivergencePart_F_P_T_Form :: S_F % DivergenceTotal )
      associate ( DT  =>  S_F % DivergenceTotal )
        call DT % Initialize ( F )
      end associate !-- DT

      RiemannSolverType = 'HLLC'
      call PROGRAM_HEADER % GetParameter &
             ( RiemannSolverType, 'RiemannSolverType' )
      if ( trim ( RiemannSolverType ) == 'HLLC' ) then
        allocate ( RiemannSolver_HLLC_P_Form :: S_F % RiemannSolver )
        associate ( RS  =>  S_F % RiemannSolver )
        call RS % Initialize ( F )
        end associate !-- RS
      end if

      !-- Combined step
      call S % Initialize &
             ( R, F, ImplicitExplicitOption = .true., &
               ComputeImplicit_CS_1D_Option = U % ApplyInteractions, &
               ComputeImplicit_CS_Option    = U % EvolveFluid, &
               ComputeExplicit_CS_1D_Option = U % ApplyStreaming, &
               ComputeExplicit_CS_Option    = U % EvolveFluid, &
               nStagesOption = 3 )

      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_R_B_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select    !-- F
      end select    !-- R
      end associate !-- S_R, S_F
      end associate !-- nR
      end select    !-- S

    end select !-- I

  end subroutine InitializeStep


end module Universe_R_B_C__Form
