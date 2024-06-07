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
    integer ( KDI ) :: &
      nRadiations = 0
    character ( LDL ) :: &
      FormalismType = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      RadiationName, &
      RadiationType
    type ( Units_R_Form ), dimension ( : ), allocatable :: &
      Units_R
    type ( Coarsening_C_RM_Form ), dimension ( : ), allocatable :: &
      Coarsening_R
    class ( Interactions_NM_G_Form ), dimension ( : ), allocatable :: &
      Interactions_NM_G
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      InitializeInteractions
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, public, pass :: &
      InitializeStep
    procedure, public, pass :: &
      InitializeIntegrator
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass :: &
      ShowDiagnostics
  end type Universe_R_CC_C_Form

    private :: &
      Compute_dT_Local, &
      InitializeSeries, &
      Analyze, &
      Write, &
      Set_T_CheckpointInterval, &
      SetSlope_F_P_DFV_N


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

    if ( U % Type  ==  '' ) &
      U % Type  =  'a Universe_R_CC_C'

    call U % Universe_H_Form % Initialize &
           ( Name, UnitsTypeOption = UnitsTypeOption )

    !-- Radiations

    U % nRadiations  =  size ( RadiationName )

    allocate ( U % RadiationName ( U % nRadiations ) )
    allocate ( U % RadiationType ( U % nRadiations ) )
    U % RadiationName  =  RadiationName
    U % RadiationType  =  RadiationType

    U % FormalismType  =  FormalismType

    !-- Initializations

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )
    call U % InitializeGravitation &
           ( GravitationType )
    call U % InitializeFluid &
           ( FluidType )
    call U % InitializeRadiation &
           ( )
    call U % InitializeInteractions &
           ( )
    call U % SetBoundaryConditions &
           ( )
    call U % InitializeStep &
           ( )
    call U % InitializeIntegrator &
           ( GravitationType, &
             FinishTimeOption = FinishTimeOption, &
             nWriteOption = nWriteOption )
    call U % InitializeDiagnostics &
           ( )

    call U % SetMeasures ( )

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % Compute_dT_Local          =>  Compute_dT_Local
    I % InitializeSeries          =>  InitializeSeries
    I % Analyze                   =>  Analyze
    I % Write                     =>  Write
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    end associate !-- I

  end subroutine Initialize_R_CC_C


  impure elemental subroutine Finalize ( U )

    type ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Interactions_NM_G ) ) &
      deallocate ( U % Interactions_NM_G )
    if ( allocated ( U % Coarsening_R ) ) &
      deallocate ( U % Coarsening_R )
    if ( allocated ( U % Units_R ) ) &
      deallocate ( U % Units_R )
    if ( allocated ( U % RadiationType ) ) &
      deallocate ( U % RadiationType )
    if ( allocated ( U % RadiationName ) ) &
      deallocate ( U % RadiationName )

  end subroutine Finalize


  subroutine AllocateIntegrator ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    select case ( trim ( U % FormalismType ) )
    case ( 'GREY' )
      allocate ( Integrator_CS_1D_C_CS_Form :: U % Integrator )
    case default
      call Show ( 'FormalismType not recognized', CONSOLE % ERROR )
      call Show ( U % FormalismType, 'FormalismType', CONSOLE % ERROR )
      call Show ( 'Universe_R_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FormalismType

  end subroutine AllocateIntegrator


  subroutine InitializeRadiation ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_C_CS_Form )

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

      allocate ( NeutrinoMoments_G_Form :: I % CurrentSet_X_1D ( nR ) )

      do iR  =  1, nR

        select case ( trim ( U % RadiationType ( iR ) ) )
        case ( 'NEUTRINOS_E', 'NEUTRINOS_E_BAR' )

        select type ( R  =>  I % CurrentSet_X_1D ( iR ) )
        class is ( NeutrinoMoments_G_Form )

        call R % Initialize &
               ( F, U % Units_R, U % RadiationType ( iR ), &
                 NameOption = U % RadiationName ( iR ) )
        if ( allocated ( U % Interactions_NM_G ) ) &
          call R % SetInteractions ( U % Interactions_NM_G ( iR ) )

        end select !-- R

        case default
          call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
          call Show ( U % RadiationType ( iR ), 'RadiationType', &
                      CONSOLE % ERROR )
          call Show ( 'Universe_R_CC_C__Form', 'module', CONSOLE % ERROR )
          call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- RadiationType

      end do !-- iR

      end associate !-- nR
      end select !-- F

    end select !-- I

  end subroutine InitializeRadiation


  subroutine InitializeInteractions ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( R  =>  I % CurrentSet_X_1D )
      class is ( RadiationMoments_BM_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_Form )

    if ( allocated ( U % Interactions_NM_G ) ) then
      do iR  =  1, U % nRadiations
        call U % Interactions_NM_G ( iR ) % Initialize &
               ( R ( iR ) , U % Units_R, F )
      end do !-- iR
    end if

    end select !-- F
    end select !-- R
    end select !-- I

  end subroutine InitializeInteractions


  subroutine SetBoundaryConditions ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U
    
    integer ( KDI ) :: &
      iR

    call U % Universe_F_CC_Form % SetBoundaryConditions ( )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_1D_C_CS_Form )

    associate &
      ( R  =>  I % CurrentSet_X_1D )

    do iR  =  1, U % nRadiations
      call R ( iR ) % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
      call R ( iR ) % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
      call R ( iR ) % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end do !-- iR

    end associate !-- R

    end select !-- I

  end subroutine SetBoundaryConditions


  subroutine InitializeStep ( U )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      iR
    character ( LDL ) :: &
      RiemannSolverType

    U % Coarsen  =  .true.
    call PROGRAM_HEADER % GetParameter ( U % Coarsen, 'Coarsen' )

    select type ( I  =>  U % Integrator )
    class is ( Integrator_CS_1D_C_CS_Form )

      allocate ( Step_RK_NM_G_1D_C_Form :: I % Step_X )
      select type ( S  =>  I % Step_X )
        class is ( Step_RK_NM_G_1D_C_Form )

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
 
        allocate ( DivergencePart_NM_G_Form :: S_R ( iR ) % DivergenceTotal )
        associate ( DT  =>  S_R ( iR ) % DivergenceTotal )
        call DT % Initialize ( R ( iR ) )
        end associate !-- DT

        allocate ( DiffusionFactor_RM_Form :: S_R ( iR ) % DiffusionFactor )
        select type ( DF  =>  S_R ( iR ) % DiffusionFactor )
        class is ( DiffusionFactor_RM_Form )
          call DF % Initialize ( U % Interactions_NM_G ( iR ) )
        end select !-- DF

      end do !-- iR

      !-- Fluid
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_HN_Form )

      allocate ( DivergencePart_F_P_HN_T_Form :: S_F % DivergenceTotal )
      associate ( DT  =>  S_F % DivergenceTotal )
        call DT % Initialize ( F )
      end associate !-- DT

      RiemannSolverType = 'HLL'
      call PROGRAM_HEADER % GetParameter &
             ( RiemannSolverType, 'RiemannSolverType' )
      if ( trim ( RiemannSolverType ) == 'HLLC' ) then
        allocate ( RiemannSolver_HLLC_P_HN_Form :: S_F % RiemannSolver )
        associate ( RS  =>  S_F % RiemannSolver )
        call RS % Initialize ( F )
        end associate !-- RS
      end if        

      S_F % SetSlopeExplicit  =>  SetSlope_F_P_DFV_N

      !-- Combined step
      call S % Initialize &
             ( R, F, ImplicitExplicitOption = .true., nStagesOption = 3 )

      !-- Coarsening
      if ( U % Coarsen ) then

        allocate ( U % Coarsening_F )
        associate &
          ( C_F  =>  U % Coarsening_F, &
            G    =>  I % Geometry_X )
          call C_F % Initialize ( F, G )
          call S_F % SetCoarsening ( C_F )
        end associate !-- C_F, etc.

        allocate ( U % Coarsening_R ( nR ) )
        associate &
          ( C_R  =>  U % Coarsening_R, &
            G    =>  I % Geometry_X )
        do iR  =  1, nR
          call C_R ( iR ) % Initialize ( R ( iR ), G )
          call S_R ( iR ) % SetCoarsening ( C_R ( iR ) )
        end do !-- iR
        end associate !-- C_R, etc.
      end if

      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_R_CC_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select    !-- F
      end select    !-- R
      end associate !-- S_R, S_F
      end associate !-- nR
      end select    !-- S

    end select !-- I

  end subroutine InitializeStep


  subroutine InitializeIntegrator &
               ( U, GravitationType, FinishTimeOption, GravityFactorOption, &
                 nWriteOption )

    class ( Universe_R_CC_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      GravityFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iR

    select type ( I => U % Integrator )
    class is ( Integrator_CS_1D_C_CS_Form )

    allocate ( I % dT_Label ( 3 ) )

    I % dT_Label ( 1 )  =  'GravitationAcceleration'
    I % dT_Label ( 2 )  =  'FluidAdvection'
    I % dT_Label ( 3 )  =  'RadiationStreaming'

    U % GravityFactor  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter &
           ( U % GravityFactor, 'GravityFactor' )

    call I % Initialize &
           ( U % nRadiations, &
             Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    do iR  =  1, U % nRadiations
      call U % Interactions_NM_G ( iR ) % SetStream ( I % Checkpoint_X )
      call U % Interactions_NM_G ( iR ) % Show ( )
    end do

    end select !-- I

  end subroutine InitializeIntegrator


  subroutine ShowParameters ( U )

    class ( Universe_R_CC_C_Form ), intent ( in ) :: &
      U

    call U % Universe_F_CC_Form % ShowParameters ( )

    call Show ( U % RadiationName,     'RadiationName',     U % IGNORABILITY )
    call Show ( U % RadiationType,     'RadiationType',     U % IGNORABILITY )
    call Show ( U % FormalismType,     'FormalismType',     U % IGNORABILITY )

  end subroutine ShowParameters


  subroutine ShowDiagnostics ( U )

      class ( Universe_R_CC_C_Form ), intent ( in ) :: &
        U

    integer ( KDI ) :: &
      iR

    if ( allocated ( U % Coarsening_F ) ) &
      call U % Coarsening_F % Show ( )
    do iR  =  1,  U % nRadiations
      if ( allocated ( U % Coarsening_R ) ) &
        call U % Coarsening_R ( iR ) % Show ( )
    end do !-- iR

    if ( allocated ( U % PositionSpace_AA ) ) then
      call U % PositionSpace_AA % Show ( )
      call U % AA_Gravitation % FieldSet_AA % Show ( )
      call U % AA_Fluid % FieldSet_AA % Show ( )
      call U % Stream_AA % Show ( )
    end if !-- allocated PositionSpace_AA

    if ( allocated ( U % PositionSpace_SA ) ) then
      call U % PositionSpace_SA % Show ( )
      call U % SA_Gravitation % FieldSet_SA % Show ( )
      call U % SA_Fluid % FieldSet_SA % Show ( )
      call U % Stream_SA % Show ( )
    end if !-- allocated PositionSpace_SA

  end subroutine ShowDiagnostics


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iR
    real ( KDR ), dimension ( : ), allocatable :: &
      dT_RS

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )
    select type ( I )
      class is ( Integrator_CS_1D_C_CS_Form )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_H_Form )
    associate &  !-- See InitializeIntegrator subroutine herein
      ( dT_1  =>  dT_Candidate ( 1 ), &
        dT_2  =>  dT_Candidate ( 2 ), &
        dT_3  =>  dT_Candidate ( 3 ) )

    !-- Gravity step

    call U % Compute_dT_G_CGS ( dT_1, iC, T_Option )
    dT_1  =  U % GravityFactor  *  dT_1    

    !-- Fluid advection step

    if ( U % Coarsen ) then
      call U % Compute_dT_CS_CGS_C &
             ( I % EigenspeedSet_X, dT_2, U % Coarsening_F, iC, T_Option )
    else !-- .not. Coarsen
      call I % Compute_dT_CS_CGS &
             ( I % EigenspeedSet_X, dT_2, iC, T_Option )
    end if !-- Coarsen
    dT_2  =  I % CourantFactor  *  dT_2
    
    !-- Radiation streaming step

    associate ( nR  =>  U % nRadiations )
    allocate ( dT_RS ( nR ) )
    do iR  =  1, nR
      if ( U % Coarsen ) then
        call U % Compute_dT_CS_CGS_C &
               ( I % EigenspeedSet_X_1D ( :, iR ), dT_RS ( iR ), &
                 U % Coarsening_R ( iR ), iC, T_Option )
      else !-- .not. Coarsen
        call I % Compute_dT_CS_CGS &
               ( I % EigenspeedSet_X_1D ( :, iR ), dT_RS ( iR ), iC, T_Option )
      end if
    end do !-- iR
    end associate !-- nR
    dT_3  =  I % CourantFactor_1D  *  minval ( dT_RS )

    end associate !-- dT_1, etc.
    end select !-- S
    end select !-- I
    end select !-- U

  end subroutine Compute_dT_Local


  subroutine InitializeSeries ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call U % InitializeSeries_F_CC ( I )

    end select !-- U

  end subroutine InitializeSeries


  subroutine Analyze ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call U % Analyze_F_CC ( I, Ignorability, T_Option )

    end select !-- U

  end subroutine Analyze


  subroutine Write ( I, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call U % Write_F_C ( I, T_Option )

    end select !-- U

  end subroutine Write


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I

    select type ( U  =>  I % System )
      class is ( Universe_R_CC_C_Form )

    call U % Set_T_CheckpointInterval_F_CC ( I )

    end select !-- U

  end subroutine Set_T_CheckpointInterval


  subroutine SetSlope_F_P_DFV_N ( S, K )

    !-- Compare SetSlope_N in Universe_F_C__Form

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    character ( 1 ) :: &
      StageNumber

!    select type ( S )
!      class is ( Step_RK_CS_Form )
!    select type ( F  =>  S % CurrentSet )
!    class is ( Fluid_P_HN_Form )
!      allocate ( Slope_DFV_N_F_P_HN_Form :: K )
!    class default
      allocate ( Slope_DFV_N_Form :: K )
!    end select !-- F
!    end select !-- S

    select type ( K )
      class is ( Slope_DFV_N_Form )
    select type ( S )
      class is ( Step_RK_CS_Form )
    select type ( F  =>  S % CurrentSet )
      class is ( Fluid_D_Form ) 

    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_B ( 1 ) )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_B ( 2 ) )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_B ( 3 ) )

    !-- Dust
    iEnergy_B  =  0

    !-- Perfect fluid
    select type ( F )
    class is ( Fluid_P_Form )
      call Search ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_B )
    end select !-- F

    if ( allocated ( S % DivergenceTotal ) ) then
      call K % Initialize &
             ( S % RiemannSolver, &
               S % DiffusionFactor, &
               S % DivergenceTotal, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B )
    else if ( allocated ( S % DivergencePart ) ) then
      call K % Initialize &
             ( S % RiemannSolver, &
               S % DiffusionFactor, &
               S % DivergencePart, &
               iVelocity_F = F % VELOCITY_U, &
               iMomentum_B = iMomentum_B, &
               iBaryonMass_F = F % BARYON_MASS, &
               iBaryonDensity_F = F % BARYON_DENSITY_B, &
               iEnergy_B = iEnergy_B )
    end if  !-- DivergenceTotal

    end select !-- F
    end select !-- S
    end select !-- K

  end subroutine SetSlope_F_P_DFV_N


end module Universe_R_CC_C__Form
