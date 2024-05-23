module Universe_F_C__Form

  !-- Universe_Fluid_Central__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_H__Form

  implicit none
  private

  type, public, extends ( Universe_H_Form ) :: Universe_F_C_Form
    real ( KDR ) :: &
      GravityFactor = 0.0_KDR
    logical ( KDL ) :: &
      Coarsen
    type ( Coarsening_C_F_Form ), allocatable :: &
      Coarsening_F
    class ( Atlas_SCG_Form ), allocatable :: &
      PositionSpace_SA, &  !-- SphericalAverage
      PositionSpace_AA     !-- AzimuthalAverage
    type ( Stream_BM_Form ), allocatable :: &
      Stream_SA, &
      Stream_AA
    type ( SphericalAverageForm ), allocatable :: &
      SA_Gravitation, &
      SA_Fluid
    type ( AzimuthalAverageForm ), allocatable :: &
      AA_Gravitation, &
      AA_Fluid
    type ( Units_F_Form ), dimension ( : ), allocatable :: &
      Units_F
  contains
    procedure, public, pass :: &
      Initialize_F_C
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator
    procedure, public, pass :: &
      InitializePositionSpace
    procedure, public, pass :: &
      InitializeGravitation
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, public, pass :: &
      InitializeStep
    procedure, public, pass :: &
      InitializeIntegrator
    procedure, private, pass :: &
      InitializeAtlas
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass :: &
      ShowDiagnostics
    procedure, public, pass :: &
      Compute_dT_CS_CGS_C
    procedure, public, nopass :: &
      Analyze_F_C
    procedure, public, nopass :: &
      Write_F_C
  end type Universe_F_C_Form

    private :: &
      SetSlope_N

    private :: &
      Compute_dT_CS_CGS_C_Kernel

    interface
    
      module subroutine Compute_dT_CS_CGS_C_Kernel &
               ( dT, ProperCell, FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 dX_1, dX_2, dX_3, Crsn_2, Crsn_3, &
                 nDimensions, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          FEP_1, FEP_2, FEP_3, &
          FEM_1, FEM_2, FEM_3, &
          dX_1, dX_2, dX_3, &
          Crsn_2, Crsn_3
        integer ( KDI ), intent ( in ) :: &
          nDimensions
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_CS_CGS_C_Kernel

    end interface


contains


  subroutine Initialize_F_C &
               ( U, FluidType, GravitationType, Name, &
                 UnitsTypeOption, FinishTimeOption, RadiusMaxOption, &
                 RadiusCoreOption, RadiusExcisionOption, RadialRatioOption, &
                 GravityFactorOption, CentralMassOption, nCellsPolarOption, &
                 nWriteOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType, &
      Name
    character ( * ), intent ( in ), optional :: &
      UnitsTypeOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption, &
      GravityFactorOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( U % Type == '' ) &
      U % Type = 'a Universe_F_C'

    call U % Universe_H_Form % Initialize &
           ( Name, UnitsTypeOption = UnitsTypeOption )

    call U % AllocateIntegrator &
           ( )
    call U % InitializePositionSpace &
           ( RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )
    call U % InitializeGravitation &
           ( GravitationType, &
             CentralMassOption = CentralMassOption )
    call U % InitializeFluid &
           ( FluidType )
    call U % SetBoundaryConditions &
           ( )
    call U % InitializeStep &
           ( )
    call U % InitializeIntegrator &
           ( GravitationType, &
             FinishTimeOption = FinishTimeOption, &
             GravityFactorOption = GravityFactorOption, &
             nWriteOption = nWriteOption )

  end subroutine Initialize_F_C


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Units_F ) ) &
      deallocate ( U % Units_F )
    if ( allocated ( U % AA_Fluid ) ) &
      deallocate ( U % AA_Fluid )
    if ( allocated ( U % AA_Gravitation ) ) &
      deallocate ( U % AA_Gravitation )
    if ( allocated ( U % SA_Fluid ) ) &
      deallocate ( U % SA_Fluid )
    if ( allocated ( U % SA_Gravitation ) ) &
      deallocate ( U % SA_Gravitation )
    if ( allocated ( U % Stream_AA ) ) &
      deallocate ( U % Stream_AA )
    if ( allocated ( U % Stream_SA ) ) &
      deallocate ( U % Stream_SA )
    if ( allocated ( U % PositionSpace_AA ) ) &
      deallocate ( U % PositionSpace_AA )
    if ( allocated ( U % PositionSpace_SA ) ) &
      deallocate ( U % PositionSpace_SA )
    if ( allocated ( U % Coarsening_F ) ) &
      deallocate ( U % Coarsening_F )
    
  end subroutine Finalize


  subroutine AllocateIntegrator ( U )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    allocate ( Integrator_CS_Form :: U % Integrator )

    if ( allocated ( U % dT_Label ) ) then
      associate ( I => U % Integrator )
      allocate ( I % dT_Label, source = U % dT_Label )
      end associate !-- I
    end if

  end subroutine AllocateIntegrator


  subroutine InitializePositionSpace &
               ( U, CommunicatorOption, RadiusMaxOption, RadiusCoreOption, &
                 RadiusExcisionOption, RadialRatioOption, nCellsPolarOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    call U % InitializeAtlas &
           ( CommunicatorOption = CommunicatorOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )

  end subroutine InitializePositionSpace


  subroutine InitializeGravitation ( U, GravitationType, CentralMassOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
    real ( KDR ), intent ( in ), optional :: &
      CentralMassOption

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAverage
    real ( KDR ) :: &
      GravitationalConstant 

    associate ( I  =>  U % Integrator )

    select case ( trim ( GravitationType ) )
    case ( 'GALILEO' )

      allocate ( Gravitation_G_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_G_Form )
      call G % Initialize &
             ( I % X, &
               DeviceMemoryOption = U % DeviceMemory, &
               PinnedMemoryOption = U % PinnedMemory, &
               DevicesCommunicateOption = U % DevicesCommunicate )

      allocate ( iaAverage ( 0 ) )

      !-- Azimuthal average
      if ( allocated ( U % PositionSpace_AA ) ) then
        allocate ( U % AA_Gravitation )
        associate &
          ( AA     =>  U % AA_Gravitation, &
             A_AA  =>  U % PositionSpace_AA )
        allocate ( Gravitation_G_Form :: AA % FieldSet_AA )
        select type ( G_AA  =>  AA % FieldSet_AA )
          type is ( Gravitation_G_Form )
        call G_AA % Initialize ( A_AA, NameOption = trim ( G % Name ) // '_AA' )
        call AA % Initialize ( G, G, A_AA, iaAverageOption = iaAverage )
        end select !-- G_AA
        end associate !-- AA, etc.
      end if !-- allocated PositionSpace_AA

      !-- Spherical average
      if ( allocated ( U % PositionSpace_SA ) ) then
        allocate ( U % SA_Gravitation )
        associate &
          ( SA     =>  U % SA_Gravitation, &
             A_SA  =>  U % PositionSpace_SA )
        allocate ( Gravitation_G_Form :: SA % FieldSet_SA )
        select type ( G_SA  =>  SA % FieldSet_SA )
          type is ( Gravitation_G_Form )
        call G_SA % Initialize ( A_SA, NameOption = trim ( G % Name ) // '_SA' )
        call SA % Initialize ( G, G, A_SA, iaAverageOption = iaAverage )
        end select !-- G_SA
        end associate !-- SA, etc.
      end if !-- allocated PositionSpace_SA

      end select !-- G

    case ( 'NEWTON_CM' )

      if ( trim ( U % UnitsType )  ==  '' ) then
        GravitationalConstant  =  1.0_KDR
      else
        GravitationalConstant  =  CONSTANT % GRAVITATIONAL
      end if

      if ( .not. present ( CentralMassOption ) ) then
        call Show ( 'CentralMassOption not present', CONSOLE % ERROR )
        call Show ( 'NEWTON_CM', 'GravitationType', CONSOLE % ERROR )
        call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      allocate ( Gravitation_N_CM_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_N_CM_Form )
      call G % Initialize &
             ( I % X, GravitationalConstant, &
               Mass = CentralMassOption, &
               DeviceMemoryOption = U % DeviceMemory, &
               PinnedMemoryOption = U % PinnedMemory, &
               DevicesCommunicateOption = U % DevicesCommunicate )

      allocate &
        ( iaAverage, source = [ G % POTENTIAL, G % POTENTIAL_GRADIENT_D ] )

      !-- Azimuthal average
      if ( allocated ( U % PositionSpace_AA ) ) then
        allocate ( U % AA_Gravitation )
        associate &
          ( AA     =>  U % AA_Gravitation, &
             A_AA  =>  U % PositionSpace_AA )
        allocate ( Gravitation_N_H_Form :: AA % FieldSet_AA )
        select type ( G_AA  =>  AA % FieldSet_AA )
          type is ( Gravitation_N_H_Form )
        call G_AA % Initialize ( A_AA, NameOption = trim ( G % Name ) // '_AA' )
        call AA % Initialize ( G, G, A_AA, iaAverageOption = iaAverage )
        end select !-- G_AA
        end associate !-- AA, etc.
      end if !-- allocated PositionSpace_AA

      !-- Spherical average
      if ( allocated ( U % PositionSpace_SA ) ) then
        allocate ( U % SA_Gravitation )
        associate &
          ( SA     =>  U % SA_Gravitation, &
             A_SA  =>  U % PositionSpace_SA )
        allocate ( Gravitation_N_H_Form :: SA % FieldSet_SA )
        select type ( G_SA  =>  SA % FieldSet_SA )
          type is ( Gravitation_N_H_Form )
        call G_SA % Initialize ( A_SA, NameOption = trim ( G % Name ) // '_SA' )
        call SA % Initialize ( G, G, A_SA, iaAverageOption = iaAverage )
        end select !-- G_SA
        end associate !-- SA, etc.
      end if !-- allocated PositionSpace_SA

      end select !-- G

    case ( 'NEWTON_SG' )

      if ( trim ( U % UnitsType )  ==  '' ) then
        GravitationalConstant  =  1.0_KDR
      else
        GravitationalConstant  =  CONSTANT % GRAVITATIONAL
      end if

      allocate ( Gravitation_N_SG_Form  ::  I % Geometry_X )
      select type ( G  =>  I % Geometry_X )
        class is ( Gravitation_N_SG_Form )
      call G % Initialize &
             ( I % X, GravitationalConstant, &
               DeviceMemoryOption = U % DeviceMemory, &
               PinnedMemoryOption = U % PinnedMemory, &
               DevicesCommunicateOption = U % DevicesCommunicate )

      allocate &
        ( iaAverage, source = [ G % POTENTIAL, G % POTENTIAL_GRADIENT_D ] )

      !-- Azimuthal average
      if ( allocated ( U % PositionSpace_AA ) ) then
        allocate ( U % AA_Gravitation )
        associate &
          ( AA     =>  U % AA_Gravitation, &
             A_AA  =>  U % PositionSpace_AA )
        allocate ( Gravitation_N_H_Form :: AA % FieldSet_AA )
        select type ( G_AA  =>  AA % FieldSet_AA )
          type is ( Gravitation_N_H_Form )
        call G_AA % Initialize ( A_AA, NameOption = trim ( G % Name ) // '_AA' )
        call AA % Initialize ( G, G, A_AA, iaAverageOption = iaAverage )
        end select !-- G_AA
        end associate !-- AA, etc.
      end if !-- allocated PositionSpace_AA

      !-- Spherical average
      if ( allocated ( U % PositionSpace_SA ) ) then
        allocate ( U % SA_Gravitation )
        associate &
          ( SA     =>  U % SA_Gravitation, &
             A_SA  =>  U % PositionSpace_SA )
        allocate ( Gravitation_N_H_Form :: SA % FieldSet_SA )
        select type ( G_SA  =>  SA % FieldSet_SA )
          type is ( Gravitation_N_H_Form )
        call G_SA % Initialize ( A_SA, NameOption = trim ( G % Name ) // '_SA' )
        call SA % Initialize ( G, G, A_SA, iaAverageOption = iaAverage )
        end select !-- G_SA
        end associate !-- SA, etc.
      end if !-- allocated PositionSpace_SA

      end select !-- G

    case default
      call Show ( 'GravitationType not recognized', CONSOLE % ERROR )
      call Show ( GravitationType, 'GravitationType', CONSOLE % ERROR )
      call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeGravitation', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GravitationType

    end associate !-- I

  end subroutine InitializeGravitation


  subroutine InitializeFluid ( U, FluidType )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType

    integer ( KDI ) :: &
      iB  !-- iBoundary

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( G  =>  I % Geometry_X )

    select type ( A  =>  I % X )
    class is ( Atlas_SCG_Form )
    associate ( C  =>  A % Chart_GS )
      allocate ( U % Units_F ( 1 ) )
      call U % Units_F ( 1 ) % Initialize &
             ( C % CoordinateUnit, TypeOption = U % UnitsType )
    end associate !-- C
    end select !-- A

    select case ( trim ( FluidType ) )
    case ( 'DUST' )

      allocate ( Fluid_D_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )

        !-- TallyInterior must be allocated before F % Initialize...
        allocate ( Tally_F_D_Form :: F % TallyInterior )
        allocate ( Tally_F_D_Form :: F % TallyTotal )
        allocate ( Tally_F_D_Form :: F % TallyChange )
        select type ( TI  =>  F % TallyInterior )
          type is ( Tally_F_D_Form )
        call TI % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TI
        select type ( TT  =>  F % TallyTotal )
          type is ( Tally_F_D_Form )
        call TT % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TT
        select type ( TC  =>  F % TallyChange )
          type is ( Tally_F_D_Form )
        call TC % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TC

        !-- ( Initialize the Fluid )
        call F % Initialize ( G, U % Units_F )

        !-- ... but TallyBoundary needs F % Initialize already called.
        allocate ( F % TallyBoundary ( F % nBoundaries ) )
        do iB  =  1,  F % nBoundaries
          allocate ( Tally_F_D_Form :: F % TallyBoundary ( iB ) % Element )
          select type ( TB  =>  F % TallyBoundary ( iB ) % Element )
            type is ( Tally_F_D_Form )
          call TB % Initialize ( G, U % Units_F ( 1 ) )
          end select !-- TB
        end do !-- iB

        !-- Boundary accumulation storage
        call F % AllocateBoundary_SCG ( nT = F % TallyInterior % nSelected )

        !-- Azimuthal average
        if ( allocated ( U % PositionSpace_AA ) ) then
          allocate ( U % AA_Fluid )
          associate &
            ( AA     =>  U % AA_Fluid, &
               A_AA  =>  U % PositionSpace_AA )
          allocate ( Fluid_D_Form :: AA % FieldSet_AA )
          select type ( F_AA  =>  AA % FieldSet_AA )
            type is ( Fluid_D_Form )
          select type ( G_AA  =>  U % AA_Gravitation % FieldSet_AA )
            class is ( Geometry_F_Form )
          call F_AA % Initialize &
                 ( G_AA, U % Units_F, NameOption = trim ( F % Name ) // '_AA' )
          call AA % Initialize ( G, F, A_AA, iaAverageOption = F % iaBalanced )
          end select !-- G_AA
          end select !-- F_AA
          end associate !-- AA, etc.
        end if !-- allocated PositionSpace_AA

        !-- Spherical average
        if ( allocated ( U % PositionSpace_SA ) ) then
          allocate ( U % SA_Fluid )
          associate &
            ( SA     =>  U % SA_Fluid, &
               A_SA  =>  U % PositionSpace_SA )
          allocate ( Fluid_D_Form :: SA % FieldSet_SA )
          select type ( F_SA  =>  SA % FieldSet_SA )
            type is ( Fluid_D_Form )
          select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
            class is ( Geometry_F_Form )
          call F_SA % Initialize &
                 ( G_SA, U % Units_F, NameOption = trim ( F % Name ) // '_SA' )
          call SA % Initialize ( G, F, A_SA, iaAverageOption = F % iaBalanced )
          end select !-- G_SA
          end select !-- F_SA
          end associate !-- SA, etc.
        end if !-- allocated PositionSpace_SA

      end select !-- F
      
    case ( 'IDEAL' )

      allocate ( Fluid_P_I_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_I_Form )

        !-- TallyInterior must be allocated before F % Initialize...
        allocate ( Tally_F_P_Form :: F % TallyInterior )
        allocate ( Tally_F_P_Form :: F % TallyTotal )
        allocate ( Tally_F_P_Form :: F % TallyChange )
        select type ( TI  =>  F % TallyInterior )
          type is ( Tally_F_P_Form )
        call TI % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TI
        select type ( TT  =>  F % TallyTotal )
          type is ( Tally_F_P_Form )
        call TT % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TT
        select type ( TC  =>  F % TallyChange )
          type is ( Tally_F_P_Form )
        call TC % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TC

        !-- ( Initialize the Fluid )
        call F % Initialize ( G, U % Units_F )

        !-- ... but TallyBoundary needs F % Initialize already called.
        allocate ( F % TallyBoundary ( F % nBoundaries ) )
        do iB  =  1,  F % nBoundaries
          allocate ( Tally_F_P_Form :: F % TallyBoundary ( iB ) % Element )
          select type ( TB  =>  F % TallyBoundary ( iB ) % Element )
            type is ( Tally_F_P_Form )
          call TB % Initialize ( G, U % Units_F ( 1 ) )
          end select !-- TB
        end do !-- iB

        !-- Boundary accumulation storage
        call F % AllocateBoundary_SCG ( nT = F % TallyInterior % nSelected )

        !-- Azimuthal average
        if ( allocated ( U % PositionSpace_AA ) ) then
          allocate ( U % AA_Fluid )
          associate &
            ( AA     =>  U % AA_Fluid, &
               A_AA  =>  U % PositionSpace_AA )
          allocate ( Fluid_P_I_Form :: AA % FieldSet_AA )
          select type ( F_AA  =>  AA % FieldSet_AA )
            type is ( Fluid_P_I_Form )
          select type ( G_AA  =>  U % AA_Gravitation % FieldSet_AA )
            class is ( Geometry_F_Form )
          call F_AA % Initialize &
                 ( G_AA, U % Units_F, NameOption = trim ( F % Name ) // '_AA' )
          call AA % Initialize ( G, F, A_AA, iaAverageOption = F % iaBalanced )
          end select !-- G_AA
          end select !-- F_AA
          end associate !-- AA, etc.
        end if !-- allocated PositionSpace_AA

        !-- Spherical average
        if ( allocated ( U % PositionSpace_SA ) ) then
          allocate ( U % SA_Fluid )
          associate &
            ( SA     =>  U % SA_Fluid, &
               A_SA  =>  U % PositionSpace_SA )
          allocate ( Fluid_P_I_Form :: SA % FieldSet_SA )
          select type ( F_SA  =>  SA % FieldSet_SA )
            type is ( Fluid_P_I_Form )
          select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
            class is ( Geometry_F_Form )
          call F_SA % Initialize &
                 ( G_SA, U % Units_F, NameOption = trim ( F % Name ) // '_SA' )
          call SA % Initialize ( G, F, A_SA, iaAverageOption = F % iaBalanced )
          end select !-- G_SA
          end select !-- F_SA
          end associate !-- SA, etc.
        end if !-- allocated PositionSpace_SA

      end select !-- F

    case ( 'HEAVY_NUCLEUS' )
      
      allocate ( Fluid_P_HN_Form  ::  I % CurrentSet_X )
      select type ( F  =>  I % CurrentSet_X )
        class is ( Fluid_P_HN_Form )

        !-- TallyInterior must be allocated before F % Initialize...
        allocate ( Tally_F_P_HN_Form :: F % TallyInterior )
        allocate ( Tally_F_P_HN_Form :: F % TallyTotal )
        allocate ( Tally_F_P_HN_Form :: F % TallyChange )
        select type ( TI  =>  F % TallyInterior )
          type is ( Tally_F_P_HN_Form )
        call TI % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TI
        select type ( TT  =>  F % TallyTotal )
          type is ( Tally_F_P_HN_Form )
        call TT % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TT
        select type ( TC  =>  F % TallyChange )
          type is ( Tally_F_P_HN_Form )
        call TC % Initialize ( G, U % Units_F ( 1 ) )
        end select !-- TC

        !-- ( Initialize the Fluid )
        call F % Initialize ( G, U % Units_F )

        !-- ... but TallyBoundary needs F % Initialize already called.
        allocate ( F % TallyBoundary ( F % nBoundaries ) )
        do iB  =  1,  F % nBoundaries
          allocate ( Tally_F_P_HN_Form :: F % TallyBoundary ( iB ) % Element )
          select type ( TB  =>  F % TallyBoundary ( iB ) % Element )
            type is ( Tally_F_P_HN_Form )
          call TB % Initialize ( G, U % Units_F ( 1 ) )
          end select !-- TB
        end do !-- iB

        !-- Boundary accumulation storage
        call F % AllocateBoundary_SCG ( nT = F % TallyInterior % nSelected )

        !-- Azimuthal average
        if ( allocated ( U % PositionSpace_AA ) ) then
          allocate ( U % AA_Fluid )
          associate &
            ( AA     =>  U % AA_Fluid, &
               A_AA  =>  U % PositionSpace_AA )
          allocate ( Fluid_P_HN_Form :: AA % FieldSet_AA )
          select type ( F_AA  =>  AA % FieldSet_AA )
            type is ( Fluid_P_HN_Form )
          select type ( G_AA  =>  U % AA_Gravitation % FieldSet_AA )
            class is ( Geometry_F_Form )
          call F_AA % Initialize &
                 ( G_AA, U % Units_F, NameOption = trim ( F % Name ) // '_AA' )
          call AA % Initialize &
                 ( G, F, A_AA, &
                   iaAverageOption = [ F % iaBalanced, F % TEMPERATURE ] )
          end select !-- G_AA
          end select !-- F_AA
          end associate !-- AA, etc.
        end if !-- allocated PositionSpace_AA

        !-- Spherical average
        if ( allocated ( U % PositionSpace_SA ) ) then
          allocate ( U % SA_Fluid )
          associate &
            ( SA     =>  U % SA_Fluid, &
               A_SA  =>  U % PositionSpace_SA )
          allocate ( Fluid_P_HN_Form :: SA % FieldSet_SA )
          select type ( F_SA  =>  SA % FieldSet_SA )
            type is ( Fluid_P_HN_Form )
          select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
            class is ( Geometry_F_Form )
          call F_SA % Initialize &
                 ( G_SA, U % Units_F, NameOption = trim ( F % Name ) // '_SA' )
          call SA % Initialize &
                 ( G, F, A_SA, &
                   iaAverageOption = [ F % iaBalanced, F % TEMPERATURE ] )
          end select !-- G_SA
          end select !-- F_SA
          end associate !-- SA, etc.
        end if !-- allocated PositionSpace_SA

      end select !-- F

    case default
      call Show ( 'FluidType not recognized', CONSOLE % ERROR )
      call Show ( FluidType, 'FluidType', CONSOLE % ERROR )
      call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeFluid', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- FluidType

    end associate !-- G
    end select !-- I

  end subroutine InitializeFluid


  subroutine SetBoundaryConditions ( U )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    call Show ( 'SetBoundaryConditions should be overridden', &
                CONSOLE % WARNING )
    call Show ( 'Universe_F_C__Form', 'module', CONSOLE % WARNING )
    call Show ( 'SetBoundaryConditions', 'subroutine', CONSOLE % WARNING )

  end subroutine SetBoundaryConditions


  subroutine InitializeStep ( U )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U

    integer ( KDI ) :: &
      nStages
    logical ( KDL ) :: &
      DivergenceParts
    character ( LDL ) :: &
      RiemannSolverType

    nStages  =  2
    call PROGRAM_HEADER % GetParameter ( nStages, 'nStages' )

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    associate &
      ( G  =>  I % Geometry_X )

    allocate ( Step_RK_CS_Form :: I % Step_X )
    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_Form )

    select type ( G )
    class is ( Gravitation_N_H_Form )
      S % SetSlopeExplicit  =>  SetSlope_N
    end select !-- G

    DivergenceParts  =  .false.
    call PROGRAM_HEADER % GetParameter ( DivergenceParts, 'DivergenceParts' )

    if ( DivergenceParts ) then
      select type ( F )
      type is ( Fluid_D_Form )
        allocate ( S % DivergencePart ( 1 ) )
        associate ( DP_1D  =>  S % DivergencePart )
          allocate ( DivergencePart_F_D_V_Form :: DP_1D ( 1 ) % Element )
          associate ( DP  =>  DP_1D ( 1 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
        end associate !-- DP_1D
      type is ( Fluid_P_I_Form )
        allocate ( S % DivergencePart ( 2 ) )
        associate ( DP_1D  =>  S % DivergencePart )
          allocate ( DivergencePart_F_P_V_Form :: DP_1D ( 1 ) % Element )
          associate ( DP  =>  DP_1D ( 1 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
          allocate ( DivergencePart_F_P_P_Form :: DP_1D ( 2 ) % Element )
          associate ( DP  =>  DP_1D ( 2 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
        end associate !-- DP_1D
      type is ( Fluid_P_HN_Form )
        allocate ( S % DivergencePart ( 2 ) )
        associate ( DP_1D  =>  S % DivergencePart )
          allocate ( DivergencePart_F_P_HN_V_Form :: DP_1D ( 1 ) % Element )
          associate ( DP  =>  DP_1D ( 1 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
          allocate ( DivergencePart_F_P_HN_P_Form :: DP_1D ( 2 ) % Element )
          associate ( DP  =>  DP_1D ( 2 ) % Element )
            call DP % Initialize ( F )
          end associate !-- DP
        end associate !-- DP_1D
      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- F
    else  !-- DivergenceTotal
      select type ( F )
      type is ( Fluid_D_Form )
        allocate ( DivergencePart_F_D_T_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
          call DT % Initialize ( F )
        end associate !-- DT
      type is ( Fluid_P_I_Form )

        allocate ( DivergencePart_F_P_T_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
          call DT % Initialize ( F )
        end associate !-- DT

        RiemannSolverType = 'HLL'
        call PROGRAM_HEADER % GetParameter &
               ( RiemannSolverType, 'RiemannSolverType' )
        if ( trim ( RiemannSolverType ) == 'HLLC' ) then
          allocate ( RiemannSolver_HLLC_P_Form :: S % RiemannSolver )
          associate ( RS  =>  S % RiemannSolver )
          call RS % Initialize ( F )
          end associate !-- RS
        end if
        
      type is ( Fluid_P_HN_Form )

        allocate ( DivergencePart_F_P_HN_T_Form :: S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
          call DT % Initialize ( F )
        end associate !-- DT

        RiemannSolverType = 'HLL'
        call PROGRAM_HEADER % GetParameter &
               ( RiemannSolverType, 'RiemannSolverType' )
        if ( trim ( RiemannSolverType ) == 'HLLC' ) then
          allocate ( RiemannSolver_HLLC_P_HN_Form :: S % RiemannSolver )
          associate ( RS  =>  S % RiemannSolver )
          call RS % Initialize ( F )
          end associate !-- RS
        end if
        
      class default
        call Show ( 'Fluid type not recognized', CONSOLE % ERROR )
        call Show ( 'Universe_F_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'InitializeStep', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- F
    end if  !-- DivergenceParts

    call S % Initialize ( F, nStagesOption = nStages )

    !-- Coarsening_F
    U % Coarsen  =  .true.
    call PROGRAM_HEADER % GetParameter ( U % Coarsen, 'Coarsen' )
    if ( U % Coarsen ) then
      allocate ( U % Coarsening_F )
      associate ( C  =>  U % Coarsening_F )
      select type ( F )
        type is ( Fluid_D_Form )
          call C % Initialize &
                 ( F, G, RadiusZeroOption = 0.0_KDR, nRadialZeroOption = 0, &
                   nPolarZeroOption = 0 )
        class default
          call C % Initialize ( F, G )
      end select !-- F
      call S % SetCoarsening ( C )
      end associate !-- C
    end if

    end select !-- S

    end associate !-- G
    end select    !-- F
    end select    !-- I

  end subroutine InitializeStep


  subroutine InitializeIntegrator &
               ( U, GravitationType, FinishTimeOption, GravityFactorOption, &
                 nWriteOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in ) :: &
      GravitationType
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      GravityFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    if ( .not. allocated ( I % dT_Label ) ) then
      if ( any ( trim ( GravitationType ) == [ 'NEWTON_SG' ] ) ) then
        allocate ( I % dT_Label ( 2 ) )
        I % dT_Label ( 1 )  =  'FluidAdvection'
        I % dT_Label ( 2 )  =  'GravitationAcceleration'
        U % GravityFactor  =  0.7_KDR
        if ( present ( GravityFactorOption ) ) &
          U % GravityFactor  =  GravityFactorOption
        call PROGRAM_HEADER % GetParameter &
               ( U % GravityFactor, 'GravityFactor' )
      else
        allocate ( I % dT_Label ( 1 ) )
        I % dT_Label ( 1 )  =  'FluidAdvection'
      end if
    end if

    call I % Initialize &
           ( CommunicatorOption = U % Communicator, &
             Unit_T_Option = U % Units_F ( 1 ) % Time, &
             T_FinishOption = FinishTimeOption, &
             nWriteOption = nWriteOption )

    !-- AzimuthalAverage Stream

    if ( allocated ( U % PositionSpace_AA ) ) then
      allocate ( U % Stream_AA )
      associate &
        (   A_AA  =>  U % PositionSpace_AA, &
            S_AA  =>  U % Stream_AA, &
          GIS     =>  I % GridImageStream, &
            S     =>  I % Checkpoint_X )

      call S_AA % Initialize &
             ( A_AA, GIS, NameOption = trim ( S % Name ) // '_AA' )

      select type ( G_AA  =>  U % AA_Gravitation % FieldSet_AA )
        class is ( Geometry_F_Form )
      select type ( F_AA  =>  U % AA_Fluid % FieldSet_AA )
        class is ( Fluid_D_Form )
      call G_AA % SetStream ( S_AA )
      call F_AA % SetStream ( S_AA )
      end select !-- F_AA
      end select !-- G_AA

      end associate !-- A_AA, etc.
    end if !-- allocated PositionSpace_AA

    !-- SphericalAverage Stream

    if ( allocated ( U % PositionSpace_SA ) ) then
      allocate ( U % Stream_SA )
      associate &
        (   A_SA  =>  U % PositionSpace_SA, &
            S_SA  =>  U % Stream_SA, &
          GIS     =>  I % GridImageStream, &
            S     =>  I % Checkpoint_X )

      call S_SA % Initialize &
             ( A_SA, GIS, NameOption = trim ( S % Name ) // '_SA' )

      select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
        class is ( Geometry_F_Form )
      select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
        class is ( Fluid_D_Form )
      call G_SA % SetStream ( S_SA )
      call F_SA % SetStream ( S_SA )
      end select !-- F_SA
      end select !-- G_SA

      end associate !-- A_SA, etc.
    end if !-- allocated PositionSpace_SA

    !-- Integrator methods

    I % Analyze  =>  Analyze_F_C
    I % Write    =>  Write_F_C

    end select !-- I

  end subroutine InitializeIntegrator


  subroutine InitializeAtlas &
               ( U, CommunicatorOption, RadiusMaxOption, RadiusCoreOption, &
                 RadiusExcisionOption, RadialRatioOption, nCellsPolarOption )

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    call Show ( 'InitializeAtlas should be overridden', CONSOLE % WARNING )
    call Show ( 'Universe_F_C__Form', 'module', CONSOLE % WARNING )
    call Show ( 'InitializeAtlas', 'subroutine', CONSOLE % WARNING )

  end subroutine InitializeAtlas


  subroutine ShowParameters ( U )

      class ( Universe_F_C_Form ), intent ( in ) :: &
        U

    call U % Universe_H_Form % ShowParameters ( )

    call Show ( U % Coarsen, 'Coarsen' )

    if ( U % GravityFactor  >  0.0_KDR ) &
      call Show ( U % GravityFactor, 'GravityFactor' )

  end subroutine ShowParameters


  subroutine ShowDiagnostics ( U )

      class ( Universe_F_C_Form ), intent ( in ) :: &
        U

    if ( allocated ( U % Coarsening_F ) ) &
      call U % Coarsening_F % Show ( )

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


  subroutine Compute_dT_CS_CGS_C ( U, dT, iC, T_Option )

    !-- Compute_dT_CurrentSet_ChartGridStructured_Central (or _Coarsened)

    class ( Universe_F_C_Form ), intent ( inout ) :: &
      U
    real ( KDR ), intent ( inout ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( ES_1  =>  I % EigenspeedSet_X ( 1 ), &
        ES_2  =>  I % EigenspeedSet_X ( 2 ), &
        ES_3  =>  I % EigenspeedSet_X ( 3 ), &
         G    =>  I % Geometry_X, &
        Crsn  =>  U % Coarsening_F )

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    call ES_1 % Compute ( iC = 1, iD = 1 )
    call ES_2 % Compute ( iC = 1, iD = 2 )
    call ES_3 % Compute ( iC = 1, iD = 3 )

    associate &
      ( EV_1  =>  ES_1 % Storage ( 1 ) % Value, &
        EV_2  =>  ES_2 % Storage ( 1 ) % Value, &
        EV_3  =>  ES_3 % Storage ( 1 ) % Value, &
        GV    =>   G   % Storage ( 1 ) % Value, &
        CV    =>  Crsn % Storage ( 1 ) % Value )

    call Compute_dT_CS_CGS_C_Kernel &
           ( dT, C % ProperCell, &
             EV_1 ( :, ES_1 % EIGENSPEED_FAST_PLUS_U ), &
             EV_2 ( :, ES_2 % EIGENSPEED_FAST_PLUS_U ), &
             EV_3 ( :, ES_3 % EIGENSPEED_FAST_PLUS_U ), &
             EV_1 ( :, ES_1 % EIGENSPEED_FAST_MINUS_U ), &
             EV_2 ( :, ES_2 % EIGENSPEED_FAST_MINUS_U ), &
             EV_3 ( :, ES_3 % EIGENSPEED_FAST_MINUS_U ), &
             GV ( :, G % WIDTH_U_1 ), &
             GV ( :, G % WIDTH_U_2 ), &
             GV ( :, G % WIDTH_U_3 ), &
             CV ( :, Crsn % COARSENING_POLAR ), &
             CV ( :, Crsn % COARSENING_AZIMUTHAL ), &
             C % nDimensions, &
             UseDeviceOption = G % DeviceMemory )

    end associate !-- EV, etc.
    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Integrator_CS_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_dT_CS_CGS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- ES_1, etc.
    end select !-- I

  end subroutine Compute_dT_CS_CGS_C


  subroutine Analyze_F_C ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( I )
      class is ( Integrator_CS_Form )
    call I % Analyze_CS ( I, Ignorability, T_Option )
    end select !-- I

    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )

    !-- Azimuthal average

    if ( allocated ( U % PositionSpace_AA ) ) then
      select type ( F_AA  =>  U % AA_Fluid % FieldSet_AA )
        class is ( Fluid_D_Form )

      call F_AA % ComputeFromInitial ( )  !-- Ensure BARYON_MASS set

      call U % AA_Gravitation % Compute ( )
      call U % AA_Fluid % Compute ( )

      call F_AA % ComputeFromBalanced ( )

      end select !-- F_AA
    end if !-- allocated PositionSpace_AA

    !-- Spherical average

    if ( allocated ( U % PositionSpace_SA ) ) then
      select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
        class is ( Fluid_D_Form )

      call F_SA % ComputeFromInitial ( )  !-- Ensure BARYON_MASS set

      call U % SA_Gravitation % Compute ( )
      call U % SA_Fluid % Compute ( )

      call F_SA % ComputeFromBalanced ( )

      end select !-- F_SA
    end if !-- allocated PositionSpace_SA

    end select !-- U

  end subroutine Analyze_F_C


  subroutine Write_F_C ( I, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_PS, &
      T_S

    !-- Reproduce and expand Write_H functionality rather than call Write_H 
    !   functionality to avoid multiple GIS % Open calls and multiple file sets.

    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )

    !-- Position space

    associate &
      ( GIS     =>  I % GridImageStream, &
          S_PS  =>  I % Checkpoint_X )
    if ( present ( T_Option ) ) then
      T_PS  =>  S_PS % TimerWrite ( Level = T_Option % Level + 1 )
    else
      T_PS  => null ( )
    end if
    if ( associated ( T_PS ) ) call T_PS % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )

    call S_PS % Write &
           ( TimeOption  =  I % T  /  I % Unit_T, &
             CycleNumberOption  =  I % iCycle )

    !-- Azimuthal average
    if ( allocated ( U % PositionSpace_AA ) ) then
      associate ( S_AA  =>  U % Stream_AA )
      call S_AA % Write &
             ( TimeOption  =  I % T  /  I % Unit_T, &
               CycleNumberOption  =  I % iCycle )
      end associate !-- S_AA
    end if !-- allocated PositionSpace_AA

    !-- Spherical average
    if ( allocated ( U % PositionSpace_SA ) ) then
      associate ( S_SA  =>  U % Stream_SA )
      call S_SA % Write &
             ( TimeOption  =  I % T  /  I % Unit_T, &
               CycleNumberOption  =  I % iCycle )
      end associate !-- S_SA
    end if !-- allocated PositionSpace_SA

    call GIS % Close ( )
    if ( associated ( T_PS ) ) call T_PS % Stop ( )
    end associate !-- GIS, etc.

    !-- Series

    if ( allocated ( I % Series ) ) then
      associate ( S  =>  I % Series )
      if ( present ( T_Option ) ) then
        T_S  =>  S % TimerWrite ( Level = T_Option % Level + 1 )
      else
        T_S  =>  null ( )
      end if
      if ( associated ( T_S ) ) call T_S % Start ( )
      call S % Write ( )
      if ( associated ( T_S ) ) call T_S % Stop ( )
      end associate !-- S
    end if

    end select !-- U

  end subroutine Write_F_C


  subroutine SetSlope_N ( S, K )

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

  end subroutine SetSlope_N


end module Universe_F_C__Form
