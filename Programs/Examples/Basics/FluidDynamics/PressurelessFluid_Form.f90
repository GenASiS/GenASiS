module PressurelessFluid_Form

  use Basics
  use DistributedMesh_Form
  use ConservedFields_Template
  use ConservationLawStep_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PRESSURELESS = 4, &
      N_CONSERVED_PRESSURELESS = 4, &
      N_FIELDS_PRESSURELESS    = 8, &
      N_VECTORS_PRESSURELESS   = 2

  type, public, extends ( ConservedFieldsTemplate ) :: &
    PressurelessFluidForm 
      integer ( KDI ) :: &
        N_PRIMITIVE_PRESSURELESS = N_PRIMITIVE_PRESSURELESS, &
        N_CONSERVED_PRESSURELESS = N_CONSERVED_PRESSURELESS, &
        N_FIELDS_PRESSURELESS    = N_FIELDS_PRESSURELESS, &
        N_VECTORS_PRESSURELESS   = N_VECTORS_PRESSURELESS, &
        COMOVING_DENSITY  = 0, &
        CONSERVED_DENSITY = 0
      integer ( KDI ), dimension ( 3 ) :: &
        VELOCITY         = 0, &
        MOMENTUM_DENSITY = 0
  contains
    procedure, public, pass :: &
      InitializeWithMesh
    generic :: &
      Initialize => InitializeWithMesh
    procedure, public, pass :: &
      ComputeConserved
    procedure, public, pass :: &
      ComputePrimitive
    procedure, public, pass :: &
      ComputeAuxiliary
    procedure, public, pass :: &
      ApplyBoundaryConditions
    procedure, public, pass :: &
      ComputeRawFluxes
    procedure, public, pass :: &
      ComputeRiemannSolverInput
    procedure, public, pass :: &
      SetOutputPressureless
    final :: &
      Finalize
  end type PressurelessFluidForm

    private :: &
      InitializeBasics, &
      ComputeConservedKernel, &
      ComputePrimitiveKernel, &
      ComputeEigenspeedsKernel, &
      ApplyBoundaryConditionsReflecting, &
      ComputeRawFluxesKernel, &
      ComputeRiemannSolverInputKernel
      
    interface
    
      module subroutine ComputeConservedKernel &
                   ( D, S_1, S_2, S_3, N, V_1, V_2, V_3, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          D, &
          S_1, S_2, S_3
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeConservedKernel

      module subroutine ComputePrimitiveKernel &
                   ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          N, &
          V_1, V_2, V_3, &
          D, &
          S_1, S_2, S_3
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputePrimitiveKernel

      module subroutine ComputeEigenspeedsKernel &
                   ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                     V_1, V_2, V_3, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          FEP_1, FEP_2, FEP_3, &
          FEM_1, FEM_2, FEM_3
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeEigenspeedsKernel


      module subroutine ApplyBoundaryConditionsReflecting &
                   ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, &
                     nB, oBE, oBI, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          N_E, &
          VI_E, VJ_E, VK_E
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          N_I, &
          VI_I, VJ_I, VK_I
        integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
          nB,  & 
          oBE, &
          oBI
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ApplyBoundaryConditionsReflecting


      module subroutine ComputeRawFluxesKernel &
                   ( F_D, F_S_1, F_S_2, F_S_3, D, S_1, S_2, S_3, V_Dim, &
                     UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          F_D, &
          F_S_1, F_S_2, F_S_3
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          S_1, S_2, S_3, &
          V_Dim
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeRawFluxesKernel


      module subroutine ComputeRiemannSolverInputKernel &
                   ( AP_I, AP_O, AM_I, AM_O, LP_I, LP_O, LM_I, LM_O, oV, iD, &
                     UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          AP_I, AP_O, &
          AM_I, AM_O
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          LP_I, LP_O, &
          LM_I, LM_O
        integer ( KDI ), intent ( in ) :: &
          oV, &
          iD
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeRiemannSolverInputKernel

    end interface
      
contains


  subroutine InitializeWithMesh &
               ( PF, DistributedMesh, VectorOption, VariableOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      PF
    type ( DistributedMeshForm ), intent ( in ), target :: &
      DistributedMesh
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption, &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector

    call InitializeBasics &
           ( PF, Vector, Variable, VariableUnit, VectorIndices, &
             VectorOption, VariableOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call PF % InitializeTemplate &
           ( DistributedMesh, VectorOption = Vector, &
             VariableOption = Variable, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeWithMesh

  
  subroutine ComputeConserved ( CF, Value, UseDeviceOption )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = OnDevice ( Value ( :, CF % CONSERVED_DENSITY ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    call ComputeConservedKernel &
      ( Value ( :, CF % CONSERVED_DENSITY ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
        Value ( :, CF % COMOVING_DENSITY ), &
        Value ( :, CF % VELOCITY ( 1 ) ), &    
        Value ( :, CF % VELOCITY ( 2 ) ), &   
        Value ( :, CF % VELOCITY ( 3 ) ), UseDevice )

  end subroutine ComputeConserved
  

  subroutine ComputePrimitive ( CF, Value, UseDeviceOption )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    logical ( KDL ) :: &
      UseDevice
          
    UseDevice = OnDevice ( Value ( :, CF % CONSERVED_DENSITY ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
   
    call ComputePrimitiveKernel &
           ( Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             Value ( :, CF % CONSERVED_DENSITY ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
             UseDevice )
  
  end subroutine ComputePrimitive
  
  
  subroutine ComputeAuxiliary ( CF, Value, UseDeviceOption )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    logical ( KDL ) :: &
      UseDevice
          
    UseDevice = OnDevice ( Value ( :, CF % VELOCITY ( 1 ) ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
   
    !-- No auxiliary variables besides eigenspeeds

    call ComputeEigenspeedsKernel &
           ( Value ( :, CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             UseDevice )
    
  end subroutine ComputeAuxiliary
  
  
  subroutine ApplyBoundaryConditions &
              ( CF, ExteriorValue, InteriorValue, iDimension, iBoundary, &
                PrimitiveOnlyOption, UseDeviceOption )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      ExteriorValue
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      InteriorValue
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iBoundary
    logical ( KDL ), intent ( in ), optional :: &
      PrimitiveOnlyOption, &
      UseDeviceOption
      
    integer ( KDI ) :: &
      jD, kD   !-- jDimension, kDimension
    integer ( KDI ), dimension ( 3 ) :: &
      oBI, &  !-- oBoundaryInterior
      oBE, &  !-- oBoundaryExterior
      nB      !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      N_I, &
      VI_I, VJ_I, VK_I, &
      N_E, &
      VI_E, VJ_E, VK_E      
    logical ( KDL ) :: &
      PrimitiveOnly, &
      UseDevice
          
    
    PrimitiveOnly = .false.
    if ( present ( PrimitiveOnlyOption ) ) PrimitiveOnly = PrimitiveOnlyOption
    
    UseDevice = OnDevice ( InteriorValue ( :, CF % COMOVING_DENSITY ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    associate &
      ( iD => iDimension, &
        iB => iBoundary, &
        DM => CF % DistributedMesh )
    
    if ( iB == -1 .and. DM % iaBrick ( iD ) /= 1 ) return
    if ( iB == +1 .and. DM % iaBrick ( iD ) /= DM % nBricks ( iD ) ) return
 
    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1

    call DM % SetVariablePointer &
           ( InteriorValue ( :, CF % COMOVING_DENSITY ), N_I )
    call DM % SetVariablePointer &
           ( InteriorValue ( :, CF % VELOCITY ( iD ) ), VI_I )
    call DM % SetVariablePointer &
           ( InteriorValue ( :, CF % VELOCITY ( jD ) ), VJ_I )
    call DM % SetVariablePointer &
           ( InteriorValue ( :, CF % VELOCITY ( kD ) ), VK_I )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % COMOVING_DENSITY ), N_E )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % VELOCITY ( iD ) ), VI_E )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % VELOCITY ( jD ) ), VJ_E )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % VELOCITY ( kD ) ), VK_E )

    !-- In setting oBI and oBE, note kernel routine does not inherit lbound

    oBI = DM % nGhostLayers
    if ( iB == +1 ) then
      oBI ( iD ) = oBI ( iD ) + DM % nCellsPerBrick ( iD ) - 1
    end if

    oBE = oBI
    oBE ( iD ) = oBE ( iD ) + iB

    nB ( iD ) = 1
    nB ( jD ) = DM % nCellsPerBrick ( jD )
    nB ( kD ) = DM % nCellsPerBrick ( kD )

    select case ( trim ( DM % BoundaryCondition ) ) 
    case ( 'PERIODIC' )
      return
    case ( 'REFLECTING' )
      call ApplyBoundaryConditionsReflecting &
             ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, nB, oBE, oBI, &
               UseDevice )
               
    case default
      call Show &
             ( 'This boundary condition is not implemented', CONSOLE % ERROR )
      call Show &
             ( DM % BoundaryCondition, 'Name', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select 

    if ( PrimitiveOnly ) return

    call CF % ComputeAuxiliary ( ExteriorValue, UseDeviceOption )
    call CF % ComputeConserved ( ExteriorValue, UseDeviceOption )

    end associate  !-- iD, etc.
  
  end subroutine ApplyBoundaryConditions


  subroutine ComputeRawFluxes &
               ( CF, RawFlux, Value, iDimension, UseDeviceOption )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    integer ( KDI ), intent ( in ) :: &
      iDimension
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
    
    integer ( KDI ) :: &
      iDensity
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = OnDevice ( Value ( :, CF % CONSERVED_DENSITY ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    call Search &
           ( CF % iaConserved, CF % CONSERVED_DENSITY, iDensity )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( 1 ), iMomentum ( 1 ) )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( 2 ), iMomentum ( 2 ) )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( 3 ), iMomentum ( 3 ) )
    
    call ComputeRawFluxesKernel &
           ( RawFlux ( :, iDensity ), &
             RawFlux ( :, iMomentum ( 1 ) ), &
             RawFlux ( :, iMomentum ( 2 ) ), &
             RawFlux ( :, iMomentum ( 3 ) ), &
             Value ( :, CF % CONSERVED_DENSITY ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
             Value ( :, CF % VELOCITY ( iDimension ) ), &
             UseDevice )
               
  end subroutine ComputeRawFluxes
  
  
  subroutine ComputeRiemannSolverInput &
               ( CF, Step, ValueInner, ValueOuter, iDimension, &
                 UseDeviceOption )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    class ( * ), intent ( inout ) :: &
      Step
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      ValueInner, &
      ValueOuter
    integer ( KDI ), intent ( in ) :: &
      iDimension
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
    
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      AP_I, AP_O, &
      AM_I, AM_O, &
      LP_I, LP_O, &
      LM_I, LM_O
    logical ( KDL ) :: &
      UseDevice
      
    select type ( S => Step )
    type is ( ConservationLawStepForm )

    call CF % DistributedMesh % SetVariablePointer &
           ( S % ModifiedSpeedsInner % Value ( :, S % ALPHA_PLUS ), AP_I )
    call CF % DistributedMesh % SetVariablePointer &
           ( S % ModifiedSpeedsOuter % Value ( :, S % ALPHA_PLUS ), AP_O )
    call CF % DistributedMesh % SetVariablePointer &
           ( S % ModifiedSpeedsInner % Value ( :, S % ALPHA_MINUS ), AM_I )
    call CF % DistributedMesh % SetVariablePointer &
           ( S % ModifiedSpeedsOuter % Value ( :, S % ALPHA_MINUS ), AM_O )
    call CF % DistributedMesh % SetVariablePointer &
           ( ValueInner ( :, CF % FAST_EIGENSPEED_PLUS ( iDimension ) ), LP_I )
    call CF % DistributedMesh % SetVariablePointer &
           ( ValueOuter ( :, CF % FAST_EIGENSPEED_PLUS ( iDimension ) ), LP_O )
    call CF % DistributedMesh % SetVariablePointer &
           ( ValueInner ( :, CF % FAST_EIGENSPEED_MINUS ( iDimension ) ), LM_I )
    call CF % DistributedMesh % SetVariablePointer &
           ( ValueOuter ( :, CF % FAST_EIGENSPEED_MINUS ( iDimension ) ), LM_O )
           
    UseDevice = OnDevice ( AP_I )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    call ComputeRiemannSolverInputKernel &
           ( AP_I, AP_O, AM_I, AM_O, LP_I, LP_O, LM_I, LM_O, &
             CF % DistributedMesh % nGhostLayers ( iDimension ), iDimension, &
             UseDevice )
             
    end select !-- S
      
  end subroutine ComputeRiemannSolverInput


  subroutine SetOutputPressureless &
               ( PF, UnitsOnlyOption, VelocityUnitOption, DensityUnitOption )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      PF
    logical ( KDL ), intent ( in ), optional :: &
      UnitsOnlyOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      VelocityUnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      DensityUnitOption

    integer ( KDI ) :: &
      iD  !-- iDimension
    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices
    type ( MeasuredValueForm ) :: &
      DensityUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit
    logical ( KDL ) :: &
      UnitsOnly

    DensityUnit = UNIT % IDENTITY
    if ( present ( DensityUnitOption ) ) DensityUnit = DensityUnitOption

    VelocityUnit = UNIT % IDENTITY
    if ( present ( VelocityUnitOption ) ) VelocityUnit = VelocityUnitOption

    UnitsOnly = .false.
    if ( present ( UnitsOnlyOption ) ) UnitsOnly = UnitsOnlyOption

    PF % Unit ( PF % COMOVING_DENSITY ) = DensityUnit
    PF % Unit ( PF % CONSERVED_DENSITY ) = DensityUnit
    do iD = 1, 3
      PF % Unit ( PF % VELOCITY ( iD ) ) = VelocityUnit ( iD )
      PF % Unit ( PF % MOMENTUM_DENSITY ( iD ) ) &
        = DensityUnit * VelocityUnit ( iD )
      PF % Unit ( PF % FAST_EIGENSPEED_PLUS ( iD ) ) = VelocityUnit ( iD )
      PF % Unit ( PF % FAST_EIGENSPEED_MINUS ( iD ) ) = VelocityUnit ( iD )
    end do

    if ( UnitsOnly ) return

    call VectorIndices ( 1 ) % Initialize ( PF % VELOCITY )

    call PF % Output ( 1 ) % Initialize &
           ( PF, iaSelectedOption = [ PF % COMOVING_DENSITY, PF % VELOCITY ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

    call PF % DistributedMesh % SetImage ( PF % Output, PROGRAM_HEADER % Name )

  end subroutine SetOutputPressureless


  subroutine Finalize ( PF )

    type ( PressurelessFluidForm ), intent ( inout ) :: &
      PF

    nullify ( PF % DistributedMesh )
    
    if ( allocated ( PF % iaConserved ) )  deallocate ( PF % iaConserved )
    if ( allocated ( PF % iaPrimitive ) ) deallocate ( PF % iaPrimitive )
 
    call Show ( 'Finalizing ' // trim ( PF % Type ), CONSOLE % INFO_3 )
   
   end subroutine Finalize


  subroutine InitializeBasics &
               ( PF, Vector, Variable, VariableUnit, VectorIndices, &
                 VectorOption, VariableOption, NameOption, VariableUnitOption, &
                 VectorIndicesOption )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      PF
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Vector, &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    !-- FIXME: intent(out) here caused ICE with Intel Compiler 15
    !          Temporarily set to intent(inout)
    !type ( Integer_1D_Form ), dimension ( : ), allocatable, &
    !  intent ( out ) :: &
    type ( Integer_1D_Form ), dimension ( : ), allocatable, &
      intent ( inout ) :: &
        VectorIndices
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption, &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      oF, &  !-- oField
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDF ) :: &
      Name

    if ( PF % Type == '' ) PF % Type = 'a PressurelessFluid'

    Name = ''
    if ( present ( NameOption ) ) Name = NameOption

    call Show ( 'Initializing ' // trim ( PF % Type ), CONSOLE % INFO_3 )
    call Show ( Name, 'Name', CONSOLE % INFO_3 )

    !-- variable indices

    oF = PF % N_FIELDS_TEMPLATE
    if ( PF % N_FIELDS == 0 ) PF % N_FIELDS = oF + PF % N_FIELDS_PRESSURELESS

    PF % COMOVING_DENSITY      = oF + 1
    PF % CONSERVED_DENSITY     = oF + 2
    PF % VELOCITY              = oF + [ 3, 4, 5 ]
    PF % MOMENTUM_DENSITY      = oF + [ 6, 7, 8 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( PF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + PF % N_FIELDS_PRESSURELESS ) &
      = [ 'ComovingDensity                ', &
          'ConservedDensity               ', &
          'Velocity_1                     ', &
          'Velocity_2                     ', &
          'Velocity_3                     ', &
          'MomentumDensity_1              ', &
          'MomentumDensity_2              ', &
          'MomentumDensity_3              ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( PF % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = PF % N_VECTORS_TEMPLATE
    if ( PF % N_VECTORS == 0 ) PF % N_VECTORS = oV + PF % N_VECTORS_PRESSURELESS

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( PF % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + PF % N_VECTORS_PRESSURELESS ) &
      = [ 'Velocity                       ', &
          'MomentumDensity                ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( PF % N_VECTORS ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( PF % VELOCITY )
    call VectorIndices ( oV + 2 ) % Initialize ( PF % MOMENTUM_DENSITY )

    !-- select primitive, conserved

    oP = PF % N_PRIMITIVE_TEMPLATE
    oC = PF % N_CONSERVED_TEMPLATE

    if ( .not. allocated ( PF % iaPrimitive ) ) then
      PF % N_PRIMITIVE = oP + PF % N_PRIMITIVE_PRESSURELESS
      allocate ( PF % iaPrimitive ( PF % N_PRIMITIVE ) )
    end if
    PF % iaPrimitive ( oP + 1 : oP + PF % N_PRIMITIVE_PRESSURELESS ) &
      = [ PF % COMOVING_DENSITY, PF % VELOCITY ]

    if ( .not. allocated ( PF % iaConserved ) ) then
      PF % N_CONSERVED = oC + PF % N_CONSERVED_PRESSURELESS
      allocate ( PF % iaConserved ( PF % N_CONSERVED ) )
    end if
    PF % iaConserved ( oC + 1 : oC + PF % N_CONSERVED_PRESSURELESS ) &
      = [ PF % CONSERVED_DENSITY, PF % MOMENTUM_DENSITY ]
    
  end subroutine InitializeBasics


end module PressurelessFluid_Form
