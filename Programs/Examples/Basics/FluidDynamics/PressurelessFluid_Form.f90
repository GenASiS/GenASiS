#include "Preprocessor"

module PressurelessFluid_Form

  use iso_c_binding
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
      ComputeConservedHost
    procedure, public, pass :: &
      ComputeConservedDevice
    procedure, public, pass :: &
      ComputePrimitiveHost
    procedure, public, pass :: &
      ComputePrimitiveDevice
    procedure, public, pass :: &
      ComputeAuxiliaryHost
    procedure, public, pass :: &
      ComputeAuxiliaryDevice
    procedure, public, pass :: &
      ApplyBoundaryConditionsHost
    procedure, public, pass :: &
      ApplyBoundaryConditionsDevice
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
      ComputeConservedKernelDevice, &
      ComputePrimitiveKernel, &
      ComputeEigenspeedsKernel, &
      ComputeEigenspeedsKernelDevice, &
      ApplyBoundaryConditionsReflecting, &
      ApplyBoundaryConditionsReflectingDevice, &
      ComputeRawFluxesKernel, &
      ComputeRiemannSolverInputKernel

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

  
  subroutine ComputeConservedHost ( CF, Value )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    
    call ComputeConservedKernel &
      ( Value ( :, CF % CONSERVED_DENSITY ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
        Value ( :, CF % COMOVING_DENSITY ), &
        Value ( :, CF % VELOCITY ( 1 ) ), &    
        Value ( :, CF % VELOCITY ( 2 ) ), &   
        Value ( :, CF % VELOCITY ( 3 ) ) )    

  end subroutine ComputeConservedHost
  

  subroutine ComputeConservedDevice ( CF, Value, D_Value )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    type  ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_Value
    
    call ComputeConservedKernelDevice &
      ( Value ( :, CF % CONSERVED_DENSITY ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
        Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
        Value ( :, CF % COMOVING_DENSITY ), &
        Value ( :, CF % VELOCITY ( 1 ) ), &    
        Value ( :, CF % VELOCITY ( 2 ) ), &   
        Value ( :, CF % VELOCITY ( 3 ) ), &
        D_Value ( CF % CONSERVED_DENSITY ), &
        D_Value ( CF % MOMENTUM_DENSITY ( 1 ) ), &
        D_Value ( CF % MOMENTUM_DENSITY ( 2 ) ), &
        D_Value ( CF % MOMENTUM_DENSITY ( 3 ) ), &
        D_Value ( CF % COMOVING_DENSITY ), &
        D_Value ( CF % VELOCITY ( 1 ) ), &    
        D_Value ( CF % VELOCITY ( 2 ) ), &   
        D_Value ( CF % VELOCITY ( 3 ) ) )

  end subroutine ComputeConservedDevice
  

  subroutine ComputePrimitiveHost ( CF, Value )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
   
    call ComputePrimitiveKernel &
           ( Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             Value ( :, CF % CONSERVED_DENSITY ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ) )
  
  end subroutine ComputePrimitiveHost
  
  
  subroutine ComputePrimitiveDevice ( CF, Value, D_Value )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    type  ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_Value
   
    call ComputePrimitiveKernelDevice &
           ( Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             Value ( :, CF % CONSERVED_DENSITY ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
             D_Value ( CF % COMOVING_DENSITY ), &
             D_Value ( CF % VELOCITY ( 1 ) ), &
             D_Value ( CF % VELOCITY ( 2 ) ), &
             D_Value ( CF % VELOCITY ( 3 ) ), &
             D_Value ( CF % CONSERVED_DENSITY ), &
             D_Value ( CF % MOMENTUM_DENSITY ( 1 ) ), &
             D_Value ( CF % MOMENTUM_DENSITY ( 2 ) ), &
             D_Value ( CF % MOMENTUM_DENSITY ( 3 ) ) )
  
  end subroutine ComputePrimitiveDevice
  
  
  subroutine ComputeAuxiliaryHost ( CF, Value )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    
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
             Value ( :, CF % VELOCITY ( 3 ) ) )
    
  end subroutine ComputeAuxiliaryHost
  
  
  subroutine ComputeAuxiliaryDevice ( CF, Value, D_Value )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_Value
    
    !-- No auxiliary variables besides eigenspeeds

    call ComputeEigenspeedsKernelDevice &
           ( Value ( :, CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             D_Value ( CF % VELOCITY ( 1 ) ), &
             D_Value ( CF % VELOCITY ( 2 ) ), &
             D_Value ( CF % VELOCITY ( 3 ) ) )
    
  end subroutine ComputeAuxiliaryDevice
  
  
  subroutine ApplyBoundaryConditionsHost &
              ( CF, ExteriorValue, InteriorValue, iDimension, iBoundary, &
                PrimitiveOnlyOption )

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
      PrimitiveOnlyOption
      
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
      PrimitiveOnly

    PrimitiveOnly = .false.
    if ( present ( PrimitiveOnlyOption ) ) PrimitiveOnly = PrimitiveOnlyOption

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
             ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, nB, oBE, oBI )
    case default
      call Show &
             ( 'This boundary condition is not implemented', CONSOLE % ERROR )
      call Show &
             ( DM % BoundaryCondition, 'Name', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select 

    if ( PrimitiveOnly ) return

    call CF % ComputeAuxiliary ( ExteriorValue )
    call CF % ComputeConserved ( ExteriorValue )

    end associate  !-- iD, etc.
  
  end subroutine ApplyBoundaryConditionsHost


  subroutine ApplyBoundaryConditionsDevice &
              ( CF, ExteriorValue, InteriorValue, iDimension, iBoundary, &
                D_ExteriorValue, D_InteriorValue, PrimitiveOnlyOption )

    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      ExteriorValue
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      InteriorValue
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iBoundary
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_ExteriorValue, &
      D_InteriorValue
    
    logical ( KDL ), intent ( in ), optional :: &
      PrimitiveOnlyOption
      
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
      PrimitiveOnly

    PrimitiveOnly = .false.
    if ( present ( PrimitiveOnlyOption ) ) PrimitiveOnly = PrimitiveOnlyOption

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
      call ApplyBoundaryConditionsReflectingDevice &
             ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, nB, oBE, oBI, &
               D_ExteriorValue ( CF % COMOVING_DENSITY ), &
               D_ExteriorValue ( CF % VELOCITY ( iD ) ), &
               D_ExteriorValue ( CF % VELOCITY ( jD ) ), &
               D_ExteriorValue ( CF % VELOCITY ( kD ) ), &
               D_InteriorValue ( CF % COMOVING_DENSITY ), &
               D_InteriorValue ( CF % VELOCITY ( iD ) ), &
               D_InteriorValue ( CF % VELOCITY ( jD ) ), &
               D_InteriorValue ( CF % VELOCITY ( kD ) ) )
    case default
      call Show &
             ( 'This boundary condition is not implemented', CONSOLE % ERROR )
      call Show &
             ( DM % BoundaryCondition, 'Name', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select 

    if ( PrimitiveOnly ) return

    call CF % ComputeAuxiliary ( ExteriorValue, D_ExteriorValue )
    call CF % ComputeConserved ( ExteriorValue, D_ExteriorValue )

    end associate  !-- iD, etc.
  
  end subroutine ApplyBoundaryConditionsDevice


  subroutine ComputeRawFluxes &
               ( CF, RawFlux, Value, iDimension, D_RawFlux, D_Value )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    integer ( KDI ), intent ( in ) :: &
      iDimension
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_RawFlux, &
      D_Value
    
    integer ( KDI ) :: &
      iDensity
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum
      
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
             Value ( :, CF % VELOCITY ( iDimension ) ) )
               
  end subroutine ComputeRawFluxes
  
  
  subroutine ComputeRiemannSolverInput &
               ( CF, Step, ValueInner, ValueOuter, iDimension, &
                 D_ValueInner, D_ValueOuter )
    
    class ( PressurelessFluidForm ), intent ( inout ) :: &
      CF
    class ( * ), intent ( inout ) :: &
      Step
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      ValueInner, &
      ValueOuter
    integer ( KDI ), intent ( in ) :: &
      iDimension
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_ValueInner, &
      D_ValueOuter
    
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      AP_I, AP_O, &
      AM_I, AM_O, &
      LP_I, LP_O, &
      LM_I, LM_O

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

    call ComputeRiemannSolverInputKernel &
           ( AP_I, AP_O, AM_I, AM_O, LP_I, LP_O, LM_I, LM_O, &
             CF % DistributedMesh % nGhostLayers ( iDimension ), iDimension, &
             S % ModifiedSpeedsInner % D_Selected ( S % ALPHA_PLUS ), &
             S % ModifiedSpeedsOuter % D_Selected ( S % ALPHA_PLUS ), &
             S % ModifiedSpeedsInner % D_Selected ( S % ALPHA_MINUS ), &
             S % ModifiedSpeedsOuter % D_Selected ( S % ALPHA_MINUS ), &
             D_ValueInner ( CF % FAST_EIGENSPEED_PLUS ( iDimension ) ), &
             D_ValueOuter ( CF % FAST_EIGENSPEED_PLUS ( iDimension ) ), &
             D_ValueInner ( CF % FAST_EIGENSPEED_MINUS ( iDimension ) ), &
             D_ValueOuter ( CF % FAST_EIGENSPEED_MINUS ( iDimension ) ) )
             
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


  subroutine ComputeConservedKernel ( D, S_1, S_2, S_3, N, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      D, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3

    D   = N
    S_1 = N * V_1
    S_2 = N * V_2
    S_3 = N * V_3

  end subroutine ComputeConservedKernel


  subroutine ComputeConservedKernelDevice &
               ( D, S_1, S_2, S_3, N, V_1, V_2, V_3, &
                 D_D, D_S_1, D_S_2, D_S_3, D_N, D_V_1, D_V_2, D_V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      D, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3
    type ( c_ptr ), intent ( in ) :: &
      D_D, &
      D_S_1, D_S_2, D_S_3, &
      D_N, &
      D_V_1, D_V_2, D_V_3
    
    integer ( KDI ) :: &
      iV
    
    call AssociateHost ( D_D, D )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )

    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( D )
      D   ( iV ) = N ( iV )
      S_1 ( iV ) = N ( iV ) * V_1 ( iV )
      S_2 ( iV ) = N ( iV ) * V_2 ( iV )
      S_3 ( iV ) = N ( iV ) * V_3 ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( D )

  end subroutine ComputeConservedKernelDevice


  subroutine ComputePrimitiveKernel ( N, V_1, V_2, V_3, D, S_1, S_2, S_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      V_1, V_2, V_3, &
      D, &
      S_1, S_2, S_3

    N = D

    where ( N > 0.0_KDR )
      V_1 = S_1 / N 
      V_2 = S_2 / N
      V_3 = S_3 / N
    elsewhere
      N   = 0.0_KDR
      V_1 = 0.0_KDR
      V_2 = 0.0_KDR
      V_3 = 0.0_KDR
      D   = 0.0_KDR
      S_1 = 0.0_KDR
      S_2 = 0.0_KDR
      S_3 = 0.0_KDR
    end where

  end subroutine ComputePrimitiveKernel


  subroutine ComputePrimitiveKernelDevice &
               ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, &
                 D_N, D_V_1, D_V_2, D_V_3, D_D, D_S_1, D_S_2, D_S_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      V_1, V_2, V_3, &
      D, &
      S_1, S_2, S_3
    type ( c_ptr ), intent ( in ) :: &
      D_N, &
      D_V_1, D_V_2, D_V_3, &
      D_D, &
      D_S_1, D_S_2, D_S_3
      
    integer ( KDI ) :: &
      iV


    call Copy ( D, D_D, D_N, N )
    
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_D, D )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( N )
      if ( N ( iV ) > 0.0_KDR ) then
        V_1 ( iV ) = S_1 ( iV ) / N ( iV )
        V_2 ( iV ) = S_2 ( iV ) / N ( iV )
        V_3 ( iV ) = S_3 ( iV ) / N ( iV )
      else
        N   ( iV )= 0.0_KDR
        V_1 ( iV ) = 0.0_KDR
        V_2 ( iV ) = 0.0_KDR
        V_3 ( iV ) = 0.0_KDR
        D   ( iV ) = 0.0_KDR
        S_1 ( iV ) = 0.0_KDR
        S_2 ( iV ) = 0.0_KDR
        S_3 ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( D )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    
  end subroutine ComputePrimitiveKernelDevice

  
  subroutine ComputeEigenspeedsKernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V_1, V_2, V_3

    FEP_1 = V_1
    FEP_2 = V_2
    FEP_3 = V_3
    FEM_1 = V_1
    FEM_2 = V_2
    FEM_3 = V_3

  end subroutine ComputeEigenspeedsKernel


  subroutine ComputeEigenspeedsKernelDevice &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3, &
                 D_FEP_1, D_FEP_2, D_FEP_3, D_FEM_1, D_FEM_2, D_FEM_3, &
                 D_V_1, D_V_2, D_V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V_1, V_2, V_3
    type ( c_ptr ), intent ( in ) :: &
      D_FEP_1, D_FEP_2, D_FEP_3, &
      D_FEM_1, D_FEM_2, D_FEM_3, &
      D_V_1, D_V_2, D_V_3
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_FEP_1, FEP_1 )
    call AssociateHost ( D_FEP_2, FEP_2 )
    call AssociateHost ( D_FEP_3, FEP_3 )
    call AssociateHost ( D_FEM_1, FEM_1 )
    call AssociateHost ( D_FEM_2, FEM_2 )
    call AssociateHost ( D_FEM_3, FEM_3 )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3,  V_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( FEP_1 )
      FEP_1 ( iV ) = V_1 ( iV )
      FEP_2 ( iV ) = V_2 ( iV )
      FEP_3 ( iV ) = V_3 ( iV )
      FEM_1 ( iV ) = V_1 ( iV )
      FEM_2 ( iV ) = V_2 ( iV )
      FEM_3 ( iV ) = V_3 ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( FEM_3 )
    call DisassociateHost ( FEM_2 )
    call DisassociateHost ( FEM_1 )
    call DisassociateHost ( FEP_3 )
    call DisassociateHost ( FEP_2 )
    call DisassociateHost ( FEP_1 )

  end subroutine ComputeEigenspeedsKernelDevice


  subroutine ApplyBoundaryConditionsReflecting &
               ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, nB, oBE, oBI )

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

    N_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
          oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
          oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = N_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
              oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
              oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    VI_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
           oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
           oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = - VI_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
                 oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
                 oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    VJ_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
           oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
           oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = VJ_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
               oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
               oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    VK_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
           oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
           oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = VK_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
               oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
               oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

  end subroutine ApplyBoundaryConditionsReflecting


  subroutine ApplyBoundaryConditionsReflectingDevice &
               ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, nB, oBE, oBI, &
                 D_N_E, D_VI_E, D_VJ_E, D_VK_E, D_N_I, D_VI_I, D_VJ_I, &
                 D_VK_I )

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
    type ( c_ptr ), intent ( in ) :: &
      D_N_E, &
      D_VI_E, D_VJ_E, D_VK_E, &
      D_N_I, &
      D_VI_I, D_VJ_I, D_VK_I
      
    integer ( KDI ) :: &
      iV, jV, kV
      
    call AssociateHost ( D_N_E, N_E )
    call AssociateHost ( D_VI_E, VI_E )
    call AssociateHost ( D_VJ_E, VJ_E )
    call AssociateHost ( D_VK_E, VK_E )
    call AssociateHost ( D_N_I, N_I )
    call AssociateHost ( D_VI_I, VI_I )
    call AssociateHost ( D_VJ_I, VJ_I )
    call AssociateHost ( D_VK_I, VK_I )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE )
    do kV = 1, nB ( 3 )
      do jV = 1, nB ( 2 )
        do iV = 1, nB ( 1 )

          N_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = N_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

          VI_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = - VI_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

          VJ_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = VJ_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

          VK_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = VK_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )

        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
               
    call DisassociateHost ( VK_I )
    call DisassociateHost ( VJ_I )
    call DisassociateHost ( VI_I )
    call DisassociateHost ( N_I )
    call DisassociateHost ( VK_E )
    call DisassociateHost ( VJ_E )
    call DisassociateHost ( VI_E )
    call DisassociateHost ( N_E )
    
  end subroutine ApplyBoundaryConditionsReflectingDevice


  subroutine ComputeRawFluxesKernel &
               ( F_D, F_S_1, F_S_2, F_S_3, D, S_1, S_2, S_3, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_D, &
      F_S_1, F_S_2, F_S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      D, &
      S_1, S_2, S_3, &
      V_Dim

    F_D   = D   * V_Dim
    F_S_1 = S_1 * V_Dim
    F_S_2 = S_2 * V_Dim
    F_S_3 = S_3 * V_Dim

  end subroutine ComputeRawFluxesKernel


  subroutine ComputeRiemannSolverInputKernel &
               ( AP_I, AP_O, AM_I, AM_O, LP_I, LP_O, LM_I, LM_O, oV, iD, &
                 D_AP_I, D_AP_O, D_AM_I, D_AM_O, D_LP_I, D_LP_O, &
                 D_LM_I, D_LM_O )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      AP_I, AP_O, &
      AM_I, AM_O
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      LP_I, LP_O, &
      LM_I, LM_O
    integer ( KDI ), intent ( in ) :: &
      oV, &
      iD
    type ( c_ptr ), intent ( in ) :: &
      D_AP_I, D_AP_O, &
      D_AM_I, D_AM_O, &
      D_LP_I, D_LP_O, &
      D_LM_I, D_LM_O
      
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS_P, iaS_M, &
      iaVS_P, iaVS_M, &
      lV, uV
      
    !AP_I = max ( 0.0_KDR, + cshift ( LP_O, shift = -1, dim = iD ), + LP_I )
    !AP_O = max ( 0.0_KDR, + LP_O, + cshift ( LP_I, shift = +1, dim = iD ) )
    !AM_I = max ( 0.0_KDR, - cshift ( LM_O, shift = -1, dim = iD ), - LM_I )
    !AM_O = max ( 0.0_KDR, - LM_O, - cshift ( LM_I, shift = +1, dim = iD ) )
    
    call AssociateHost ( D_AP_I, AP_I )
    call AssociateHost ( D_AP_O, AP_O )
    call AssociateHost ( D_AM_I, AM_I )
    call AssociateHost ( D_AM_O, AM_O )
    call AssociateHost ( D_LP_I, LP_I )
    call AssociateHost ( D_LP_O, LP_O )
    call AssociateHost ( D_LM_I, LM_I )
    call AssociateHost ( D_LM_O, LM_O )

    lV = 1
    where ( shape ( LP_O ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV
    
    uV = 1
    where ( shape ( LP_O ) > 1 )
      uV = shape ( LP_O ) - oV
    end where
    uV ( iD ) = size ( LP_O, dim = iD ) - 1
    
    iaS_M = 0
    iaS_M ( iD ) = - 1
    iaS_P = 0
    iaS_P ( iD ) = + 1
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( iaVS_M, iaVS_P )
    do kV = lV ( 3 ), uV ( 3 )    
      do jV = lV ( 2 ), uV ( 2 )  
        do iV = lV ( 1 ), uV ( 1 )
        
          iaVS_M = [ iV, jV, kV ] + iaS_M
          iaVS_P = [ iV, jV, kV ] + iaS_P
          
          AP_I ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    + LP_O ( iaVS_M ( 1 ), iaVS_M ( 2 ), iaVS_M ( 3 ) ), &
                    + LP_I ( iV, jV, kV ) )
          AP_O ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    + LP_O ( iV, jV, kV ), &
                    + LP_I ( iaVS_P ( 1 ), iaVS_P ( 2 ), iaVS_P ( 3 ) ) )
          AM_I ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    - LM_O ( iaVS_M ( 1 ), iaVS_M ( 2 ), iaVS_M ( 3 ) ), &
                    - LM_I ( iV, jV, kV ) )
          AM_O ( iV, jV, KV ) &
            = max ( 0.0_KDR, &
                    - LM_O ( iV, jV, kV ), &
                    - LM_I ( iaVS_P ( 1 ), iaVS_P ( 2 ), iaVS_P ( 3 ) ) )
        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( LM_O )
    call DisassociateHost ( LM_I )
    call DisassociateHost ( LP_O )
    call DisassociateHost ( LP_I )
    call DisassociateHost ( AM_O )
    call DisassociateHost ( AM_I )
    call DisassociateHost ( AP_O )
    call DisassociateHost ( AP_I )
    
  end subroutine ComputeRiemannSolverInputKernel


end module PressurelessFluid_Form
