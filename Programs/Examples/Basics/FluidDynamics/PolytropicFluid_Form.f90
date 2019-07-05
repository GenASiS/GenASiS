module PolytropicFluid_Form

  use iso_c_binding
  use Basics
  use DistributedMesh_Form
  use PressurelessFluid_Form
  use PolytropicFluid_Kernel

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_POLYTROPIC = 2, &
      N_CONSERVED_POLYTROPIC = 1, &
      N_FIELDS_POLYTROPIC    = 6, &
      N_VECTORS_POLYTROPIC   = 0

  type, public, extends ( PressurelessFluidForm ) :: PolytropicFluidForm
    integer ( KDI ) :: &
      N_PRIMITIVE_POLYTROPIC = N_PRIMITIVE_POLYTROPIC, &
      N_CONSERVED_POLYTROPIC = N_CONSERVED_POLYTROPIC, &
      N_FIELDS_POLYTROPIC    = N_FIELDS_POLYTROPIC, &
      N_VECTORS_POLYTROPIC   = N_VECTORS_POLYTROPIC, &
      INTERNAL_ENERGY      = 0, &
      CONSERVED_ENERGY     = 0, &
      PRESSURE             = 0, &
      ADIABATIC_INDEX      = 0, &
      SOUND_SPEED          = 0, &
      POLYTROPIC_PARAMETER = 0
  contains
    procedure, public, pass :: &
      InitializeWithMesh
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
      ComputeAuxiliaryFromPressure
    procedure, public, pass :: &
      ApplyBoundaryConditionsHost
    procedure, public, pass :: &
      ApplyBoundaryConditionsDevice
    procedure, public, pass :: &
      ComputeRawFluxes
    procedure, public, pass :: &
      SetOutputPolytropic
    final :: &
      Finalize
  end type PolytropicFluidForm

    private :: &
      InitializeBasics
      
contains


  subroutine InitializeWithMesh &
               ( PF, DistributedMesh, VectorOption, VariableOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
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

    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit

    call InitializeBasics &
           ( PF, Variable, VariableUnit, VariableOption, UnitOption )

    call PF % PressurelessFluidForm % InitializeWithMesh &
           ( DistributedMesh, VectorOption = VectorOption, &
             VariableOption = Variable, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeWithMesh


  subroutine ComputeConservedHost ( CF, Value )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    
    call CF % PressurelessFluidForm % ComputeConserved ( Value )

    call ComputeConservedKernel &
      ( Value ( :, CF % CONSERVED_ENERGY ), &
        Value ( :, CF % INTERNAL_ENERGY ), &    
        Value ( :, CF % COMOVING_DENSITY ), &
        Value ( :, CF % VELOCITY ( 1 ) ), &    
        Value ( :, CF % VELOCITY ( 2 ) ), &   
        Value ( :, CF % VELOCITY ( 3 ) ) )    

  end subroutine ComputeConservedHost


  subroutine ComputeConservedDevice ( CF, Value, D_Value )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_Value
    
    call CF % PressurelessFluidForm % ComputeConserved ( Value, D_Value )

    call ComputeConservedKernelDevice &
      ( Value ( :, CF % CONSERVED_ENERGY ), &
        Value ( :, CF % INTERNAL_ENERGY ), &    
        Value ( :, CF % COMOVING_DENSITY ), &
        Value ( :, CF % VELOCITY ( 1 ) ), &    
        Value ( :, CF % VELOCITY ( 2 ) ), &   
        Value ( :, CF % VELOCITY ( 3 ) ), &
        D_Value ( CF % CONSERVED_ENERGY ), &
        D_Value ( CF % INTERNAL_ENERGY ), &    
        D_Value ( CF % COMOVING_DENSITY ), &
        D_Value ( CF % VELOCITY ( 1 ) ), &    
        D_Value ( CF % VELOCITY ( 2 ) ), &   
        D_Value ( CF % VELOCITY ( 3 ) ) )

  end subroutine ComputeConservedDevice


  subroutine ComputePrimitiveHost ( CF, Value )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    
    call CF % PressurelessFluidForm % ComputePrimitive ( Value )

    call ComputePrimitiveKernel &
           ( Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % CONSERVED_ENERGY ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ) )
  
  end subroutine ComputePrimitiveHost
  
  
  subroutine ComputePrimitiveDevice ( CF, Value, D_Value )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_Value
    
    call CF % PressurelessFluidForm % ComputePrimitive ( Value, D_Value )

    call ComputePrimitiveKernelDevice &
           ( Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % CONSERVED_ENERGY ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             D_Value ( CF % INTERNAL_ENERGY ), &
             D_Value ( CF % CONSERVED_ENERGY ), &
             D_Value ( CF % COMOVING_DENSITY ), &
             D_Value ( CF % VELOCITY ( 1 ) ), &
             D_Value ( CF % VELOCITY ( 2 ) ), &
             D_Value ( CF % VELOCITY ( 3 ) ) )
  
  end subroutine ComputePrimitiveDevice
  
  
  subroutine ComputeAuxiliaryHost ( CF, Value )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    
    call ComputeAuxiliaryKernel &
           ( Value ( :, CF % PRESSURE ), &
             Value ( :, CF % POLYTROPIC_PARAMETER ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % ADIABATIC_INDEX ) )
  
    call ComputeEigenspeedsKernel &
           ( Value ( :, CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             Value ( :, CF % SOUND_SPEED ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % ADIABATIC_INDEX ) )
    
  end subroutine ComputeAuxiliaryHost
  
  
  subroutine ComputeAuxiliaryDevice ( CF, Value, D_Value )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    type ( c_ptr ), dimension ( : ), intent ( in ) :: &
      D_Value
    
    call ComputeAuxiliaryKernelDevice &
           ( Value ( :, CF % PRESSURE ), &
             Value ( :, CF % POLYTROPIC_PARAMETER ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % ADIABATIC_INDEX ), &
             D_Value ( CF % PRESSURE ), &
             D_Value ( CF % POLYTROPIC_PARAMETER ), &
             D_Value ( CF % COMOVING_DENSITY ), &
             D_Value ( CF % INTERNAL_ENERGY ), &
             D_Value ( CF % ADIABATIC_INDEX ) )
  
    call ComputeEigenspeedsKernelDevice &
           ( Value ( :, CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             Value ( :, CF % SOUND_SPEED ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % ADIABATIC_INDEX ), &
             D_Value ( CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             D_Value ( CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             D_Value ( CF % SOUND_SPEED ), &
             D_Value ( CF % COMOVING_DENSITY ), &
             D_Value ( CF % VELOCITY ( 1 ) ), &
             D_Value ( CF % VELOCITY ( 2 ) ), &
             D_Value ( CF % VELOCITY ( 3 ) ), &
             D_Value ( CF % PRESSURE ), &
             D_Value ( CF % ADIABATIC_INDEX ) )
    
  end subroutine ComputeAuxiliaryDevice
  
  
  subroutine ComputeAuxiliaryFromPressure ( CF, Value )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    
    call ComputeAuxiliaryFromPressureKernel &
           ( Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % POLYTROPIC_PARAMETER ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % ADIABATIC_INDEX ) )
  
    call ComputeEigenspeedsKernel &
           ( Value ( :, CF % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             Value ( :, CF % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             Value ( :, CF % SOUND_SPEED ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % ADIABATIC_INDEX ) )
    
  end subroutine ComputeAuxiliaryFromPressure
  
  
  subroutine ApplyBoundaryConditionsHost &
              ( CF, ExteriorValue, InteriorValue, iDimension, iBoundary, &
                PrimitiveOnlyOption )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
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
      E_I, &
      Gamma_I, &
      E_E, &
      Gamma_E
    logical ( KDL ) :: &
      PrimitiveOnly

    call CF % PressurelessFluidForm % ApplyBoundaryConditions &
           ( ExteriorValue, InteriorValue, iDimension, iBoundary, &
             PrimitiveOnlyOption = .true. )

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
           ( InteriorValue ( :, CF % INTERNAL_ENERGY ), E_I )
    call DM % SetVariablePointer &
           ( InteriorValue ( :, CF % ADIABATIC_INDEX ), Gamma_I )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % INTERNAL_ENERGY ), E_E )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % ADIABATIC_INDEX ), Gamma_E )

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
             ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI )
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

    class ( PolytropicFluidForm ), intent ( inout ) :: &
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
      E_I, &
      Gamma_I, &
      E_E, &
      Gamma_E
    logical ( KDL ) :: &
      PrimitiveOnly

    call CF % PressurelessFluidForm % ApplyBoundaryConditions &
           ( ExteriorValue, InteriorValue, iDimension, iBoundary, &
             D_ExteriorValue, D_InteriorValue, PrimitiveOnlyOption = .true. )

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
           ( InteriorValue ( :, CF % INTERNAL_ENERGY ), E_I )
    call DM % SetVariablePointer &
           ( InteriorValue ( :, CF % ADIABATIC_INDEX ), Gamma_I )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % INTERNAL_ENERGY ), E_E )
    call DM % SetVariablePointer &
           ( ExteriorValue ( :, CF % ADIABATIC_INDEX ), Gamma_E )

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
             ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI, &
               D_ExteriorValue ( CF % INTERNAL_ENERGY ), &
               D_ExteriorValue ( CF % ADIABATIC_INDEX ), &
               D_InteriorValue ( CF % INTERNAL_ENERGY ), &
               D_InteriorValue ( CF % ADIABATIC_INDEX ) )
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
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
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
      iDensity, &
      iEnergy, &
      iMomentumDim
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    !call CF % PressurelessFluidForm % ComputeRawFluxes &
    !       ( RawFlux, Value, iDimension, D_RawFlux, D_Value ) 

    call Search &
           ( CF % iaConserved, CF % CONSERVED_DENSITY, iDensity )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( 1 ), iMomentum ( 1 ) )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( 2 ), iMomentum ( 2 ) )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( 3 ), iMomentum ( 3 ) )
    call Search &
           ( CF % iaConserved, CF % CONSERVED_ENERGY, iEnergy )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( iDimension ), &
             iMomentumDim )
    
    call ComputeRawFluxesKernel &
           ( RawFlux ( :, iDensity ), &
             RawFlux ( :, iMomentum ( 1 ) ), &
             RawFlux ( :, iMomentum ( 2 ) ), &
             RawFlux ( :, iMomentum ( 3 ) ), &
             RawFlux ( :, iMomentumDim ), &
             RawFlux ( :, iEnergy ), &
             Value ( :, CF % CONSERVED_DENSITY ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
             Value ( :, CF % CONSERVED_ENERGY ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % VELOCITY ( iDimension ) ), &
             D_RawFlux ( iDensity ), &
             D_RawFlux ( iMomentum ( 1 ) ), &
             D_RawFlux ( iMomentum ( 2 ) ), &
             D_RawFlux ( iMomentum ( 3 ) ), &
             D_RawFlux ( iMomentumDim ), &
             D_RawFlux ( iEnergy ), &
             D_Value ( CF % CONSERVED_DENSITY ), &
             D_Value ( CF % MOMENTUM_DENSITY ( 1 ) ), &
             D_Value ( CF % MOMENTUM_DENSITY ( 2 ) ), &
             D_Value ( CF % MOMENTUM_DENSITY ( 3 ) ), &
             D_Value ( CF % CONSERVED_ENERGY ), &
             D_Value ( CF % PRESSURE ), &
             D_Value ( CF % VELOCITY ( iDimension ) ) )
               
  end subroutine ComputeRawFluxes
  
  
  subroutine SetOutputPolytropic &
               ( PF, UnitsOnlyOption, VelocityUnitOption, DensityUnitOption, &
                 EnergyUnitOption  )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
      PF
    logical ( KDL ), intent ( in ), optional :: &
      UnitsOnlyOption
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ), optional :: &
      VelocityUnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      DensityUnitOption, &
      EnergyUnitOption

    integer ( KDI ) :: &
      iD  !-- iDimension
    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices
    real ( KDR ) :: &
      Gamma
    type ( MeasuredValueForm ) :: &
      DensityUnit, &
      EnergyUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      VelocityUnit
    logical ( KDL ) :: &
      UnitsOnly

    call PF % SetOutputPressureless &
           ( UnitsOnlyOption = .true., VelocityUnitOption = VelocityUnitOption, &
             DensityUnitOption = DensityUnitOption )

    DensityUnit = UNIT % IDENTITY
    if ( present ( DensityUnitOption ) ) DensityUnit = DensityUnitOption

    EnergyUnit = UNIT % IDENTITY
    if ( present ( EnergyUnitOption ) ) EnergyUnit = EnergyUnitOption

    VelocityUnit = UNIT % IDENTITY
    if ( present ( VelocityUnitOption ) ) VelocityUnit = VelocityUnitOption

    UnitsOnly = .false.
    if ( present ( UnitsOnlyOption ) ) UnitsOnly = UnitsOnlyOption

    Gamma = maxval ( PF % Value ( :, PF % ADIABATIC_INDEX ) )

    PF % Unit ( PF % INTERNAL_ENERGY ) = EnergyUnit
    PF % Unit ( PF % CONSERVED_ENERGY ) = EnergyUnit
    PF % Unit ( PF % PRESSURE ) = EnergyUnit
    PF % Unit ( PF % ADIABATIC_INDEX ) = UNIT % IDENTITY
    PF % Unit ( PF % SOUND_SPEED ) = VelocityUnit ( 1 )
!    PF % Unit ( PF % POLYTROPIC_PARAMETER ) = EnergyUnit / DensityUnit ** Gamma

    if ( UnitsOnly ) return

    call VectorIndices ( 1 ) % Initialize ( PF % VELOCITY )

    call PF % Output ( 1 ) % Initialize &
           ( PF, &
             iaSelectedOption &
               = [ PF % COMOVING_DENSITY, PF % VELOCITY, &
                   PF % PRESSURE, PF % POLYTROPIC_PARAMETER ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

    call PF % DistributedMesh % SetImage ( PF % Output, PROGRAM_HEADER % Name )

  end subroutine SetOutputPolytropic


  subroutine Finalize ( PF )

    type ( PolytropicFluidForm ), intent ( inout ) :: &
      PF

    !-- Presence of this routine triggers finalization of parent

  end subroutine Finalize


  subroutine InitializeBasics &
               ( PF, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
      PF
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    type ( MeasuredValueForm ), dimension ( : ), optional, intent ( in ) :: &
      VariableUnitOption

    integer ( KDI ) :: &
      oF, &  !-- oField
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oC     !-- oConserved

    if ( PF % Type == '' ) PF % Type = 'a PolytropicFluid'

    !-- variable indices

    oF = PF % N_FIELDS_TEMPLATE + PF % N_FIELDS_PRESSURELESS
    if ( PF % N_FIELDS == 0 ) PF % N_FIELDS = oF + PF % N_FIELDS_POLYTROPIC

    PF % INTERNAL_ENERGY      = oF + 1
    PF % CONSERVED_ENERGY     = oF + 2
    PF % PRESSURE             = oF + 3
    PF % ADIABATIC_INDEX      = oF + 4
    PF % SOUND_SPEED          = oF + 5
    PF % POLYTROPIC_PARAMETER = oF + 6

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( PF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + PF % N_FIELDS_POLYTROPIC ) &
      = [ 'InternalEnergy                 ', &
          'ConservedEnergy                ', &
          'Pressure                       ', &
          'AdiabaticIndex                 ', &
          'SoundSpeed                     ', &
          'PolytropicParameter            ' ]
    
    !-- units

    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( PF % N_FIELDS ) )
      VariableUnit ( oF + 1 : oF + PF % N_FIELDS_POLYTROPIC ) = UNIT % IDENTITY
    end if

    !-- vectors

    oV = PF % N_VECTORS_TEMPLATE + PF % N_VECTORS_PRESSURELESS
    if ( PF % N_VECTORS == 0 ) PF % N_VECTORS = oV + PF % N_VECTORS_POLYTROPIC

    !-- select primitive, conserved

    oP = PF % N_PRIMITIVE_TEMPLATE + PF % N_PRIMITIVE_PRESSURELESS
    oC = PF % N_CONSERVED_TEMPLATE + PF % N_CONSERVED_PRESSURELESS

    if ( .not. allocated ( PF % iaPrimitive ) ) then
      PF % N_PRIMITIVE = oP + PF % N_PRIMITIVE_POLYTROPIC
      allocate ( PF % iaPrimitive ( PF % N_PRIMITIVE ) )
    end if
    PF % iaPrimitive ( oP + 1 : oP + PF % N_PRIMITIVE_POLYTROPIC ) &
      = [ PF % INTERNAL_ENERGY, PF % ADIABATIC_INDEX ]

    if ( .not. allocated ( PF % iaConserved ) ) then
      PF % N_CONSERVED = oC + PF % N_CONSERVED_POLYTROPIC
      allocate ( PF % iaConserved ( PF % N_CONSERVED ) )
    end if
    PF % iaConserved ( oC + 1 : oC + PF % N_CONSERVED_POLYTROPIC ) &
      = [ PF % CONSERVED_ENERGY ]
    
  end subroutine InitializeBasics
  
  
end module PolytropicFluid_Form
