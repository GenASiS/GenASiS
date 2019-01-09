#include "Preprocessor"

module PolytropicFluid_Form

  use iso_c_binding
  use Basics
  use DistributedMesh_Form
  use PressurelessFluid_Form

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
      InitializeBasics, &
      ComputeConservedKernel, &
      ComputeConservedKernelDevice, &
      ComputePrimitiveKernel, &
      ComputeAuxiliaryKernel, &
      ComputeAuxiliaryKernelDevice, &
      ComputeAuxiliaryFromPressureKernel, &
      ComputeEigenspeedsKernel, &
      ComputeEigenspeedsKernelDevice, &
      ApplyBoundaryConditionsReflecting, &
      ApplyBoundaryConditionsReflectingDevice, &
      ComputeRawFluxesKernel

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
  
  
  subroutine ComputeConservedKernel ( G, E, N, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      E, &
      N, &
      V_1, V_2, V_3

    G = E  +  0.5_KDR * N * ( V_1 **2  +  V_2 ** 2  +  V_3 ** 2 )

  end subroutine ComputeConservedKernel


  subroutine ComputeConservedKernelDevice &
               ( G, E, N, V_1, V_2, V_3, &
                 D_G, D_E, D_N, D_V_1, D_V_2, D_V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      E, &
      N, &
      V_1, V_2, V_3
    type ( c_ptr ), intent ( in ) :: &
      D_G, &
      D_E, &
      D_N, &
      D_V_1, D_V_2, D_V_3
      
    integer ( KDI ) :: &
      iV
    
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( G )
      G ( iV ) = E ( iV ) + 0.5_KDR * N ( iV ) &
                 * ( V_1 ( iV ) * V_1 ( iV ) &
                     + V_2 ( iV ) * V_2 ( iV ) &
                     + V_3 ( iV ) * V_3 ( iV ) )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    
  end subroutine ComputeConservedKernelDevice


  subroutine ComputePrimitiveKernel ( E, G, N, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3
      
    real ( KDR ), dimension ( size ( N ) ) :: &
      KE
      
    !-- FIXME: 'associate' construct causes segfault with IBM XL
    !associate ( KE => 0.5_KDR * N * ( V_1 ** 2 + V_2 ** 2 + V_3 ** 2 ) )
    
    KE = 0.5_KDR * N * ( V_1 ** 2 + V_2 ** 2 + V_3 ** 2 )
    E  = G - KE
    where ( E < 0.0_KDR )
      E = 0.0_KDR
      G = KE
    end where

    !end associate !-- KE

  end subroutine ComputePrimitiveKernel


  subroutine ComputePrimitiveKernelDevice &
               ( E, G, N, V_1, V_2, V_3, &
                 D_E, D_G, D_N, D_V_1, D_V_2, D_V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3
    type ( c_ptr ), intent ( in ) :: &
      D_E, &
      D_G, &
      D_N, &
      D_V_1, D_V_2, D_V_3
      
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      KE
    
    
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE ) private ( KE )
    do iV = 1, size ( E )
    
      KE = 0.5_KDR * N ( iV ) &
             * ( V_1 ( iV ) * V_1 ( iV )  +  V_2 ( iV ) * V_2 ( iV ) &
                 + V_3 ( iV ) * V_3 ( iV ) )
      E ( iV )  = G ( iV ) - KE

      if ( E ( iV ) < 0.0_KDR ) then
        E ( iV ) = 0.0_KDR
        G ( iV ) = KE
      end if

    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( G )
    call DisassociateHost ( E )


  end subroutine ComputePrimitiveKernelDevice


  subroutine ComputeAuxiliaryKernel ( P, K, N, E, Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      E, &
      Gamma

    P = E * ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end subroutine ComputeAuxiliaryKernel


  subroutine ComputeAuxiliaryKernelDevice &
               ( P, K, N, E, Gamma, &
                 D_P, D_K, D_N, D_E, D_Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      P, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      E, &
      Gamma
    type ( c_ptr ), intent ( in ) :: &
      D_P, &
      D_K, &
      D_N, &
      D_E, &
      D_Gamma
    
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_K, K )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_E, E )
    call AssociateHost ( D_Gamma, Gamma )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( P )
      P ( iV ) = E ( iV ) * ( Gamma ( iV ) - 1.0_KDR )
      if ( N ( iV ) > 0.0_KDR ) then
        K ( iV ) = P ( iV ) / ( N ( iV ) ** Gamma ( iV ) )
      else
        K ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    call DisassociateHost ( Gamma )
    call DisassociateHost ( E )
    call DisassociateHost ( N )
    call DisassociateHost ( K )
    call DisassociateHost ( P )

  end subroutine ComputeAuxiliaryKernelDevice


  subroutine ComputeAuxiliaryFromPressureKernel ( E, K, N, P, Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      K
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      P, &
      Gamma

    E = P / ( Gamma - 1.0_KDR ) 

    where ( N ** Gamma > 0.0_KDR )
      K = P / ( N ** Gamma )
    elsewhere
      K = 0.0_KDR
    end where

  end subroutine ComputeAuxiliaryFromPressureKernel


  subroutine ComputeEigenspeedsKernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, N, &
                 V_1, V_2, V_3, P, Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3, &
      P, &
      Gamma

    where ( N > 0.0_KDR .and. P > 0.0_KDR )
      CS = sqrt ( Gamma * P / N )
    elsewhere
      CS = 0.0_KDR
    end where 
    
    FEP_1 = V_1 + CS
    FEP_2 = V_2 + CS
    FEP_3 = V_3 + CS
    FEM_1 = V_1 - CS
    FEM_2 = V_2 - CS
    FEM_3 = V_3 - CS

  end subroutine ComputeEigenspeedsKernel


  subroutine ComputeEigenspeedsKernelDevice &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, N, &
                 V_1, V_2, V_3, P, Gamma, &
                 D_FEP_1, D_FEP_2, D_FEP_3, D_FEM_1, D_FEM_2, D_FEM_3, &
                 D_CS, D_N, D_V_1, D_V_2, D_V_3, D_P, D_Gamma )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      CS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3, &
      P, &
      Gamma
    type ( c_ptr ), intent ( in ) :: &
      D_FEP_1, D_FEP_2, D_FEP_3, &
      D_FEM_1, D_FEM_2, D_FEM_3, &
      D_CS, &
      D_N, &
      D_V_1, D_V_2, D_V_3, &
      D_P, &
      D_Gamma
      
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_FEP_1, FEP_1 )
    call AssociateHost ( D_FEP_2, FEP_2 )
    call AssociateHost ( D_FEP_3, FEP_3 )
    call AssociateHost ( D_FEM_1, FEM_1 )
    call AssociateHost ( D_FEM_2, FEM_2 )
    call AssociateHost ( D_FEM_3, FEM_3 )
    call AssociateHost ( D_CS, CS )
    call AssociateHost ( D_N, N )
    call AssociateHost ( D_V_1, V_1 )
    call AssociateHost ( D_V_2, V_2 )
    call AssociateHost ( D_V_3, V_3 )
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_Gamma, Gamma )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( N )
      if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
        CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / N ( iV ) )
      else
        CS ( iV ) = 0.0_KDR
      end if
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( N )
      FEP_1 ( iV ) = V_1 ( iV ) + CS ( iV )
      FEP_2 ( iV ) = V_2 ( iV ) + CS ( iV )
      FEP_3 ( iV ) = V_3 ( iV ) + CS ( iV )
      FEM_1 ( iV ) = V_1 ( iV ) - CS ( iV )
      FEM_2 ( iV ) = V_2 ( iV ) - CS ( iV )
      FEM_3 ( iV ) = V_3 ( iV ) - CS ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( Gamma )
    call DisassociateHost ( P )
    call DisassociateHost ( V_3 )
    call DisassociateHost ( V_2 )
    call DisassociateHost ( V_1 )
    call DisassociateHost ( N )
    call DisassociateHost ( CS )
    call DisassociateHost ( FEM_3 )
    call DisassociateHost ( FEM_2 )
    call DisassociateHost ( FEM_1 )
    call DisassociateHost ( FEP_3 )
    call DisassociateHost ( FEP_2 )
    call DisassociateHost ( FEP_1 )

  end subroutine ComputeEigenspeedsKernelDevice
  
  
  subroutine ApplyBoundaryConditionsReflecting &
               ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      E_E, &
      Gamma_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      E_I, &
      Gamma_I
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB,  & 
      oBE, &
      oBI

    E_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
          oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
          oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = E_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
              oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
              oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

    Gamma_E ( oBE ( 1 ) + 1 : oBE ( 1 ) + nB ( 1 ), &
          oBE ( 2 ) + 1 : oBE ( 2 ) + nB ( 2 ), &
          oBE ( 3 ) + 1 : oBE ( 3 ) + nB ( 3 ) ) &
      = Gamma_I ( oBI ( 1 ) + 1 : oBI ( 1 ) + nB ( 1 ), &
              oBI ( 2 ) + 1 : oBI ( 2 ) + nB ( 2 ), &
              oBI ( 3 ) + 1 : oBI ( 3 ) + nB ( 3 ) )

  end subroutine ApplyBoundaryConditionsReflecting


  subroutine ApplyBoundaryConditionsReflectingDevice &
               ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI, &
                 D_E_E, D_Gamma_E, D_E_I, D_Gamma_I )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      E_E, &
      Gamma_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      E_I, &
      Gamma_I
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nB,  & 
      oBE, &
      oBI
    type ( c_ptr ), intent ( in ) :: &
      D_E_E, &
      D_Gamma_E, &
      D_E_I, &
      D_Gamma_I
      
    integer ( KDI ) :: &
      iV, jV, kV
      
    call AssociateHost ( D_E_E, E_E )
    call AssociateHost ( D_Gamma_E, Gamma_E )
    call AssociateHost ( D_E_I, E_I )
    call AssociateHost ( D_Gamma_I, Gamma_I )
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
    !$OMP& schedule ( OMP_SCHEDULE )
    do kV = 1, nB ( 3 )
      do jV = 1, nB ( 2 )
        do iV = 1, nB ( 1 )
        
          E_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = E_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )
          
          Gamma_E ( oBE ( 1 ) + iV, oBE ( 2 ) + jV, oBE ( 3 ) + kV ) &
            = Gamma_I ( oBI ( 1 ) + iV, oBI ( 2 ) + jV, oBI ( 3 ) + kV )
            
        end do
      end do
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    call DisassociateHost ( Gamma_I )
    call DisassociateHost ( E_I )
    call DisassociateHost ( Gamma_E )
    call DisassociateHost ( E_E )

  end subroutine ApplyBoundaryConditionsReflectingDevice


  subroutine ComputeRawFluxesKernel &
               ( F_D, F_S_1, F_S_2, F_S_3, F_S_Dim, F_G, D, S_1, S_2, S_3, &
                 G, P, V_Dim, D_F_D, D_F_S_1, D_F_S_2, D_F_S_3, D_F_S_Dim, &
                 D_F_G, D_D, D_S_1, D_S_2, D_S_3, D_G, D_P, D_V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_D, &
      F_S_1, F_S_2, F_S_3, &
      F_S_Dim, &
      F_G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      D, &
      S_1, S_2, S_3, &
      G, &
      P, &
      V_Dim
    type ( c_ptr ), intent ( in ) :: &
      D_F_D, &
      D_F_S_1, D_F_S_2, D_F_S_3, &
      D_F_S_Dim, &
      D_F_G, &
      D_D, &
      D_S_1, D_S_2, D_S_3, &
      D_G, &
      D_P, &
      D_V_Dim
      
    integer ( KDI ) :: &
      iV
      
    call AssociateHost ( D_F_D, F_D )
    call AssociateHost ( D_F_S_1, F_S_1 )
    call AssociateHost ( D_F_S_2, F_S_2 )
    call AssociateHost ( D_F_S_3, F_S_3 )
    call AssociateHost ( D_F_S_Dim, F_S_Dim )
    call AssociateHost ( D_F_G, F_G )
    call AssociateHost ( D_D, D )
    call AssociateHost ( D_S_1, S_1 )
    call AssociateHost ( D_S_2, S_2 )
    call AssociateHost ( D_S_3, S_3 )
    call AssociateHost ( D_G, G )
    call AssociateHost ( D_P, P )
    call AssociateHost ( D_V_Dim, V_Dim )
    
    !F_D   = D   * V_Dim
    !F_S_1 = S_1 * V_Dim
    !F_S_2 = S_2 * V_Dim
    !F_S_3 = S_3 * V_Dim
    !F_S_Dim = F_S_Dim + P
    !F_G = ( G + P ) * V_Dim
    
    !$OMP  OMP_TARGET_DIRECTIVE parallel do &
    !$OMP& schedule ( OMP_SCHEDULE )
    do iV = 1, size ( F_D )
      F_D ( iV )     = D ( iV )   * V_Dim ( iV ) 
      F_S_1 ( iV )   = S_1 ( iV ) * V_Dim ( iV ) 
      F_S_2 ( iV )   = S_2 ( iV ) * V_Dim ( iV ) 
      F_S_3 ( iV )   = S_3 ( iV ) * V_Dim ( iV ) 
      F_S_Dim ( iV ) = F_S_Dim ( iV ) + P ( iV ) 
      F_G ( iV )     = ( G ( iV ) + P ( iV ) ) * V_Dim ( iV )
    end do
    !$OMP end OMP_TARGET_DIRECTIVE parallel do

    call DisassociateHost ( V_Dim )
    call DisassociateHost ( P )
    call DisassociateHost ( G )
    call DisassociateHost ( S_3 )
    call DisassociateHost ( S_2 )
    call DisassociateHost ( S_1 )
    call DisassociateHost ( D )
    call DisassociateHost ( F_G )
    call DisassociateHost ( F_S_Dim )
    call DisassociateHost ( F_S_3 )
    call DisassociateHost ( F_S_2 )
    call DisassociateHost ( F_S_1 )
    call DisassociateHost ( F_D )
    
  end subroutine ComputeRawFluxesKernel


end module PolytropicFluid_Form
