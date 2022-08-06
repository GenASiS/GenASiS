module PolytropicFluid_Form

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
      ComputeConserved
    procedure, public, pass :: &
      ComputePrimitive
    procedure, public, pass :: &
      ComputeAuxiliary
    procedure, public, pass :: &
      ComputeAuxiliaryFromPressure
    procedure, public, pass :: &
      ApplyBoundaryConditions
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
      ComputePrimitiveKernel, &
      ComputeAuxiliaryKernel, &
      ComputeAuxiliaryFromPressureKernel, &
      ComputeEigenspeedsKernel, &
      ApplyBoundaryConditionsReflecting, &
      ComputeRawFluxesKernel
      

    interface
    
      module subroutine ComputeConservedKernel &
                   ( G, E, N, V_1, V_2, V_3, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          G
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          E, &
          N, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeConservedKernel
      
      module subroutine ComputePrimitiveKernel &
                   ( E, G, N, V_1, V_2, V_3, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          E, &
          G  
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N, &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputePrimitiveKernel

      module subroutine ComputeAuxiliaryKernel &
                   ( P, K, N, E, Gamma, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          P, &
          K
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N, &
          E, &
          Gamma
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeAuxiliaryKernel

      module subroutine ComputeAuxiliaryFromPressureKernel ( E, K, N, P, Gamma )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          E, &
          K
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N, &
          P, &
          Gamma
      end subroutine ComputeAuxiliaryFromPressureKernel

      module subroutine ComputeEigenspeedsKernel &
                   ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, N, &
                     V_1, V_2, V_3, P, Gamma, UseDevice )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          FEP_1, FEP_2, FEP_3, &
          FEM_1, FEM_2, FEM_3, &
          CS
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N, &
          V_1, V_2, V_3, &
          P, &
          Gamma
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ComputeEigenspeedsKernel
      
      module subroutine ApplyBoundaryConditionsReflecting &
                   ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI, UseDevice )
        use Basics
        implicit none
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
        logical ( KDL ), intent ( in ) :: &
          UseDevice
      end subroutine ApplyBoundaryConditionsReflecting

      module subroutine ComputeRawFluxesKernel &
                   ( F_D, F_S_1, F_S_2, F_S_3, F_G, &
                     D, S_1, S_2, S_3, G, P, V_Dim, UseDevice, iDim )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          F_D, &
          F_S_1, &
          F_S_2, &
          F_S_3, &
          F_G
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          S_1, S_2, S_3, &
          G, &
          P, &
          V_Dim
        logical ( KDL ), intent ( in ) :: &
          UseDevice
        integer ( KDI ), intent ( in ) :: &
          iDim
      end subroutine ComputeRawFluxesKernel

    end interface
    
      
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
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      VariableUnit

    call InitializeBasics &
           ( PF, Variable, VariableUnit, VariableOption, UnitOption )

    call PF % PressurelessFluidForm % InitializeWithMesh &
           ( DistributedMesh, VectorOption = VectorOption, &
             VariableOption = Variable, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeWithMesh


  subroutine ComputeConserved ( CF, Value, UseDeviceOption )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
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
    
    call CF % PressurelessFluidForm &
            % ComputeConserved ( Value, UseDeviceOption )

    call ComputeConservedKernel &
      ( Value ( :, CF % CONSERVED_ENERGY ), &
        Value ( :, CF % INTERNAL_ENERGY ), &    
        Value ( :, CF % COMOVING_DENSITY ), &
        Value ( :, CF % VELOCITY ( 1 ) ), &    
        Value ( :, CF % VELOCITY ( 2 ) ), &   
        Value ( :, CF % VELOCITY ( 3 ) ), &
        UseDevice )

  end subroutine ComputeConserved


  subroutine ComputePrimitive ( CF, Value, UseDeviceOption )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
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
    
    call CF % PressurelessFluidForm &
            % ComputePrimitive ( Value, UseDeviceOption )

    call ComputePrimitiveKernel &
           ( Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % CONSERVED_ENERGY ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % VELOCITY ( 1 ) ), &
             Value ( :, CF % VELOCITY ( 2 ) ), &
             Value ( :, CF % VELOCITY ( 3 ) ), &
             UseDevice )
               
  end subroutine ComputePrimitive
  
  
  subroutine ComputeAuxiliary ( CF, Value, UseDeviceOption )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      Value
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = OnDevice ( Value ( :, CF % PRESSURE ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    call ComputeAuxiliaryKernel &
           ( Value ( :, CF % PRESSURE ), &
             Value ( :, CF % POLYTROPIC_PARAMETER ), &
             Value ( :, CF % COMOVING_DENSITY ), &
             Value ( :, CF % INTERNAL_ENERGY ), &
             Value ( :, CF % ADIABATIC_INDEX ), &
             UseDevice )
  
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
             Value ( :, CF % ADIABATIC_INDEX ), &
             UseDevice )
    
  end subroutine ComputeAuxiliary
  
  
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
             Value ( :, CF % ADIABATIC_INDEX ), &
             UseDevice = .false. )  !-- Only called at initialization
    
  end subroutine ComputeAuxiliaryFromPressure
  
  
  subroutine ApplyBoundaryConditions &
              ( CF, ExteriorValue, InteriorValue, iDimension, iBoundary, &
                PrimitiveOnlyOption, UseDeviceOption )

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
      PrimitiveOnlyOption, &
      UseDeviceOption

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
      PrimitiveOnly, &
      UseDevice

    call CF % PressurelessFluidForm % ApplyBoundaryConditions &
           ( ExteriorValue, InteriorValue, iDimension, iBoundary, &
             PrimitiveOnlyOption = .true., UseDeviceOption = UseDeviceOption )

    PrimitiveOnly = .false.
    if ( present ( PrimitiveOnlyOption ) ) PrimitiveOnly = PrimitiveOnlyOption
    
    UseDevice = OnDevice ( InteriorValue ( :, CF % INTERNAL_ENERGY ) )
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
             ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI, UseDevice )
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
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
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
      iDensity, &
      iEnergy, &
      iMomentumDim
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = OnDevice ( Value ( :, CF % CONSERVED_DENSITY ) )
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    !call CF % PressurelessFluidForm % ComputeRawFluxes &
    !       ( RawFlux, Value, iDimension, UseDeviceOption = UseDeviceOption ) 

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
    
    call ComputeRawFluxesKernel &
           ( RawFlux ( :, iDensity ), &
             RawFlux ( :, iMomentum ( 1 ) ), &
             RawFlux ( :, iMomentum ( 2 ) ), &
             RawFlux ( :, iMomentum ( 3 ) ), &
             RawFlux ( :, iEnergy ), &
             Value ( :, CF % CONSERVED_DENSITY ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 1 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 2 ) ), &
             Value ( :, CF % MOMENTUM_DENSITY ( 3 ) ), &
             Value ( :, CF % CONSERVED_ENERGY ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % VELOCITY ( iDimension ) ), &
             UseDevice, iDimension )

  end subroutine ComputeRawFluxes
  
  
  subroutine SetOutputPolytropic &
               ( PF, UnitsOnlyOption, VelocityUnitOption, DensityUnitOption, &
                 EnergyUnitOption  )

    class ( PolytropicFluidForm ), intent ( inout ) :: &
      PF
    logical ( KDL ), intent ( in ), optional :: &
      UnitsOnlyOption
    type ( QuantityForm ), dimension ( 3 ), intent ( in ), optional :: &
      VelocityUnitOption
    type ( QuantityForm ), intent ( in ), optional :: &
      DensityUnitOption, &
      EnergyUnitOption

    integer ( KDI ) :: &
      iD  !-- iDimension
    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices
    real ( KDR ) :: &
      Gamma
    type ( QuantityForm ) :: &
      DensityUnit, &
      EnergyUnit
    type ( QuantityForm ), dimension ( 3 ) :: &
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
                   PF % INTERNAL_ENERGY, PF % PRESSURE, &
                   PF % POLYTROPIC_PARAMETER ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

    call PF % DistributedMesh % SetImage &
           ( Output = PF % Output, Name = PROGRAM_HEADER % Name )

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
    type ( QuantityForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    type ( QuantityForm ), dimension ( : ), optional, intent ( in ) :: &
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
