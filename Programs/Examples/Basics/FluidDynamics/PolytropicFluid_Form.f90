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


  subroutine ComputeConserved ( CF, Value )

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

  end subroutine ComputeConserved


  subroutine ComputePrimitive ( CF, Value )
    
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
  
  end subroutine ComputePrimitive
  
  
  subroutine ComputeAuxiliary ( CF, Value )
    
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
             Value ( :, CF % ADIABATIC_INDEX ) )
    
  end subroutine ComputeAuxiliaryFromPressure
  
  
  subroutine ApplyBoundaryConditions &
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
  
  end subroutine ApplyBoundaryConditions


  subroutine ComputeRawFluxes ( CF, RawFlux, Value, iDimension )
    
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      CF
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value
    integer ( KDI ), intent ( in ) :: &
      iDimension
    
    integer ( KDI ) :: &
      iEnergy, &
      iMomentum

    call CF % PressurelessFluidForm % ComputeRawFluxes &
           ( RawFlux, Value, iDimension ) 

    call Search &
           ( CF % iaConserved, CF % CONSERVED_ENERGY, iEnergy )
    call Search &
           ( CF % iaConserved, CF % MOMENTUM_DENSITY ( iDimension ), &
             iMomentum )
    
    call ComputeRawFluxesKernel &
           ( RawFlux ( :, iMomentum ), &
             RawFlux ( :, iEnergy ), &
             Value ( :, CF % CONSERVED_ENERGY ), &
             Value ( :, CF % PRESSURE ), &
             Value ( :, CF % VELOCITY ( iDimension ) ) )
               
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
    PF % Unit ( PF % POLYTROPIC_PARAMETER ) = EnergyUnit / DensityUnit ** Gamma

    if ( UnitsOnly ) return

    call VectorIndices ( 1 ) % Initialize ( PF % VELOCITY )

    call PF % Output ( 1 ) % Initialize &
           ( PF, &
             iaSelectedOption &
               = [ PF % COMOVING_DENSITY, PF % VELOCITY, &
                   PF % PRESSURE ], &
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


  subroutine ComputePrimitiveKernel ( E, G, N, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      G  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      V_1, V_2, V_3

    associate ( KE => 0.5_KDR * N * ( V_1 ** 2 + V_2 ** 2 + V_3 ** 2 ) )

    E = G - KE
    where ( E < 0.0_KDR )
      E = 0.0_KDR
      G = KE
    end where

    end associate !-- KE

  end subroutine ComputePrimitiveKernel


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


  subroutine ComputeRawFluxesKernel ( F_S_Dim, F_G, G, P, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_S_Dim, &
      F_G
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      G, &
      P, &
      V_Dim

    F_S_Dim = F_S_Dim + P
    F_G = ( G + P ) * V_Dim

  end subroutine ComputeRawFluxesKernel


end module PolytropicFluid_Form
