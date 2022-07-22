module ConservedFields_Template

  use Basics
  use DistributedMesh_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_TEMPLATE = 0, &
      N_CONSERVED_TEMPLATE = 0, &
      N_FIELDS_TEMPLATE    = 6, &
      N_VECTORS_TEMPLATE   = 0

  type, public, extends ( StorageForm ), abstract :: &
    ConservedFieldsTemplate
      integer ( KDI ) :: &
        N_PRIMITIVE = 0, &
        N_CONSERVED = 0, &
        N_FIELDS    = 0, &
        N_VECTORS   = 0, &
        N_PRIMITIVE_TEMPLATE = N_PRIMITIVE_TEMPLATE, &
        N_CONSERVED_TEMPLATE = N_CONSERVED_TEMPLATE, &
        N_FIELDS_TEMPLATE    = N_FIELDS_TEMPLATE, &
        N_VECTORS_TEMPLATE   = N_VECTORS_TEMPLATE
      integer ( KDI ), dimension ( 3 ) :: &
        FAST_EIGENSPEED_PLUS  = 0, &
        FAST_EIGENSPEED_MINUS = 0
      integer ( KDI ), dimension ( : ), allocatable :: &
        iaPrimitive, &
        iaConserved
      character ( LDL ) :: &
        Type = ''
      type ( StorageForm ), dimension ( 1 ) :: &
        Output
      type ( DistributedMeshForm ), pointer :: &
        DistributedMesh => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      SetOutputTemplate
    procedure ( ComputeInterface ), public, pass, deferred :: &
      ComputeConserved
    procedure ( ComputeInterface ), public, pass, deferred :: &
      ComputePrimitive
    procedure ( ComputeInterface ), public, pass, deferred :: &
      ComputeAuxiliary
    procedure ( ApplyBoundaryConditionsInterface ), public, pass, deferred :: &
      ApplyBoundaryConditions
    procedure ( ComputeRawFluxesInterface ), public, pass, deferred :: &
      ComputeRawFluxes
    procedure ( ComputeRiemannSolverInputInterface ), public, pass, &
      deferred :: &
        ComputeRiemannSolverInput
  end type ConservedFieldsTemplate

  abstract interface 

    subroutine ComputeInterface ( CF, Value, UseDeviceOption )
      use Basics
      import ConservedFieldsTemplate
      class ( ConservedFieldsTemplate ), intent ( inout ) :: &
        CF
      real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
        Value
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeInterface
    
    subroutine ApplyBoundaryConditionsInterface &
                 ( CF, ExteriorValue, InteriorValue, iDimension, iBoundary, &
                   PrimitiveOnlyOption, UseDeviceOption )
      use Basics
      import ConservedFieldsTemplate
      class ( ConservedFieldsTemplate ), intent ( inout ) :: &
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
    end subroutine ApplyBoundaryConditionsInterface
    
    subroutine ComputeRawFluxesInterface &
                 ( CF, RawFlux, Value, iDimension, UseDeviceOption )
      use Basics
      import ConservedFieldsTemplate
      class ( ConservedFieldsTemplate ), intent ( inout ) :: &
        CF
      real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
        RawFlux
      real ( KDR ), dimension ( :, : ), intent ( in ) :: &
        Value
      integer ( KDI ), intent ( in ) :: &
        iDimension
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeRawFluxesInterface
    
    subroutine ComputeRiemannSolverInputInterface &
                 ( CF, Step, ValueInner, ValueOuter, iDimension, &
                   UseDeviceOption )
      use Basics
      import ConservedFieldsTemplate    
      class ( ConservedFieldsTemplate ), intent ( inout ) :: &
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
    end subroutine ComputeRiemannSolverInputInterface
    
  end interface

contains


  subroutine InitializeTemplate &
               ( CF, DistributedMesh, VectorOption, VariableOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( ConservedFieldsTemplate ), intent ( inout ) :: &
      CF
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

    type ( QuantityForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable

    call InitializeBasics &
           ( CF, Variable, VariableUnit, VariableOption, NameOption, &
             UnitOption )

    associate ( DM => DistributedMesh )
    call CF % Initialize &
           ( [ DM % nProperCells + DM % nGhostCells, CF % N_FIELDS ], &
             VectorIndicesOption = VectorIndicesOption, &
             UnitOption = VariableUnit, VectorOption = VectorOption, &
             VariableOption = Variable, NameOption = NameOption, &
             ClearOption = ClearOption, PinnedOption = .true. )
    end associate !-- DM
    
    CF % DistributedMesh => DistributedMesh

  end subroutine InitializeTemplate


  subroutine SetOutputTemplate ( CF )

    class ( ConservedFieldsTemplate ), intent ( inout ) :: &
      CF

  end subroutine SetOutputTemplate
  
  
  subroutine InitializeBasics &
               ( CF, Variable, VariableUnit, VariableOption, NameOption, &
                 VariableUnitOption )

    class ( ConservedFieldsTemplate ), intent ( inout ) :: &
      CF
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable
    type ( QuantityForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption

    integer ( KDI ) :: &
      iF  !-- iField
    character ( LDF ) :: &
      Name
    character ( LDF ), dimension ( : ), allocatable :: &
      PrimitiveName, &
      ConservedName

    if ( CF % Type == '' ) CF % Type = 'a ConservedFields'

    Name = ''
    if ( present ( NameOption ) ) Name = NameOption

    call Show ( 'Initializing ' // trim ( CF % Type ), CONSOLE % INFO_3 )
    call Show ( Name, 'Name', CONSOLE % INFO_3 )

    !-- variable indices

    if ( CF % N_FIELDS == 0 ) CF % N_FIELDS = CF % N_FIELDS_TEMPLATE

    CF % FAST_EIGENSPEED_PLUS  = [ 1, 2, 3 ]
    CF % FAST_EIGENSPEED_MINUS = [ 4, 5, 6 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( CF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : CF % N_FIELDS_TEMPLATE ) &
      = [ 'FastEigenspeedPlus_1           ', &
          'FastEigenspeedPlus_2           ', &
          'FastEigenspeedPlus_3           ', &
          'FastEigenspeedMinus_1          ', &
          'FastEigenspeedMinus_2          ', &
          'FastEigenspeedMinus_3          ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( CF % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( CF % N_VECTORS == 0 ) CF % N_VECTORS = CF % N_VECTORS_TEMPLATE

    !-- show primitive, conserved

    allocate ( PrimitiveName ( CF % N_PRIMITIVE ) )
    do iF = 1, CF % N_PRIMITIVE
      PrimitiveName ( iF ) = Variable ( CF % iaPrimitive ( iF ) )
    end do
    call Show ( PrimitiveName, 'Primitive', CONSOLE % INFO_3 )
 
    allocate ( ConservedName ( CF % N_CONSERVED ) )
    do iF = 1, CF % N_CONSERVED
      ConservedName ( iF ) = Variable ( CF % iaConserved ( iF ) )
    end do
    call Show ( ConservedName, 'Conserved', CONSOLE % INFO_3 )
 
  end subroutine InitializeBasics


end module ConservedFields_Template
