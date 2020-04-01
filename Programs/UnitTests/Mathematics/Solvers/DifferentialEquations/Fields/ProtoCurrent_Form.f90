module ProtoCurrent_Form

  use Basics
  use Manifolds
  use Current_Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_PC = 1, &
      N_CONSERVED_PC = 1, &
      N_FIELDS_PC    = 5, &
      N_VECTORS_PC   = 1

  type, public, extends ( CurrentTemplate ) :: ProtoCurrentForm
    integer ( KDI ) :: &
      N_PRIMITIVE_PC = N_PRIMITIVE_PC, &
      N_CONSERVED_PC = N_CONSERVED_PC, &
      N_FIELDS_PC    = N_FIELDS_PC, &
      N_VECTORS_PC   = N_VECTORS_PC, &
      COMOVING_DENSITY  = 0, &
      CONSERVED_DENSITY = 0
    integer ( KDI ), dimension ( 3 ) :: &
      VELOCITY = 0
    real ( KDR ) :: &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber
  contains
    procedure, public, pass :: &
      InitializeAllocate_PC
    generic, public :: &
      Initialize => InitializeAllocate_PC
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, nopass :: &
      ComputeConservedDensityKernel
    procedure, public, nopass :: &
      ComputeComovingDensityKernel
    procedure, public, nopass :: &
      ComputeEigenspeedsKernel
  end type ProtoCurrentForm

    private :: &
      InitializeBasics, &
      ComputeVelocityKernel, &
      ComputeRawFluxesKernel
    
    interface
    
      module subroutine ComputeConservedDensityKernel &
               ( D, N, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          D
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeConservedDensityKernel


      module subroutine ComputeComovingDensityKernel ( N, D, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeComovingDensityKernel


      module subroutine ComputeEigenspeedsKernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3, &
                 UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          FEP_1, FEP_2, FEP_3, &
          FEM_1, FEM_2, FEM_3
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          V_1, V_2, V_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeEigenspeedsKernel


      module subroutine ComputeVelocityKernel &
               ( V_1, V_2, V_3, K, V, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          V_1, V_2, V_3
        real ( KDR ), dimension ( 3 ), intent ( in ) :: &
          K
        real ( KDR ), intent ( in ) :: &
          V
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeVelocityKernel


      module subroutine ComputeRawFluxesKernel &
                ( F_D, D, V_Dim, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          F_D
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          V_Dim
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeRawFluxesKernel
    
    end interface

contains


  subroutine InitializeAllocate_PC &
               ( PC, Velocity_U_Unit, nValues, VariableOption, VectorOption, &
                 NameOption, ClearOption, UnitOption, VectorIndicesOption )

    class ( ProtoCurrentForm ), intent ( inout ) :: &
      PC
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector

    call InitializeBasics &
           ( PC, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call PC % InitializeTemplate &
           ( RiemannSolverType = 'HLL', ReconstructedType = 'PRIMITIVE', &
             UseLimiter = .true., Velocity_U_Unit = Velocity_U_Unit, &
             LimiterParameter = 2.0_KDR, nValues = nValues, &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = ClearOption, &
             UnitOption = VariableUnit, VectorIndicesOption = VectorIndices )

    call PC % SetPrimitiveConserved ( )

  end subroutine InitializeAllocate_PC
  

  subroutine SetPrimitiveConserved ( C )

    class ( ProtoCurrentForm ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC, &  !-- oConserved
      Ignorability
    character ( LDL ), dimension ( C % N_PRIMITIVE_PC ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_PC ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE
    oC = C % N_CONSERVED_TEMPLATE

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_PC
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_PC ) &
      = [ C % COMOVING_DENSITY ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_PC
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_PC ) &
      = [ C % CONSERVED_DENSITY ]

    do iF = 1, C % N_PRIMITIVE_PC
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_PC
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )
    
  end subroutine SetPrimitiveConserved


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Storage_C, Storage_G, iDimension, &
                 nValuesOption, oValueOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      RawFlux
    class ( ProtoCurrentForm ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( inout ) :: &
      Storage_C, &
      Storage_G
    integer ( KDI ), intent ( in ) :: &
      iDimension
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption
    
    integer ( KDI ) :: &
      iDensity
    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( Storage_C % Value, dim = 1 )
    end if
    
    call Search ( C % iaConserved, C % CONSERVED_DENSITY, iDensity )
    
    associate &
      ( F_D   => RawFlux % Value ( oV + 1 : oV + nV, iDensity ), &
        D     => Storage_C % Value ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ), &
        V_Dim => Storage_C % Value ( oV + 1 : oV + nV, C % VELOCITY ( iDimension ) ) )

    call ComputeRawFluxesKernel &
           ( F_D, D, V_Dim, UseDeviceOption = Storage_C % AllocatedDevice )

    end associate !-- F_D, etc.

  end subroutine ComputeRawFluxes
  
  
  subroutine SetOutput ( PC, Output )

    class ( ProtoCurrentForm ), intent ( inout ) :: &
      PC
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( PC % VELOCITY )
    call Output % Initialize &
           ( PC, iaSelectedOption = [ PC % COMOVING_DENSITY, PC % VELOCITY ], &
             VectorOption = [ 'Velocity' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  impure elemental subroutine Finalize ( PC )

    type ( ProtoCurrentForm ), intent ( inout ) :: &
      PC

    call PC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeFromPrimitiveCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( ProtoCurrentForm ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => Storage_C % Value, &
        GV => Storage_G % Value )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( FV, dim = 1 )
    end if

    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ) )
    associate &
      ( FEP_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ) )

    call ComputeVelocityKernel &
           ( V_1, V_2, V_3, C % Wavenumber, C % Speed, &
             UseDeviceOption = Storage_C % AllocatedDevice )
    call C % ComputeConservedDensityKernel &
           ( D, N, UseDeviceOption = Storage_C % AllocatedDevice )
    call C % ComputeEigenspeedsKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3, &
             UseDeviceOption = Storage_C % AllocatedDevice )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ), target :: &
      Storage_C
    class ( ProtoCurrentForm ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => Storage_C % Value, &
        GV => Storage_G % Value )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( FV, dim = 1 )
    end if

    associate &
      ( M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => FV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ) )

    call ComputeVelocityKernel &
           ( V_1, V_2, V_3, C % Wavenumber, C % Speed, &
             UseDeviceOption = Storage_C % AllocatedDevice )
    call C % ComputeComovingDensityKernel &
           ( N, D, UseDeviceOption = Storage_C % AllocatedDevice )
    call C % ComputeEigenspeedsKernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3, &
             UseDeviceOption = Storage_C % AllocatedDevice )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.

  end subroutine ComputeFromConservedCommon


  subroutine InitializeBasics &
               ( PC, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( ProtoCurrentForm ), intent ( inout ) :: &
      PC
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    character ( LDF ), intent ( out ) :: &
      Name
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
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      oF, &  !-- oField
      oV, &  !-- oVector
      iV     !-- iVector

    if ( PC % Type == '' ) &
      PC % Type = 'a ProtoCurrent'

    Name = 'ProtoCurrent'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = PC % N_FIELDS_TEMPLATE
    if ( PC % N_FIELDS == 0 ) PC % N_FIELDS = oF + PC % N_FIELDS_PC

    PC % COMOVING_DENSITY  = oF + 1
    PC % CONSERVED_DENSITY = oF + 2
    PC % VELOCITY          = oF + [ 3, 4, 5 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( PC % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + PC % N_FIELDS_PC ) &
      = [ 'ComovingDensity ', &
          'ConservedDensity', &
          'Velocity_1      ', &
          'Velocity_2      ', &
          'Velocity_3      ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( PC % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = PC % N_VECTORS_TEMPLATE
    if ( PC % N_VECTORS == 0 ) &
      PC % N_VECTORS = oV + PC % N_VECTORS_PC

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( PC % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + PC % N_VECTORS_PC ) &
      = [ 'Velocity' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + PC % N_VECTORS_PC + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( PC % N_VECTORS ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( PC % VELOCITY )

  end subroutine InitializeBasics


end module ProtoCurrent_Form
