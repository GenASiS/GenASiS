module CurrentSet_C__Form

  !-- CurrentSet_Chart_Form

  use Basics
  use FieldSets
  use Streams
  use Geometries

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_CS    = 0, &
      N_VECTORS_CS   = 0, &
      N_PRIMITIVE_CS = 0, &
      N_BALANCED_CS  = 0

  type, public, extends ( FieldSet_C_Form ) :: CurrentSet_C_Form
    !-- Fields and vectors
    integer ( KDI ) :: &
      N_FIELDS_CS    = N_FIELDS_CS, &
      N_VECTORS_CS   = N_VECTORS_CS
    !-- Defaults not generally used
    integer ( KDI ) :: &
      DENSITY_DEFAULT      = 0, &
      VELOCITY_DEFAULT_U_1 = 0, &
      VELOCITY_DEFAULT_U_2 = 0, &
      VELOCITY_DEFAULT_U_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      VELOCITY_DEFAULT_U = 0
    !-- Primtive and Balanced
    integer ( KDI ) :: &
      nPrimitive = 0, &
      nBalanced  = 0
    integer ( KDI ) :: &
      N_PRIMITIVE_CS = N_PRIMITIVE_CS, &
      N_BALANCED_CS  = N_BALANCED_CS
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
    character ( LDL ), dimension ( : ), allocatable :: &
      Primitive, &
      Balanced
    class ( Geometry_F_C_Form ), pointer :: &
      Geometry_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
    procedure, public, pass :: &
      SetVelocityConstant
    procedure, public, pass :: &
      SetVelocityLinear
    procedure, public, pass ( CSC ) :: &
      SetStream
    procedure, private, pass :: &
      Show_FS
    procedure, public, pass :: &
      ComputeFromInitial
    procedure, public, pass :: &
      ComputeFromConserved
    procedure, public, pass ( CSC ) :: &
      ComputeFluxes
    procedure, public, pass ( CSC ) :: &
      ComputeEigenspeeds
    final :: &
      Finalize
  end type CurrentSet_C_Form

    private :: &
      ComputeFluxesKernel, &
      ComputeEigenspeedsKernel

    interface

      module subroutine ComputeFluxesKernel ( D, V_Dim, F_D, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          V_Dim
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          F_D
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeFluxesKernel

      module subroutine ComputeEigenspeedsKernel &
               ( V_Dim, EF_P, EF_M, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          V_Dim
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          EF_P, EF_M
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeEigenspeedsKernel

    end interface

contains


  subroutine InitializeAllocate_CS &
               ( CSC, GC, FieldOption, VectorOption, NameOption, UnitOption, &
                 DensityUnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC
    class ( Geometry_F_C_Form ), intent ( in ), target :: &
      GC
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      DensityUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iP, &  !-- iPrimitive
      iB, &  !-- iBalanced
      iF, &  !-- iField
      iV, &  !-- iVector
      nFields, &
      nVectors
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( CSC % Type  ==  '' ) &
      CSC % Type  =  'a CurrentSet_C' 
    
    Name  =  'CurrentSet'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    CSC % Geometry_C  =>  GC

    associate &
      ( DeviceMemory  =>  GC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  GC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  GC % GhostExchange_FSC % DevicesCommunicate ) 

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else

      CSC % DENSITY_DEFAULT       =  1
      CSC % VELOCITY_DEFAULT_U_1  =  2
      CSC % VELOCITY_DEFAULT_U_2  =  3
      CSC % VELOCITY_DEFAULT_U_3  =  4

      nFields  =  CSC % N_FIELDS_CS  +  4

      CSC % VELOCITY_DEFAULT_U  =  [ CSC % VELOCITY_DEFAULT_U_1, &
                                     CSC % VELOCITY_DEFAULT_U_2, &
                                     CSC % VELOCITY_DEFAULT_U_3 ]

    end if

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
      Field ( CSC % N_FIELDS_CS + 1 )  =  'Density'
      Field ( CSC % N_FIELDS_CS + 2 )  =  'Velocity_U_1'
      Field ( CSC % N_FIELDS_CS + 3 )  =  'Velocity_U_2'
      Field ( CSC % N_FIELDS_CS + 4 )  =  'Velocity_U_3'
    end if !-- FieldOption

    !-- Units

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields ) )
      if ( present ( DensityUnitOption ) ) &
        Unit ( CSC % DENSITY_DEFAULT )  =  DensityUnitOption
    end if !-- UnitOption

    !-- Vector indices

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  CSC % N_VECTORS_CS + 1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  CSC % N_VECTORS_CS + 1
      allocate ( VectorIndices ( nVectors ) )
      call VectorIndices ( CSC % N_VECTORS_CS + 1 ) % Initialize &
             ( CSC % VELOCITY_DEFAULT_U )
    end if

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
      Vector ( CSC % N_VECTORS_CS + 1 )  =  'Velocity'
    end if !-- FieldOption

    !-- Primitive fields

    if ( present ( iaPrimitiveOption ) ) then
      CSC % nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( CSC % iaPrimitive, source = iaPrimitiveOption )
    else
      CSC % nPrimitive  =  CSC % N_PRIMITIVE_CS + 1
      allocate ( CSC % iaPrimitive ( CSC % nPrimitive ) )
      CSC % iaPrimitive  =  [ CSC % DENSITY_DEFAULT ]
    end if !-- iaPrimitiveOption

    associate ( nP  =>  CSC % nPrimitive )
    allocate ( CSC % Primitive ( nP ) )
    do iP  =  1, nP
      iF  =  CSC % iaPrimitive ( iP )
      CSC % Primitive ( iP )  =  Field ( iF )
    end do !-- iP
    end associate !-- nP

    !-- Balanced fields

    if ( present ( iaBalancedOption ) ) then
      CSC % nBalanced  =  size ( iaBalancedOption )
      allocate ( CSC % iaBalanced, source = iaBalancedOption )
    else
      CSC % nBalanced  =  CSC % N_BALANCED_CS + 1
      allocate ( CSC % iaBalanced ( CSC % nBalanced ) )
      CSC % iaBalanced  =  [ CSC % DENSITY_DEFAULT ]
    end if !-- iaBalancedOption

    associate ( nB  =>  CSC % nBalanced )
    allocate ( CSC % Balanced ( nB ) )
    do iB  =  1, nB
      iF  =  CSC % iaBalanced ( iB )
      CSC % Balanced ( iB )  =  Field ( iF )
    end do !-- iB
    end associate !-- nP

    !-- FieldSet

    call CSC % FieldSet_C_Form % Initialize &
           ( GC % Chart, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndices, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    end associate !-- DeviceMemory, etc.

  end subroutine InitializeAllocate_CS


  subroutine SetVelocityConstant ( CSC, Direction, Speed )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Direction
    real ( KDR ), intent ( in ) :: &
      Speed

    associate &
      ( CSV  =>  CSC % Storage_FSC % Storage % Value )
    associate &
      (     V_1  =>  CSV ( :, CSC % VELOCITY_DEFAULT_U_1 ), &
            V_2  =>  CSV ( :, CSC % VELOCITY_DEFAULT_U_2 ), &
            V_3  =>  CSV ( :, CSC % VELOCITY_DEFAULT_U_3 ), &
            K    =>  Direction, &
        Abs_K    =>  sqrt ( dot_product ( Direction, Direction ) ) )

    V_1  =  Speed  *  K ( 1 )  /  Abs_K
    V_2  =  Speed  *  K ( 2 )  /  Abs_K
    V_3  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- V_1, etc.
    end associate !-- CSV

  end subroutine SetVelocityConstant


  subroutine SetVelocityLinear ( CSC, Speed, Length )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC
    real ( KDR ), intent ( in ) :: &
      Speed, &
      Length

    associate &
      ( GC  =>  CSC % Geometry_C )
    associate &
      ( CSV  =>  CSC % Storage_FSC % Storage % Value, &
         GV  =>   GC % Storage_FSC % Storage % Value )
    associate &
      (    V_1  =>  CSV ( :, CSC % VELOCITY_DEFAULT_U_1 ), &
           V_2  =>  CSV ( :, CSC % VELOCITY_DEFAULT_U_2 ), &
           V_3  =>  CSV ( :, CSC % VELOCITY_DEFAULT_U_3 ), &
           X_1  =>   GV ( :,  GC % CENTER_U_1 ) )

    V_1  =  Speed  *  ( X_1 / Length )
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR
    
    end associate !-- V_1, etc.
    end associate !-- CSV, etc.
    end associate !-- GC

  end subroutine SetVelocityLinear


  subroutine SetStream ( SC, CSC )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( CurrentSet_C_Form ), intent ( in ) :: &
      CSC

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    if ( CSC % DENSITY_DEFAULT  >  0 ) then
      allocate ( iaSelected ( 1 ) )
      iaSelected ( 1 )  =  CSC % DENSITY_DEFAULT
    else
      allocate ( iaSelected ( 0 ) )
    end if

    call SC % AddFieldSet &
           ( CSC, &
             iaSelectedOption  =  iaSelected )

  end subroutine SetStream


  subroutine Show_FS ( FSC )

    class ( CurrentSet_C_Form ), intent ( in ) :: &
      FSC

    call FSC % FieldSet_C_Form % Show ( )

    call Show ( FSC %  nPrimitive,  'nPrimitive', FSC % IGNORABILITY )
    call Show ( FSC % iaPrimitive, 'iaPrimitive', FSC % IGNORABILITY )
    call Show ( FSC %   Primitive,   'Primitive', FSC % IGNORABILITY )

    call Show ( FSC %  nBalanced,  'nBalanced', FSC % IGNORABILITY )
    call Show ( FSC % iaBalanced, 'iaBalanced', FSC % IGNORABILITY )
    call Show ( FSC %   Balanced,   'Balanced', FSC % IGNORABILITY )

  end subroutine Show_FS


  subroutine ComputeFromInitial ( CSC )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC

  end subroutine ComputeFromInitial


  subroutine ComputeFromConserved ( CSC )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC

  end subroutine ComputeFromConserved


  subroutine ComputeFluxes ( FSC, CSC, iD )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( CurrentSet_C_Form ), intent ( in ) :: &
      CSC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension
    
    integer ( KDI ) :: &
      iDensity

    if ( CSC % DENSITY_DEFAULT > 0 ) then

      call Search ( CSC % iaBalanced, CSC % DENSITY_DEFAULT, iDensity )

      associate &
        ( FSS  =>  FSC % Storage_FSC % Storage, &
          CSS  =>  CSC % Storage_FSC % Storage, &
          DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory )
      associate &
        ( F_D      =>  FSS % Value ( :, iDensity ), &
            D      =>  CSS % Value ( :, CSC % DENSITY_DEFAULT ), & 
            V_Dim  =>  CSS % Value ( :, CSC % VELOCITY_DEFAULT_U ( iD ) ) ) 
 
      call ComputeFluxesKernel &
             ( D, V_Dim, F_D, UseDeviceOption = DeviceMemory )
  
      end associate !-- F_D, etc.
      end associate !-- FSS, etc.

    end if !-- Density default

  end subroutine ComputeFluxes


  subroutine ComputeEigenspeeds ( FSC, CSC, iaEigenspeeds, iD )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( CurrentSet_C_Form ), intent ( in ) :: &
      CSC
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaEigenspeeds
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension

    if ( CSC % DENSITY_DEFAULT > 0 ) then

      associate &
        ( FSS  =>  FSC % Storage_FSC % Storage, &
          CSS  =>  CSC % Storage_FSC % Storage, &
          DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory )
      associate &
        ( EF_P    =>  FSS % Value ( :, iaEigenspeeds ( 1 ) ), &
          EF_M    =>  FSS % Value ( :, iaEigenspeeds ( 2 ) ), &
           V_Dim  =>  CSS % Value ( :, CSC % VELOCITY_DEFAULT_U ( iD ) ) ) 
 
      call ComputeEigenspeedsKernel &
             ( V_Dim, EF_P, EF_M, UseDeviceOption = DeviceMemory )
  
      end associate !-- EF_P, etc.
      end associate !-- FSS, etc.

    end if !-- Density default

  end subroutine ComputeEigenspeeds


  impure elemental subroutine Finalize ( CSC )

    type ( CurrentSet_C_Form ), intent ( inout ) :: &
      CSC

    nullify ( CSC % Geometry_C )

    if ( allocated ( CSC % Balanced ) ) &
      deallocate ( CSC % Balanced )
    if ( allocated ( CSC % Primitive ) ) &
      deallocate ( CSC % Primitive )
    if ( allocated ( CSC % iaBalanced ) ) &
      deallocate ( CSC % iaBalanced )
    if ( allocated ( CSC % iaPrimitive ) ) &
      deallocate ( CSC % iaPrimitive )

  end subroutine Finalize

  
end module CurrentSet_C__Form
