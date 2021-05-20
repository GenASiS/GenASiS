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
      DENSITY_DEFAULT = 0
    real ( KDR ), dimension ( 3 ) :: &
      VelocityDefault_U = 0.0_KDR
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
    procedure, public, pass ( CSC ) :: &
      SetStream
    procedure, private, pass :: &
      Show_FSC
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
          D
        real ( KDR ), intent ( in ) :: &
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
        real ( KDR ), intent ( in ) :: &
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

    ! CSC % FAST_EIGENSPEED_PLUS_U_1   =  1
    ! CSC % FAST_EIGENSPEED_PLUS_U_2   =  2
    ! CSC % FAST_EIGENSPEED_PLUS_U_3   =  3
    ! CSC % FAST_EIGENSPEED_MINUS_U_1  =  4
    ! CSC % FAST_EIGENSPEED_MINUS_U_2  =  5
    ! CSC % FAST_EIGENSPEED_MINUS_U_3  =  6

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      CSC % DENSITY_DEFAULT  =  1
      nFields  =  CSC % N_FIELDS_CS  +  1
    end if

   ! CSC % FAST_EIGENSPEED_PLUS_U  &
   !   =  [ CSC % FAST_EIGENSPEED_PLUS_U_1, &
   !        CSC % FAST_EIGENSPEED_PLUS_U_2, &
   !        CSC % FAST_EIGENSPEED_PLUS_U_3 ]
   ! CSC % FAST_EIGENSPEED_MINUS_U  &
   !   =  [ CSC % FAST_EIGENSPEED_MINUS_U_1, &
   !        CSC % FAST_EIGENSPEED_MINUS_U_2, &
   !        CSC % FAST_EIGENSPEED_MINUS_U_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
      Field ( CSC % N_FIELDS_CS + 1 )  =  'Density'
    end if !-- FieldOption

    ! Field ( 1 : CSC % N_FIELDS_CS ) &
    !   =  [ 'FastEigenspeedPlus_U_1 ', &
    !        'FastEigenspeedPlus_U_2 ', &
    !        'FastEigenspeedPlus_U_3 ', &
    !        'FastEigenspeedMinus_U_1', &
    !        'FastEigenspeedMinus_U_2', &
    !        'FastEigenspeedMinus_U_3' ]
          
    !-- Units

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields ) )
      if ( present ( DensityUnitOption ) ) &
        Unit ( CSC % DENSITY_DEFAULT )  =  DensityUnitOption
    end if !-- UnitOption

!    Unit ( CSC % FAST_EIGENSPEED_PLUS_U_1 : CSC % FAST_EIGENSPEED_PLUS_U_3 ) &
!      =  Velocity_U_Unit
!    Unit ( CSC % FAST_EIGENSPEED_MINUS_U_1 : CSC % FAST_EIGENSPEED_MINUS_U_3 ) &
!      =  Velocity_U_Unit

    !-- Vector indices

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  CSC % N_VECTORS_CS + 1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  CSC % N_VECTORS_CS
      allocate ( VectorIndices ( nVectors ) )
    end if

!    call VectorIndices ( 1 ) % Initialize ( CSC % FAST_EIGENSPEED_PLUS_U )
!    call VectorIndices ( 2 ) % Initialize ( CSC % FAST_EIGENSPEED_MINUS_U )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    ! Vector ( 1 : CSC % N_VECTORS_CS ) &
    !   = [ 'FastEigenspeedPlus ', &
    !       'FastEigenspeedMinus' ]

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

  end subroutine InitializeAllocate_CS


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


  subroutine Show_FSC ( FSC )

    class ( CurrentSet_C_Form ), intent ( in ) :: &
      FSC

    integer ( KDI ) :: &
      iF, &  !-- iField
      iS     !-- iSelected

    call FSC % FieldSet_C_Form % Show ( )

    call Show ( FSC %  nPrimitive,  'nPrimitive', FSC % IGNORABILITY )
    call Show ( FSC % iaPrimitive, 'iaPrimitive', FSC % IGNORABILITY )
    call Show ( FSC %   Primitive,   'Primitive', FSC % IGNORABILITY )

    call Show ( FSC %  nBalanced,  'nBalanced', FSC % IGNORABILITY )
    call Show ( FSC % iaBalanced, 'iaBalanced', FSC % IGNORABILITY )
    call Show ( FSC %   Balanced,   'Balanced', FSC % IGNORABILITY )

  end subroutine Show_FSC


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
    real ( KDR ) :: &
      V_Dim

    if ( CSC % DENSITY_DEFAULT > 0 ) then

      call Search ( CSC % iaBalanced, CSC % DENSITY_DEFAULT, iDensity )

      V_Dim  =  CSC % VelocityDefault_U ( iD )

      associate &
        ( FSS  =>  FSC % Storage_FSC % Storage, &
          CSS  =>  CSC % Storage_FSC % Storage, &
          DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory )
      associate &
        ( F_D  =>  FSS % Value ( :, iDensity ), &
            D  =>  CSS % Value ( :, CSC % DENSITY_DEFAULT ) ) 
 
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
    real ( KDR ) :: &
      V_Dim

    if ( CSC % DENSITY_DEFAULT > 0 ) then

      V_Dim  =  CSC % VelocityDefault_U ( iD )

      associate &
        ( FSS  =>  FSC % Storage_FSC % Storage, &
          DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory )
      associate &
        ( EF_P  =>  FSS % Value ( :, iaEigenspeeds ( 1 ) ), &
          EF_M  =>  FSS % Value ( :, iaEigenspeeds ( 2 ) ) ) 
 
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
