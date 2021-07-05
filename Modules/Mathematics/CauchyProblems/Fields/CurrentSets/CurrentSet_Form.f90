module CurrentSet_Form

  use Basics
  use FieldSets
  use Geometries

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_CS    = 0, &
      N_VECTORS_CS   = 0, &
      N_PRIMITIVE_CS = 0, &
      N_BALANCED_CS  = 0

  type, public, extends ( FieldSetForm ) :: CurrentSetForm
    !-- Fields and vectors
    integer ( KDI ) :: &
      N_FIELDS_CS    = N_FIELDS_CS, &
      N_VECTORS_CS   = N_VECTORS_CS
    !-- Defaults not generally used
    integer ( KDI ) :: &
      DENSITY_CS      = 0, &
      VELOCITY_CS_U_1 = 0, &
      VELOCITY_CS_U_2 = 0, &
      VELOCITY_CS_U_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      VELOCITY_CS_U = 0
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
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
!     procedure, public, pass :: &
!       SetVelocityLinear
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_CS
    procedure, public, pass :: &
      ComputeFromInitial
    procedure, public, pass :: &
      ComputeFromBalanced
    procedure, public, pass ( CS ) :: &
      ComputeFluxes
    procedure, public, pass ( CS ) :: &
      ComputeEigenspeeds
    final :: &
      Finalize
  end type CurrentSetForm

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
               ( CS, G, FieldOption, VectorOption, NameOption, UnitOption, &
                 VectorIndicesOption, iaPrimitiveOption, iaBalancedOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
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
    type ( MeasuredValueForm ), dimension ( :, : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( CS % Type  ==  '' ) &
      CS % Type  =  'a CurrentSet' 
    
    Name  =  'Currents'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    CS % Geometry  =>  G

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else

      CS % DENSITY_CS       =  1
      CS % VELOCITY_CS_U_1  =  2
      CS % VELOCITY_CS_U_2  =  3
      CS % VELOCITY_CS_U_3  =  4

      nFields  =  CS % N_FIELDS_CS  +  4

      CS % VELOCITY_CS_U  =  [ CS % VELOCITY_CS_U_1, &
                               CS % VELOCITY_CS_U_2, &
                               CS % VELOCITY_CS_U_3 ]

    end if

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
      Field ( CS % N_FIELDS_CS + 1 )  =  'Density'
      Field ( CS % N_FIELDS_CS + 2 )  =  'Velocity_U_1'
      Field ( CS % N_FIELDS_CS + 3 )  =  'Velocity_U_2'
      Field ( CS % N_FIELDS_CS + 4 )  =  'Velocity_U_3'
    end if !-- FieldOption

    !-- Units

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields, G % Atlas % nCharts ) )
    end if !-- UnitOption

    !-- Vector indices

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  CS % N_VECTORS_CS + 1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  CS % N_VECTORS_CS + 1
      allocate ( VectorIndices ( nVectors ) )
      call VectorIndices ( CS % N_VECTORS_CS + 1 ) % Initialize &
             ( CS % VELOCITY_CS_U )
    end if

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
      Vector ( CS % N_VECTORS_CS + 1 )  =  'Velocity'
    end if !-- FieldOption

    !-- Primitive fields

    if ( present ( iaPrimitiveOption ) ) then
      CS % nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( CS % iaPrimitive, source = iaPrimitiveOption )
    else
      CS % nPrimitive  =  CS % N_PRIMITIVE_CS + 1
      allocate ( CS % iaPrimitive ( CS % nPrimitive ) )
      CS % iaPrimitive  =  [ CS % DENSITY_CS ]
    end if !-- iaPrimitiveOption

    associate ( nP  =>  CS % nPrimitive )
    allocate ( CS % Primitive ( nP ) )
    do iP  =  1, nP
      iF  =  CS % iaPrimitive ( iP )
      CS % Primitive ( iP )  =  Field ( iF )
    end do !-- iP
    end associate !-- nP

    !-- Balanced fields

    if ( present ( iaBalancedOption ) ) then
      CS % nBalanced  =  size ( iaBalancedOption )
      allocate ( CS % iaBalanced, source = iaBalancedOption )
    else
      CS % nBalanced  =  CS % N_BALANCED_CS + 1
      allocate ( CS % iaBalanced ( CS % nBalanced ) )
      CS % iaBalanced  =  [ CS % DENSITY_CS ]
    end if !-- iaBalancedOption

    associate ( nB  =>  CS % nBalanced )
    allocate ( CS % Balanced ( nB ) )
    do iB  =  1, nB
      iF  =  CS % iaBalanced ( iB )
      CS % Balanced ( iB )  =  Field ( iF )
    end do !-- iB
    end associate !-- nP

    !-- FieldSet

    call CS % FieldSetForm % Initialize &
           ( G % Atlas, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             DeviceMemoryOption = G % DeviceMemory, &
             PinnedMemoryOption = G % PinnedMemory, &
             DevicesCommunicateOption = G % DevicesCommunicate, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndices, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_CS


!   subroutine SetVelocityLinear ( CS, Speed, Length )

!     class ( CurrentSetForm ), intent ( inout ) :: &
!       C
!     real ( KDR ), intent ( in ) :: &
!       Speed, &
!       Length

!     associate &
!       ( G  =>  CS % Geometry )
!     associate &
!       ( CV  =>  CS % Storage_FS % Storage % Value, &
!          GV  =>   G % Storage_FS % Storage % Value )
!     associate &
!       (    V_1  =>  CV ( :, CS % VELOCITY_CS_U_1 ), &
!            V_2  =>  CV ( :, CS % VELOCITY_CS_U_2 ), &
!            V_3  =>  CV ( :, CS % VELOCITY_CS_U_3 ), &
!            X_1  =>   GV ( :,  G % CENTER_U_1 ) )

!     V_1  =  Speed  *  ( X_1 / Length )
!     V_2  =  0.0_KDR
!     V_3  =  0.0_KDR
    
!     end associate !-- V_1, etc.
!     end associate !-- CV, etc.
!     end associate !-- G

!   end subroutine SetVelocityLinear


  subroutine SetStream ( S, CS )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ) :: &
      CS

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    if ( CS % DENSITY_CS  >  0 ) then
      allocate ( iaSelected ( 4 ) )
      iaSelected ( 1 )  =  CS % DENSITY_CS
      iaSelected ( 2 )  =  CS % VELOCITY_CS_U_1
      iaSelected ( 3 )  =  CS % VELOCITY_CS_U_2
      iaSelected ( 4 )  =  CS % VELOCITY_CS_U_3
    else
      allocate ( iaSelected ( 0 ) )
    end if

    call S % AddFieldSet ( CS, iaSelectedOption = iaSelected )

  end subroutine SetStream


  subroutine Show_CS ( FS )

    class ( CurrentSetForm ), intent ( in ) :: &
      FS

    call FS % FieldSetForm % Show ( )

    call Show ( FS %  nPrimitive,  'nPrimitive', FS % IGNORABILITY )
    call Show ( FS % iaPrimitive, 'iaPrimitive', FS % IGNORABILITY )
    call Show ( FS %   Primitive,   'Primitive', FS % IGNORABILITY )

    call Show ( FS %  nBalanced,  'nBalanced', FS % IGNORABILITY )
    call Show ( FS % iaBalanced, 'iaBalanced', FS % IGNORABILITY )
    call Show ( FS %   Balanced,   'Balanced', FS % IGNORABILITY )

  end subroutine Show_CS


  subroutine ComputeFromInitial ( CS )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS

  end subroutine ComputeFromInitial


  subroutine ComputeFromBalanced ( CS )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS

  end subroutine ComputeFromBalanced


  subroutine ComputeFluxes ( FS, CS, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( CurrentSetForm ), intent ( in ) :: &
      CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD  !-- iDimension
    
    integer ( KDI ) :: &
      iDensity

    if ( CS % DENSITY_CS > 0 ) then

      call Search ( CS % iaBalanced, CS % DENSITY_CS, iDensity )

      associate &
        ( FSS  =>  FS % Storage ( iC ), &
          CSS  =>  CS % Storage ( iC ) )
      associate &
        ( F_D      =>  FSS % Value ( :, iDensity ), &
            D      =>  CSS % Value ( :, CS % DENSITY_CS ), & 
            V_Dim  =>  CSS % Value ( :, CS % VELOCITY_CS_U ( iD ) ) ) 
 
      call ComputeFluxesKernel &
             ( D, V_Dim, F_D, UseDeviceOption = CS % DeviceMemory )
  
      end associate !-- F_D, etc.
      end associate !-- FSS, etc.
  
    end if !-- Density default

  end subroutine ComputeFluxes


  subroutine ComputeEigenspeeds ( FS, CS, iaEigenspeeds, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( CurrentSetForm ), intent ( in ) :: &
      CS
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaEigenspeeds
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    if ( CS % DENSITY_CS > 0 ) then

      associate &
        ( FSS  =>  FS % Storage ( iC ), &
          CSS  =>  CS % Storage ( iC ) )
      associate &
        ( EF_P    =>  FSS % Value ( :, iaEigenspeeds ( 1 ) ), &
          EF_M    =>  FSS % Value ( :, iaEigenspeeds ( 2 ) ), &
           V_Dim  =>  CSS % Value ( :, CS % VELOCITY_CS_U ( iD ) ) ) 
 
      call ComputeEigenspeedsKernel &
             ( V_Dim, EF_P, EF_M, UseDeviceOption = CS % DeviceMemory )
  
      end associate !-- EF_P, etc.
      end associate !-- FSS, etc.
  
    end if !-- Density default

  end subroutine ComputeEigenspeeds


  impure elemental subroutine Finalize ( CS )

    type ( CurrentSetForm ), intent ( inout ) :: &
      CS

    nullify ( CS % Geometry )

    if ( allocated ( CS % Balanced ) ) &
      deallocate ( CS % Balanced )
    if ( allocated ( CS % Primitive ) ) &
      deallocate ( CS % Primitive )
    if ( allocated ( CS % iaBalanced ) ) &
      deallocate ( CS % iaBalanced )
    if ( allocated ( CS % iaPrimitive ) ) &
      deallocate ( CS % iaPrimitive )

  end subroutine Finalize

  
end module CurrentSet_Form
