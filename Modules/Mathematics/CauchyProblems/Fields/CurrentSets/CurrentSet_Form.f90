module CurrentSet_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use Tally_CS__Form

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
    !-- Geometry
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
    !-- Tally
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      BoundaryFlux_SCG, &
      BoundaryFluence_SCG
    class ( Tally_CS_Form ), allocatable :: &
      TallyInterior, &
      TallyTotal, &
      TallyChange
    class ( Tally_CS_Element ), dimension ( : ), allocatable :: &
      TallyBoundary
  contains
    procedure, private, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_CS
    procedure, public, pass :: &
      ComputeFromInitial
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    procedure, public, pass :: &
      ComputeFromBalanced
    procedure, public, pass ( CS ) :: &
      ComputeEigenspeeds
    procedure, private, pass :: &
      AccumulateBoundaryFluence_SCG
    generic, public :: &
      AccumulateBoundaryFluence &
        => AccumulateBoundaryFluence_SCG
    procedure, public, pass :: &
      ComputeTally
    final :: &
      Finalize
    procedure, public, pass :: &
      AllocateBoundary_SCG
  end type CurrentSetForm

    private :: &
      ComputeEigenspeedsKernel

    interface

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
               ( CS, G, FieldOption, VectorOption, TallyVariableOption, &
                 NameOption, AllocateTallyOption, UnitOption, TallyUnitOption, &
                 VectorIndicesOption, iaPrimitiveOption, iaBalancedOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      TallyVariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      AllocateTallyOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      TallyUnitOption
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
      iB, &  !-- iBalanced, iBoundary
      iF, &  !-- iField
      iV, &  !-- iVector
      nFields, &
      nVectors
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      Unit
    logical ( KDL ) :: &
      AllocateTally
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
      CS % nPrimitive  =  CS % N_PRIMITIVE_CS + 4
      allocate ( CS % iaPrimitive ( CS % nPrimitive ) )
      CS % iaPrimitive  =  [ CS % DENSITY_CS, CS % VELOCITY_CS_U ]
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

    !-- Tally

    AllocateTally = .true.
    if ( present ( AllocateTallyOption ) ) &
      AllocateTally = AllocateTallyOption

    if ( .not. allocated ( CS % TallyInterior ) .and. AllocateTally ) then

      allocate ( CS % TallyInterior )
      allocate ( CS % TallyTotal )
      allocate ( CS % TallyChange )
      allocate ( CS % TallyBoundary ( CS % nBoundaries ) )
      do iB  =  1,  CS % nBoundaries 
        allocate ( CS % TallyBoundary ( iB ) % Element )
      end do !-- iB

      associate &
        ( TI   =>  CS % TallyInterior, &
          TT   =>  CS % TallyTotal, &
          TC   =>  CS % TallyChange, &
          TB  =>  CS % TallyBoundary ( : ) )
      call TI % Initialize &
             ( CS, G, CS % iaBalanced, VariableOption = TallyVariableOption, &
               UnitOption = TallyUnitOption )
      call TC % Initialize &
             ( CS, G, CS % iaBalanced, VariableOption = TallyVariableOption, &
               UnitOption = TallyUnitOption )
      call TT % Initialize &
             ( CS, G, CS % iaBalanced, VariableOption = TallyVariableOption, &
               UnitOption = TallyUnitOption )
      do iB  =  1,  CS % nBoundaries
        call TB ( iB ) % Element % Initialize &
               ( CS, G, CS % iaBalanced, VariableOption = TallyVariableOption, &
                 UnitOption = TallyUnitOption )
      end do !-- iB
      end associate !-- TI, etc.

      call CS % AllocateBoundary_SCG ( nT = CS % nBalanced )

    end if !-- AllocateTally

  end subroutine InitializeAllocate_CS


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

    integer ( KDI ) :: &
      iFS

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


  subroutine ComputeFromPrimitive ( FS_CS, CS )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS_CS
    class ( CurrentSetForm ), intent ( in ) :: &
      CS

  end subroutine ComputeFromPrimitive


  subroutine ComputeFromBalanced ( CS, T_Option )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

  end subroutine ComputeFromBalanced


  subroutine ComputeEigenspeeds ( ES, CS, FS_CS, iaEigenspeeds, iC, iD )

    class ( FieldSetForm ), intent ( inout ) :: &
      ES
    class ( CurrentSetForm ), intent ( in ) :: &
      CS
    class ( FieldSetForm ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaEigenspeeds
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    if ( CS % DENSITY_CS > 0 ) then

      associate &
        ( ESS  =>  ES    % Storage ( iC ), &
          CSS  =>  FS_CS % Storage ( iC ) )
      associate &
        ( EF_P    =>  ESS % Value ( :, iaEigenspeeds ( 1 ) ), &
          EF_M    =>  ESS % Value ( :, iaEigenspeeds ( 2 ) ), &
           V_Dim  =>  CSS % Value ( :, CS % VELOCITY_CS_U ( iD ) ) ) 
 
      call ComputeEigenspeedsKernel &
             ( V_Dim, EF_P, EF_M, UseDeviceOption = CS % DeviceMemory )
  
      end associate !-- EF_P, etc.
      end associate !-- FSS, etc.
  
    end if !-- Density default

  end subroutine ComputeEigenspeeds


  subroutine AccumulateBoundaryFluence_SCG ( CS, Factor )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    real ( KDR ), intent ( in ) :: &
      Factor  !-- dT * Weight_RK

    associate &
      ( BFc  =>  CS % BoundaryFluence_SCG, &
        BFx  =>  CS % BoundaryFlux_SCG )
        
    call BFc % MultiplyAdd ( BFx, Factor )

    end associate !-- BFc, etc.

  end subroutine AccumulateBoundaryFluence_SCG


  subroutine ComputeTally ( CS, ChangeOption, IgnorabilityOption )
    
    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    logical ( KDL ), intent ( in ), optional :: &
      ChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    
    integer ( KDI ) :: &
      iB  !-- iBoundary
    real ( KDR ), dimension ( : ), allocatable :: &
      OldTotal
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      OldBoundary
    logical ( KDL ) :: &
      Change

    call Show ( 'Computing Tally', IgnorabilityOption )

    Change = .true.
    if ( present ( ChangeOption ) ) &
      Change  =  ChangeOption

    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( nI  =>  CS % TallyInterior % nIntegrals )

    allocate ( OldTotal ( nI ) )
    OldTotal  =  CS % TallyTotal % Value

    allocate ( OldBoundary ( CS % nBoundaries ) )
    do iB  =  1,  CS % nBoundaries
      call OldBoundary ( iB ) % Initialize ( nI )
      OldBoundary ( iB ) % Value  &
        =  CS % TallyBoundary ( iB ) % Element % Value
    end do !-- iB

    !-- Interior

    call CS % TallyInterior % ComputeInterior ( CS )

    !-- Boundary

    call CS % BoundaryFluence_SCG % UpdateHost ( )

    associate ( iExtent => 1 )  !-- only boundary for Atlas_SCG
    call CS % TallyBoundary ( iExtent ) % Element &
           % ComputeBoundary ( CS, CS % BoundaryFluence_SCG )
    end associate !-- iExtent

    call CS % BoundaryFluence_SCG % Clear ( )

    !-- Total

    CS % TallyTotal % Value  =  CS % TallyInterior % Value

    do iB  =  1,  CS % nBoundaries

      CS % TallyTotal % Value &
        =  CS % TallyTotal % Value  &
           +  CS % TallyBoundary ( iB ) % Element % Value

    end do !-- iB

    !-- Change

    if ( Change ) then
      CS % TallyChange % Value &
        =  CS % TallyChange % Value &
           +  ( CS % TallyTotal % Value - OldTotal )
    end if
  
    !-- Display

    call CS % TallyInterior % Show &
           ( 'Interior Tally ' // trim ( CS % Name ), IgnorabilityOption )

    do iB  =  1,  CS % nBoundaries
      call CS % TallyBoundary ( iB ) % Element % Show &
             ( 'Boundary ' // trim ( CS % Boundaries ( 1 ) % Boundary ( iB ) ) &
               // ' Tally ' // trim ( CS % Name ), IgnorabilityOption )
    end do

    call CS % TallyTotal % Show &
           ( 'Total Tally ' // trim ( CS % Name ), IgnorabilityOption )

    if ( Change ) then
      call CS % TallyChange % Show &
             ( 'Change in Total Tally ' // trim ( CS % Name ), &
               IgnorabilityOption )
    end if

    end associate !-- nI

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'CurrentSet_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeTally', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine ComputeTally
  
  
  impure elemental subroutine Finalize ( CS )

    type ( CurrentSetForm ), intent ( inout ) :: &
      CS

    nullify ( CS % Geometry )

    if ( allocated ( CS % TallyBoundary ) ) &
      deallocate ( CS % TallyBoundary )
    if ( allocated ( CS % TallyChange ) ) &
      deallocate ( CS % TallyChange )
    if ( allocated ( CS % TallyTotal ) ) &
      deallocate ( CS % TallyTotal )
    if ( allocated ( CS % TallyInterior ) ) &
      deallocate ( CS % TallyInterior )
    if ( allocated ( CS % BoundaryFluence_SCG ) ) &
      deallocate ( CS % BoundaryFluence_SCG )
    if ( allocated ( CS % BoundaryFlux_SCG ) ) &
      deallocate ( CS % BoundaryFlux_SCG )
    if ( allocated ( CS % Balanced ) ) &
      deallocate ( CS % Balanced )
    if ( allocated ( CS % Primitive ) ) &
      deallocate ( CS % Primitive )
    if ( allocated ( CS % iaBalanced ) ) &
      deallocate ( CS % iaBalanced )
    if ( allocated ( CS % iaPrimitive ) ) &
      deallocate ( CS % iaPrimitive )

  end subroutine Finalize

  
  subroutine AllocateBoundary_SCG ( CS, nT )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    integer ( KDI ), intent ( in ) :: &
      nT  !-- nTally

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- iDimension, etc.
      iT             !-- iTally
    integer ( KDI ), dimension ( 3 ) :: &
      nSurface

    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS )
    associate &
      ( nD   =>  C % nDimensions, &
        nF   =>  C % Connectivity % nFaces, &
        iaI  =>  C % Connectivity % iaInner ( : ), &
        iaO  =>  C % Connectivity % iaOuter ( : ) )

    allocate &
      ( CS % BoundaryFluence_SCG ( nT, nF ), &
        CS % BoundaryFlux_SCG ( nT, nF ) )
    associate &
      ( BFc  =>  CS % BoundaryFluence_SCG, &
        BFx  =>  CS % BoundaryFlux_SCG )
        
    do iD  =  1, nD
      jD  =  mod ( iD, 3 ) + 1
      kD  =  mod ( jD, 3 ) + 1
      nSurface ( iD )  =  1
      nSurface ( jD )  =  C % nCellsBrick ( jD ) 
      nSurface ( kD )  =  C % nCellsBrick ( kD )
      do iT  =  1, nT
        call BFc ( iT, iaI ( iD ) ) &
               % Initialize ( nSurface, ClearOption = .true. )
        call BFc ( iT, iaO ( iD ) ) &
               % Initialize ( nSurface, ClearOption = .true. )
        call BFx ( iT, iaI ( iD ) ) &
               % Initialize ( nSurface, ClearOption = .true. )
        call BFx ( iT, iaO ( iD ) ) &
               % Initialize ( nSurface, ClearOption = .true. )
        if ( CS % DeviceMemory ) then
          call BFc ( iT, iaI ( iD ) ) % AllocateDevice ( )
          call BFc ( iT, iaO ( iD ) ) % AllocateDevice ( )
          call BFx ( iT, iaI ( iD ) ) % AllocateDevice ( )
          call BFx ( iT, iaO ( iD ) ) % AllocateDevice ( )
        end if
      end do !-- iT
    end do !-- iD

    end associate !-- BFc, etc.
    end associate !-- nD, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine AllocateBoundary_SCG


end module CurrentSet_Form
