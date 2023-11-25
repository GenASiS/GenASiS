module Fluid_D__Form

  !-- Fluid_Dust__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_D    = 9, &
      N_VECTORS_D   = 2, &
      N_PRIMITIVE_D = 4, &
      N_BALANCED_D  = 4

  type, public, extends ( CurrentSetForm ) :: Fluid_D_Form
    integer ( KDI ) :: &
      N_FIELDS_D    = N_FIELDS_D, &
      N_VECTORS_D   = N_VECTORS_D, &
      N_PRIMITIVE_D = N_PRIMITIVE_D, &
      N_BALANCED_D  = N_BALANCED_D
    integer ( KDI ) :: &
      BARYON_DENSITY_C = 0, &  !-- Comoving
      BARYON_DENSITY_B = 0, &  !-- Balanced
      BARYON_MASS      = 0
    integer ( KDI ) :: &
      VELOCITY_U_1 = 0, &
      VELOCITY_U_2 = 0, &
      VELOCITY_U_3 = 0, &
      MOMENTUM_DENSITY_D_1 = 0, &
      MOMENTUM_DENSITY_D_2 = 0, &
      MOMENTUM_DENSITY_D_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      VELOCITY_U, &
      MOMENTUM_DENSITY_D
    integer ( KDI ) :: &
      iTimer_CFB = 0
    real ( KDR ) :: &
      BaryonMass, &
      BaryonDensityMin
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    procedure, private, pass :: &
      SetBaryonDensityMinValue
    procedure, private, pass :: &
      SetBaryonDensityMinFind
    generic, public :: &
      SetBaryonDensityMin => SetBaryonDensityMinValue, SetBaryonDensityMinFind
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      ComputeFromInitial
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    procedure, public, pass :: &
      ComputeFromBalanced
    procedure, public, pass ( CS ) :: &
      ComputeEigenspeeds
    final :: &
      Finalize
  end type Fluid_D_Form

    private :: &
      Compute_M_Kernel, &
      Compute_D_S_G_Kernel, &
      Compute_N_V_G_Kernel, &
      Compute_ES_G_Kernel

  interface
  
    module subroutine Compute_M_Kernel &
             ( M_Ref, M, UseDeviceOption )
      !-- Compute_BaryonMass_Kernel
      use Basics
      implicit none
      real ( KDR ), intent ( in ) :: &
        M_Ref
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        M
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_M_Kernel

    module subroutine Compute_D_S_G_Kernel & 	 	 
             ( N, V_1, V_2, V_3, M, M_DD_11, M_DD_22, M_DD_33, N_Min, &
               D, S_1, S_2, S_3, UseDeviceOption )
      !-- Compute_DensityB_Momentum_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
        N, & 	 	 
        V_1, V_2, V_3
      real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
        M, & 	 	 
        M_DD_11, M_DD_22, M_DD_33
      real ( KDR ), intent ( in ) :: &
        N_Min
      real ( KDR ), dimension ( : ), intent ( out ) :: & 	 	 
        D, & 	 	 
        S_1, S_2, S_3 	 	 
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_D_S_G_Kernel 	 	 

    module subroutine Compute_N_V_G_Kernel &
             ( D, S_1, S_2, S_3, M, M_UU_11, M_UU_22, M_UU_33, N_Min, &
               N, V_1, V_2, V_3, UseDeviceOption )
      !-- Compute_DensityC_Velocity_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        D, &
        S_1, S_2, S_3
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        M_UU_11, M_UU_22, M_UU_33
      real ( KDR ), intent ( in ) :: &
        N_Min
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        N, &
        V_1, V_2, V_3
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_N_V_G_Kernel

    module subroutine Compute_ES_G_Kernel &
             ( V_Dim, EF_P, EF_M, UseDeviceOption )
      !-- Compute_EigenspeedSet_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        V_Dim
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        EF_P, EF_M
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Compute_ES_G_Kernel

  end interface


contains


  subroutine InitializeAllocate_F &
               ( F, G, Units_F, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Units_F_Form ), dimension ( : ), intent ( in ) :: &
      Units_F
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
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
      iV, &  !-- iVector
      iP, &  !-- iPrimitive
      iB, &  !-- iBalanced
      iC, &  !-- iChart
      oF, &  !-- oField
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oB, &  !-- oBalanced
      nFields, &
      nVectors, &
      nPrimitive, &
      nBalanced
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( F % Type  ==  '' ) &
      F % Type  =  'a Fluid_D' 
    
    Name  =  'Fluid'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  F % N_FIELDS_CS

    F % BARYON_DENSITY_C      =  oF + 1
    F % BARYON_DENSITY_B      =  oF + 2
    F % BARYON_MASS           =  oF + 3
    F % VELOCITY_U_1          =  oF + 4
    F % VELOCITY_U_2          =  oF + 5
    F % VELOCITY_U_3          =  oF + 6
    F % MOMENTUM_DENSITY_D_1  =  oF + 7
    F % MOMENTUM_DENSITY_D_2  =  oF + 8
    F % MOMENTUM_DENSITY_D_3  =  oF + 9

    nFields  =  oF  +  F % N_FIELDS_D
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    F % VELOCITY_U          =  [ F % VELOCITY_U_1, &
                                 F % VELOCITY_U_2, &
                                 F % VELOCITY_U_3 ]
    F % MOMENTUM_DENSITY_D  =  [ F % MOMENTUM_DENSITY_D_1, &
                                 F % MOMENTUM_DENSITY_D_2, &
                                 F % MOMENTUM_DENSITY_D_3 ]
 
    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + F % N_FIELDS_D ) &
      = [ 'BaryonDensity_C    ', &
          'BaryonDensity_B    ', &
          'BaryonMass         ', &
          'Velocity_U_1       ', &
          'Velocity_U_2       ', &
          'Velocity_U_3       ', &
          'MomentumDensity_D_1', &
          'MomentumDensity_D_2', &
          'MomentumDensity_D_3' ]
          
    !-- Units

    associate ( nC  =>  G % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( F % BARYON_DENSITY_C, iC ) &
        =  Units_F ( iC ) % NumberDensity
      FieldUnit ( F % BARYON_DENSITY_B, iC ) &
        =  Units_F ( iC ) % SqrtDet_M  *  Units_F ( iC ) % NumberDensity
      FieldUnit ( F % BARYON_MASS, iC ) &
        =  Units_F ( iC ) % BaryonMass
      FieldUnit ( F % VELOCITY_U_1, iC ) &
        =  Units_F ( iC ) % Velocity_U ( 1 )
      FieldUnit ( F % VELOCITY_U_2, iC ) &
        =  Units_F ( iC ) % Velocity_U ( 2 )
      FieldUnit ( F % VELOCITY_U_3, iC ) &
        =  Units_F ( iC ) % Velocity_U ( 3 )
      FieldUnit ( F % MOMENTUM_DENSITY_D_1, iC ) &
        =  Units_F ( iC ) % MomentumDensity_D ( 1 )
      FieldUnit ( F % MOMENTUM_DENSITY_D_2, iC ) &
        =  Units_F ( iC ) % MomentumDensity_D ( 2 )
      FieldUnit ( F % MOMENTUM_DENSITY_D_3, iC ) &
        =  Units_F ( iC ) % MomentumDensity_D ( 3 )
    end do !-- iC

    end associate !-- nC

    !-- Vector indices

    oV  =  F % N_VECTORS_CS

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  oV  +  F % N_VECTORS_D  +  1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  oV  +  F % N_VECTORS_D
      allocate ( VectorIndices ( nVectors ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( F % VELOCITY_U )
    call VectorIndices ( oV + 2 ) % Initialize ( F % MOMENTUM_DENSITY_D )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    Vector ( oV  +  1 : oV  +  F % N_VECTORS_D ) &
      = [ 'Velocity_U       ', &
          'MomentumDensity_D' ]

    !-- Primitive fields

    oP  =  F % N_PRIMITIVE_CS

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  F % N_PRIMITIVE_D
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  F % N_PRIMITIVE_D )  &
      =  [ F % BARYON_DENSITY_C, F % VELOCITY_U ]

    !-- Balanced fields

    oB  =  F % N_BALANCED_CS

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oB  +  F % N_BALANCED_D
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaPrimitiveOption

    iaBalanced ( oB  +  1 : oB  +  F % N_BALANCED_D )  &
      =  [ F % BARYON_DENSITY_B, F % MOMENTUM_DENSITY_D ]

    !-- CurrentSet

    call F % CurrentSetForm % Initialize &
           ( G, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             UnitOption = FieldUnit, &
             VectorIndicesOption = VectorIndices, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    !-- Parameters

    if ( Units_F ( 1 ) % BaryonMass  ==  UNIT % IDENTITY ) then
      F % BaryonMass  =  1.0_KDR
    else
      F % BaryonMass  =  CONSTANT % ATOMIC_MASS_UNIT
    end if

    F % BaryonDensityMin  =  1.0e2_KDR * sqrt ( tiny ( 0.0_KDR ) )
    call PROGRAM_HEADER % GetParameter &
           ( F % BaryonDensityMin, 'BaryonDensityMin' )

  end subroutine InitializeAllocate_F


  subroutine SetBaryonDensityMinValue ( F, BaryonDensityMin )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      BaryonDensityMin

    F % BaryonDensityMin  =  BaryonDensityMin

    call Show ( 'Setting BaryonDensityMin of a Fluid', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % BaryonDensityMin, &
                F % Unit ( F % BARYON_DENSITY_C, 1 ), 'BaryonDensityMin', &
                F % IGNORABILITY + 1 )

  end subroutine SetBaryonDensityMinValue


  subroutine SetBaryonDensityMinFind ( F )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
 
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        FV  =>  F % Storage_GS % Value )

    call CO % Initialize &
           ( C % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )

    associate &
      ( My_N_Min => CO % Outgoing % Value ( 1 ), &
           N_Min => CO % Incoming % Value ( 1 ) )
 
    My_N_Min  =  minval ( FV ( :, F % BARYON_DENSITY_C ) )

    call CO % Reduce ( REDUCTION % MIN )

    F % BaryonDensityMin  =  N_Min

    call Show ( 'Setting BaryonDensityMin of a Fluid', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % BaryonDensityMin, &
                F % Unit ( F % BARYON_DENSITY_C, 1 ), 'BaryonDensityMin', &
                F % IGNORABILITY + 1 )

    end associate !-- My_N_Min, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetBaryonDensityMinFind


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Fluid_D_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % BARYON_DENSITY_C, CS % VELOCITY_U ] )

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( Fluid_D_Form ), intent ( in ) :: &
      FS

    call FS % CurrentSetForm % Show ( )

    call Show ( FS % BaryonMass, &
                FS % Unit ( FS % BARYON_MASS, 1 ), 'BaryonMass', &
                FS % IGNORABILITY )
    call Show ( FS % BaryonDensityMin, &
                FS % Unit ( FS % BARYON_DENSITY_C, 1 ), 'BaryonDensityMin', &
                FS % IGNORABILITY )

  end subroutine Show_FS


  subroutine ComputeFromInitial ( CS )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      CS

    call Show ( 'ComputeFromInitial', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    call CS % ComputeFromPrimitive ( CS )

    select type ( G  =>  CS % Geometry )
      class is ( Gravitation_N_H_Form )
    call G % Solve &
           ( CS, &
             iBaryonMass = CS % BARYON_MASS, &
             iBaryonDensity = CS % BARYON_DENSITY_B )
    end select !-- G

  end subroutine ComputeFromInitial


  subroutine ComputeFromPrimitive ( FS_CS, CS )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_CS
    class ( Fluid_D_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromPrimitive', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, CS % Atlas % nCharts

      associate &
        (   CSV  =>  FS_CS % Storage ( iC ) % Value, &
          M_Ref  =>  CS    % BaryonMass, &
          N_Min  =>  CS    % BaryonDensityMin )
      associate &
        ( M    =>  CSV ( :, CS % BARYON_MASS ), &
          N    =>  CSV ( :, CS % BARYON_DENSITY_C ), &
          V_1  =>  CSV ( :, CS % VELOCITY_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_U_3 ), &
          D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
          S_1  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ) )
   
      call Compute_M_Kernel &
             ( M_Ref, M, UseDeviceOption = CS % DeviceMemory )

      select type ( G  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ) )

        call Compute_D_S_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, M, M_DD_11, M_DD_22, M_DD_33, N_Min, &
                 D, S_1, S_2, S_3, UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_D__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromPrimitive', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      end associate !-- M, etc.
      end associate !-- CSV, etc.

    end do !-- iC

  end subroutine ComputeFromPrimitive


  subroutine ComputeFromBalanced ( CS, T_Option )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC
    type ( TimerForm ), pointer :: &
      T_K

    call Show ( 'ComputeFromBalanced', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    if ( present ( T_Option ) ) then
      T_K  =>  PROGRAM_HEADER % Timer &
                 ( Handle = CS % iTimer_CFB, &
                   Name = trim ( T_Option % Name ) // '_Krnl', &
                   Level = T_Option % Level + 1 )
    else
      T_K  =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )
    do iC  =  1, CS % Atlas % nCharts

      associate &
        (   CSV  =>  CS % Storage ( iC ) % Value, &
          M_Ref  =>  CS % BaryonMass, &
          N_Min  =>  CS % BaryonDensityMin )
      associate &
        ( M    =>  CSV ( :, CS % BARYON_MASS ), &
          N    =>  CSV ( :, CS % BARYON_DENSITY_C ), &
          V_1  =>  CSV ( :, CS % VELOCITY_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_U_3 ), &
          D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
          S_1  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ) )
   
      select type ( G  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_UU_11  =>  GSV ( :, G % METRIC_F_UU_11 ), &
            M_UU_22  =>  GSV ( :, G % METRIC_F_UU_22 ), &
            M_UU_33  =>  GSV ( :, G % METRIC_F_UU_33 ) )

        call Compute_N_V_G_Kernel &
               ( D, S_1, S_2, S_3, M, M_UU_11, M_UU_22, M_UU_33, N_Min, &
                 N, V_1, V_2, V_3, UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_UU_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_D__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromBalanced', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      call Compute_M_Kernel &
             ( M_Ref, M, UseDeviceOption = CS % DeviceMemory )

      end associate !-- M, etc.
      end associate !-- CSV, etc.

    end do !-- iC
    if ( associated ( T_K ) ) call T_K % Stop ( )

  end subroutine ComputeFromBalanced


  subroutine ComputeEigenspeeds ( ES, CS, FS_CS, iaEigenspeeds, iC, iD )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      ES
    class ( Fluid_D_Form ), intent ( in ) :: &
      CS
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaEigenspeeds
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    associate &
      ( ESV  =>  ES    % Storage ( iC ) % Value, &
        CSV  =>  FS_CS % Storage ( iC ) % Value )
    associate &
      ( EF_P    =>  ESV ( :, iaEigenspeeds ( 1 ) ), &
        EF_M    =>  ESV ( :, iaEigenspeeds ( 2 ) ), & 
         V_Dim  =>  CSV ( :, CS % VELOCITY_U ( iD ) ) )
 
    call Compute_ES_G_Kernel &
           ( V_Dim, EF_P, EF_M, UseDeviceOption = CS % DeviceMemory )
  
    end associate !-- EF_P, etc.
    end associate !-- FSV, etc.

  end subroutine ComputeEigenspeeds


  impure elemental subroutine Finalize ( F )

    type ( Fluid_D_Form ), intent ( inout ) :: &
      F

  end subroutine Finalize


end module Fluid_D__Form
