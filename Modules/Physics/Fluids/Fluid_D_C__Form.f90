module Fluid_D_C__Form

  !-- Fluid_Dust_Chart_Form

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

  type, public, extends ( CurrentSet_C_Form ) :: Fluid_D_C_Form
    integer ( KDI ), private :: &
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
    real ( KDR ) :: &
      BaryonMassReference, &
      BaryonDensityMin
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
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
  end type Fluid_D_C_Form


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
      !-- Compute_ConservedDensity_Momentum_Galileo_Kernel
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
      !-- Compute_ComovingBaryonDensity_Velocity_Galileo_Kernel
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

    module subroutine ComputeFluxes_G_Kernel &
             ( D, S_1, S_2, S_3, V_Dim, F_D, F_S_1, F_S_2, F_S_3, &
               UseDeviceOption )
      !-- ComputeFluxes_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        D, &
        S_1, S_2, S_3, &
        V_Dim
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        F_D, &
        F_S_1, F_S_2, F_S_3
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeFluxes_G_Kernel

    module subroutine ComputeEigenspeeds_G_Kernel &
             ( V_Dim, EF_P, EF_M, UseDeviceOption )
      !-- Compute_Eigenspeeds_Galileo_Kernel
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        V_Dim
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        EF_P, EF_M
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeEigenspeeds_G_Kernel
    
  end interface


contains


  subroutine InitializeAllocate_F &
               ( FC, GC, Units_F, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( Fluid_D_C_Form ), intent ( inout ) :: &
      FC
    class ( Geometry_F_C_Form ), intent ( in ) :: &
      GC
    class ( Units_F_Form ), intent ( in ) :: &
      Units_F
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
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
      iF, &  !-- iField
      iV, &  !-- iVector
      iP, &  !-- iPrimitive
      iB, &  !-- iBalanced
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
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( FC % Type  ==  '' ) &
      FC % Type  =  'a Fluid_D_C' 
    
    Name  =  'Fluid'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  FC % N_FIELDS_CS

    FC % BARYON_DENSITY_C      =  oF + 1
    FC % BARYON_DENSITY_B      =  oF + 2
    FC % BARYON_MASS           =  oF + 3
    FC % VELOCITY_U_1          =  oF + 4
    FC % VELOCITY_U_2          =  oF + 5
    FC % VELOCITY_U_3          =  oF + 6
    FC % MOMENTUM_DENSITY_D_1  =  oF + 7
    FC % MOMENTUM_DENSITY_D_2  =  oF + 8
    FC % MOMENTUM_DENSITY_D_3  =  oF + 9

    nFields  =  oF  +  FC % N_FIELDS_D
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    FC % VELOCITY_U          =  [ FC % VELOCITY_U_1, &
                                  FC % VELOCITY_U_2, &
                                  FC % VELOCITY_U_3 ]
    FC % MOMENTUM_DENSITY_D  =  [ FC % MOMENTUM_DENSITY_D_1, &
                                  FC % MOMENTUM_DENSITY_D_2, &
                                  FC % MOMENTUM_DENSITY_D_3 ]
 
    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + FC % N_FIELDS_D ) &
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

    if ( present ( UnitOption ) ) then
      allocate ( Unit, source = UnitOption )
    else
      allocate ( Unit ( nFields ) )
    end if !-- FieldOption

    Unit ( FC % BARYON_DENSITY_C )      =  Units_F % NumberDensity
    Unit ( FC % BARYON_DENSITY_B )      =  Units_F % SqrtDet_M  &
                                                *  Units_F % NumberDensity
    Unit ( FC % BARYON_MASS )           =  Units_F % BaryonMass
    Unit ( FC % VELOCITY_U_1 )          =  Units_F % Velocity_U ( 1 )
    Unit ( FC % VELOCITY_U_2 )          =  Units_F % Velocity_U ( 2 )
    Unit ( FC % VELOCITY_U_3 )          =  Units_F % Velocity_U ( 3 )
    Unit ( FC % MOMENTUM_DENSITY_D_1 )  =  Units_F % MomentumDensity_D ( 1 )
    Unit ( FC % MOMENTUM_DENSITY_D_2 )  =  Units_F % MomentumDensity_D ( 2 )
    Unit ( FC % MOMENTUM_DENSITY_D_3 )  =  Units_F % MomentumDensity_D ( 3 )

    !-- Vector indices

    oV  =  FC % N_VECTORS_CS

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  oV  +  FC % N_VECTORS_D  +  1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  oV  +  FC % N_VECTORS_D
      allocate ( VectorIndices ( nVectors ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( FC % VELOCITY_U )
    call VectorIndices ( oV + 2 ) % Initialize ( FC % MOMENTUM_DENSITY_D )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    Vector ( oV  +  1 : oV  +  FC % N_VECTORS_D ) &
      = [ 'Velocity_U       ', &
          'MomentumDensity_D' ]

    !-- Primitive fields

    oP  =  FC % N_PRIMITIVE_CS

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  FC % N_PRIMITIVE_D
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  FC % N_PRIMITIVE_D )  &
      =  [ FC % BARYON_DENSITY_C, FC % VELOCITY_U ]

    !-- Balanced fields

    oB  =  FC % N_BALANCED_CS

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oB  +  FC % N_BALANCED_D
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaPrimitiveOption

    iaBalanced ( oB  +  1 : oB  +  FC % N_BALANCED_D )  &
      =  [ FC % BARYON_DENSITY_B, FC % MOMENTUM_DENSITY_D ]

    !-- CurrentSet

    call FC % CurrentSet_C_Form % Initialize &
           ( GC, &
             FieldOption = Field, &
             VectorOption = Vector, &
             NameOption = Name, &
             UnitOption = Unit, &
             VectorIndicesOption = VectorIndices, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    !-- Parameters

    if ( Units_F % BaryonMass % Number  ==  1.0_KDR ) then
      FC % BaryonMassReference  =  1.0_KDR
    else
      FC % BaryonMassReference  =  CONSTANT % ATOMIC_MASS_UNIT
    end if

  end subroutine InitializeAllocate_F


  subroutine ComputeFromInitial ( CSC )

    class ( Fluid_D_C_Form ), intent ( inout ) :: &
      CSC

    associate &
      ( CSS  =>  CSC % Storage_FSC % Storage, &
        M_Ref  =>  CSC % BaryonMassReference, &
        N_Min  =>  CSC % BaryonDensityMin, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory )
    associate &
      ( M    =>  CSS % Value ( :, CSC % BARYON_MASS ), &
        N    =>  CSS % Value ( :, CSC % BARYON_DENSITY_C ), &
        V_1  =>  CSS % Value ( :, CSC % VELOCITY_U_1 ), &
        V_2  =>  CSS % Value ( :, CSC % VELOCITY_U_2 ), &
        V_3  =>  CSS % Value ( :, CSC % VELOCITY_U_3 ), &
        D    =>  CSS % Value ( :, CSC % BARYON_DENSITY_B ), &
        S_1  =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D_1 ), &
        S_2  =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D_2 ), &
        S_3  =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D_3 ) )
 
    call Compute_M_Kernel &
           ( M_Ref, M, UseDeviceOption = DeviceMemory )

    select type ( GC  =>  CSC % Geometry_C )
    class is ( Gravitation_G_C_Form )

      associate &
        ( GS  =>  GC % Storage_FSC % Storage )
      associate &
        ( M_DD_11  =>  GS % Value ( :, GC % METRIC_F_DD_11 ), &
          M_DD_22  =>  GS % Value ( :, GC % METRIC_F_DD_22 ), &
          M_DD_33  =>  GS % Value ( :, GC % METRIC_F_DD_33 ) )

      call Compute_D_S_G_Kernel & 	 	 
             ( N, V_1, V_2, V_3, M, M_DD_11, M_DD_22, M_DD_33, N_Min, &
               D, S_1, S_2, S_3, UseDeviceOption = DeviceMemory )

      end associate !-- M_DD_11, etc.
      end associate !-- GC

    class default
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_D_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromInitial', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GC

    end associate !-- EF_P, etc.
    end associate !-- CSS, etc.

  end subroutine ComputeFromInitial


  subroutine ComputeFromConserved ( CSC )

    class ( Fluid_D_C_Form ), intent ( inout ) :: &
      CSC

    associate &
      ( CSS  =>  CSC % Storage_FSC % Storage, &
        M_Ref  =>  CSC % BaryonMassReference, &
        N_Min  =>  CSC % BaryonDensityMin, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory )
    associate &
      ( M    =>  CSS % Value ( :, CSC % BARYON_MASS ), &
        N    =>  CSS % Value ( :, CSC % BARYON_DENSITY_C ), &
        V_1  =>  CSS % Value ( :, CSC % VELOCITY_U_1 ), &
        V_2  =>  CSS % Value ( :, CSC % VELOCITY_U_2 ), &
        V_3  =>  CSS % Value ( :, CSC % VELOCITY_U_3 ), &
        D    =>  CSS % Value ( :, CSC % BARYON_DENSITY_B ), &
        S_1  =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D_1 ), &
        S_2  =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D_2 ), &
        S_3  =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D_3 ) )
 
    select type ( GC  =>  CSC % Geometry_C )
    class is ( Gravitation_G_C_Form )

      associate &
        ( GS  =>  GC % Storage_FSC % Storage )
      associate &
        ( M_UU_11  =>  GS % Value ( :, GC % METRIC_F_UU_11 ), &
          M_UU_22  =>  GS % Value ( :, GC % METRIC_F_UU_22 ), &
          M_UU_33  =>  GS % Value ( :, GC % METRIC_F_UU_33 ) )

      call Compute_N_V_G_Kernel &
             ( D, S_1, S_2, S_3, M, M_UU_11, M_UU_22, M_UU_33, N_Min, &
               N, V_1, V_2, V_3, UseDeviceOption = DeviceMemory )

      end associate !-- M_DD_11, etc.
      end associate !-- GC

    class default
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_D_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromConserved', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GC

    call Compute_M_Kernel &
           ( M_Ref, M, UseDeviceOption = DeviceMemory )

    end associate !-- EF_P, etc.
    end associate !-- CSS, etc.

  end subroutine ComputeFromConserved


  subroutine ComputeFluxes ( FSC, CSC, iD )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( Fluid_D_C_Form ), intent ( in ) :: &
      CSC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension
    
    integer ( KDI ) :: &
      iDensity
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    call Search &
           ( CSC % iaBalanced, CSC % BARYON_DENSITY_B, iDensity )
    call Search &
           ( CSC % iaBalanced, CSC % MOMENTUM_DENSITY_D ( 1 ), iMomentum ( 1 ) )
    call Search &
           ( CSC % iaBalanced, CSC % MOMENTUM_DENSITY_D ( 2 ), iMomentum ( 2 ) )
    call Search &
           ( CSC % iaBalanced, CSC % MOMENTUM_DENSITY_D ( 3 ), iMomentum ( 3 ) )

    associate &
      ( FSS  =>  FSC % Storage_FSC % Storage, &
        CSS  =>  CSC % Storage_FSC % Storage, &
        DeviceMemory  =>  FSC % Storage_FSC % DeviceMemory )
    associate &
      ( F_D      =>  FSS % Value ( :, iDensity ), &
        F_S_1    =>  FSS % Value ( :, iMomentum ( 1 ) ), &
        F_S_2    =>  FSS % Value ( :, iMomentum ( 2 ) ), &
        F_S_3    =>  FSS % Value ( :, iMomentum ( 3 ) ), &
          D      =>  CSS % Value ( :, CSC % DENSITY_DEFAULT ), &
          S_1    =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D ( 1 ) ), &
          S_2    =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D ( 2 ) ), &
          S_3    =>  CSS % Value ( :, CSC % MOMENTUM_DENSITY_D ( 3 ) ), &
          V_Dim  =>  CSS % Value ( :, CSC % VELOCITY_U ( iD ) ) )
 
    call ComputeFluxes_G_Kernel &
           ( D, S_1, S_2, S_3, V_Dim, F_D, F_S_1, F_S_2, F_S_3, &
             UseDeviceOption = DeviceMemory )
  
    end associate !-- F_D, etc.
    end associate !-- FSS, etc.

  end subroutine ComputeFluxes


  subroutine ComputeEigenspeeds ( FSC, CSC, iaEigenspeeds, iD )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( Fluid_D_C_Form ), intent ( in ) :: &
      CSC
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaEigenspeeds
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimension

    associate &
      ( FSS  =>  FSC % Storage_FSC % Storage, &
        CSS  =>  CSC % Storage_FSC % Storage, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory )
    associate &
      ( EF_P    =>  FSS % Value ( :, iaEigenspeeds ( 1 ) ), &
        EF_M    =>  FSS % Value ( :, iaEigenspeeds ( 2 ) ), & 
         V_Dim  =>  CSS % Value ( :, CSC % VELOCITY_U ( iD ) ) )
 
    call ComputeEigenspeeds_G_Kernel &
           ( V_Dim, EF_P, EF_M, UseDeviceOption = DeviceMemory )
  
    end associate !-- EF_P, etc.
    end associate !-- FSS, etc.

  end subroutine ComputeEigenspeeds


  impure elemental subroutine Finalize ( FC )

    type ( Fluid_D_C_Form ), intent ( inout ) :: &
      FC

  end subroutine Finalize


end module Fluid_D_C__Form
