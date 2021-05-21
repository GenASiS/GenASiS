module Fluid_D_C__Form

  !-- Fluid_Dust_Chart_Form

  use Basics
  use Mathematics
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
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    final :: &
      Finalize
  end type Fluid_D_C_Form


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
    FC % MOMENTUM_DENSITY_D_2  =  oF + 7
    FC % MOMENTUM_DENSITY_D_3  =  oF + 7

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

  end subroutine InitializeAllocate_F


  impure elemental subroutine Finalize ( FC )

    type ( Fluid_D_C_Form ), intent ( inout ) :: &
      FC

  end subroutine Finalize


end module Fluid_D_C__Form
