module Fluid_P_I__Form

  !-- Fluid_Perfect_Ideal__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_I = 0, &
      N_CONSERVED_I = 0, &
      N_FIELDS_I    = 0, &
      N_VECTORS_I   = 0

  type, public, extends ( Fluid_P_Form ) :: Fluid_P_I_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_I = N_PRIMITIVE_I, &
      N_CONSERVED_I = N_CONSERVED_I, &
      N_FIELDS_I    = N_FIELDS_I, &
      N_VECTORS_I   = N_VECTORS_I
    real ( KDR ) :: &
      AdiabaticIndex!, &
    !   MeanMolecularWeight, &
    !   SpecificHeatVolume, &  !-- per baryon
    !   FiducialBaryonDensity, &
    !   FiducialPressure
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      SetAdiabaticIndex
  !   procedure, public, pass :: &
  !     SetMeanMolecularWeight
  !   procedure, public, pass :: &
  !     SetSpecificHeatVolume
  !   procedure, public, pass :: &
  !     SetFiducialParameters
  !   procedure, public, pass :: &
  !     SetOutput
    procedure, public, pass :: &
      Show => Show_FS
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromTemperature
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromPrimitiveCommon
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromConservedCommon
  !   procedure, public, pass ( C ) :: &
  !     ComputeRawFluxes
  !   procedure, public, nopass :: &
  !     Apply_EOS_I_T_Kernel
  !   procedure, public, nopass :: &
  !     Apply_EOS_I_SB_E_Kernel
  !   procedure, public, nopass :: &
  !     Apply_EOS_I_E_Kernel
    final :: &
      Finalize
  end type Fluid_P_I_Form


contains


  subroutine InitializeAllocate_F &
               ( F, G, Units_F, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
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

    if ( F % Type  ==  '' ) &
      F % Type  =  'a Fluid_P_I' 
    
    !-- Field indices: no additional fields
    
    !-- Field names: no additional fields

    !-- Units: no additional fields

    !-- Vector indices: no additional vectors

    !-- Vector names: no additional vectors

    !-- Primitive fields: no additional fields

    !-- Balanced fields: no additional fields

    !-- Fluid_P

    call F % Fluid_P_Form % Initialize &
           ( G, Units_F, &
             FieldOption = FieldOption, &
             VectorOption = VectorOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitiveOption, &
             iaBalancedOption = iaBalancedOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

    !-- Parameters

    F % AdiabaticIndex  =  1.4_KDR
    call PROGRAM_HEADER % GetParameter &
           ( F % AdiabaticIndex, 'AdiabaticIndex' )

  end subroutine InitializeAllocate_F

  
  subroutine SetAdiabaticIndex ( F, AdiabaticIndex )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      AdiabaticIndex

    F % AdiabaticIndex  =  AdiabaticIndex

    call Show ( 'Setting AdiabaticIndex of a Fluid_P_I', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % AdiabaticIndex, 'AdiabaticIndex', F % IGNORABILITY + 1 )

  end subroutine SetAdiabaticIndex


  subroutine Show_FS ( FS )

    class ( Fluid_P_I_Form ), intent ( in ) :: &
      FS

    call FS % Fluid_D_Form % Show ( )

    call Show ( FS % AdiabaticIndex, 'AdiabaticIndex', FS % IGNORABILITY )

  end subroutine Show_FS


  impure elemental subroutine Finalize ( F )

    type ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

  end subroutine Finalize


end module Fluid_P_I__Form
