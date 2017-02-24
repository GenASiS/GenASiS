module Fluid_D__Form

  !-- Fluid_Dust__Form

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_DUST = 4, &
      N_CONSERVED_DUST = 4, &
      N_FIELDS_DUST    = 9, &
      N_VECTORS_DUST   = 2

  type, public, extends ( CurrentTemplate ) :: Fluid_D_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_DUST  = N_PRIMITIVE_DUST, &
      N_CONSERVED_DUST  = N_CONSERVED_DUST, &
      N_FIELDS_DUST     = N_FIELDS_DUST, &
      N_VECTORS_DUST    = N_VECTORS_DUST, &
      COMOVING_DENSITY  = 0, &
      CONSERVED_DENSITY = 0, &
      BARYON_MASS       = 0
    integer ( KDI ), dimension ( 3 ) :: &
      VELOCITY_U         = 0, &
      MOMENTUM_DENSITY_D = 0
  contains
    procedure, public, pass :: &
      InitializeAllocate_D
    generic, public :: &
      Initialize => InitializeAllocate_D
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
      ComputeBaryonMassKernel
    procedure, public, nopass :: &
      ComputeDensityMomentumKernel
    procedure, public, nopass :: &
      ComputeDensityVelocityKernel
    procedure, public, nopass :: &
      ComputeEigenspeedsKernel_D
  end type Fluid_D_Form

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeRawFluxesKernel

  public :: &
    ApplySourcesCurvilinear_Fluid_D

    private :: &
      ApplySourcesCurvilinearKernel
  
contains


  subroutine InitializeAllocate_D &
               ( F, VelocityUnit, MassDensityUnit, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit
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
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
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
           ( F, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits ( VariableUnit, F, VelocityUnit, MassDensityUnit )

    call F % InitializeTemplate &
           ( VelocityUnit, nValues, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_D
  

  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                 nValuesOption, oValueOption )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    class ( Fluid_D_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_C, &
      Value_G
    integer ( KDI ), intent ( in ) :: &
      iDimension
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      iDensity
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum
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
      nV = size ( Value_C, dim = 1 )
    end if
    
    call Search &
           ( C % iaConserved, C % CONSERVED_DENSITY, iDensity )
    call Search &
           ( C % iaConserved, C % MOMENTUM_DENSITY_D ( 1 ), iMomentum ( 1 ) )
    call Search &
           ( C % iaConserved, C % MOMENTUM_DENSITY_D ( 2 ), iMomentum ( 2 ) )
    call Search &
           ( C % iaConserved, C % MOMENTUM_DENSITY_D ( 3 ), iMomentum ( 3 ) )
    
    associate &
      ( F_D   => RawFlux ( oV + 1 : oV + nV, iDensity ), &
        F_S_1 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 1 ) ), &
        F_S_2 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 2 ) ), &
        F_S_3 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 3 ) ), &
        D     => Value_C ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ), &
        S_1   => Value_C ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => Value_C ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => Value_C ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        V_Dim => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxesKernel &
           ( F_D, F_S_1, F_S_2, F_S_3, D, S_1, S_2, S_3, V_Dim )

    end associate !-- F_D, etc.

  end subroutine ComputeRawFluxes
  
  
  subroutine SetOutput ( F, Output )

    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption = [ F % COMOVING_DENSITY, F % VELOCITY_U ], &
             VectorOption = [ 'Velocity                       ' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  impure elemental subroutine Finalize ( F )

    type ( Fluid_D_Form ), intent ( inout ) :: &
      F

    call F % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeFromPrimitiveCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( Fluid_D_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => Value_C, &
        GV => Value_G )

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
        M     => FV ( oV + 1 : oV + nV, C % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ) )

    call C % ComputeBaryonMassKernel &
           ( M )
    call C % ComputeDensityMomentumKernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % ComputeEigenspeedsKernel_D &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3 )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- FV, etc.

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( Fluid_D_Form ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Value_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( FV => Value_C, &
        GV => Value_G )

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
        M     => FV ( oV + 1 : oV + nV, C % BARYON_MASS ), &
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ) )

    call C % ComputeBaryonMassKernel &
           ( M )
    call C % ComputeDensityVelocityKernel &
           ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )
    call C % ComputeEigenspeedsKernel_D &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.

  end subroutine ComputeFromConservedCommon


  subroutine ComputeBaryonMassKernel ( M )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( M )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      M ( iV ) = 1.0_KDR
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeBaryonMassKernel


  subroutine ComputeDensityMomentumKernel & 	 	 
	       ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
 	 
    real ( KDR ), dimension ( : ), intent ( inout ) :: & 	 	 
      D, & 	 	 
      S_1, S_2, S_3, & 	 	 
      N 	 	 
    real ( KDR ), dimension ( : ), intent ( in ) :: & 	 	 
      M, & 	 	 
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33
 	 	 
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( D )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( N ( iV )  <  0.0_KDR ) &
        N ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      D ( iV ) = N ( iV ) 	 	 
 	 	 
      S_1 ( iV )  =  M ( iV ) * N ( iV ) * V_1 ( iV )	 	 
      S_2 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_22 ( iV ) * V_2 ( iV )
      S_3 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_33 ( iV ) * V_3 ( iV )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeDensityMomentumKernel 	 	 


  subroutine ComputeDensityVelocityKernel &
               ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      V_1, V_2, V_3, &
      D, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      M_UU_22, M_UU_33

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( N )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( D ( iV )  >  0.0_KDR ) then
        N ( iV )    =  D ( iV )
        V_1 ( iV )  =  S_1 ( iV ) / ( M ( iV ) * D ( iV ) )
        V_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV ) / ( M ( iV ) * D ( iV ) )
        V_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV ) / ( M ( iV ) * D ( iV ) )
      else
        N   ( iV ) = 0.0_KDR
        V_1 ( iV ) = 0.0_KDR
        V_2 ( iV ) = 0.0_KDR
        V_3 ( iV ) = 0.0_KDR
        D   ( iV ) = 0.0_KDR
        S_1 ( iV ) = 0.0_KDR
        S_2 ( iV ) = 0.0_KDR
        S_3 ( iV ) = 0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeDensityVelocityKernel


  subroutine ComputeEigenspeedsKernel_D &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( FEP_1 )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      FEP_1 ( iV ) = V_1 ( iV ) 
      FEP_2 ( iV ) = V_2 ( iV ) 
      FEP_3 ( iV ) = V_3 ( iV ) 
      FEM_1 ( iV ) = V_1 ( iV ) 
      FEM_2 ( iV ) = V_2 ( iV ) 
      FEM_3 ( iV ) = V_3 ( iV ) 
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeEigenspeedsKernel_D


  subroutine InitializeBasics &
               ( F, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
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
      iV, &  !-- iVector
      oF, &  !-- oField
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oC     !-- oConserved

    if ( F % Type == '' ) &
      F % Type = 'a Fluid_D'

    Name = 'Fluid'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_DUST

    F % COMOVING_DENSITY   = oF + 1
    F % CONSERVED_DENSITY  = oF + 2
    F % BARYON_MASS        = oF + 3
    F % VELOCITY_U         = oF + [ 4, 5, 6 ]
    F % MOMENTUM_DENSITY_D = oF + [ 7, 8, 9 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_DUST ) &
      = [ 'ComovingDensity                ', &
          'ConservedDensity               ', &
          'BaryonMass                     ', &
          'Velocity_1                     ', &
          'Velocity_2                     ', &
          'Velocity_3                     ', &
          'MomentumDensity_1              ', &
          'MomentumDensity_2              ', &
          'MomentumDensity_3              ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( F % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = F % N_VECTORS_TEMPLATE
    if ( F % N_VECTORS == 0 ) F % N_VECTORS = oV + F % N_VECTORS_DUST

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( F % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + F % N_VECTORS_DUST ) &
      = [ 'Velocity                       ', &
          'MomentumDensity                ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + F % N_VECTORS_DUST + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( F % N_VECTORS ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( F % VELOCITY_U )
    call VectorIndices ( oV + 2 ) % Initialize ( F % MOMENTUM_DENSITY_D )

    !-- select primitive, conserved

    oP = F % N_PRIMITIVE_TEMPLATE
    oC = F % N_CONSERVED_TEMPLATE

    if ( .not. allocated ( F % iaPrimitive ) ) then
      F % N_PRIMITIVE = oP + F % N_PRIMITIVE_DUST
      allocate ( F % iaPrimitive ( F % N_PRIMITIVE ) )
    end if
    F % iaPrimitive ( oP + 1 : oP + F % N_PRIMITIVE_DUST ) &
      = [ F % COMOVING_DENSITY, F % VELOCITY_U ]

    if ( .not. allocated ( F % iaConserved ) ) then
      F % N_CONSERVED = oC + F % N_CONSERVED_DUST
      allocate ( F % iaConserved ( F % N_CONSERVED ) )
    end if
    F % iaConserved ( oC + 1 : oC + F % N_CONSERVED_DUST ) &
      = [ F % CONSERVED_DENSITY, F % MOMENTUM_DENSITY_D ]
    
  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, F, VelocityUnit, MassDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      MassDensityUnit

    integer ( KDI ) :: &
      iD

    VariableUnit ( F % COMOVING_DENSITY )  = MassDensityUnit
    VariableUnit ( F % CONSERVED_DENSITY ) = MassDensityUnit

    do iD = 1, 3
      VariableUnit ( F % VELOCITY_U ( iD ) ) = VelocityUnit ( iD )
      VariableUnit ( F % MOMENTUM_DENSITY_D ( iD ) ) &
        = MassDensityUnit * VelocityUnit ( iD )
    end do

  end subroutine SetUnits


  subroutine ComputeRawFluxesKernel &
               ( F_D, F_S_1, F_S_2, F_S_3, D, S_1, S_2, S_3, V_Dim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_D, &
      F_S_1, F_S_2, F_S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      D, &
      S_1, S_2, S_3, &
      V_Dim
    
    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( F_D )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      F_D   ( iV ) = D   ( iV ) * V_Dim ( iV ) 
      F_S_1 ( iV ) = S_1 ( iV ) * V_Dim ( iV ) 
      F_S_2 ( iV ) = S_2 ( iV ) * V_Dim ( iV ) 
      F_S_3 ( iV ) = S_3 ( iV ) * V_Dim ( iV ) 
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRawFluxesKernel


  subroutine ApplySourcesCurvilinear_Fluid_D ( S, Increment, Fluid, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iMomentum_1, &
      iMomentum_2

    select type ( F => Fluid )
    class is ( Fluid_D_Form )

    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )

    if ( trim ( Grid % CoordinateSystem ) == 'CARTESIAN' ) &
      return

    call ApplySourcesCurvilinearKernel &
           ( Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             Grid % CoordinateSystem, Grid % IsProperCell, &
             F % Value ( :, F % MOMENTUM_DENSITY_D ( 2 ) ), &
             F % Value ( :, F % MOMENTUM_DENSITY_D ( 3 ) ), &
             F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             F % Value ( :, F % VELOCITY_U ( 3 ) ), &
             S % dLogVolumeJacobian_dX ( 1 ) % Value, &
             S % dLogVolumeJacobian_dX ( 2 ) % Value, &
             TimeStep, Grid % nDimensions )

    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'Fluid_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ApplySourcesCurvilinear', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end select !-- F
    
  end subroutine ApplySourcesCurvilinear_Fluid_D


  subroutine ApplySourcesCurvilinearKernel &
               ( KVM_1, KVM_2, CoordinateSystem, IsProperCell, &
                 S_2, S_3, V_2, V_3, dLVJ_dX1, dLVJ_dX2, dT, nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVM_1, KVM_2
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      S_2, S_3, &
      V_2, V_3, &
      dLVJ_dX1, dLVJ_dX2
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( KVM_1 )

    select case ( trim ( CoordinateSystem ) )
    case ( 'CYLINDRICAL' )

      !$OMP parallel do private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        KVM_1 ( iV ) &
          =  KVM_1 ( iV )  +  ( V_3 ( iV ) * S_3 ( iV ) ) &
                              * dLVJ_dX1 ( iV ) * dT
      end do
      !$OMP end parallel do

    case ( 'SPHERICAL' )

      !$OMP parallel do private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        KVM_1 ( iV ) &
          =  KVM_1 ( iV )  +  ( V_2 ( iV ) * S_2 ( iV )  &
                                +  V_3 ( iV ) * S_3 ( iV ) ) &
                              * 0.5_KDR * dLVJ_dX1 ( iV ) * dT
      end do
      !$OMP end parallel do

      if ( nDimensions > 1 ) then
        !$OMP parallel do private ( iV )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          KVM_2 ( iV )  &
            =  KVM_2 ( iV )  +  ( V_3 ( iV ) * S_3 ( iV ) ) &
                                * dLVJ_dX2 ( iV ) * dT
        end do
        !$OMP end parallel do
      end if

    end select !-- CoordinateSystem

  end subroutine ApplySourcesCurvilinearKernel


end module Fluid_D__Form
