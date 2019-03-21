module Fluid_D__Form

  !-- Fluid_Dust__Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergyBasics
  use FluidFeatures_Template

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
      N_VECTORS_DUST    = N_VECTORS_DUST
    integer ( KDI ) :: &
      COMOVING_BARYON_DENSITY  = 0, &
      CONSERVED_BARYON_DENSITY = 0, &
      BARYON_MASS              = 0
    integer ( KDI ), dimension ( 3 ) :: &
      VELOCITY_U         = 0, &
      MOMENTUM_DENSITY_D = 0
    real ( KDR ) :: &
      BaryonMassReference = 1.0_KDR, &
      BaryonDensityMin = 0.0_KDR
    class ( FluidFeaturesTemplate ), pointer :: &
      Features => null ( )
  contains
    procedure, public, pass :: &
      Initialize_D
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetBaryonDensityMin
    procedure, public, pass :: &
      SetOutput
    procedure, public, pass :: &
      SetFeatures
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, nopass :: &
      Compute_M_Kernel
    procedure, public, nopass :: &
      Compute_D_S_G_Kernel
    procedure, public, nopass :: &
      Compute_N_V_G_Kernel
    procedure, public, nopass :: &
      Compute_FE_D_G_Kernel
  end type Fluid_D_Form

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeRawFluxes_G_Kernel, &
      ComputeRawFluxes_N_S_Kernel

contains


  subroutine Initialize_D &
               ( F, RiemannSolverType, ReconstructedType, UseLimiter, Units, &
                 BaryonMassReference, LimiterParameter, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    character ( * ), intent ( in ) :: &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    real ( KDR ), intent ( in ) :: &
      BaryonMassReference, &
      LimiterParameter
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

    call SetUnits ( VariableUnit, F, Units )

    call F % InitializeTemplate &
           ( RiemannSolverType, ReconstructedType, &
             UseLimiter, Units % Velocity_U, LimiterParameter, nValues, &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = ClearOption, &
             UnitOption = VariableUnit, VectorIndicesOption = VectorIndices )

    F % BaryonMassReference = BaryonMassReference
    call Show ( F % BaryonMassReference, F % Unit ( F % BARYON_MASS ), &
                'BaryonMassReference', F % IGNORABILITY )

  end subroutine Initialize_D
  

  subroutine SetPrimitiveConserved ( C )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_DUST ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_DUST ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE
    oC = C % N_CONSERVED_TEMPLATE

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_DUST
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_DUST ) &
      = [ C % COMOVING_BARYON_DENSITY, C % VELOCITY_U ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_DUST
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_DUST ) &
      = [ C % CONSERVED_BARYON_DENSITY, C % MOMENTUM_DENSITY_D ]

    do iF = 1, C % N_PRIMITIVE_DUST
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_DUST
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )

  end subroutine SetPrimitiveConserved


  subroutine SetBaryonDensityMin ( F, BaryonDensityMin )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      BaryonDensityMin

    F % BaryonDensityMin = BaryonDensityMin

    call Show ( 'Setting BaryonDensityMin of a Fluid', F % IGNORABILITY )
    call Show ( F % Name, 'Name', F % IGNORABILITY )
    call Show ( F % BaryonDensityMin, &
                F % Unit ( F % COMOVING_BARYON_DENSITY ), 'BaryonDensityMin', &
                F % IGNORABILITY )

  end subroutine SetBaryonDensityMin


  subroutine SetOutput ( F, Output )

    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    type ( StorageForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( F % VELOCITY_U )
    call Output % Initialize &
           ( F, iaSelectedOption = [ F % COMOVING_BARYON_DENSITY, &
                                     F % VELOCITY_U ], &
             VectorOption = [ 'Velocity' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  subroutine SetFeatures ( F, Features )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    class ( FluidFeaturesTemplate ), intent ( in ), target :: &
      Features

    call Show ( 'Setting Features', F % IGNORABILITY )
    call Show ( F % Name, 'Name', F % IGNORABILITY )
    call Show ( Features % Name, 'Features', F % IGNORABILITY )

    F % Features => Features

  end subroutine SetFeatures


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
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_BARYON_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_BARYON_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ) )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference )
    call C % Compute_D_S_G_Kernel &
           ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
    call C % Compute_FE_D_G_Kernel &
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
        N     => FV ( oV + 1 : oV + nV, C % COMOVING_BARYON_DENSITY ), &
        V_1   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 1 ) ), &
        V_2   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 2 ) ), &
        V_3   => FV ( oV + 1 : oV + nV, C % VELOCITY_U ( 3 ) ), &
        D     => FV ( oV + 1 : oV + nV, C % CONSERVED_BARYON_DENSITY ), &
        S_1   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => FV ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ) )

    call C % Compute_M_Kernel &
           ( M, C % BaryonMassReference )
    call C % Compute_N_V_G_Kernel &
           ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33, &
             C % BaryonDensityMin )
    call C % Compute_FE_D_G_Kernel &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3 )

    end associate !-- FEP_1, etc.
    end associate !-- M_UU_22, etc.
    end associate !-- FV, etc.

  end subroutine ComputeFromConservedCommon


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
           ( C % iaConserved, C % CONSERVED_BARYON_DENSITY, iDensity )
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
        D     => Value_C ( oV + 1 : oV + nV, C % CONSERVED_BARYON_DENSITY ), &
        S_1   => Value_C ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2   => Value_C ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3   => Value_C ( oV + 1 : oV + nV, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        V_Dim => Value_C ( oV + 1 : oV + nV, C % VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxes_G_Kernel &
           ( F_D, F_S_1, F_S_2, F_S_3, D, S_1, S_2, S_3, V_Dim )

    select type ( G )
    class is ( Geometry_N_S_Form )
      associate &
        ( GradPhi_1 &
            => Value_G ( oV + 1 : oV + nV, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
          GradPhi_2 &
            => Value_G ( oV + 1 : oV + nV, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
          GradPhi_3 &
            => Value_G ( oV + 1 : oV + nV, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
          M_UU_22 &
            => Value_G ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
          M_UU_33 &
            => Value_G ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
      call ComputeRawFluxes_N_S_Kernel &
             ( F_S_1, F_S_2, F_S_3, GradPhi_1, GradPhi_2, GradPhi_3, &
               M_UU_22, M_UU_33, iDimension )
      end associate !-- GradPhi_1, etc.
    end select

    end associate !-- F_D, etc.

  end subroutine ComputeRawFluxes
  
  
  subroutine Compute_M_Kernel ( M, BaryonMassReference )

    !-- Compute_BaryonMass_Kernel

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      M
    real ( KDR ), intent ( in ) :: &
      BaryonMassReference

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( M )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      M ( iV )  =  BaryonMassReference
    end do !-- iV
    !$OMP end parallel do

  end subroutine Compute_M_Kernel


  subroutine Compute_D_S_G_Kernel & 	 	 
	       ( D, S_1, S_2, S_3, N, M, V_1, V_2, V_3, M_DD_22, M_DD_33 )
 	 
    !-- Compute_ConservedDensity_Momentum_Galilean_Kernel

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

  end subroutine Compute_D_S_G_Kernel 	 	 


  subroutine Compute_N_V_G_Kernel &
               ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, M, M_UU_22, M_UU_33, &
                 N_Min )

    !-- Compute_ComovingBaryonDensity_Velocity_Galilean

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      V_1, V_2, V_3, &
      D, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      M_UU_22, M_UU_33
    real ( KDR ), intent ( in ) :: &
      N_Min

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( N )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( D ( iV )  >  N_Min ) then
        N ( iV )    =  D ( iV )
        V_1 ( iV )  =  S_1 ( iV ) / ( M ( iV ) * D ( iV ) )
        V_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV ) / ( M ( iV ) * D ( iV ) )
        V_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV ) / ( M ( iV ) * D ( iV ) )
      else
        N   ( iV )  =  N_Min
        V_1 ( iV )  =  0.0_KDR
        V_2 ( iV )  =  0.0_KDR
        V_3 ( iV )  =  0.0_KDR
        D   ( iV )  =  N_Min
        S_1 ( iV )  =  0.0_KDR
        S_2 ( iV )  =  0.0_KDR
        S_3 ( iV )  =  0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine Compute_N_V_G_Kernel


  subroutine Compute_FE_D_G_Kernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, V_1, V_2, V_3 )

    !-- Compute_FastEigenspeeds_Dust_Galilean_Kernel

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
      FEP_1 ( iV )  =  V_1 ( iV ) 
      FEP_2 ( iV )  =  V_2 ( iV ) 
      FEP_3 ( iV )  =  V_3 ( iV ) 
      FEM_1 ( iV )  =  V_1 ( iV ) 
      FEM_2 ( iV )  =  V_2 ( iV ) 
      FEM_3 ( iV )  =  V_3 ( iV ) 
    end do !-- iV
    !$OMP end parallel do

  end subroutine Compute_FE_D_G_Kernel


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
      oV     !-- oVector

    if ( F % Type == '' ) &
      F % Type = 'a Fluid_D'

    Name = 'Fluid'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = F % N_FIELDS_TEMPLATE
    if ( F % N_FIELDS == 0 ) &
      F % N_FIELDS = oF + F % N_FIELDS_DUST

    F % COMOVING_BARYON_DENSITY   =  oF + 1
    F % CONSERVED_BARYON_DENSITY  =  oF + 2
    F % BARYON_MASS               =  oF + 3
    F % VELOCITY_U                =  oF + [ 4, 5, 6 ]
    F % MOMENTUM_DENSITY_D        =  oF + [ 7, 8, 9 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( F % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + F % N_FIELDS_DUST ) &
      = [ 'ComovingBaryonDensity ', &
          'ConservedBaryonDensity', &
          'BaryonMass            ', &
          'Velocity_1            ', &
          'Velocity_2            ', &
          'Velocity_3            ', &
          'MomentumDensity_1     ', &
          'MomentumDensity_2     ', &
          'MomentumDensity_3     ' ]
          
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
    if ( F % N_VECTORS == 0 ) &
      F % N_VECTORS = oV + F % N_VECTORS_DUST

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( F % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + F % N_VECTORS_DUST ) &
      = [ 'Velocity       ', &
          'MomentumDensity' ]

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

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, F, Units )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Fluid_D_Form ), intent ( in ) :: &
      F
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units

    integer ( KDI ) :: &
      iD

    VariableUnit ( F % BARYON_MASS )              = Units % BaryonMass
    VariableUnit ( F % COMOVING_BARYON_DENSITY )  = Units % NumberDensity
    VariableUnit ( F % CONSERVED_BARYON_DENSITY ) = Units % NumberDensity

    do iD = 1, 3
      VariableUnit ( F % VELOCITY_U ( iD ) ) &
        = Units % Velocity_U ( iD )
      VariableUnit ( F % MOMENTUM_DENSITY_D ( iD ) ) &
        = Units % MomentumDensity_D ( iD )
    end do

  end subroutine SetUnits


  subroutine ComputeRawFluxes_G_Kernel &
               ( F_D, F_S_1, F_S_2, F_S_3, D, S_1, S_2, S_3, V_Dim )

    !-- Galilean

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
      F_D   ( iV )  =  D   ( iV )  *  V_Dim ( iV ) 
      F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV ) 
      F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV ) 
      F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV ) 
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRawFluxes_G_Kernel


  subroutine ComputeRawFluxes_N_S_Kernel &
               ( F_S_1, F_S_2, F_S_3, GradPhi_1, GradPhi_2, GradPhi_3, &
                 M_UU_22, M_UU_33, iDimension )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_S_1, F_S_2, F_S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      GradPhi_1, GradPhi_2, GradPhi_3, &
      M_UU_22, M_UU_33
    integer ( KDI ) :: &
      iDimension

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      One_FourPi_G

    One_FourPi_G  =  1.0_KDR  /  ( 4.0_KDR  *  CONSTANT % PI  &
                                   *  CONSTANT % GRAVITATIONAL )

    nValues = size ( F_S_1 )

    select case ( iDimension )
    case ( 1 )
      !$OMP parallel do private ( iV ) 
      do iV = 1, nValues
        F_S_1 ( iV )  &
          =  F_S_1 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_1 ( iV )  *  GradPhi_1 ( iV ) &
             -  One_FourPi_G  *  0.5_KDR  &
                *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                     +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                     +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
        F_S_2 ( iV )  &
          =  F_S_2 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_2 ( iV )  *  GradPhi_1 ( iV ) 
        F_S_3 ( iV )  &
          =  F_S_3 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_3 ( iV )  *  GradPhi_1 ( iV ) 
      end do !-- iV
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV ) 
      do iV = 1, nValues
        F_S_1 ( iV )  &
          =  F_S_1 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_1 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV ) 
        F_S_2 ( iV )  &
          =  F_S_2 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_2 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV ) &
             -  One_FourPi_G  *  0.5_KDR  &
                *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                     +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                     +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
        F_S_3 ( iV )  &
          =  F_S_3 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_3 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV )
      end do !-- iV
      !$OMP end parallel do
    case ( 3 )
      !$OMP parallel do private ( iV ) 
      do iV = 1, nValues
        F_S_1 ( iV )  &
          =  F_S_1 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_1 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV ) 
        F_S_2 ( iV )  &
          =  F_S_2 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_2 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV )
        F_S_3 ( iV )  &
          =  F_S_3 ( iV )  &
             +  One_FourPi_G  &
                *  GradPhi_3 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV ) &
             -  One_FourPi_G  *  0.5_KDR  &
                *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                     +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                     +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
      end do !-- iV
      !$OMP end parallel do
    end select !-- iDimension

  end subroutine ComputeRawFluxes_N_S_Kernel


end module Fluid_D__Form
