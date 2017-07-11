module RadiationMoments_Form

  use Basics
  use Mathematics
  use Interactions_Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_RM =  7, &
      N_CONSERVED_RM =  4, &
      N_FIELDS_RM    = 15, &
      N_VECTORS_RM   =  3

  type, public, extends ( CurrentTemplate ) :: RadiationMomentsForm
    integer ( KDI ) :: &
      N_PRIMITIVE_RM     = N_PRIMITIVE_RM, &
      N_CONSERVED_RM     = N_CONSERVED_RM, &
      N_FIELDS_RM        = N_FIELDS_RM, &
      N_VECTORS_RM       = N_VECTORS_RM, &
      COMOVING_ENERGY    = 0, &
      CONSERVED_ENERGY   = 0, &
      FLUX_FACTOR        = 0, &
      STRESS_FACTOR      = 0, &
      COMOVING_ENERGY_EQ = 0, &
      DIFFUSION_FACTOR_E = 0
    integer ( KDI ), dimension ( 3 ) :: &
      COMOVING_MOMENTUM_U  = 0, &
      CONSERVED_MOMENTUM_D = 0, &
      FLUID_VELOCITY_U     = 0
    class ( InteractionsTemplate ), pointer :: &
      Interactions => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    generic, public :: &
      Initialize => InitializeAllocate_RM
    procedure, public, pass :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetInteractions
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass ( C ) :: &
      ComputeDiffusionFactor_HLL
  end type RadiationMomentsForm

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeConservedEnergyMomentum, &
      ComputeComovingEnergyMomentum, &
      ComputeEigenspeeds, &
      ComputeRawFluxesKernel, &
      ComputeDiffusionFactor_HLL_CSL

      private :: &
        ComputeMomentFactors, &
        ComputeComovingStress_D, &
        ComputeComovingNonlinearSolve

contains


  subroutine InitializeAllocate_RM &
               ( RM, RiemannSolverType, UseLimiter, Velocity_U_Unit, &
                 MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
                 EnergyDensityUnit, LimiterParameter, &
                 nValues, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    character ( * ), intent ( in ) :: &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit
    real ( KDR ), intent ( in ) :: &
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
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
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
           ( RM, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits &
           ( VariableUnit, RM, Velocity_U_Unit, MomentumDensity_U_Unit, &
             MomentumDensity_D_Unit, EnergyDensityUnit )

    call RM % InitializeTemplate &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             LimiterParameter, nValues, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )
    
  end subroutine InitializeAllocate_RM


  subroutine SetPrimitiveConserved ( C )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iF, &  !-- iField
      oP, &  !-- oPrimitive
      oC     !-- oConserved
    character ( LDL ), dimension ( C % N_PRIMITIVE_RM ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED_RM ) :: &
      ConservedName

    oP = C % N_PRIMITIVE_TEMPLATE
    oC = C % N_CONSERVED_TEMPLATE

    if ( .not. allocated ( C % iaPrimitive ) ) then
      C % N_PRIMITIVE = oP + C % N_PRIMITIVE_RM
      allocate ( C % iaPrimitive ( C % N_PRIMITIVE ) )
    end if
    C % iaPrimitive ( oP + 1 : oP + C % N_PRIMITIVE_RM ) &
      = [ C % COMOVING_ENERGY, C % COMOVING_MOMENTUM_U, C % FLUID_VELOCITY_U ]

    if ( .not. allocated ( C % iaConserved ) ) then
      C % N_CONSERVED = oC + C % N_CONSERVED_RM
      allocate ( C % iaConserved ( C % N_CONSERVED ) )
    end if
    C % iaConserved ( oC + 1 : oC + C % N_CONSERVED_RM ) &
      = [ C % CONSERVED_ENERGY, C % CONSERVED_MOMENTUM_D ]
    
    do iF = 1, C % N_PRIMITIVE_RM
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( oP + iF ) )
    end do
    do iF = 1, C % N_CONSERVED_RM
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( oC + iF ) )
    end do
    call Show ( PrimitiveName, 'Adding primitive variables', &
                C % IGNORABILITY, oIndexOption = oP )
    call Show ( ConservedName, 'Adding conserved variables', &
                C % IGNORABILITY, oIndexOption = oC )
    
  end subroutine SetPrimitiveConserved


  subroutine SetInteractions ( RM, Interactions )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    class ( InteractionsTemplate ), intent ( in ), target :: &
      Interactions

    RM % Interactions => Interactions

  end subroutine SetInteractions


  subroutine SetOutput ( RM, Output )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( RM % COMOVING_MOMENTUM_U )

    call Output % Initialize &
           ( RM, iaSelectedOption = [ RM % COMOVING_ENERGY, &
                                      RM % COMOVING_MOMENTUM_U, &
                                      RM % FLUID_VELOCITY_U, &
                                      RM % FLUX_FACTOR, &
                                      RM % STRESS_FACTOR, &
                                      RM % COMOVING_ENERGY_EQ ], &
             VectorOption = [ 'ComovingMomentum_U' ], &
             VectorIndicesOption = VectorIndices )

  end subroutine SetOutput


  impure elemental subroutine Finalize ( RM )

    type ( RadiationMomentsForm ), intent ( inout ) :: &
      RM

    nullify ( RM % Interactions )

    call RM % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeFromPrimitiveCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( RadiationMomentsForm ), intent ( in ) :: &
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
    real ( KDR ), dimension ( :, : ), pointer :: &
      RMV
      
    RMV => Value_C
    associate ( GV => Value_G )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( RMV, dim = 1 )
    end if
    
    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ), &
        M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        J     => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        H_1   => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 1 ) ), &
        H_2   => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 2 ) ), &
        H_3   => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 3 ) ), &
        E     => RMV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        S_1   => RMV ( oV + 1 : oV + nV, C % CONSERVED_MOMENTUM_D ( 1 ) ), &
        S_2   => RMV ( oV + 1 : oV + nV, C % CONSERVED_MOMENTUM_D ( 2 ) ), &
        S_3   => RMV ( oV + 1 : oV + nV, C % CONSERVED_MOMENTUM_D ( 3 ) ), &
        FF    => RMV ( oV + 1 : oV + nV, C % FLUX_FACTOR ), &
        SF    => RMV ( oV + 1 : oV + nV, C % STRESS_FACTOR ), &
        V_1   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 1 ) ), &
        V_2   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ), &
        V_3   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 3 ) ) )

    call ComputeConservedEnergyMomentum &
           ( E, S_1, S_2, S_3, J, H_1, H_2, H_3, FF, SF, M_DD_22, M_DD_33, &
             V_1, V_2, V_3 )
    call ComputeEigenspeeds &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, J, H_1, H_2, H_3, &
             M_UU_22, M_UU_33, CONSTANT % SPEED_OF_LIGHT )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- GV
    nullify ( RMV )

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( Value_C, C, G, Value_G, nValuesOption, oValueOption )

    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value_C
    class ( RadiationMomentsForm ), intent ( in ) :: &
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
    real ( KDR ), dimension ( :, : ), pointer :: &
      RMV
      
    RMV => Value_C
    associate ( GV => Value_G )

    if ( present ( oValueOption ) ) then
      oV = oValueOption
    else
      oV = 0
    end if

    if ( present ( nValuesOption ) ) then
      nV = nValuesOption
    else
      nV = size ( RMV, dim = 1 )
    end if
    
    associate &
      ( M_DD_22 => GV ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => GV ( oV + 1 : oV + nV, G % METRIC_DD_33 ), &
        M_UU_22 => GV ( oV + 1 : oV + nV, G % METRIC_UU_22 ), &
        M_UU_33 => GV ( oV + 1 : oV + nV, G % METRIC_UU_33 ) )
    associate &
      ( FEP_1 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
        FEP_2 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
        FEP_3 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
        FEM_1 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
        FEM_2 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
        FEM_3 => RMV ( oV + 1 : oV + nV, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
        J     => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        H_1   => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 1 ) ), &
        H_2   => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 2 ) ), &
        H_3   => RMV ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 3 ) ), &
        E     => RMV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY ), &
        S_1   => RMV ( oV + 1 : oV + nV, C % CONSERVED_MOMENTUM_D ( 1 ) ), &
        S_2   => RMV ( oV + 1 : oV + nV, C % CONSERVED_MOMENTUM_D ( 2 ) ), &
        S_3   => RMV ( oV + 1 : oV + nV, C % CONSERVED_MOMENTUM_D ( 3 ) ), &
        FF    => RMV ( oV + 1 : oV + nV, C % FLUX_FACTOR ), &
        SF    => RMV ( oV + 1 : oV + nV, C % STRESS_FACTOR ), &
        V_1   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 1 ) ), &
        V_2   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ), &
        V_3   => RMV ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 3 ) ) )

    call ComputeComovingEnergyMomentum &
           ( J, H_1, H_2, H_3, E, S_1, S_2, S_3, FF, SF, C, M_DD_22, M_DD_33, &
             M_UU_22, M_UU_33, V_1, V_2, V_3 )
    call ComputeEigenspeeds &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, J, H_1, H_2, H_3, &
             M_UU_22, M_UU_33, CONSTANT % SPEED_OF_LIGHT )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- GV
    nullify ( RMV )

  end subroutine ComputeFromConservedCommon


  subroutine ComputeRawFluxes &
               ( RawFlux, C, G, Value_C, Value_G, iDimension, &
                 nValuesOption, oValueOption )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      RawFlux
    class ( RadiationMomentsForm ), intent ( in ) :: &
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
      iEnergy
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
           ( C % iaConserved, C % CONSERVED_ENERGY, iEnergy )
    call Search &
           ( C % iaConserved, C % CONSERVED_MOMENTUM_D ( 1 ), iMomentum ( 1 ) )
    call Search &
           ( C % iaConserved, C % CONSERVED_MOMENTUM_D ( 2 ), iMomentum ( 2 ) )
    call Search &
           ( C % iaConserved, C % CONSERVED_MOMENTUM_D ( 3 ), iMomentum ( 3 ) )
    
    associate &
      ( M_DD_22 => Value_G ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => Value_G ( oV + 1 : oV + nV, G % METRIC_DD_33 ) )
    associate &
      ( F_E   => RawFlux ( oV + 1 : oV + nV, iEnergy ), &
        F_S_1 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 1 ) ), &
        F_S_2 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 2 ) ), &
        F_S_3 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 3 ) ), &
        F_S_Dim => RawFlux ( oV + 1 : oV + nV, iMomentum ( iDimension ) ), & 
        J   => Value_C ( oV + 1 : oV + nV, C % COMOVING_ENERGY ), &
        H_1 => Value_C ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 1 ) ), &
        H_2 => Value_C ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 2 ) ), &
        H_3 => Value_C ( oV + 1 : oV + nV, C % COMOVING_MOMENTUM_U ( 3 ) ), &
        H_Dim => Value_C ( oV + 1 : oV + nV, &
                           C % COMOVING_MOMENTUM_U ( iDimension ) ), &
        FF    => Value_C ( oV + 1 : oV + nV, C % FLUX_FACTOR ), &
        SF    => Value_C ( oV + 1 : oV + nV, C % STRESS_FACTOR ), &
        V_1   => Value_C ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 1 ) ), &
        V_2   => Value_C ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 2 ) ), &
        V_3   => Value_C ( oV + 1 : oV + nV, C % FLUID_VELOCITY_U ( 3 ) ), &
        V_Dim => Value_C ( oV + 1 : oV + nV, &
                           C % FLUID_VELOCITY_U ( iDimension ) ) )

    call ComputeRawFluxesKernel &
           ( F_E, F_S_1, F_S_2, F_S_3, F_S_Dim, J, H_1, H_2, H_3, H_Dim, &
             FF, SF, V_1, V_2, V_3, V_Dim, M_DD_22, M_DD_33, iDimension )

    end associate !-- F_E, etc.
    end associate !-- M_DD_33, etc.

  end subroutine ComputeRawFluxes
  

  subroutine ComputeDiffusionFactor_HLL ( DF_I, Grid, C, iDimension )

    type ( VariableGroupForm ), intent ( inout ) :: &
      DF_I
    class ( * ), intent ( in ), target :: &
      Grid
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iV, &
      iD, &
      iEnergy, &
      iMomentum
    real ( KDR ), dimension ( : ), allocatable :: &
      M_DD_11
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dX, &
      M_DD, &
      TO, &
      SF, &
      DF_I_E
    class ( GeometryFlatForm ), pointer :: &
      G

    call C % SetDiffusionFactorUnity ( DF_I % Value )

    if ( .not. associated ( C % Interactions ) ) &
      return

    select type ( Grid )
    class is ( Chart_SL_Template )

      G => Grid % Geometry ( )

      call Grid % SetVariablePointer &
             ( G % Value ( :, G % WIDTH ( iDimension ) ), dX )

      select case ( iDimension )
      case ( 1 )
        allocate ( M_DD_11 ( C % nValues ) )
        !$OMP parallel do private ( iV )
        do iV = 1, C % nValues
          M_DD_11 ( iV )  =  1.0_KDR
        end do !-- iV
        !$OMP end parallel do
        call Grid % SetVariablePointer ( M_DD_11, M_DD )
      case ( 2 )
        call Grid % SetVariablePointer &
               ( G % Value ( :, G % METRIC_DD_22 ), M_DD )
      case ( 3 )
        call Grid % SetVariablePointer &
               ( G % Value ( :, G % METRIC_DD_33 ), M_DD )
      end select !-- iDimension

      associate ( I => C % Interactions )
      call Grid % SetVariablePointer &
             ( I % Value ( :, I % OPACITY_H ), TO )
      end associate !-- I

      call Grid % SetVariablePointer &
             ( C % Value ( :, C % STRESS_FACTOR ), SF )

      call Search ( C % iaConserved, C % CONSERVED_ENERGY, iEnergy )
      call Grid % SetVariablePointer &
             ( DF_I % Value ( :, iEnergy ), DF_I_E )

      call ComputeDiffusionFactor_HLL_CSL &
             ( DF_I_E, SF, TO, M_DD, dX, iDimension, &
               Grid % nGhostLayers ( iDimension ) )

      do iD = 1, 3
        call Search &
               ( C % iaConserved, C % CONSERVED_MOMENTUM_D ( iD ), iMomentum )
        call Copy ( DF_I % Value ( :, iEnergy ), &
                    DF_I % Value ( :, iMomentum ) )
      end do !-- iD

      call Copy ( DF_I % Value ( :, iEnergy ), &
                  C % Value ( :, C % DIFFUSION_FACTOR_E ) )

      select case ( iDimension )
      case ( 1 )
        deallocate ( M_DD_11 ) 
      end select !-- iDimension

    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'RadiationMoments_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeDiffusionFactor_HLL', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    nullify ( G )
    nullify ( dX, M_DD, TO, SF, DF_I_E )

  end subroutine ComputeDiffusionFactor_HLL


  subroutine InitializeBasics &
               ( RM, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
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

    if ( RM % Type == '' ) &
      RM % Type = 'RadiationMoments'

    Name = 'Radiation'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = RM % N_FIELDS_TEMPLATE
    if ( RM % N_FIELDS == 0 ) &
      RM % N_FIELDS = oF + RM % N_FIELDS_RM

    RM % COMOVING_ENERGY       =  oF +  1
    RM % CONSERVED_ENERGY      =  oF +  2
    RM % COMOVING_MOMENTUM_U   =  oF + [ 3,  4,  5 ]
    RM % CONSERVED_MOMENTUM_D  =  oF + [ 6,  7,  8 ]
    RM % FLUID_VELOCITY_U      =  oF + [ 9, 10, 11 ]
    RM % FLUX_FACTOR           =  oF + 12
    RM % STRESS_FACTOR         =  oF + 13
    RM % COMOVING_ENERGY_EQ    =  oF + 14
    RM % DIFFUSION_FACTOR_E    =  oF + 15

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( RM % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + RM % N_FIELDS_RM ) &
      = [ 'ComovingEnergy       ', &
          'ConservedEnergy      ', &
          'ComovingMomentum_U_1 ', &
          'ComovingMomentum_U_2 ', &
          'ComovingMomentum_U_3 ', &
          'ConservedMomentum_D_1', &
          'ConservedMomentum_D_2', &
          'ConservedMomentum_D_3', &
          'FluidVelocity_U_1    ', &
          'FluidVelocity_U_2    ', &
          'FluidVelocity_U_3    ', &
          'FluxFactor           ', &
          'StressFactor         ', &
          'ComovingEnergy_EQ    ', &
          'DiffusionFactor_E    ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( RM % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = RM % N_VECTORS_TEMPLATE
    if ( RM % N_VECTORS == 0 ) &
      RM % N_VECTORS = oV + RM % N_VECTORS_RM

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( RM % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + RM % N_VECTORS_RM ) &
      = [ 'ComovingMomentum_U ', &
          'ConservedMomentum_D', &
          'FluidVelocity_U    ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + RM % N_VECTORS_RM + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( RM % N_VECTORS ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( RM % COMOVING_MOMENTUM_U )
    call VectorIndices ( oV + 2 ) % Initialize ( RM % CONSERVED_MOMENTUM_D )
    call VectorIndices ( oV + 3 ) % Initialize ( RM % FLUID_VELOCITY_U )

  end subroutine InitializeBasics


  subroutine SetUnits &
               ( VariableUnit, RM, Velocity_U_Unit, MomentumDensity_U_Unit, &
                 MomentumDensity_D_Unit, EnergyDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit

    integer ( KDI ) :: &
      iD

    VariableUnit ( RM % COMOVING_ENERGY )  = EnergyDensityUnit
    VariableUnit ( RM % CONSERVED_ENERGY ) = EnergyDensityUnit

    do iD = 1, 3
      VariableUnit ( RM % COMOVING_MOMENTUM_U ( iD ) ) &
        = MomentumDensity_U_Unit ( iD )
      VariableUnit ( RM % CONSERVED_MOMENTUM_D ( iD ) ) &
        = MomentumDensity_D_Unit ( iD )      
      VariableUnit ( RM % FLUID_VELOCITY_U ( iD ) ) &
        = Velocity_U_Unit ( iD )
    end do

    VariableUnit ( RM % COMOVING_ENERGY_EQ ) = EnergyDensityUnit

  end subroutine SetUnits


  subroutine ComputeConservedEnergyMomentum &
               ( E, S_1, S_2, S_3, J, H_1, H_2, H_3, FF, SF, &
                 M_DD_22, M_DD_33, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3, &
      FF, SF
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33, &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ), dimension ( 3 ) :: &
      K_U_Dim_D
    real ( KDR ), dimension ( size ( E ) ) :: &
      H

    nValues = size ( E )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( J ( iV )  <  0.0_KDR ) &
        J ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      H ( iV )  =  sqrt ( H_1 ( iV ) ** 2  &
                          +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                          +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( H ( iV )  >  J ( iV ) ) then
        H_1 ( iV )  =  ( H_1 ( iV )  /  H ( iV ) )  *  J ( iV )
        H_2 ( iV )  =  ( H_2 ( iV )  /  H ( iV ) )  *  J ( iV )
        H_3 ( iV )  =  ( H_3 ( iV )  /  H ( iV ) )  *  J ( iV )
      end if
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      E ( iV )  &
        =  J ( iV )  &
           +  2.0_KDR * (                       V_1 ( iV )  *  H_1 ( iV )  &
                          +  M_DD_22 ( iV )  *  V_2 ( iV )  *  H_2 ( iV )  &
                          +  M_DD_33 ( iV )  *  V_3 ( iV )  *  H_3 ( iV ) )
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      call ComputeMomentFactors &
             ( SF ( iV ), FF ( iV ), J ( iV ), H_1 ( iV ), H_2 ( iV ), &
               H_3 ( iV ), M_DD_22 ( iV ), M_DD_33 ( iV ) )

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( 1 ), J ( iV ), H_1 ( iV ), H_2 ( iV ), &
               H_3 ( iV ), H_1 ( iV ), FF ( iV ), SF ( iV ), &
               M_DD_22 ( iV ), M_DD_33 ( iV ) )
      S_1 ( iV )  =  H_1 ( iV )  +  J ( iV ) * V_1 ( iV ) &
                     +  K_U_Dim_D ( 1 )  *  V_1 ( iV )  &
                     +  K_U_Dim_D ( 2 )  *  V_2 ( iV )  &
                     +  K_U_Dim_D ( 3 )  *  V_3 ( iV )

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( 2 ), J ( iV ), H_1 ( iV ), H_2 ( iV ), &
               H_3 ( iV ), H_2 ( iV ), FF ( iV ), SF ( iV ), &
               M_DD_22 ( iV ), M_DD_33 ( iV ) )
      S_2 ( iV )  =  M_DD_22 ( iV )  &
                     *  ( H_2 ( iV )  +  J ( iV ) * V_2 ( iV )  &
                          +  K_U_Dim_D ( 1 )  *  V_1 ( iV )  &
                          +  K_U_Dim_D ( 2 )  *  V_2 ( iV )  &
                          +  K_U_Dim_D ( 3 )  *  V_3 ( iV ) )

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( 3 ), J ( iV ), H_1 ( iV ), H_2 ( iV ), &
               H_3 ( iV ), H_3 ( iV ), FF ( iV ), SF ( iV ), &
               M_DD_22 ( iV ), M_DD_33 ( iV ) )
      S_3 ( iV )  =  M_DD_33 ( iV )  &
                     *  ( H_3 ( iV )  +  J ( iV ) * V_3 ( iV )  &
                          +  K_U_Dim_D ( 1 )  *  V_1 ( iV )  &
                          +  K_U_Dim_D ( 2 )  *  V_2 ( iV )  &
                          +  K_U_Dim_D ( 3 )  *  V_3 ( iV ) )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedEnergyMomentum


  subroutine ComputeComovingEnergyMomentum &
               ( J, H_1, H_2, H_3, E, S_1, S_2, S_3, FF, SF, RM, &
                 M_DD_22, M_DD_33, M_UU_22, M_UU_33, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      S_1, S_2, S_3, &
      FF, SF
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33, &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Delta_J_J, Delta_H_H
    real ( KDR ), dimension ( size ( E ) ) :: &
      H
    logical ( KDL ) :: &
      Success

    nValues  =  size ( J )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( E ( iV )  >  0.0_KDR ) then
        call ComputeComovingNonlinearSolve &
               ( J ( iV ), H_1 ( iV ), H_2 ( iV ), H_3 ( iV ), FF ( iV ), &
                 SF ( iV ), E ( iV ), S_1 ( iV ), S_2 ( iV ), S_3 ( iV ), &
                 M_DD_22 ( iV ), M_DD_33 ( iV ), M_UU_22 ( iV ), &
                 M_UU_33 ( iV ), V_1 ( iV ), V_2 ( iV ), V_3 ( iV ), &
                 Success, Delta_J_J, Delta_H_H )
        if ( .not. Success ) then
          call Show ( '>>> ComputeComoving fail', CONSOLE % ERROR )
          call Show ( RM % Name, '>>> Species', CONSOLE % ERROR )
          call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
                      CONSOLE % ERROR )
          call Show ( iV, '>>> iV', CONSOLE % ERROR )
          call Show ( J ( iV ), '>>> J', CONSOLE % ERROR )
          call Show ( H_1 ( iV ), '>>> H_1', CONSOLE % ERROR )
          call Show ( H_2 ( iV ), '>>> H_2', CONSOLE % ERROR )
          call Show ( H_3 ( iV ), '>>> H_3', CONSOLE % ERROR )
          call Show ( Delta_J_J, '>>> Delta_J_J', CONSOLE % ERROR )
          call Show ( Delta_H_H, '>>> Delta_H_H', CONSOLE % ERROR )
        end if
      else
        J   ( iV ) = 0.0_KDR
        H_1 ( iV ) = 0.0_KDR
        H_2 ( iV ) = 0.0_KDR
        H_3 ( iV ) = 0.0_KDR
        E   ( iV ) = 0.0_KDR
        S_1 ( iV ) = 0.0_KDR
        S_2 ( iV ) = 0.0_KDR
        S_3 ( iV ) = 0.0_KDR
        FF  ( iV ) = 0.0_KDR
        SF  ( iV ) = 0.0_KDR
      end if
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      H ( iV )  =  sqrt ( H_1 ( iV ) ** 2  &
                          +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                          +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( H ( iV )  >  J ( iV ) ) then
        H_1 ( iV )  =  ( H_1 ( iV )  /  H ( iV ) )  *  J ( iV )
        H_2 ( iV )  =  ( H_2 ( iV )  /  H ( iV ) )  *  J ( iV )
        H_3 ( iV )  =  ( H_3 ( iV )  /  H ( iV ) )  *  J ( iV )
        S_1 ( iV )  =  H_1 ( iV )
        S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
        S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )
        call ComputeMomentFactors &
               ( SF ( iV ), FF ( iV ), J ( iV ), H_1 ( iV ), H_2 ( iV ), &
                 H_3 ( iV ), M_DD_22 ( iV ), M_DD_33 ( iV ) )
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeComovingEnergyMomentum


  subroutine ComputeEigenspeeds &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 J, H_1, H_2, H_3, M_UU_22, M_UU_33, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, &
      M_UU_22, M_UU_33
    real ( KDR ) :: &
      c  !-- speed of light

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( FEP_1 )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      FEP_1 ( iV )  =  + c
      FEP_2 ( iV )  =  + sqrt ( M_UU_22 ( iV ) )  *  c
      FEP_3 ( iV )  =  + sqrt ( M_UU_33 ( iV ) )  *  c
      FEM_1 ( iV )  =  - c
      FEM_2 ( iV )  =  - sqrt ( M_UU_22 ( iV ) )  *  c
      FEM_3 ( iV )  =  - sqrt ( M_UU_33 ( iV ) )  *  c
    end do
    !$OMP end parallel do

  end subroutine ComputeEigenspeeds


  subroutine ComputeRawFluxesKernel &
               ( F_E, F_S_1, F_S_2, F_S_3, F_S_Dim, J, H_1, H_2, H_3, H_Dim, &
                 FF, SF, V_1, V_2, V_3, V_Dim, M_DD_22, M_DD_33, iDim )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_E, &
      F_S_1, F_S_2, F_S_3, F_S_Dim
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, H_Dim, &
      FF, SF, &
      V_1, V_2, V_3, V_Dim, &
      M_DD_22, M_DD_33
    integer ( KDI ), intent ( in ) :: &
      iDim

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ), dimension ( 3 ) :: &
      K_U_Dim_D

    nValues = size ( J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( iDim ), J ( iV ), H_1 ( iV ), H_2 ( iV ), &
               H_3 ( iV ), H_Dim ( iV ), FF ( iV ), SF ( iV ), &
               M_DD_22 ( iV ), M_DD_33 ( iV ) )
      
      F_E   ( iV )  =  H_Dim ( iV )  +  J ( iV ) * V_Dim ( iV )  &
                       +  K_U_Dim_D ( 1 )  *  V_1 ( iV )  &
                       +  K_U_Dim_D ( 2 )  *  V_2 ( iV )  &
                       +  K_U_Dim_D ( 3 )  *  V_3 ( iV )

      F_S_1 ( iV )  =  K_U_Dim_D ( 1 )  &
                       +  H_Dim ( iV )  *  V_1 ( iV )  &
                       +  V_Dim ( iV )  *  H_1 ( iV )

      F_S_2 ( iV )  =  K_U_Dim_D ( 2 )  &
                       +  M_DD_22 ( iV )  *  ( H_Dim ( iV )  *  V_2 ( iV )  &
                                               +  V_Dim ( iV )  *  H_2 ( iV ) )

      F_S_3 ( iV )  =  K_U_Dim_D ( 3 )  &
                       +  M_DD_33 ( iV )  *  ( H_Dim ( iV )  *  V_3 ( iV )  &
                                               +  V_Dim ( iV )  *  H_3 ( iV ) )

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRawFluxesKernel


  subroutine ComputeDiffusionFactor_HLL_CSL &
               ( DF_I_E, SF, TO, M_DD, dX, iD, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      DF_I_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      SF, &
      TO, &
      M_DD, &
      dX
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV   

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    real ( KDR ) :: &
      DF_L, DF_R

    lV = 1
    where ( shape ( DF_I_E ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( DF_I_E ) > 1 )
      uV = shape ( DF_I_E ) - oV
    end where
    uV ( iD ) = size ( DF_I_E, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = -1
          
    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          DF_L  &
            =  SF ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
               /  max (  sqrt ( M_DD ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                         *  dX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )   &
                         *  TO ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ),  &
                         tiny ( 0.0_KDR ) )

          DF_R  &
            =  SF ( iV, jV, kV )  &
               /  max (  sqrt ( M_DD ( iV, jV, kV ) ) &
                         *  dX ( iV, jV, kV )   &
                         *  TO ( iV, jV, kV ),  &
                         tiny ( 0.0_KDR ) )

          DF_I_E ( iV, jV, kV )  =  min ( 1.0_KDR, max ( DF_L, DF_R ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
  end subroutine ComputeDiffusionFactor_HLL_CSL


  subroutine ComputeMomentFactors &
               ( SF, FF, J, H_1, H_2, H_3, M_DD_22, M_DD_33 )

    real ( KDR ), intent ( inout ) :: &
      SF, &
      FF
    real ( KDR ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, &
      M_DD_22, M_DD_33

    real ( KDR ) :: &
      H

    H  =  sqrt ( H_1 ** 2  &
                 +  M_DD_22  *  H_2 ** 2  &
                 +  M_DD_33  *  H_3 ** 2 )
    
    if ( H  <=  J ) then
      FF  =  H  /  max ( J, tiny ( 0.0_KDR ) )
    else
      FF  =  1.0_KDR
    end if

    SF  =  1.0_KDR / 3.0_KDR &
           + 2.0_KDR / 3.0_KDR &
             * ( FF ** 2  /  5.0_KDR  &
                 * ( 3.0_KDR  -  FF  +  3.0_KDR  *  FF ** 2 ) )

  end subroutine ComputeMomentFactors


  subroutine ComputeComovingStress_D &
               ( K_1, K_2, K_3, K_Dim, J, H_1, H_2, H_3, H_Dim, FF, SF, &
                 M_DD_22, M_DD_33 )

    real ( KDR ), intent ( inout ) :: &
      K_1, K_2, K_3, K_Dim
    real ( KDR ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, H_Dim, &
      FF, SF, &
      M_DD_22, M_DD_33

    real ( KDR ) :: &
      H_Sq

    H_Sq  =  max (                H_1 ** 2  &
                   +  M_DD_22  *  H_2 ** 2  &
                   +  M_DD_33  *  H_3 ** 2, &
                   tiny ( 0.0_KDR ) )   

    K_1  =  0.5_KDR  *  ( 3.0_KDR * SF  -  1.0_KDR )  &
            *  H_Dim  *            H_1  /  H_Sq   *  J

    K_2  =  0.5_KDR  *  ( 3.0_KDR * SF  -  1.0_KDR )  &
            *  H_Dim  *  M_DD_22 * H_2  /  H_Sq   *  J

    K_3  =  0.5_KDR  *  ( 3.0_KDR * SF  -  1.0_KDR )  &
            *  H_Dim  *  M_DD_33 * H_3  /  H_Sq   *  J

    K_Dim  =  K_Dim  +  0.5_KDR * ( 1.0_KDR  -  SF )  *  J

  end subroutine ComputeComovingStress_D


  subroutine ComputeComovingNonlinearSolve &
               ( J, H_1, H_2, H_3, FF, SF, E, S_1, S_2, S_3, &
                 M_DD_22, M_DD_33, M_UU_22, M_UU_33, V_1, V_2, V_3, Success, &
                 Delta_J_J, Delta_H_H )

    real ( KDR ), intent ( inout ) :: &
      J, H_1, H_2, H_3, &
      FF, SF
    real ( KDR ), intent ( in ) :: &
      E, S_1, S_2, S_3, &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33, &
      V_1, V_2, V_3
    logical ( KDL ), intent ( out ) :: &
      Success
    real ( KDR ), intent ( out ) :: &
      Delta_J_J, Delta_H_H

    integer ( KDI ) :: &
      iS, &  !-- iSolve
      MaxIterations
    real ( KDR ) :: &
      J_Old, H_Old_1, H_Old_2, H_Old_3, &
      Norm_H, NormDelta_H
    real ( KDR ), dimension ( 3 ) :: &
      K_U_Dim_D

    MaxIterations  =  20

    Success  =  .false.

    do iS = 1, MaxIterations

      J_Old    =  J
      H_Old_1  =  H_1
      H_Old_2  =  H_2
      H_Old_3  =  H_3 
                            
      J  =  E  -  2.0_KDR  * (                V_1  *  H_Old_1  &
                               +  M_DD_22  *  V_2  *  H_Old_2  &
                               +  M_DD_33  *  V_3  *  H_Old_3 )

      call ComputeMomentFactors &
             ( SF, FF, J_Old, H_Old_1, H_Old_2, H_Old_3, M_DD_22, M_DD_33  )

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( 1 ), J_Old, H_Old_1, H_Old_2, H_Old_3, H_Old_1, &
               FF, SF, M_DD_22, M_DD_33 )
      H_1  =  S_1  -  J_Old * V_1   &
                   -  K_U_Dim_D ( 1 )  *  V_1   &
                   -  K_U_Dim_D ( 2 )  *  V_2   &
                   -  K_U_Dim_D ( 3 )  *  V_3 

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( 2 ), J_Old, H_Old_1, H_Old_2, H_Old_3, H_Old_2, &
               FF, SF, M_DD_22, M_DD_33 )
      H_2  =  M_UU_22 * S_2  -  J_Old  * V_2   &
                             -  K_U_Dim_D ( 1 )  *  V_1   &
                             -  K_U_Dim_D ( 2 )  *  V_2   &
                             -  K_U_Dim_D ( 3 )  *  V_3 

      call ComputeComovingStress_D &
             ( K_U_Dim_D ( 1 ), K_U_Dim_D ( 2 ), K_U_Dim_D ( 3 ), &
               K_U_Dim_D ( 3 ), J_Old, H_Old_1, H_Old_2, H_Old_3, H_Old_3, &
               FF, SF, M_DD_22, M_DD_33 )
      H_3  =  M_UU_33 * S_3  -  J_Old * V_3   &
                             -  K_U_Dim_D ( 1 )  *  V_1   &
                             -  K_U_Dim_D ( 2 )  *  V_2   &
                             -  K_U_Dim_D ( 3 )  *  V_3 

      Norm_H       =  sqrt ( (                H_1 ** 2  &
                               +  M_DD_22  *  H_2 ** 2  &
                               +  M_DD_33  *  H_3 ** 2 ) )

      NormDelta_H  =  sqrt ( (                ( H_1 - H_Old_1 ) ** 2  &
                               +  M_DD_22  *  ( H_2 - H_Old_2 ) ** 2  &
                               +  M_DD_33  *  ( H_3 - H_Old_3 ) ** 2 ) )

      Delta_J_J  =  abs ( J - J_Old )  /  max ( abs ( J ), tiny ( 0.0_KDR ) )
      Delta_H_H  =  NormDelta_H  /  max ( Norm_H, tiny ( 0.0_KDR ) )

      if ( Delta_J_J  <  1.0e-10_KDR &
           .and. Delta_H_H  <  1.0e-10_KDR ) &
      then
        Success  =  .true.
        exit
      end if

    end do !-- iS

  end subroutine ComputeComovingNonlinearSolve


end module RadiationMoments_Form
