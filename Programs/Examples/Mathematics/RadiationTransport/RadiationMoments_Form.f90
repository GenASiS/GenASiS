module RadiationMoments_Form

  use Basics
  use Mathematics
  use Interactions_Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_RM = 4, &
      N_CONSERVED_RM = 4, &
      N_FIELDS_RM    = 10, &
      N_VECTORS_RM   = 2

  type, public, extends ( CurrentTemplate ) :: RadiationMomentsForm
    integer ( KDI ) :: &
      N_PRIMITIVE_RM            = N_PRIMITIVE_RM, &
      N_CONSERVED_RM            = N_CONSERVED_RM, &
      N_FIELDS_RM               = N_FIELDS_RM, &
      N_VECTORS_RM              = N_VECTORS_RM, &
      COMOVING_ENERGY_DENSITY   = 0, &
      CONSERVED_ENERGY_DENSITY  = 0, &
      FLUX_FACTOR               = 0, &
      VARIABLE_EDDINGTON_FACTOR = 0
    integer ( KDI ), dimension ( 3 ) :: &
      COMOVING_MOMENTUM_DENSITY_U  = 0, &
      CONSERVED_MOMENTUM_DENSITY_D = 0
    class ( InteractionsTemplate ), pointer :: &
      Interactions => null ( )
  contains
    procedure, public, pass :: &
      InitializeAllocate_RM
    generic, public :: &
      Initialize => InitializeAllocate_RM
    procedure, public, pass :: &
      SetInteractions
    procedure, public, pass ( C ) :: &
      ComputeRawFluxes
    procedure, public, pass ( C ) :: &
      ComputeDiffusionFactor
    procedure, public, pass :: &
      SetOutput
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      ComputeFromPrimitiveCommon
    procedure, public, pass ( C ) :: &
      ComputeFromConservedCommon
  end type RadiationMomentsForm

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeConservedEnergyMomentum, &
      ComputePrimitiveEnergyMomentum, &
      ComputeEigenspeeds, &
      ComputeVariableEddingtonFactor, &
      ComputeRawFluxesKernel, &
      ComputeDiffusionFactor_CSL

  public :: &
    ApplyRelaxation_Interactions, &
    ApplySourcesCurvilinear_RadiationMoments

    private :: &
      ApplyRelaxationKernel

contains


  subroutine InitializeAllocate_RM &
               ( RM, VelocityUnit, EnergyDensityUnit, nValues, &
                 VariableOption, VectorOption, NameOption, ClearOption, &
                 UnitOption, VectorIndicesOption )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit
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

    call SetUnits ( VariableUnit, RM, VelocityUnit, EnergyDensityUnit )

    call RM % InitializeTemplate &
           ( VelocityUnit, nValues, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )
    
  end subroutine InitializeAllocate_RM


  subroutine SetInteractions ( RM, Interactions )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    class ( InteractionsTemplate ), intent ( in ), target :: &
      Interactions

    RM % Interactions => Interactions

  end subroutine SetInteractions


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
           ( C % iaConserved, C % CONSERVED_ENERGY_DENSITY, iEnergy )
    call Search &
           ( C % iaConserved, C % CONSERVED_MOMENTUM_DENSITY_D ( 1 ), &
             iMomentum ( 1 ) )
    call Search &
           ( C % iaConserved, C % CONSERVED_MOMENTUM_DENSITY_D ( 2 ), &
             iMomentum ( 2 ) )
    call Search &
           ( C % iaConserved, C % CONSERVED_MOMENTUM_DENSITY_D ( 3 ), &
             iMomentum ( 3 ) )
    
    associate &
      ( M_DD_22 => Value_G ( oV + 1 : oV + nV, G % METRIC_DD_22 ), &
        M_DD_33 => Value_G ( oV + 1 : oV + nV, G % METRIC_DD_33 ) )
    associate &
      ( F_E   => RawFlux ( oV + 1 : oV + nV, iEnergy ), &
        F_S_1 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 1 ) ), &
        F_S_2 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 2 ) ), &
        F_S_3 => RawFlux ( oV + 1 : oV + nV, iMomentum ( 3 ) ), &
        F_S_Dim => RawFlux ( oV + 1 : oV + nV, iMomentum ( iDimension ) ), & 
        J     => Value_C ( oV + 1 : oV + nV, C % COMOVING_ENERGY_DENSITY ), &
        H_1   => Value_C ( oV + 1 : oV + nV, &
                         C % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
        H_2   => Value_C ( oV + 1 : oV + nV, &
                         C % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
        H_3   => Value_C ( oV + 1 : oV + nV, &
                         C % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
        H_Dim => Value_C ( oV + 1 : oV + nV, &
                         C % COMOVING_MOMENTUM_DENSITY_U ( iDimension ) ), &
        VEF => Value_C ( oV + 1 : oV + nV, C % VARIABLE_EDDINGTON_FACTOR ) )

    call ComputeRawFluxesKernel &
           ( F_E, F_S_1, F_S_2, F_S_3, F_S_Dim, J, H_1, H_2, H_3, H_Dim, VEF, &
             M_DD_22, M_DD_33 )

    end associate !-- F_E, etc.
    end associate !-- M_DD_33, etc.

  end subroutine ComputeRawFluxes
  

  subroutine ComputeDiffusionFactor ( DF_I, Grid, C, iDimension )

    type ( VariableGroupForm ), intent ( inout ) :: &
      DF_I
    class ( * ), intent ( in ), target :: &
      Grid
    class ( RadiationMomentsForm ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iV
    real ( KDR ), dimension ( : ), allocatable :: &
      M_DD_11
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dX, &
      M_DD, &
      TO, &
      VEF, &
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
             ( I % Value ( :, I % TRANSPORT_OPACITY ), TO )
      end associate !-- I

      call Grid % SetVariablePointer &
             ( C % Value ( :, C % VARIABLE_EDDINGTON_FACTOR ), VEF )
      call Grid % SetVariablePointer &
             ( DF_I % Value ( :, 1 ), DF_I_E )

      call ComputeDiffusionFactor_CSL &
             ( DF_I_E, VEF, TO, M_DD, dX, iDimension, &
               Grid % nGhostLayers ( iDimension ) )

      select case ( iDimension )
      case ( 1 )
        deallocate ( M_DD_11 ) 
      end select !-- iDimension

    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'RadiationMoments_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeDiffusionFactor', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    nullify ( G )
    nullify ( dX, M_DD, TO, VEF, DF_I_E )

  end subroutine ComputeDiffusionFactor


  subroutine SetOutput ( RM, Output )

    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices

    call VectorIndices ( 1 ) % Initialize ( RM % COMOVING_MOMENTUM_DENSITY_U )
    call Output % Initialize &
           ( RM, iaSelectedOption = [ RM % COMOVING_ENERGY_DENSITY, &
                                      RM % COMOVING_MOMENTUM_DENSITY_U, &
                                      RM % FLUX_FACTOR, &
                                      RM % VARIABLE_EDDINGTON_FACTOR ], &
             VectorOption = [ 'ComovingMomentumDensity        ' ], &
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
      
    associate &
      ( RMV => Value_C, &
        GV => Value_G )

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
        J     => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY_DENSITY ), &
        H_1 => RMV ( oV + 1 : oV + nV, &
                     C % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
        H_2 => RMV ( oV + 1 : oV + nV, &
                     C % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
        H_3 => RMV ( oV + 1 : oV + nV, &
                     C % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
        E   => RMV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY_DENSITY ), &
        S_1 => RMV ( oV + 1 : oV + nV, &
                     C % CONSERVED_MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2 => RMV ( oV + 1 : oV + nV, &
                     C % CONSERVED_MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3 => RMV ( oV + 1 : oV + nV, &
                     C % CONSERVED_MOMENTUM_DENSITY_D ( 3 ) ), &
        FF  => RMV ( oV + 1 : oV + nV, C % FLUX_FACTOR ), &
        VEF => RMV ( oV + 1 : oV + nV, C % VARIABLE_EDDINGTON_FACTOR ) )

    call ComputeConservedEnergyMomentum &
           ( E, S_1, S_2, S_3, J, H_1, H_2, H_3, M_DD_22, M_DD_33 )
    call ComputeEigenspeeds &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, J, H_1, H_2, H_3, &
             M_UU_22, M_UU_33, CONSTANT % SPEED_OF_LIGHT )
    call ComputeVariableEddingtonFactor &
           ( VEF, FF, J, H_1, H_2, H_3, M_DD_22, M_DD_33 )

    if ( associated ( C % Value, Value_C ) &
         .and. associated ( C % Interactions ) ) &
      call C % Interactions % Compute ( )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- RMV, etc.

  end subroutine ComputeFromPrimitiveCommon


  subroutine ComputeFromConservedCommon &
               ( C, G, Value, nValuesOption, oValueOption )

    class ( RadiationMomentsForm ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( :, : ), intent ( inout ), target :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    integer ( KDI ) :: &
      oV, &  !-- oValue
      nV     !-- nValues
      
    associate &
      ( RMV => Value, &
        GV  => G % Value )

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
        J   => RMV ( oV + 1 : oV + nV, C % COMOVING_ENERGY_DENSITY ), &
        H_1 => RMV ( oV + 1 : oV + nV, &
                     C % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
        H_2 => RMV ( oV + 1 : oV + nV, &
                     C % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
        H_3 => RMV ( oV + 1 : oV + nV, &
                     C % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
        E   => RMV ( oV + 1 : oV + nV, C % CONSERVED_ENERGY_DENSITY ), &
        S_1 => RMV ( oV + 1 : oV + nV, &
                     C % CONSERVED_MOMENTUM_DENSITY_D ( 1 ) ), &
        S_2 => RMV ( oV + 1 : oV + nV, &
                     C % CONSERVED_MOMENTUM_DENSITY_D ( 2 ) ), &
        S_3 => RMV ( oV + 1 : oV + nV, &
                     C % CONSERVED_MOMENTUM_DENSITY_D ( 3 ) ), &
        FF  => RMV ( oV + 1 : oV + nV, C % FLUX_FACTOR ), &
        VEF => RMV ( oV + 1 : oV + nV, C % VARIABLE_EDDINGTON_FACTOR ) )

    call ComputePrimitiveEnergyMomentum &
           ( J, H_1, H_2, H_3, E, S_1, S_2, S_3, M_DD_22, M_DD_33, &
             M_UU_22, M_UU_33 )
    call ComputeEigenspeeds &
           ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, J, H_1, H_2, H_3, &
             M_UU_22, M_UU_33, CONSTANT % SPEED_OF_LIGHT )
    call ComputeVariableEddingtonFactor &
           ( VEF, FF, J, H_1, H_2, H_3, M_DD_22, M_DD_33 )

    if ( associated ( C % Value, Value ) &
         .and. associated ( C % Interactions ) ) &
      call C % Interactions % Compute ( )

    end associate !-- FEP_1, etc.
    end associate !-- M_DD_22, etc.
    end associate !-- RMV, etc.

  end subroutine ComputeFromConservedCommon


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
      oV, &  !-- oVector
      oP, &  !-- oPrimitive
      oC     !-- oConserved

    if ( RM % Type == '' ) RM % Type = 'RadiationMoments'

    Name = 'RadiationMoments'
    if ( present ( NameOption ) ) Name = NameOption

    !-- variable indices

    oF = RM % N_FIELDS_TEMPLATE
    if ( RM % N_FIELDS == 0 ) &
      RM % N_FIELDS = oF + RM % N_FIELDS_RM

    RM % COMOVING_ENERGY_DENSITY      = oF + 1
    RM % CONSERVED_ENERGY_DENSITY     = oF + 2
    RM % COMOVING_MOMENTUM_DENSITY_U  = oF + [ 3, 4, 5 ]
    RM % CONSERVED_MOMENTUM_DENSITY_D = oF + [ 6, 7, 8 ]
    RM % FLUX_FACTOR                  = oF + 9
    RM % VARIABLE_EDDINGTON_FACTOR    = oF + 10

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( RM % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + RM % N_FIELDS_RM ) &
      = [ 'ComovingEnergyDensity          ', &
          'ConservedEnergyDensity         ', &
          'ComovingMomentumDensity_1      ', &
          'ComovingMomentumDensity_2      ', &
          'ComovingMomentumDensity_3      ', &
          'ConservedMomentumDensity_1     ', &
          'ConservedMomentumDensity_2     ', &
          'ConservedMomentumDensity_3     ', &
          'FluxFactor                     ', &
          'VariableEddingtonFactor        ']
          
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
      = [ 'ComovingMomentumDensity        ', &
          'ConservedMomentumDensity       ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + RM % N_VECTORS_RM + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( RM % N_VECTORS ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize &
           ( RM % COMOVING_MOMENTUM_DENSITY_U )
    call VectorIndices ( oV + 2 ) % Initialize &
           ( RM % CONSERVED_MOMENTUM_DENSITY_D )

    !-- select primitive, conserved

    oP = RM % N_PRIMITIVE_TEMPLATE
    oC = RM % N_CONSERVED_TEMPLATE

    if ( .not. allocated ( RM % iaPrimitive ) ) then
      RM % N_PRIMITIVE = oP + RM % N_PRIMITIVE_RM
      allocate ( RM % iaPrimitive ( RM % N_PRIMITIVE ) )
    end if
    RM % iaPrimitive ( oP + 1 : oP + RM % N_PRIMITIVE_RM ) &
      = [ RM % COMOVING_ENERGY_DENSITY, RM % COMOVING_MOMENTUM_DENSITY_U ]

    if ( .not. allocated ( RM % iaConserved ) ) then
      RM % N_CONSERVED = oC + RM % N_CONSERVED_RM
      allocate ( RM % iaConserved ( RM % N_CONSERVED ) )
    end if
    RM % iaConserved ( oC + 1 : oC + RM % N_CONSERVED_RM ) &
      = [ RM % CONSERVED_ENERGY_DENSITY, RM % CONSERVED_MOMENTUM_DENSITY_D ]
    
  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, RM, VelocityUnit, EnergyDensityUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      VelocityUnit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit

    integer ( KDI ) :: &
      iD

    VariableUnit ( RM % COMOVING_ENERGY_DENSITY )  = EnergyDensityUnit
    VariableUnit ( RM % CONSERVED_ENERGY_DENSITY ) = EnergyDensityUnit

    do iD = 1, 3
      VariableUnit ( RM % COMOVING_MOMENTUM_DENSITY_U ( iD ) ) &
        = EnergyDensityUnit / VelocityUnit ( iD )
      VariableUnit ( RM % CONSERVED_MOMENTUM_DENSITY_D ( iD ) ) &
        = EnergyDensityUnit / VelocityUnit ( iD )      
    end do

    VariableUnit ( RM % VARIABLE_EDDINGTON_FACTOR ) = UNIT % IDENTITY

  end subroutine SetUnits


  subroutine ComputeConservedEnergyMomentum &
               ( E, S_1, S_2, S_3, J, H_1, H_2, H_3, &
                 M_DD_22, M_DD_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33

    integer ( KDI ) :: &
      iV, &
      nValues
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
      E ( iV )  =  J ( iV )
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
      S_1 ( iV )  =  H_1 ( iV )
      S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
      S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeConservedEnergyMomentum


  subroutine ComputePrimitiveEnergyMomentum &
               ( J, H_1, H_2, H_3, E, S_1, S_2, S_3, M_DD_22, M_DD_33, &
                 M_UU_22, M_UU_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      J, &
      H_1, H_2, H_3
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      E, &
      S_1, S_2, S_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_DD_22, M_DD_33, &
      M_UU_22, M_UU_33

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ), dimension ( size ( E ) ) :: &
      H

    nValues = size ( J )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( E ( iV )  >  0.0_KDR ) then
        J ( iV )    =  E ( iV )
        H_1 ( iV )  =  S_1 ( iV )
        H_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV )
        H_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV )
      else
        J   ( iV ) = 0.0_KDR
        H_1 ( iV ) = 0.0_KDR
        H_2 ( iV ) = 0.0_KDR
        H_3 ( iV ) = 0.0_KDR
        E   ( iV ) = 0.0_KDR
        S_1 ( iV ) = 0.0_KDR
        S_2 ( iV ) = 0.0_KDR
        S_3 ( iV ) = 0.0_KDR
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
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputePrimitiveEnergyMomentum


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


  subroutine ComputeVariableEddingtonFactor &
               ( VEF, FF, J, H_1, H_2, H_3, M_DD_22, M_DD_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      VEF, &
      FF
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, &
      M_DD_22, M_DD_33

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ), dimension ( size ( VEF ) ) :: &
      H

    nValues = size ( VEF )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      H ( iV )  =  sqrt ( H_1 ( iV ) ** 2  &
                          +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                          +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
    end do
    !$OMP end parallel do
    
    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      if ( H ( iV )  <=  J ( iV ) ) then
        FF ( iV )  =  H ( iV )  /  max ( J ( iV ), tiny ( 0.0_KDR ) )
      else
        FF ( iV )  =  1.0_KDR
      end if
    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      VEF ( iV )  &
        =  1.0_KDR / 3.0_KDR &
           + 2.0_KDR / 3.0_KDR &
             * ( FF ( iV ) ** 2  /  5.0_KDR  &
                 * ( 3.0_KDR  -  FF ( iV )  +  3.0_KDR  *  FF ( iV ) ** 2 ) )
    end do
    !$OMP end parallel do

  end subroutine ComputeVariableEddingtonFactor


  subroutine ComputeRawFluxesKernel &
               ( F_E, F_S_1, F_S_2, F_S_3, F_S_Dim, J, H_1, H_2, H_3, H_Dim, &
                 VEF, M_DD_22, M_DD_33 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      F_E, &
      F_S_1, F_S_2, F_S_3, &
      F_S_Dim
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, &
      H_Dim, &
      VEF, &
      M_DD_22, M_DD_33

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ), dimension ( size ( J ) ) :: &
      H_Sq

    nValues = size ( J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      H_Sq ( iV )  =  max ( H_1 ( iV ) ** 2  &
                            +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                            +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2, &
                            tiny ( 0.0_KDR ) )   
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues

      F_E   ( iV )  =   H_Dim ( iV )

      F_S_1 ( iV )  =   0.5_KDR  *  ( 3.0_KDR * VEF ( iV )  -  1.0_KDR )  &
                        *  H_1 ( iV ) / H_Sq ( iV )  &
                        *  J ( iV )  *  H_Dim ( iV )

      F_S_2 ( iV )  =   0.5_KDR  *  ( 3.0_KDR * VEF ( iV )  -  1.0_KDR )  &
                        *  M_DD_22 ( iV )  *  H_2 ( iV )  /  H_Sq ( iV )  &
                        *  J ( iV )  *  H_Dim ( iV )

      F_S_3 ( iV )  =   0.5_KDR  *  ( 3.0_KDR * VEF ( iV )  -  1.0_KDR )  &
                        *  M_DD_33 ( iV )  *  H_3 ( iV )  /  H_Sq ( iV )  &
                        *  J ( iV )  *  H_Dim ( iV )

      F_S_Dim ( iV )  =  F_S_Dim ( iV )  &
                         +  0.5_KDR * ( 1.0_KDR  -  VEF ( iV ) ) * J ( iV )
  
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeRawFluxesKernel


  ! subroutine ComputeDiffusionFactorKernel &
  !      ( DF_E, DF_F_1, DF_F_2, DF_F_3, dX, &
  !        VEF_L, VEF_R, Chi_L, Chi_R, Sigma_L, Sigma_R )
    
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     DF_E, &
  !     DF_F_1, DF_F_2, DF_F_3
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     dX, &
  !     VEF_L, VEF_R, &
  !     Chi_L, CHI_R, &
  !     Sigma_L, Sigma_R

  !   real ( KDR ), dimension ( size ( DF_E ) ) :: &
  !     DF_E_L, DF_E_R
    
  !   DF_E_L   =  VEF_L / max ( dX * ( Sigma_L + Chi_L ), tiny ( 0.0_KDR ) ) 
    
  !   DF_E_R   =  VEF_R / max ( dX * ( Sigma_R + Chi_R ), tiny ( 0.0_KDR ) )
    
  !   DF_E   = min ( DF_E_L, DF_E_R )
  !        ! = max ( DF_E_L, DF_E_R )
  !        ! = ( DF_E_L + DF_E_R ) / 2
  !   DF_E   = min ( 1.0_KDR, DF_E )
  !   DF_F_1 = 1.0_KDR
  !   DF_F_2 = 1.0_KDR
  !   DF_F_3 = 1.0_KDR

  ! end subroutine ComputeDiffusionFactorKernel
  

  subroutine ComputeDiffusionFactor_CSL ( DF_I_E, VEF, TO, M_DD, dX, iD, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      DF_I_E
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      VEF, &
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
            =  VEF ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
               /  max (  sqrt ( M_DD ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                         *  dX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )   &
                         *  TO ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ),  &
                         tiny ( 0.0_KDR ) )

          DF_R  &
            =  VEF ( iV, jV, kV )  &
               /  max (  sqrt ( M_DD ( iV, jV, kV ) ) &
                         *  dX ( iV, jV, kV )   &
                         *  TO ( iV, jV, kV ),  &
                         tiny ( 0.0_KDR ) )

          DF_I_E ( iV, jV, kV )  =  min ( 1.0_KDR, max ( DF_L, DF_R ) )

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
  end subroutine ComputeDiffusionFactor_CSL


  subroutine ApplyRelaxation_Interactions &
               ( S, IncrementExplicit, DampingCoefficient, Current, &
                 TimeStep )

    class ( Step_RK_C_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      Current
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iEnergy, &
      iMomentum_1, &
      iMomentum_2, &
      iMomentum_3

    select type ( RM => Current )
    class is ( RadiationMomentsForm )

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY_DENSITY, &
                  iEnergy )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_DENSITY_D ( 1 ), &
                  iMomentum_1 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_DENSITY_D ( 2 ), &
                  iMomentum_2 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_DENSITY_D ( 3 ), &
                  iMomentum_3 )

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )

    associate ( I => RM % Interactions )

    call ApplyRelaxationKernel &
               ( IncrementExplicit % Value ( :, iEnergy ), &
                 DampingCoefficient % Value ( :, iEnergy ), &
                 DampingCoefficient % Value ( :, iMomentum_1 ), &
                 DampingCoefficient % Value ( :, iMomentum_2 ), &
                 DampingCoefficient % Value ( :, iMomentum_3 ), &
                 Grid % IsProperCell, &
                 I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
                 I % Value ( :, I % EFFECTIVE_OPACITY ), &
                 I % Value ( :, I % TRANSPORT_OPACITY ), &
                 TimeStep )

    end associate !-- I
    end select !-- Grid
    end select !-- RM
    
  end subroutine ApplyRelaxation_Interactions
    

  subroutine ApplySourcesCurvilinear_RadiationMoments &
               ( S, Increment, Current, TimeStep )

    class ( Step_RK_C_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ) :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Current
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iMomentum_1, &
      iMomentum_2
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( RM => Current )
    class is ( RadiationMomentsForm )

    call Search &
           ( RM % iaConserved, &
             RM % CONSERVED_MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search &
           ( RM % iaConserved, &
             RM % CONSERVED_MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )

    select type ( Grid => S % Grid )
    class is ( Chart_SL_Template )

    if ( trim ( Grid % CoordinateSystem ) == 'CARTESIAN' ) &
      return

    G => Grid % Geometry ( )
    
    associate &
      ( M_DD_22 => G % Value ( :, G % METRIC_DD_22 ), &
        M_DD_33 => G % Value ( :, G % METRIC_DD_33 ) )
    
    call ApplySourcesCurvilinearKernel &
           ( Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             Grid % CoordinateSystem, Grid % IsProperCell, &
             RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
             RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
             RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
             RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
             RM % Value ( :, RM % VARIABLE_EDDINGTON_FACTOR ), &
             M_DD_22, M_DD_33, &
             S % dLogVolumeJacobian_dX ( 1 ) % Value, &
             S % dLogVolumeJacobian_dX ( 2 ) % Value, &
             TimeStep, Grid % nDimensions, G % Value ( :, G % CENTER ( 1 ) ) )

    end associate !-- M_DD_22, etc.
    class default
      call Show ( 'Grid type not found', CONSOLE % ERROR )
      call Show ( 'RadiationMoments_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ApplySourcesCurvilinear', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end select !-- RM
    
  end subroutine ApplySourcesCurvilinear_RadiationMoments


  subroutine ApplyRelaxationKernel &
               ( KV_E, DCV_E, DCV_S_1, DCV_S_2, DCV_S_3, IsProperCell, &
                 ED, EO, TO, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_E, &
      DCV_E, &
      DCV_S_1, &
      DCV_S_2, &
      DCV_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      ED, &
      EO, &
      TO
    real ( KDR ), intent ( in ) :: &
      dT

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( KV_E )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      KV_E    ( iV )  =  KV_E ( iV )  +  EO ( iV ) * ED ( iV ) * dT
      DCV_E   ( iV )  =  EO ( iV )
      DCV_S_1 ( iV )  =  TO ( iV )
      DCV_S_2 ( iV )  =  TO ( iV )
      DCV_S_3 ( iV )  =  TO ( iV )
    end do
    !$OMP end parallel do

  end subroutine ApplyRelaxationKernel
  

  subroutine ApplySourcesCurvilinearKernel &
               ( KVM_1, KVM_2, CoordinateSystem, IsProperCell, &
                 J, H_1, H_2, H_3, VEF, M_DD_22, M_DD_33, &
                 dLVJ_dX1, dLVJ_dX2, dT, nDimensions, R )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVM_1, KVM_2
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, &
      VEF, & 
      M_DD_22, M_DD_33, &
      dLVJ_dX1, dLVJ_dX2, &
      R
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      K_22, &
      K_33
    real ( KDR ), dimension ( : ), allocatable :: &
      H

    nV = size ( KVM_1 )

    allocate ( H ( size ( KVM_1 ) ) )
    
    !$OMP parallel do private ( iV )
    do iV = 1, nV
      H ( iV )  =  max ( sqrt ( H_1 ( iV ) ** 2  &
                                +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                                +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 ), &
                    tiny ( 0.0_KDR ) ) 
    end do
    !$OMP end parallel do
   
    select case ( trim ( CoordinateSystem ) )
    case ( 'CYLINDRICAL' )

      !$OMP parallel do private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) cycle

        K_33 &
          = 0.5_KDR * ( 1.0_KDR - VEF ( iV ) ) * J ( iV )  + 0.5_KDR &
            * ( 3 * VEF ( iV ) - 1.0_KDR ) * M_DD_33 ( iV ) &
            *  H_3 ( iV ) ** 2 / H ( iV )

        KVM_1 ( iV ) &
          = KVM_1 ( iV ) + K_33 * dLVJ_dX1 ( iV ) * dT
      end do
      !$OMP end parallel do

    case ( 'SPHERICAL' )

      !$OMP parallel do private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) cycle
          
        K_22 &
          = 0.5_KDR * ( 1.0_KDR - VEF ( iV ) ) * J ( iV ) + 0.5_KDR &
            * ( 3 * VEF ( iV ) - 1.0_KDR ) * M_DD_22 ( iV ) &
            *  H_2 ( iV ) ** 2 / H ( iV )
        K_33 &
          = 0.5_KDR * ( 1.0_KDR - VEF ( iV ) ) * J ( iV ) + 0.5_KDR &
            * ( 3 * VEF ( iV ) - 1.0_KDR ) * M_DD_33 ( iV ) &
            *  H_3 ( iV ) ** 2 / H ( iV )

        KVM_1 ( iV ) &
          =  KVM_1 ( iV )  + 0.5_KDR * ( K_22 + K_33 ) * dLVJ_dX1 ( iV ) * dT
        if ( nDimensions > 1 ) then
          KVM_2 ( iV ) &
            =  KVM_2 ( iV )  + K_33 * dLVJ_dX2 ( iV ) * dT
        end if
      end do
      !$OMP end parallel do
    end select !-- CoordinateSystem

  end subroutine ApplySourcesCurvilinearKernel


end module RadiationMoments_Form
