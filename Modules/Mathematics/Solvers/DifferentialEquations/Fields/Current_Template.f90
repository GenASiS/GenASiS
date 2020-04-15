!-- Current is a template for a 4-current, of the type appearing in
!   conservation laws.

module Current_Template

  use Basics
  use Manifolds
  use Sources_C__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_TEMPLATE     = 0, &
      N_CONSERVED_TEMPLATE     = 0, &
      N_FIELDS_TEMPLATE        = 6, &
      N_VECTORS_TEMPLATE       = 2

  type, public, extends ( StorageForm ), abstract :: CurrentTemplate
    integer ( KDI ) :: &
      IGNORABILITY    = 0, &
      N_PRIMITIVE     = 0, &
      N_CONSERVED     = 0, &
      N_FIELDS        = 0, &
      N_VECTORS       = 0, &
      N_RECONSTRUCTED = 0, &
      N_PRIMITIVE_TEMPLATE = N_PRIMITIVE_TEMPLATE, &
      N_CONSERVED_TEMPLATE = N_CONSERVED_TEMPLATE, &
      N_FIELDS_TEMPLATE    = N_FIELDS_TEMPLATE, &
      N_VECTORS_TEMPLATE   = N_VECTORS_TEMPLATE
    integer ( KDI ), dimension ( 3 ) :: &
      FAST_EIGENSPEED_PLUS  = 0, &
      FAST_EIGENSPEED_MINUS = 0
    integer ( KDI ) :: &  !-- Indices in SolverSpeed storage
      ALPHA_PLUS      = 1, &
      ALPHA_MINUS     = 2, &
      ALPHA_CENTER    = 3, &
      N_SOLVER_SPEEDS = 3
    integer ( KDI ) :: &
      iTimerComputeFromPrimitive, &
      iTimerComputeFromConserved, &
      iTimerComputeEigenspeed, &
      iTimerApply_EOS, &
      iTimerComputeSolverSpeeds_HLL, &
      iTimerComputeFluxes_HLL, &
      iTimerComputeRawFluxes, &
      iTimerComputeDiffusionFactor_HLL
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaConserved, &
      iaReconstructed
    real ( KDR ) :: &
      LimiterParameter
    logical ( KDL ) :: &
      UseLimiter
    character ( LDL ) :: &
      Type = '', &
      ReconstructedType = '', &
      RiemannSolverType = ''
    class ( Sources_C_Form ), pointer :: &
      Sources => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure ( SPC ), public, pass, deferred :: &
      SetPrimitiveConserved
    procedure, public, pass :: &
      SetReconstructed
    procedure, public, pass :: &
      ShowPrimitiveConserved
    procedure, public, pass :: &
      SetSources
    procedure, private, pass :: &
      ComputeFromPrimitiveSelf
    procedure, private, pass ( C ) :: &
      ComputeFromPrimitiveOther
    procedure, private, pass :: &
      ComputeFromPrimitiveSelectGeometry
    generic, public :: &
      ComputeFromPrimitive &
        => ComputeFromPrimitiveSelf, ComputeFromPrimitiveOther, &
           ComputeFromPrimitiveSelectGeometry
    procedure, private, pass :: &
      ComputeFromConservedSelf
    procedure, private, pass ( C ) :: &
      ComputeFromConservedOther
    !procedure, private, pass :: &
    !  ComputeFromConservedSelectGeometry
    generic, public :: &
      ComputeFromConserved &
        => ComputeFromConservedSelf, ComputeFromConservedOther !, &
           !ComputeFromConservedSelectGeometry
    procedure, public, pass ( C ) :: &
      ComputeFluxes
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( CFPC ), public, pass ( C ), deferred :: &
      ComputeFromPrimitiveCommon
    procedure ( CFCC ), public, pass ( C ), deferred :: &
      ComputeFromConservedCommon
    procedure, public, pass ( C ) :: &
      ComputeFluxes_HLL
    procedure ( CRF ), public, pass ( C ), deferred :: &
      ComputeRawFluxes
    procedure, public, pass ( C ) :: &
      ComputeDiffusionFactor_HLL
    procedure, public, nopass :: &
      SetDiffusionFactorUnity
  end type CurrentTemplate

  abstract interface

    subroutine SPC ( C )
      import CurrentTemplate
      class ( CurrentTemplate ), intent ( inout ) :: &
        C
    end subroutine SPC
    
    subroutine CFPC ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )
      use Basics
      use Manifolds
      import CurrentTemplate
      class ( StorageForm ), intent ( inout ), target :: &
        Storage_C
      class ( CurrentTemplate ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      class ( StorageForm ), intent ( in ) :: &
        Storage_G
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CFPC

    subroutine CFCC ( Storage_C, C, G, Storage_G, DetectFeaturesOption, &
                      nValuesOption, oValueOption )
      use Basics
      use Manifolds
      import CurrentTemplate
      class ( StorageForm ), intent ( inout ), target :: &
        Storage_C
      class ( CurrentTemplate ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      class ( StorageForm ), intent ( in ) :: &
        Storage_G
      logical ( KDL ), intent ( in ), optional :: &
        DetectFeaturesOption
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CFCC

    subroutine CRF ( RawFlux, C, G, Storage_C, Storage_G, iDimension, &
                     nValuesOption, oValueOption )
      use Basics
      use Manifolds
      import CurrentTemplate
      class ( StorageForm ), intent ( inout ) :: &
        RawFlux
      class ( CurrentTemplate ), intent ( in ) :: &
        C
      class ( GeometryFlatForm ), intent ( in ) :: &
        G
      class ( StorageForm ), intent ( in ) :: &
        Storage_C, &
        Storage_G
      integer ( KDI ), intent ( in ) :: &
        iDimension
      integer ( KDI ), intent ( in ), optional :: &
        nValuesOption, &
        oValueOption
    end subroutine CRF
    
  end interface

  type, public :: CurrentPointerForm
    class ( CurrentTemplate ), pointer :: &
      Pointer => null ( )
  end type CurrentPointerForm

    private :: &
      InitializeBasics, &
      SetUnits, &
      ComputeSolverSpeeds_HLL_Kernel, &
      ComputeFluxes_HLL_Kernel
      
  interface
  
    module subroutine SetDiffusionFactorUnity ( DFV_I, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        DFV_I
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine SetDiffusionFactorUnity
  
    module subroutine ComputeSolverSpeeds_HLL_Kernel &
                 ( AP_I, AM_I, LP_IL, LP_IR, LM_IL, LM_IR, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        AP_I, &
        AM_I
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        LP_IL, LP_IR, &
        LM_IL, LM_IR
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeSolverSpeeds_HLL_Kernel

    module subroutine ComputeFluxes_HLL_Kernel &
                 ( F_IL, F_IR, U_IL, U_IR, AP_I, AM_I, DF_I, F_I, &
                   UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        F_IL, F_IR, &
        U_IL, U_IR, &
        AP_I, &
        AM_I, &
        DF_I
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        F_I
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeFluxes_HLL_Kernel
  
  end interface

contains


  subroutine InitializeTemplate &
               ( C, RiemannSolverType, ReconstructedType, UseLimiter, &
                 Velocity_U_Unit, LimiterParameter, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, PinnedOption, &
                 UnitOption, VectorIndicesOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    character ( * ), intent ( in ) :: &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit
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
      ClearOption, &
      PinnedOption
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
    logical ( KDL ) :: &
      Clear

    call InitializeBasics &
           ( C, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits ( VariableUnit, C, Velocity_U_Unit )

    Clear = .true.
    if ( present ( ClearOption ) ) &
      Clear = ClearOption

    call C % StorageForm % Initialize &
           ( [ nValues, C % N_FIELDS ], &
             VariableOption = Variable, VectorOption = Vector, &
             NameOption = Name, ClearOption = Clear, &
             PinnedOption = PinnedOption, &
             UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

    C % RiemannSolverType = RiemannSolverType
    C % ReconstructedType = ReconstructedType
    C % UseLimiter        = UseLimiter
    C % LimiterParameter  = LimiterParameter
    call Show ( C % RiemannSolverType, 'RiemannSolverType', C % IGNORABILITY )
    call Show ( C % ReconstructedType, 'ReconstructedType', C % IGNORABILITY )
    call Show ( C % UseLimiter, 'UseLimiter', C % IGNORABILITY )
    if ( C % UseLimiter ) &
      call Show ( C % LimiterParameter, 'LimiterParameter', C % IGNORABILITY )
      
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeFromPrimitive', C % iTimerComputeFromPrimitive, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeFromConserved', C % iTimerComputeFromConserved, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeEigenspeed', C % iTimerComputeEigenspeed, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'Apply_EOS', C % iTimerApply_EOS, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeSolverSpeeds_HLL', C % iTimerComputeSolverSpeeds_HLL, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeFluxes_HLL', C % iTimerComputeFluxes_HLL, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeRawFluxes', C % iTimerComputeRawFluxes, &
             Level = 6 )
    call PROGRAM_HEADER % AddTimer &
           ( 'ComputeDiffusionFactor_HLL', &
             C % iTimerComputeDiffusionFactor_HLL, Level = 6 )

  end subroutine InitializeTemplate


  subroutine SetReconstructed ( C )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C

    select case ( trim ( C % ReconstructedType ) )
    case ( 'PRIMITIVE' )
      C % N_RECONSTRUCTED = C % N_PRIMITIVE
      allocate ( C % iaReconstructed, source = C % iaPrimitive )
    case ( 'CONSERVED' )
      C % N_RECONSTRUCTED = C % N_CONSERVED
      allocate ( C % iaReconstructed, source = C % iaConserved )
    case ( 'ALL' )
      C % N_RECONSTRUCTED = C % N_FIELDS
      allocate ( C % iaReconstructed, source = C % iaSelected )
    case default
      call Show ( 'ReconstructedType not recognized', CONSOLE % ERROR )
      call Show ( C % ReconstructedType, 'ReconstructedType', &
                  CONSOLE % ERROR )
      call Show ( 'SetReconstructed', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Current_Template', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

  end subroutine SetReconstructed


  subroutine ShowPrimitiveConserved ( C, IgnorabilityOption )

    class ( CurrentTemplate ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iF, &  !-- iField
      Ignorability
    character ( LDL ), dimension ( C % N_PRIMITIVE ) :: &
      PrimitiveName
    character ( LDL ), dimension ( C % N_CONSERVED ) :: &
      ConservedName
    character ( LDL ), dimension ( C % N_RECONSTRUCTED ) :: &
      ReconstructedName

    Ignorability = C % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      Ignorability = IgnorabilityOption

    do iF = 1, C % N_PRIMITIVE
      PrimitiveName ( iF )  =  C % Variable ( C % iaPrimitive ( iF ) )
    end do
    do iF = 1, C % N_CONSERVED
      ConservedName ( iF )  =  C % Variable ( C % iaConserved ( iF ) )
    end do
    do iF = 1, C % N_RECONSTRUCTED
      ReconstructedName ( iF )  =  C % Variable ( C % iaReconstructed ( iF ) )
    end do
    call Show ( PrimitiveName,     'Primitive variables',     Ignorability )
    call Show ( ConservedName,     'Conserved variables',     Ignorability )
    call Show ( ReconstructedName, 'Reconstructed variables', Ignorability )
    
  end subroutine ShowPrimitiveConserved


  subroutine SetSources ( C, Sources )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( Sources_C_Form ), intent ( in ), target :: &
      Sources

    call Show ( 'Setting Sources', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )
    call Show ( Sources % Name, 'Sources', C % IGNORABILITY )

    C % Sources => Sources

  end subroutine SetSources


  subroutine ComputeFromPrimitiveSelf ( C, G, nValuesOption, oValueOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromPrimitiveCommon &
           ( C, G, G, nValuesOption, oValueOption )
    
  end subroutine ComputeFromPrimitiveSelf


  subroutine ComputeFromPrimitiveOther &
               ( Storage_C, C, G, Storage_G, nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ) :: &
      Storage_C
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromPrimitiveCommon &
           ( Storage_C, G, Storage_G, nValuesOption, oValueOption )
    
  end subroutine ComputeFromPrimitiveOther


  subroutine ComputeFromPrimitiveSelectGeometry &
               ( C, iStrgeometryValue, G, nValuesOption, oValueOption )

    !-- Violating argument position for iStrgeometryValue to distinguish
    !   overloaded interface

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iStrgeometryValue
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption
      
    type ( StorageForm ) :: &
      Storage_G

    associate &
      ( nV => C % nValues, &
        iV => iStrgeometryValue, &
        iD => 1 )
    
    call Storage_G % Initialize ( [ nV, G % nVariables ] )
    Storage_G % Value = spread ( G % Value ( iV, : ), iD, nV )
    if ( C % AllocatedDevice ) then
      call Storage_G % AllocateDevice ( )
      call Storage_G % UpdateDevice ( )
    end if
    
    call C % ComputeFromPrimitiveCommon &
           ( C, G, Storage_G, nValuesOption, oValueOption )
    
    end associate !-- nV, etc.

  end subroutine ComputeFromPrimitiveSelectGeometry


  subroutine ComputeFromConservedSelf &
               ( C, G, DetectFeaturesOption, nValuesOption, oValueOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    logical ( KDL ), intent ( in ), optional :: &
      DetectFeaturesOption
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromConservedCommon &
           ( C, G, G, DetectFeaturesOption, nValuesOption, oValueOption )
    
  end subroutine ComputeFromConservedSelf


  subroutine ComputeFromConservedOther &
               ( Storage_C, C, G, Storage_G, DetectFeaturesOption, &
                 nValuesOption, oValueOption )

    class ( StorageForm ), intent ( inout ) :: &
      Storage_C
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    class ( StorageForm ), intent ( in ) :: &
      Storage_G
    logical ( KDL ), intent ( in ), optional :: &
      DetectFeaturesOption
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption

    call C % ComputeFromConservedCommon &
           ( Storage_C, G, Storage_G, DetectFeaturesOption, &
             nValuesOption, oValueOption )
    
  end subroutine ComputeFromConservedOther


  subroutine ComputeFromConservedSelectGeometry &
               ( C, iStrgeometryValue, G, DetectFeaturesOption, &
                 nValuesOption, oValueOption )

    !-- Violating argument position for iStrgeometryValue to distinguish
    !   overloaded interface

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iStrgeometryValue
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    logical ( KDL ), intent ( in ), optional :: &
      DetectFeaturesOption
    integer ( KDI ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption
      
    type ( StorageForm ) :: &
      Storage_G

    associate &
      ( nV => C % nValues, &
        iV => iStrgeometryValue, &
        iD => 1 )
        
    call Storage_G % Initialize ( [ nV, G % nVariables ] )
    Storage_G % Value = spread ( G % Value ( iV, : ), iD, nV )
    if ( C % AllocatedDevice ) then
      call Storage_G % AllocateDevice ( )
      call Storage_G % UpdateDevice ( )
    end if
    
    call C % ComputeFromConservedCommon &
           ( C, G, Storage_G, DetectFeaturesOption, nValuesOption, &
             oValueOption )
    
    end associate !-- nV, etc.

  end subroutine ComputeFromConservedSelectGeometry


  subroutine ComputeFluxes &
               ( Increment, F_I, F_IL, F_IR, SS_I, DF_I, C, Grid, G, &
                 C_IL, C_IR, G_I, iDimension )

    class ( * ), intent ( inout ) :: &
      Increment
    type ( StorageForm ), intent ( inout ) :: &
      F_I, &
      F_IL, F_IR, &
      SS_I, &
      DF_I
    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( * ), intent ( in ), target :: &
      Grid
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( StorageForm ), intent ( in ) :: &
      C_IL, C_IR, &
      G_I
    integer ( KDI ), intent ( in ) :: &
      iDimension

    call C % ComputeFluxes_HLL &
           ( Increment, F_I, F_IL, F_IR, SS_I, DF_I, Grid, G, C_IL, C_IR, &
             G_I, iDimension )

  end subroutine ComputeFluxes


  impure elemental subroutine FinalizeTemplate ( C )
    
    class ( CurrentTemplate ), intent ( inout ) :: &
      C

    nullify ( C % Sources )

    if ( allocated ( C % iaReconstructed ) ) &
      deallocate ( C % iaReconstructed )
    if ( allocated ( C % iaConserved ) ) &
      deallocate ( C % iaConserved )
    if ( allocated ( C % iaPrimitive ) ) &
      deallocate ( C % iaPrimitive )
 
    call Show ( 'Finalizing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )
   
  end subroutine FinalizeTemplate


  subroutine ComputeFluxes_HLL &
               ( Increment, F_I, F_IL, F_IR, SS_I, DF_I, C, Grid, G, &
                 C_IL, C_IR, G_I, iDimension )

    class ( * ), intent ( inout ) :: &
      Increment
    type ( StorageForm ), intent ( inout ) :: &
      F_I, &
      F_IL, F_IR, &
      SS_I, &
      DF_I
    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    class ( * ), intent ( in ), target :: &
      Grid
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( StorageForm ), intent ( in ) :: &
      C_IL, C_IR, &
      G_I
    integer ( KDI ), intent ( in ) :: &
      iDimension

    integer ( KDI ) :: &
      iF  !-- iField

    call C % ComputeRawFluxes &
           ( F_IL, G, C_IL, G_I, iDimension )
    call C % ComputeRawFluxes &
           ( F_IR, G, C_IR, G_I, iDimension )
           
    associate &
      ( T_SS => PROGRAM_HEADER % Timer ( C % iTimerComputeSolverSpeeds_HLL ), &
        T_F => PROGRAM_HEADER % Timer ( C % iTimerComputeFluxes_HLL ) )
      
    associate &
      ( SS_I_V => SS_I % Value, &
        C_IL_V => C_IL % Value, &
        C_IR_V => C_IR % Value )
    
    call T_SS % Start ( )
    
    call ComputeSolverSpeeds_HLL_Kernel &
           ( SS_I_V ( :, C % ALPHA_PLUS ), &
             SS_I_V ( :, C % ALPHA_MINUS ), &
             C_IL_V ( :, C % FAST_EIGENSPEED_PLUS ( iDimension ) ), &
             C_IR_V ( :, C % FAST_EIGENSPEED_PLUS ( iDimension ) ), &
             C_IL_V ( :, C % FAST_EIGENSPEED_MINUS ( iDimension ) ), &
             C_IR_V ( :, C % FAST_EIGENSPEED_MINUS ( iDimension ) ), &
             UseDeviceOption = C % AllocatedDevice )
    
    call T_SS % Stop ( )
    
    end associate   !-- SS_I_V
    
    call C % ComputeDiffusionFactor_HLL ( DF_I, Grid, iDimension )

    associate ( iaC => C % iaConserved )
    
    call T_F % Start ( )
    
    do iF = 1, C % N_CONSERVED
      call ComputeFluxes_HLL_Kernel &
             ( F_IL % Value ( :, iF ), &
               F_IR % Value ( :, iF ), &
               C_IL % Value ( :, iaC ( iF ) ), &
               C_IR % Value ( :, iaC ( iF ) ), &
               SS_I % Value ( :, C % ALPHA_PLUS ), &
               SS_I % Value ( :, C % ALPHA_MINUS ), &
               DF_I % Value ( :, iF ), &
               F_I % Value ( :, iF ), &
               UseDeviceOption = C % AllocatedDevice )
    end do !-- iF
    
    call T_F % Stop ( )
    
    end associate !-- iaC
    
    end associate !-- T_SS, T_F

  end subroutine ComputeFluxes_HLL


  subroutine ComputeDiffusionFactor_HLL ( DF_I, Grid, C, iDimension )

    type ( StorageForm ), intent ( inout ) :: &
      DF_I
    class ( * ), intent ( in ), target :: &
      Grid
    class ( CurrentTemplate ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iDimension
      
    integer ( KDI ) :: &
      iVariable
    
    associate &
      ( T_CDF => PROGRAM_HEADER % Timer &
                  ( C % iTimerComputeDiffusionFactor_HLL ) )
    
    call T_CDF % Start ( )
    
    do iVariable = 1, DF_I % nVariables
      call C % SetDiffusionFactorUnity &
             ( DF_I % Value ( :, iVariable ), &
               UseDeviceOption = DF_I % AllocatedDevice )
    end do
    
    call T_CDF % Stop ( )
    
    end associate !-- T_CDF

  end subroutine ComputeDiffusionFactor_HLL


  subroutine InitializeBasics &
               ( C, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( CurrentTemplate ), intent ( inout ) :: &
      C
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
      iV  !-- iVector

    if ( C % Type == '' ) &
      C % Type = 'a Current'

    Name = 'Current'
    if ( present ( NameOption ) ) &
      Name = NameOption

    C % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( Name, 'Name', C % IGNORABILITY )

    !-- variable indices

    if ( C % N_FIELDS == 0 ) &
      C % N_FIELDS = C % N_FIELDS_TEMPLATE

    C % FAST_EIGENSPEED_PLUS  = [ 1, 2, 3 ]
    C % FAST_EIGENSPEED_MINUS = [ 4, 5, 6 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( C % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : C % N_FIELDS_TEMPLATE ) &
      = [ 'FastEigenspeedPlus_1           ', &
          'FastEigenspeedPlus_2           ', &
          'FastEigenspeedPlus_3           ', &
          'FastEigenspeedMinus_1          ', &
          'FastEigenspeedMinus_2          ', &
          'FastEigenspeedMinus_3          ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( C % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( C % N_VECTORS == 0 ) &
      C % N_VECTORS = C % N_VECTORS_TEMPLATE

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( C % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( 1 : C % N_VECTORS_TEMPLATE ) &
      = [ 'FastEigenspeedPlus             ', &
          'FastEigenspeedMinus            ' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = C % N_VECTORS_TEMPLATE + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( C % N_VECTORS ) )
    end if

    call VectorIndices ( 1 ) % Initialize ( C % FAST_EIGENSPEED_PLUS )
    call VectorIndices ( 2 ) % Initialize ( C % FAST_EIGENSPEED_MINUS )

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, C, Velocity_U_Unit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit

    integer ( KDI ) :: &
      iD

    do iD = 1, 3
      VariableUnit ( C % FAST_EIGENSPEED_PLUS ( iD ) )  &
        = Velocity_U_Unit ( iD )
      VariableUnit ( C % FAST_EIGENSPEED_MINUS ( iD ) )  &
        = Velocity_U_Unit ( iD )
    end do

  end subroutine SetUnits


end module Current_Template
