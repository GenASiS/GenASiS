module FluidFeatures_P__Form

  !-- FluidFeatures_Perfect__Form

  use Basics
  use Mathematics
  use FluidFeatures_Template
  use Fluid_P__Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_PERFECT  = 5, &
      N_VECTORS_PERFECT = 0

  type, public, extends ( FluidFeaturesTemplate ) :: FluidFeatures_P_Form
    integer ( KDI ) :: &
      N_FIELDS_PERFECT  = N_FIELDS_PERFECT, &
      N_VECTORS_PERFECT = N_VECTORS_PERFECT, &
      EOS_ERROR = 0, &
      SHOCK     = 0
    integer ( KDI ), dimension ( 3 ) :: &
      SHOCK_I = 0
    real ( KDR ) :: &
      ShockThreshold
  contains
    procedure, public, pass :: &
      InitializeAllocate_P
    generic, public :: &
      Initialize => InitializeAllocate_P
    procedure, public, pass :: &
      Detect
    final :: &
      Finalize
    procedure, public, pass :: &
      SetOutput
  end type FluidFeatures_P_Form

    private :: &
      InitializeBasics, &
      DetectShocks_CSL

      private :: &
        DetectShocks_CSL_Kernel, &
        ClearBoundary_CSL_Kernel
  
  interface
  
    module subroutine DetectShocks_CSL_Kernel &
                 ( S, S_I_iD, DF_I_jD, DF_I_kD, P, V_iD, ST, &
                   iD, jD, kD, oV, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        S, &
        S_I_iD, &
        DF_I_jD, &
        DF_I_kD
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        P, &
        V_iD
      real ( KDR ), intent ( in ) :: &
        ST
      integer ( KDI ), intent ( in ) :: &
        iD, jD, kD, &
        oV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine DetectShocks_CSL_Kernel
    
    module subroutine ClearBoundary_CSL_Kernel &
               ( S, S_I_iD, DF_I_iD, DF_I_jD, DF_I_kD, &
                 InnerBoundary, OuterBoundary, iD, jD, kD, oV, &
                 UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        S, &
        S_I_iD, &
        DF_I_iD, &
        DF_I_jD, &
        DF_I_kD
      logical ( KDL ), intent ( in ) :: &
        InnerBoundary, &
        OuterBoundary
      integer ( KDI ), intent ( in ) :: &
        iD, jD, kD, &
        oV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ClearBoundary_CSL_Kernel

  end interface

contains


  subroutine InitializeAllocate_P &
               ( FF, Fluid, Grid, ShockThreshold, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    class ( StorageForm ), intent ( in ) :: &
      Fluid
    class ( * ), intent ( in ) :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      ShockThreshold
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

    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable

    call InitializeBasics &
           ( FF, Variable, VariableUnit, VariableOption, UnitOption )

    call FF % InitializeTemplate &
           ( Fluid, Grid, nValues, VariableOption = Variable, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndicesOption )

    FF % ShockThreshold = ShockThreshold
    call Show ( FF % ShockThreshold, 'ShockThreshold', FF % IGNORABILITY )
     
  end subroutine InitializeAllocate_P


  subroutine Detect ( FF )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF

    call Show ( 'Detecting Fluid features', CONSOLE % INFO_3 )
    
    call Clear ( FF % Value, UseDeviceOption = FF % AllocatedDevice )

    select type ( F => FF % Fluid )
    class is ( Fluid_P_Template )

    select type ( Grid => FF % Grid )
    class is ( Chart_SL_Template )

      call DetectShocks_CSL ( FF, F, Grid )

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_P__Template', 'module', CONSOLE % ERROR )
      call Show ( 'DetectFeaturesTemplate_P', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end select !-- F
    
  end subroutine Detect


  impure elemental subroutine Finalize ( FF )

    type ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF

    call FF % FinalizeTemplate ( )
   
  end subroutine Finalize


  subroutine SetOutput ( FF, Output )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    type ( StorageForm ), intent ( inout ) :: &
      Output

    call Output % Initialize &
           ( FF, iaSelectedOption = [ FF % SHOCK, FF % DIFFUSIVE_FLUX_I ] )

  end subroutine SetOutput


  subroutine InitializeBasics &
               ( FF, Variable, VariableUnit, VariableOption, &
                 VariableUnitOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      oF, &  !-- oField
      oV     !-- oVector

    if ( FF % Type == '' ) &
      FF % Type = 'a FluidFeatures_P'

    !-- variable indices

    oF = FF % N_FIELDS_TEMPLATE
    if ( FF % N_FIELDS == 0 ) &
      FF % N_FIELDS = oF + FF % N_FIELDS_PERFECT

    FF % EOS_ERROR  =  oF + 1
    FF % SHOCK      =  oF + 2
    FF % SHOCK_I    =  oF + [ 3, 4, 5 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( FF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + FF % N_FIELDS_PERFECT ) &
      = [ 'EOS_Error          ', &
          'Shock              ', &
          'Shock_I_1          ', &
          'Shock_I_2          ', &
          'Shock_I_3          ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( FF % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = FF % N_VECTORS_TEMPLATE
    if ( FF % N_VECTORS == 0 ) &
      FF % N_VECTORS = oV + FF % N_VECTORS_PERFECT

  end subroutine InitializeBasics


  subroutine DetectShocks_CSL ( FF, F, CSL )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    class ( Fluid_P_Template ), intent ( in ) :: &
      F
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL

    integer ( KDI ) :: &
      iD, jD, kD
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      DF_I_iD, &
      DF_I_jD, &
      DF_I_kD, &
      S, &
      S_I_iD, &
      P, &
      V_iD

    call CSL % SetVariablePointer &
           ( FF % Value ( :, FF % SHOCK ), S )

    call CSL % SetVariablePointer &
           ( F % Value ( :, F % PRESSURE ), P )

    do iD = 1, CSL % nDimensions

      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1

      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % DIFFUSIVE_FLUX_I ( jD ) ), DF_I_jD )
      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % DIFFUSIVE_FLUX_I ( kD ) ), DF_I_kD )
      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % SHOCK_I ( iD ) ), S_I_iD )

      call CSL % SetVariablePointer &
             ( F % Value ( :, F % VELOCITY_U ( iD ) ), V_iD )

      call DetectShocks_CSL_Kernel &
             ( S, S_I_iD, DF_I_jD, DF_I_kD, P, V_iD, &
               FF % ShockThreshold, iD, jD, kD, &
               CSL % nGhostLayers ( iD ), &
               UseDeviceOption = FF % AllocatedDevice )

    end do !-- iD

    !-- Separate dimension loop needed to avoid transverse effects
    !   (i.e. setting previously cleared cells) 
    do iD = 1, CSL % nDimensions

      jD = mod ( iD, 3 ) + 1
      kD = mod ( jD, 3 ) + 1

      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % DIFFUSIVE_FLUX_I ( iD ) ), DF_I_iD )
      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % DIFFUSIVE_FLUX_I ( jD ) ), DF_I_jD )
      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % DIFFUSIVE_FLUX_I ( kD ) ), DF_I_kD )
      call CSL % SetVariablePointer &
             ( FF % Value ( :, FF % SHOCK_I ( iD ) ), S_I_iD )

      call ClearBoundary_CSL_Kernel &
             ( S, S_I_iD, DF_I_iD, DF_I_jD, DF_I_kD, &
               InnerBoundary = ( CSL % iaBrick ( iD ) == 1 ), &
               OuterBoundary = ( CSL % iaBrick ( iD ) &
                                   == CSL % nBricks ( iD ) ), &
               iD = iD, jD = jD, kD = kD, oV = CSL % nGhostLayers ( iD ), &
               UseDeviceOption = FF % AllocatedDevice )

    end do !-- iD

    nullify ( DF_I_jD, DF_I_kD, S, S_I_iD, P, V_iD )

  end subroutine DetectShocks_CSL


end module FluidFeatures_P__Form
