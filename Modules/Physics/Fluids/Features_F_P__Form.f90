module Features_F_P__Form

  !-- Features_Fluid_Perfect__Form

  use Basics
  use Mathematics
  use Fluid_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_P  = 8, &
      N_VECTORS_P = 0

  type, public, extends ( Features_CS_Form ) :: Features_F_P_Form
    integer ( KDI ) :: &
      N_FIELDS_P  = N_FIELDS_P, &
      N_VECTORS_P = N_VECTORS_P, &
      EOS_ERROR        = 0, &
      SHOCK            = 0, &
      PHASE_TRANSITION = 0, &
      JAGGED_ENTROPY   = 0, &
      HEAVY_NUCLEUS    = 0
    integer ( KDI ), dimension ( 3 ) :: &
      SHOCK_I = 0
    real ( KDR ) :: &
      ShockThreshold
    type ( FieldSet_BM_Form ), allocatable :: &
      FeaturesExchange
  contains
    procedure, public, pass :: &
      InitializeAllocate_P
    generic, public :: &
      Initialize => InitializeAllocate_P
    procedure, public, pass ( F ) :: &
      SetStream
    procedure, public, pass :: &
      Detect
    final :: &
      Finalize
  !   procedure, public, pass :: &
  !     SetOutput
  end type Features_F_P_Form

    private :: &
      DetectShocksKernel, &
      DetectPhaseTransitionKernel, &
      DetectJaggedEntropyKernel, &
      DetectHeavyNucleusKernel

  interface
  
    module subroutine DetectShocksKernel &
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
    end subroutine DetectShocksKernel

    module subroutine DetectPhaseTransitionKernel &
             ( PT, DF_I_iD, DF_I_jD, DF_I_kD, Gamma, PTT, iD, jD, kD, oV, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        PT, &
        DF_I_iD, &
        DF_I_jD, &
        DF_I_kD
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        Gamma
      real ( KDR ), intent ( in ) :: &
        PTT
      integer ( KDI ), intent ( in ) :: &
        iD, jD, kD, &
        oV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine DetectPhaseTransitionKernel

    module subroutine DetectJaggedEntropyKernel &
             ( JE, DF_I_iD, DF_I_jD, DF_I_kD, SB, JET, iD, jD, kD, oV, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        JE, &
        DF_I_iD, &
        DF_I_jD, &
        DF_I_kD
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        SB
      real ( KDR ), intent ( in ) :: &
        JET
      integer ( KDI ), intent ( in ) :: &
        iD, jD, kD, &
        oV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine DetectJaggedEntropyKernel

    module subroutine DetectHeavyNucleusKernel &
             ( HN, DF_I_iD, DF_I_jD, DF_I_kD, XA, HNT, iD, jD, kD, oV, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        HN, &
        DF_I_iD, &
        DF_I_jD, &
        DF_I_kD
      real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
        XA
      real ( KDR ), intent ( in ) :: &
        HNT
      integer ( KDI ), intent ( in ) :: &
        iD, jD, kD, &
        oV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine DetectHeavyNucleusKernel

    module subroutine ClearBoundaryKernel &
               ( S, PT, S_I_iD, DF_I_iD, DF_I_jD, DF_I_kD, &
                 InnerBoundary, OuterBoundary, iD, jD, kD, oV, &
                 UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
        S, &
        PT, &
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
    end subroutine ClearBoundaryKernel

  end interface

contains


  subroutine InitializeAllocate_P &
               ( F, FP, ShockThreshold, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( Features_F_P_Form ), intent ( inout ) :: &
      F
    class ( Fluid_P_Form ), intent ( in ), target :: &
      FP
    real ( KDR ), intent ( in ) :: &
      ShockThreshold
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iV, &  !-- iVector
      oF, &  !-- oField
      oV, &  !-- oVector
      nFields
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( F % Type  ==  '' ) &
      F % Type  =  'a Features_F_P' 

    !-- Field indices

    oF = F % N_FIELDS_CS

    F % EOS_ERROR        =  oF + 1
    F % SHOCK            =  oF + 2
    F % PHASE_TRANSITION =  oF + 3
    F % JAGGED_ENTROPY   =  oF + 4
    F % HEAVY_NUCLEUS    =  oF + 5
    F % SHOCK_I          =  oF + [ 6, 7, 8 ]

    nFields  =  oF  +  F % N_FIELDS_P
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + F % N_FIELDS_P ) &
      = [ 'EOS_Error      ', &
          'Shock          ', &
          'PhaseTransition', &
          'JaggedEntropy  ', &
          'HeavyNucleus   ', &
          'Shock_I_1      ', &
          'Shock_I_2      ', &
          'Shock_I_3      ' ]

    !-- Units: none

    !-- Vector indices: no additional vectors

    !-- Vector names: no additional vectors

    !-- Features_CS

    call F % Features_CS_Form % Initialize &
           ( FP % Geometry, &
             FP, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )
     
    !-- FeaturesExchange

    allocate ( F % FeaturesExchange )
    call F % FeaturesExchange % Initialize ( F, [ F % SHOCK ] )

    !-- Parameters

    F % ShockThreshold  =  ShockThreshold
    call Show ( F % ShockThreshold, 'ShockThreshold', F % IGNORABILITY )
     
  end subroutine InitializeAllocate_P


  subroutine Detect ( F )

    class ( Features_F_P_Form ), intent ( inout ) :: &
      F

    integer ( KDI ) :: &
      iC, &
      iD, jD, kD
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      DF_I_iD, &
      DF_I_jD, &
      DF_I_kD, &
      S, &
      S_I_iD, &
      PT, &
      JE, &
      HN, &
      P, &
      Gamma, &
      SB, &
      XA, &
      V_iD
integer ( KDI ) :: &
  iHeavy

    call Show ( 'Detecting Fluid features', CONSOLE % INFO_6 )
    call Show ( F % Name, 'Features', CONSOLE % INFO_6 )
    
    call F % Clear ( )

    select type ( FP  =>  F % CurrentSet )
    class is ( Fluid_P_Form )

    do iC  =  1, F % Atlas % nCharts

      select type ( C  =>  F % Atlas % Chart ( iC ) % Element )
      class is ( Chart_GS_Form )

      associate &
        ( FV   =>  F  % Storage ( iC ) % Value, &
          FPV  =>  FP % Storage ( iC ) % Value )

iHeavy = 22
      call C % SetFieldPointer ( FV  ( :, F  % SHOCK ),              S )
      call C % SetFieldPointer ( FV  ( :, F  % PHASE_TRANSITION ),   PT )
      call C % SetFieldPointer ( FV  ( :, F  % JAGGED_ENTROPY ),     JE )
      call C % SetFieldPointer ( FV  ( :, F  % HEAVY_NUCLEUS ),      HN )
      call C % SetFieldPointer ( FPV ( :, FP % PRESSURE ),           P )
      call C % SetFieldPointer ( FPV ( :, FP % ADIABATIC_INDEX ),    Gamma )
      call C % SetFieldPointer ( FPV ( :, FP % ENTROPY_PER_BARYON ), SB )
      call C % SetFieldPointer ( FPV ( :, FP % ENTROPY_PER_BARYON ), SB )
      call C % SetFieldPointer ( FPV ( :, iHeavy ), XA )

      do iD = 1, C % nDimensions

        jD = mod ( iD, 3 ) + 1
        kD = mod ( jD, 3 ) + 1

        call C % SetFieldPointer &
               ( FV ( :, F % DIFFUSIVE_FLUX_I ( iD ) ), DF_I_iD )
        call C % SetFieldPointer &
               ( FV ( :, F % DIFFUSIVE_FLUX_I ( jD ) ), DF_I_jD )
        call C % SetFieldPointer &
               ( FV ( :, F % DIFFUSIVE_FLUX_I ( kD ) ), DF_I_kD )
        call C % SetFieldPointer &
               ( FV ( :, F % SHOCK_I ( iD ) ), S_I_iD )

        call C % SetFieldPointer &
               ( FPV ( :, FP % VELOCITY_U ( iD ) ), V_iD )

        call DetectShocksKernel &
               ( S, S_I_iD, DF_I_jD, DF_I_kD, P, V_iD, &
                 F % ShockThreshold, iD, jD, kD, C % nGhostLayers ( iD ), &
                 UseDeviceOption = F % DeviceMemory )
        call DetectPhaseTransitionKernel &
               ( PT, DF_I_iD, DF_I_jD, DF_I_kD, Gamma, &
                 0.01_KDR, iD, jD, kD, C % nGhostLayers ( iD ), &
                 UseDeviceOption = F % DeviceMemory )
        call DetectJaggedEntropyKernel &
               ( JE, DF_I_iD, DF_I_jD, DF_I_kD, SB, &
                 0.01_KDR, iD, jD, kD, C % nGhostLayers ( iD ), &
                 UseDeviceOption = F % DeviceMemory )
        call DetectHeavyNucleusKernel &
               ( HN, DF_I_iD, DF_I_jD, DF_I_kD, XA, &
                 0.01_KDR, iD, jD, kD, C % nGhostLayers ( iD ), &
                 UseDeviceOption = F % DeviceMemory )

      end do !-- iD

      !-- Separate dimension loop needed to avoid transverse effects
      !   (i.e. setting previously cleared cells) 
      do iD = 1, C % nDimensions

        jD = mod ( iD, 3 ) + 1
        kD = mod ( jD, 3 ) + 1

        call C % SetFieldPointer &
               ( FV ( :, F % DIFFUSIVE_FLUX_I ( iD ) ), DF_I_iD )
        call C % SetFieldPointer &
               ( FV ( :, F % DIFFUSIVE_FLUX_I ( jD ) ), DF_I_jD )
        call C % SetFieldPointer &
               ( FV ( :, F % DIFFUSIVE_FLUX_I ( kD ) ), DF_I_kD )
        call C % SetFieldPointer &
               ( FV ( :, F % SHOCK_I ( iD ) ), S_I_iD )

        call ClearBoundaryKernel &
               ( S, PT, S_I_iD, DF_I_iD, DF_I_jD, DF_I_kD, &
                 InnerBoundary = ( C % iaBrick ( iD ) == 1 ), &
                 OuterBoundary = ( C % iaBrick ( iD ) &
                                     == C % nBricks ( iD ) ), &
                 iD = iD, jD = jD, kD = kD, oV = C % nGhostLayers ( iD ), &
                 UseDeviceOption = F % DeviceMemory )

      end do !-- iD

      end associate !-- FV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Features_F_P__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Detect', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    
    end do !-- iC

    end select !-- FP

!    call F % FeaturesExchange % ExchangeGhostData ( )
    call F % ExchangeGhostData ( )

    nullify ( DF_I_jD, DF_I_kD, S, S_I_iD, P, V_iD )

  end subroutine Detect


  subroutine SetStream ( S, F )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Features_F_P_Form ), intent ( in ) :: &
      F

    call S % AddFieldSet &
           ( F, iaSelectedOption &
                  = [ F % DIFFUSIVE_FLUX_I, &
                      F % SHOCK, &
                      F % PHASE_TRANSITION, &
                      F % JAGGED_ENTROPY, &
                      F % HEAVY_NUCLEUS ] )

  end subroutine SetStream


  impure elemental subroutine Finalize ( F )

    type ( Features_F_P_Form ), intent ( inout ) :: &
      F

    if ( allocated ( F % FeaturesExchange ) ) &
      deallocate ( F % FeaturesExchange )

  end subroutine Finalize


end module Features_F_P__Form
