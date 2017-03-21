module FluidFeatures_P__Form

  !-- FluidFeatures_Perfect__Form

  use Basics
  use Mathematics
  use FluidFeatures_Template
  use Fluid_P__Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_PERFECT  = 8, &
      N_VECTORS_PERFECT = 0

  type, public, extends ( FluidFeaturesTemplate ) :: FluidFeatures_P_Form
    integer ( KDI ) :: &
      N_FIELDS_PERFECT  = N_FIELDS_PERFECT, &
      N_VECTORS_PERFECT = N_VECTORS_PERFECT, &
      SHOCK             = 0, &
      PHASE_TRANSITION  = 0
    integer ( KDI ), dimension ( 3 ) :: &
      SHOCK_I            = 0, &
      PHASE_TRANSITION_I = 0
    real ( KDR ) :: &
      ShockThreshold, &
      PhaseTransitionThreshold
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
      DetectShocks_CSL, &
      DetectPhaseTransition_CSL

      private :: &
        DetectShocks_CSL_Kernel, &
        DetectPhaseTransition_CSL_Kernel

contains


  subroutine InitializeAllocate_P &
               ( FF, Fluid, Grid, ShockThreshold, PhaseTransitionThreshold, &
                 nValues, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    class ( VariableGroupForm ), intent ( in ) :: &
      Fluid
    class ( * ), intent ( in ) :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      ShockThreshold, &
      PhaseTransitionThreshold
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

    FF % ShockThreshold           = ShockThreshold
    FF % PhaseTransitionThreshold = PhaseTransitionThreshold

  end subroutine InitializeAllocate_P


  subroutine Detect ( FF )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF

    call Clear ( FF % Value )

    select type ( F => FF % Fluid )
    class is ( Fluid_P_Template )

    select type ( Grid => FF % Grid )
    class is ( Chart_SL_Template )

      call DetectShocks_CSL ( FF, F, Grid )
      call DetectPhaseTransition_CSL ( FF, F, Grid )

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
    type ( VariableGroupForm ), intent ( inout ) :: &
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

    FF % SHOCK               =  oF + 1
    FF % PHASE_TRANSITION    =  oF + 2
    FF % SHOCK_I             =  oF + [ 3, 4, 5 ]
    FF % PHASE_TRANSITION_I  =  oF + [ 6, 7, 8 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( FF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + FF % N_FIELDS_PERFECT ) &
      = [ 'Shock              ', &
          'PhaseTransition    ', &
          'Shock_I_1          ', &
          'Shock_I_2          ', &
          'Shock_I_3          ', &
          'PhaseTransition_I_1', &
          'PhaseTransition_I_2', &
          'PhaseTransition_I_3' ]
          
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
             ( S, S_I_iD, DF_I_jD, DF_I_kD, CSL, P, V_iD, &
               FF % ShockThreshold, iD, jD, kD, &
               CSL % nGhostLayers ( iD ) )

    end do !-- iD

    nullify ( S_I_iD, P, V_iD )

  end subroutine DetectShocks_CSL


  subroutine DetectPhaseTransition_CSL ( FF, F, CSL )

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
      PT, &
      PT_I_iD, &
      Gamma

    if ( FF % PhaseTransitionThreshold == 0.0_KDR ) &
      return

    call CSL % SetVariablePointer &
           ( FF % Value ( :, FF % PHASE_TRANSITION ), PT )
    call Clear ( PT )
    
    call CSL % SetVariablePointer &
           ( F % Value ( :, F % ADIABATIC_INDEX ), Gamma )

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
             ( FF % Value ( :, FF % PHASE_TRANSITION_I ( iD ) ), PT_I_iD )

      call DetectPhaseTransition_CSL_Kernel &
             ( PT, PT_I_iD, DF_I_iD, DF_I_jD, DF_I_kD, CSL, Gamma, &
               FF % PhaseTransitionThreshold, iD, jD, kD, &
               CSL % nGhostLayers ( iD ) )

    end do !-- iD

  end subroutine DetectPhaseTransition_CSL


  subroutine DetectShocks_CSL_Kernel &
               ( S, S_I_iD, DF_I_jD, DF_I_kD, CSL, P, V_iD, ST, &
                 iD, jD, kD, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      S, &
      S_I_iD, &
      DF_I_jD, &
      DF_I_kD
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      P, &
      V_iD
    real ( KDR ), intent ( in ) :: &
      ST
    integer ( KDI ), intent ( in ) :: &
      iD, jD, kD, &
      oV

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS_i, iaS_j, iaS_k, &
      iaV_i, iaV_j, iaV_ij, iaV_k, iaV_ik, &
      lV, uV
    real ( KDR ) :: &
      dP, &
      P_Min, &
      dLnP, &
      dV_iD, &
      SqrtTiny

    lV = 1
    where ( shape ( S ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( S ) > 1 )
      uV = shape ( S ) - oV
    end where
    uV ( iD ) = size ( S, dim = iD ) - oV + 1 
      
    iaS_i = 0
    iaS_i ( iD ) = -1
      
    iaS_j = 0
    if ( size ( S, dim = jD ) > 1 ) &
      iaS_j ( jD ) = +1

    iaS_k = 0
    if ( size ( S, dim = kD ) > 1 ) &
      iaS_k ( kD ) = +1

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    !$OMP parallel do &
    !$OMP   private ( iV, jV, kV, iaV_i, iaV_j, iaV_ij, iaV_k, iaV_ik )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaV_i  = [ iV, jV, kV ] + iaS_i
          iaV_j  = [ iV, jV, kV ] + iaS_j
          iaV_ij = [ iV, jV, kV ] + iaS_i + iaS_j
          iaV_k  = [ iV, jV, kV ] + iaS_k
          iaV_ik = [ iV, jV, kV ] + iaS_i + iaS_k

          dP  =  abs ( P ( iV, jV, kV )  &
                       -  P ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) )
          P_Min  =  max ( min ( P ( iV, jV, kV ), &
                                P ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )), &
                          SqrtTiny )
          dLnP  =  dP / P_Min

          dV_iD  =  V_iD ( iV, jV, kV ) &
                    -  V_iD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )
 
          if ( dLnP > ST .and. dV_iD <= 0.0_KDR ) then

            S_I_iD ( iV, jV, kV )  &
              =  1.0_KDR
            S ( iV, jV, kV )  &
              =  S ( iV, jV, kV )  +  1.0_KDR
            S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  &
              =  S ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  +  1.0_KDR

            !-- Use diffuse flux in transverse directions, on both sides of 
            !   the shock

            DF_I_jD ( iV, jV, kV ) &
              =  1.0_KDR
            DF_I_jD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
              =  1.0_KDR
            DF_I_jD ( iaV_j ( 1 ), iaV_j ( 2 ), iaV_j ( 3 ) ) &
              =  1.0_KDR
            DF_I_jD ( iaV_ij ( 1 ), iaV_ij ( 2 ), iaV_ij ( 3 ) ) &
              =  1.0_KDR

            DF_I_kD ( iV, jV, kV ) &
              =  1.0_KDR
            DF_I_kD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
              =  1.0_KDR
            DF_I_kD ( iaV_k ( 1 ), iaV_k ( 2 ), iaV_k ( 3 ) ) &
              =  1.0_KDR
            DF_I_kD ( iaV_ik ( 1 ), iaV_ik ( 2 ), iaV_ik ( 3 ) ) &
              =  1.0_KDR

          end if

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
  end subroutine DetectShocks_CSL_Kernel


  subroutine DetectPhaseTransition_CSL_Kernel &
               ( PT, PT_I_iD, DF_I_iD, DF_I_jD, DF_I_kD, CSL, Gamma, PTT, &
                 iD, jD, kD, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      PT, &
      PT_I_iD, &
      DF_I_iD, &
      DF_I_jD, &
      DF_I_kD
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      Gamma
    real ( KDR ), intent ( in ) :: &
      PTT
    integer ( KDI ), intent ( in ) :: &
      iD, jD, kD, &
      oV

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS_i, iaS_j, iaS_k, &
      iaV_i, iaV_j, iaV_ij, iaV_k, iaV_ik, &
      lV, uV
    real ( KDR ) :: &
      dGamma, &
      Gamma_Min, &
      dLnGamma, &
      SqrtTiny

    lV = 1
    where ( shape ( PT ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( PT ) > 1 )
      uV = shape ( PT ) - oV
    end where
    uV ( iD ) = size ( PT, dim = iD ) - oV + 1 

    iaS_i = 0
    iaS_i ( iD ) = -1
      
    iaS_j = 0
    if ( size ( PT, dim = jD ) > 1 ) &
      iaS_j ( jD ) = +1

    iaS_k = 0
    if ( size ( PT, dim = kD ) > 1 ) &
      iaS_k ( kD ) = +1

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    !$OMP parallel do private ( iV, jV, kV )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaV_i  = [ iV, jV, kV ] + iaS_i
          iaV_j  = [ iV, jV, kV ] + iaS_j
          iaV_ij = [ iV, jV, kV ] + iaS_i + iaS_j
          iaV_k  = [ iV, jV, kV ] + iaS_k
          iaV_ik = [ iV, jV, kV ] + iaS_i + iaS_k

          dGamma  =  abs ( Gamma ( iV, jV, kV )  &
                           -  Gamma ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) )
          Gamma_Min  =  max ( min ( Gamma ( iV, jV, kV ), &
                                    Gamma ( iaV_i ( 1 ), iaV_i ( 2 ), &
                                            iaV_i ( 3 ) ) ), &
                              SqrtTiny )
          dLnGamma  =  dGamma / Gamma_Min

          if ( dLnGamma > PTT ) then

            PT_I_iD ( iV, jV, kV )  &
              =  1.0_KDR
            PT ( iV, jV, kV )  &
              =  PT ( iV, jV, kV )  +  1.0_KDR
            PT ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  &
              =  PT ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) )  +  1.0_KDR

            !-- Use diffuse flux all directions, on both sides of 
            !   the interface shock

            DF_I_iD ( iV, jV, kV ) &
              =  1.0_KDR
 
            DF_I_jD ( iV, jV, kV ) &
              =  1.0_KDR
            DF_I_jD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
              =  1.0_KDR
            DF_I_jD ( iaV_j ( 1 ), iaV_j ( 2 ), iaV_j ( 3 ) ) &
              =  1.0_KDR
            DF_I_jD ( iaV_ij ( 1 ), iaV_ij ( 2 ), iaV_ij ( 3 ) ) &
              =  1.0_KDR

            DF_I_kD ( iV, jV, kV ) &
              =  1.0_KDR
            DF_I_kD ( iaV_i ( 1 ), iaV_i ( 2 ), iaV_i ( 3 ) ) &
              =  1.0_KDR
            DF_I_kD ( iaV_k ( 1 ), iaV_k ( 2 ), iaV_k ( 3 ) ) &
              =  1.0_KDR
            DF_I_kD ( iaV_ik ( 1 ), iaV_ik ( 2 ), iaV_ik ( 3 ) ) &
              =  1.0_KDR

          end if

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
  end subroutine DetectPhaseTransition_CSL_Kernel


end module FluidFeatures_P__Form
