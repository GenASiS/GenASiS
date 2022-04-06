module RiemannSolver_HLLC_P__Form

  !-- RiemannSolver_HartenLaxVanLeerContact_Perfect__Form

  use Basics
  use Mathematics
  use Fluid_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_SOLVER_SPEEDS_HLLC  =  1

  type, public, extends ( RiemannSolver_HLL_Form ) :: RiemannSolver_HLLC_P_Form
    integer ( KDI ) :: &
      N_SOLVER_SPEEDS_HLLC = N_SOLVER_SPEEDS_HLLC
    integer ( KDI ) :: &
      ALPHA_CENTER_U = 0
    type ( FieldSetForm ), allocatable :: &
      Metric_I, &
      CurrentSet_ICL, CurrentSet_ICR
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type RiemannSolver_HLLC_P_Form

    private :: &
      ComputeCenterSpeedKernel

    interface 

      module subroutine ComputeCenterSpeedKernel &
                   ( AC_I, F_D_IL, F_D_IR, F_S_IL, F_S_IR, M_IL, M_IR, &
                     D_IL, D_IR, S_IL, S_IR, AP_I, AM_I, M_UU, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          AC_I
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          F_D_IL, F_D_IR, &
          F_S_IL, F_S_IR, &
          M_IL, M_IR, &
          D_IL, D_IR, &
          S_IL, S_IR, &
          AP_I, &
          AM_I, &
          M_UU
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeCenterSpeedKernel

    end interface


contains


  subroutine InitializeAllocate_RS &
               ( RS, CS, FieldOption, ReconstructedSetOption, PrefixOption, &
                 nFieldsOption )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      ReconstructedSetOption, &
      PrefixOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      oF, &  !-- oField
      nFields
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RS % Type  ==  '' ) &
      RS % Type  =  'a RiemannSolver_HLLC_P' 
    
    associate ( nB  =>  CS % nBalanced )

    !-- Field indices

    oF  =  nB  +  RS % N_SOLVER_SPEEDS_HLL

    nFields  =  oF  +  RS % N_SOLVER_SPEEDS_HLLC
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    RS % ALPHA_CENTER_U   =  oF  +  1

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RS % N_SOLVER_SPEEDS_HLLC ) &
      =  [ 'AlphaCenter_U' ]
          
    !-- FieldSet

    call RS % RiemannSolver_HLL_Form % Initialize &
           ( CS, &
             FieldOption = Field, &
             ReconstructedSetOption = ReconstructedSetOption, &
             PrefixOption = PrefixOption, &
             nFieldsOption = nFields )

    end associate !-- nB

    !-- Reconstructed metric

    allocate ( RS % Metric_I )
    associate ( M_I  =>  RS % Metric_I )
    call M_I % Initialize &
           ( CS % Atlas, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = 6 )
    end associate !-- M_I

    !-- Center states

    allocate &
      ( RS % CurrentSet_ICL, &
        RS % CurrentSet_ICR )
    associate &
      ( CS_ICL  =>  RS % CurrentSet_ICL, &
        CS_ICR  =>  RS % CurrentSet_ICR )
    call CS_ICL % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Field, &
             NameOption = trim ( CS % Name ) // '_ICL', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             UnitOption = CS % Unit, &
             nFieldsOption = CS % nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    call CS_ICR % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Field, &
             NameOption = trim ( CS % Name ) // '_ICR', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             UnitOption = CS % Unit, &
             nFieldsOption = CS % nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- CS_ICL, etc.

  end subroutine InitializeAllocate_RS


  subroutine Compute ( RS, DP, iC, iD, T_Option )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ), target :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iDensity, &
      iMomentum
    real ( KDR ), dimension ( : ), pointer :: &
      M_UU

    call RS % RiemannSolver_HLL_Form % Compute ( DP, iC, iD, T_Option )

    select type ( CS  =>  RS % CurrentSet )
      class is ( Fluid_P_Form )
    associate &
      ( G    =>  CS % Geometry, &
        M_I  =>  RS % Metric_I)
    
    call G % ComputeReconstruction ( M_I, iC, iD )

    associate &
      ( RSV      =>  RS % Storage ( iC ) % Value, &
        FS_IL_V  =>  RS % FluxSet_IL % Storage ( iC ) % Value, &
        FS_IR_V  =>  RS % FluxSet_IR % Storage ( iC ) % Value, &
        CS_IL_V  =>  RS % CurrentSet_IL % Storage ( iC ) % Value, &
        CS_IR_V  =>  RS % CurrentSet_IR % Storage ( iC ) % Value )

    associate &
      ( M_DD_11  =>  M_I % Storage ( iC ) % Value ( :, 1 ), &
        M_DD_22  =>  M_I % Storage ( iC ) % Value ( :, 2 ), &
        M_DD_33  =>  M_I % Storage ( iC ) % Value ( :, 3 ), &
        M_UU_11  =>  M_I % Storage ( iC ) % Value ( :, 4 ), &
        M_UU_22  =>  M_I % Storage ( iC ) % Value ( :, 5 ), &
        M_UU_33  =>  M_I % Storage ( iC ) % Value ( :, 6 ) )

    select case ( iD )
    case ( 1 )
      M_UU  =>  M_UU_11
    case ( 2 ) 
      M_UU  =>  M_UU_22
    case ( 3 )
      M_UU  =>  M_UU_33
    end select

    call Search ( CS % iaBalanced, CS % BARYON_DENSITY_C,          iDensity )
    call Search ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D ( iD ), iMomentum )

    call ComputeCenterSpeedKernel &
           (  AC_I  = RSV     ( :, RS % ALPHA_CENTER_U ), &
             F_D_IL = FS_IL_V ( :, iDensity ), &
             F_D_IR = FS_IR_V ( :, iDensity ), &
             F_S_IL = FS_IL_V ( :, iMomentum ), &
             F_S_IR = FS_IR_V ( :, iMomentum ), &
               M_IL = CS_IL_V ( :, CS % BARYON_MASS ), &
               M_IR = CS_IR_V ( :, CS % BARYON_MASS ), &
               D_IL = CS_IL_V ( :, CS % BARYON_DENSITY_C ), &
               D_IR = CS_IR_V ( :, CS % BARYON_DENSITY_C ), &
               S_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
               S_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
              AP_I  = RSV     ( :, RS % ALPHA_PLUS_U ), &
              AM_I  = RSV     ( :, RS % ALPHA_MINUS_U ), &
               M_UU = M_UU, &
             UseDeviceOption = RS % DeviceMemory )

    end associate !-- M_DD_11, etc.
    end associate !-- RSV, etc.
    end associate !-- G, etc.
    end select !-- CS

  end subroutine Compute


  impure elemental subroutine Finalize ( RS )

    type ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS

    if ( allocated ( RS % CurrentSet_ICR ) ) &
      deallocate ( RS % CurrentSet_ICR )
    if ( allocated ( RS % CurrentSet_ICL ) ) &
      deallocate ( RS % CurrentSet_ICL )
    if ( allocated ( RS % Metric_I ) ) &
      deallocate ( RS % Metric_I )

  end subroutine Finalize


end module RiemannSolver_HLLC_P__Form
