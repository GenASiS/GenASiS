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
      ALPHA_CENTER = 0
    type ( FieldSetForm ), allocatable :: &
      Metric_I
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type RiemannSolver_HLLC_P_Form


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

    RS % ALPHA_CENTER   =  nB  +  1

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RS % N_SOLVER_SPEEDS_HLLC ) &
      =  [ 'AlphaCenter' ]
          
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

  end subroutine InitializeAllocate_RS


  subroutine Compute ( RS, DP, iC, iD, T_Option )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    call RS % RiemannSolver_HLL_Form % Compute ( DP, iC, iD, T_Option )

    associate &
      ( G    =>  RS % CurrentSet % Geometry, &
        M_I  =>  RS % Metric_I)
    
    call G % ComputeReconstruction ( M_I, iC, iD )

    end associate !-- G, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( RS )

    type ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS

    if ( allocated ( RS % Metric_I ) ) &
      deallocate ( RS % Metric_I )

  end subroutine Finalize


end module RiemannSolver_HLLC_P__Form
