module RiemannSolver_HLL_C__Form

  !-- RiemannSolver_HartenLaxVanLeer_Chart_Form

  use Basics
  use Fields
  use Reconstruction_C__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_SOLVER_SPEEDS_HLL  =  2

  type, public, extends ( FieldSet_C_Form ) :: RiemannSolver_HLL_C_Form
    integer ( KDI ) :: &
      iTimer       = 0, &
      iTimerKernel = 0
    integer ( KDI ) :: &
      N_SOLVER_SPEEDS_HLL = N_SOLVER_SPEEDS_HLL
    integer ( KDI ) :: &
      ALPHA_PLUS_U    = 0, &
      ALPHA_MINUS_U   = 0, &
      N_SOLVER_SPEEDS = 0
    class ( CurrentSet_C_Form ), pointer :: &
      CurrentSet_C => null ( )
    class ( FluxSet_C_Form ), pointer :: &
      FluxSet_C => null ( )
    class ( Eigenspeeds_F_C_Form ), pointer :: &
      Eigenspeeds_C => null ( )
    class ( Reconstruction_C_Form ), pointer :: &
      Reconstruction_B_C => null ( ), &
      Reconstruction_F_C => null ( ), &
      Reconstruction_E_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    generic, public :: &
      Initialize => InitializeAllocate_RS
    procedure, public, pass :: &
      Compute
    procedure, private, pass :: &
      Show_FSC
    final :: &
      Finalize
  end type RiemannSolver_HLL_C_Form

    private :: &
      ComputeKernel

    interface
      
      module subroutine ComputeKernel &
               ( F_IL, F_IR, U_IL, U_IR, EP_IL, EP_IR, EM_IL, EM_IR, &
                 F_I, AP_I, AM_I, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          F_IL, F_IR, &
          U_IL, U_IR
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          EP_IL, EP_IR, &
          EM_IL, EM_IR
        real ( KDR ), dimension ( :, : ), intent ( out ) :: &
          F_I
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          AP_I, AM_I
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface

contains


  subroutine InitializeAllocate_RS &
               ( RSC, RBC, RFC, REC, EC, FSC, CSC, FieldOption, NameOption, &
                 nFieldsOption )

    class ( RiemannSolver_HLL_C_Form ), intent ( inout ) :: &
      RSC
    class ( Reconstruction_C_Form ), intent ( in ), target :: &
      RFC, &
      REC, &
      RBC
    class ( Eigenspeeds_F_C_Form ), intent ( in ), target :: &
      EC
    class ( FluxSet_C_Form ), intent ( in ), target :: &
      FSC
    class ( CurrentSet_C_Form ), intent ( in ), target :: &
      CSC
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RSC % Type  ==  '' ) &
      RSC % Type  =  'a RiemannSolver_HLL_C' 
    
    Name  =  'RS_' // trim ( CSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    RSC % CurrentSet_C        =>  CSC
    RSC % FluxSet_C           =>  FSC
    RSC % Eigenspeeds_C       =>  EC
    RSC % Reconstruction_B_C  =>  RBC
    RSC % Reconstruction_F_C  =>  RFC
    RSC % Reconstruction_E_C  =>  REC

    associate &
      ( nB  =>  CSC % nBalanced, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  CSC % GhostExchange_FSC % DevicesCommunicate ) 

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      RSC % N_SOLVER_SPEEDS  =  RSC % N_SOLVER_SPEEDS_HLL
      nFields  =  nB  +  RSC % N_SOLVER_SPEEDS
    end if

    RSC % ALPHA_PLUS_U   =  nB  +  1
    RSC % ALPHA_MINUS_U  =  nB  +  2

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( : nB )  =  CSC % Balanced

    Field ( nB + 1 : nB + RSC % N_SOLVER_SPEEDS ) &
      =  [ 'AlphaPlus_U ', &
           'AlphaMinus_U' ]
          
    !-- FieldSet

    call RSC % FieldSet_C_Form % Initialize &
           ( CSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = CSC % IGNORABILITY )

    end associate !-- nB, etc.

  end subroutine InitializeAllocate_RS


  subroutine Compute ( RSC, iD, TimerLevelOption )

    class ( RiemannSolver_HLL_C_Form ), intent ( inout ) :: &
      RSC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T, &
      T_Kernel

    associate ( iT  =>  RSC % iTimer )
    if ( iT == 0 ) then
      TimerName  =  RSC % Name
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( RSC % iTimer )
    call T % Start ( )

    call Show ( 'Computing ' // trim ( RSC % Type ), RSC % IGNORABILITY + 4 )
    call Show ( RSC % Name, 'Name', RSC % IGNORABILITY + 4 )

    associate &
      ( CSC  =>  RSC % CurrentSet_C, &
        FSC  =>  RSC % FluxSet_C, &
         EC  =>  RSC % Eigenspeeds_C, &
        RBC  =>  RSC % Reconstruction_B_C, &
        RFC  =>  RSC % Reconstruction_F_C, &
        REC  =>  RSC % Reconstruction_E_C )

    call FSC % Compute ( iD, TimerLevelOption = T % Level + 1 )
    call  EC % Compute ( iD, TimerLevelOption = T % Level + 1 )
    call RBC % Compute ( iD, TimerLevelOption = T % Level + 1 )
    call RFC % Compute ( iD, TimerLevelOption = T % Level + 1 )
    call REC % Compute ( iD, TimerLevelOption = T % Level + 1 )

    associate ( iT_K  =>  RSC % iTimerKernel )
    if ( iT_K == 0 ) then
      TimerName  =  trim ( T % Name ) // '_Kernel' 
      call PROGRAM_HEADER % AddTimer ( TimerName, iT_K, Level = T % Level + 1 )
    end if
    end associate !-- iT_K

    T_Kernel  =>  PROGRAM_HEADER % TimerPointer ( RSC % iTimerKernel )
    call T_Kernel % Start ( )

    associate &
      ( RSS     =>  RSC % Storage_FSC % Storage, &
        RBS_IL  =>  RBC % Output_IL_C % Storage_FSC % Storage, &
        RBS_IR  =>  RBC % Output_IR_C % Storage_FSC % Storage, &
        RFS_IL  =>  RFC % Output_IL_C % Storage_FSC % Storage, &
        RFS_IR  =>  RFC % Output_IR_C % Storage_FSC % Storage, &
        RES_IL  =>  REC % Output_IL_C % Storage_FSC % Storage, &
        RES_IR  =>  REC % Output_IR_C % Storage_FSC % Storage, &
        DeviceMemory  =>  RSC % Storage_FSC % DeviceMemory )
    associate &
      (  F_I   =>  RSS % Value ( :, 1 : CSC % nBalanced ), &
        AP_I   =>  RSS % Value ( :, RSC % ALPHA_PLUS_U ), &
        AM_I   =>  RSS % Value ( :, RSC % ALPHA_MINUS_U ), &
         F_IL  =>  RFS_IL % Value ( :, : ), &
         F_IR  =>  RFS_IR % Value ( :, : ), &
         U_IL  =>  RBS_IL % Value ( :, : ), &
         U_IR  =>  RBS_IR % Value ( :, : ), &
        EP_IL  =>  RES_IL % Value ( :, EC % EIGENSPEED_FAST_PLUS_U ), &
        EP_IR  =>  RES_IR % Value ( :, EC % EIGENSPEED_FAST_PLUS_U ), &
        EM_IL  =>  RES_IL % Value ( :, EC % EIGENSPEED_FAST_MINUS_U ), &
        EM_IR  =>  RES_IR % Value ( :, EC % EIGENSPEED_FAST_MINUS_U ) )

    call ComputeKernel &
           ( F_IL, F_IR, U_IL, U_IR, EP_IL, EP_IR, EM_IL, EM_IR, &
             F_I, AP_I, AM_I, UseDeviceOption = DeviceMemory )

    end associate !-- F_I, etc.
    end associate !-- RSS, etc.

    call T_Kernel % Stop

    end associate !-- CSC, etc.

    call T % Stop ( )

  end subroutine Compute


  subroutine Show_FSC ( FSC )

    class ( RiemannSolver_HLL_C_Form ), intent ( in ) :: &
      FSC

    call FSC % FieldSet_C_Form % Show ( )
    call FSC % FluxSet_C % Show ( )
    call FSC % Eigenspeeds_C % Show ( )
    call FSC % Reconstruction_B_C % Show ( )
    call FSC % Reconstruction_F_C % Show ( )
    call FSC % Reconstruction_E_C % Show ( )

  end subroutine Show_FSC


  impure elemental subroutine Finalize ( RSC )

    type ( RiemannSolver_HLL_C_Form ), intent ( inout ) :: &
      RSC

    nullify ( RSC % Reconstruction_E_C )
    nullify ( RSC % Reconstruction_F_C )
    nullify ( RSC % Reconstruction_B_C )
    nullify ( RSC % Eigenspeeds_C )
    nullify ( RSC % FluxSet_C )
    nullify ( RSC % CurrentSet_C )

  end subroutine Finalize


end module RiemannSolver_HLL_C__Form
