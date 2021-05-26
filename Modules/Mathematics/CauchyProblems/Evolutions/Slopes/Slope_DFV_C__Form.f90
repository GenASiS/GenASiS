module Slope_DFV_C__Form

  !-- Slope_DivergenceFiniteVolume_Chart_Form

  use Basics
  use Manifolds
  use Fields
  use RiemannSolver_HLL_C__Form
  use Slope_H_C__Form

  implicit none
  private

  type, public, extends ( Slope_H_C_Form ) :: Slope_DFV_C_Form
    integer ( KDI ) :: &
      iTimer       = 0, &
      iTimerKernel = 0
    class ( CurrentSet_C_Form ), pointer :: &
      CurrentSet_C => null ( )
    class ( RiemannSolver_HLL_C_Form ), pointer :: &
      RiemannSolver_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_S
    generic, public :: &
      Initialize => InitializeAllocate_S
    procedure, private, pass :: &
      Show_FSC
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_C_Form

    private :: &
      ComputeKernel

    interface

      module subroutine ComputeKernel &
               ( S, F_I, A_I, V, iD, oV, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          S
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          F_I
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          A_I, &
          V
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_S ( SC, RSC, NameOption )

    class ( Slope_DFV_C_Form ), intent ( inout ) :: &
      SC
    class ( RiemannSolver_HLL_C_Form ), intent ( in ), target :: &
      RSC
    character ( * ), intent ( in ), optional :: &
      NameOption    

    character ( LDL ) :: &
      Name

    if ( SC % Type  ==  '' ) &
      SC % Type  =  'a Slope_DFV_C' 
    
    Name  =  'S_DFV_' // trim ( RSC % CurrentSet_C % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    SC % CurrentSet_C     =>  RSC % CurrentSet_C
    SC % RiemannSolver_C  =>  RSC

    associate &
      ( CSC  =>  RSC % CurrentSet_C )
    associate &
      ( nB  =>  CSC % nBalanced, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  CSC % GhostExchange_FSC % DevicesCommunicate ) 

    call SC % Slope_H_C_Form % Initialize &
           ( CSC % Chart, &
             FieldOption = CSC % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nB, &
             IgnorabilityOption = CSC % IGNORABILITY )

    end associate !-- nB, etc.
    end associate !-- CSC

  end subroutine InitializeAllocate_S


  subroutine Show_FSC ( FSC )

    class ( Slope_DFV_C_Form ), intent ( in ) :: &
      FSC

    call FSC % FieldSet_C_Form % Show ( )
    call FSC % RiemannSolver_C % Show ( )

  end subroutine Show_FSC


  subroutine Compute ( SC, TimerLevelOption, iS_Option )

    class ( Slope_DFV_C_Form ), intent ( inout ) :: &
      SC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption, &
      iS_Option

    integer ( KDI ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I, &
      V
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      S, &
      F_I
    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T, &
      T_Kernel

    associate ( iT  =>  SC % iTimer )
    if ( iT == 0 ) then
      TimerName  =  SC % Name
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( SC % iTimer )
    call T % Start ( )

    call Show ( 'Computing ' // trim ( SC % Type ), SC % IGNORABILITY + 4 )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY + 4 )

    associate &
      ( RSC  =>  SC % RiemannSolver_C, &
         GC  =>  SC % CurrentSet_C % Geometry_C, &
          C  =>  SC % Chart )

    do iD  =  1, C % nDimensions

      call RSC % Compute &
             ( iD, &
               TimerLevelOption = T % Level + 1, &
               iS_Option = iS_Option )

      associate ( iT_K  =>  SC % iTimerKernel )
      if ( iT_K == 0 ) then
        TimerName  =  trim ( T % Name ) // '_Kernel' 
        call PROGRAM_HEADER % AddTimer &
               ( TimerName, iT_K, Level = T % Level + 1 )
      end if
      end associate !-- iT_K

      T_Kernel  =>  PROGRAM_HEADER % TimerPointer ( SC % iTimerKernel )
      call T_Kernel % Start ( )

      associate &
        (  SV  =>   SC % Storage_FSC % Storage % Value, &
          RSV  =>  RSC % Storage_FSC % Storage % Value, &
           GV  =>   GC % Storage_FSC % Storage % Value, &
          DeviceMemory  =>  SC % Storage_FSC % DeviceMemory )

      select type ( C )
      class is ( Chart_GS_Form )

        call C % SetFieldPointer (  SV ( :, : ), S )
        call C % SetFieldPointer ( RSV ( :, : ), F_I )
        call C % SetFieldPointer (  GV ( :, GC % AREA_I_D ( iD ) ), A_I )
        call C % SetFieldPointer (  GV ( :, GC % VOLUME ), V )

        call ComputeKernel &
               ( S, F_I, A_I, V, iD, C % nGhostLayers ( iD ), &
                 UseDeviceOption = DeviceMemory )

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_DFV_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C

      end associate !-- SV, etc.

      call T_Kernel % Stop

    end do !-- iD

    end associate !-- RSC, etc.

    call T % Stop ( )

  end subroutine Compute


  impure elemental subroutine Finalize ( SC )

    type ( Slope_DFV_C_Form ), intent ( inout ) :: &
      SC

    nullify ( SC % RiemannSolver_C )

  end subroutine Finalize


end module Slope_DFV_C__Form
