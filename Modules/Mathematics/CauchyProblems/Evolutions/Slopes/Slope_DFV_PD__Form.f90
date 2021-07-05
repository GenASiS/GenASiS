module Slope_DFV_PD__Form

  !-- Slope_DivergenceFiniteVolume_PartialDerivative__Form

  use Basics
  use Manifolds
  use Fields
  use RiemannSolver_HLL__Form
  use Slope_H__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_PD_Form
    integer ( KDI ) :: &
      iTimerKernel = 0
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
    class ( RiemannSolver_HLL_Form ), pointer :: &
      RiemannSolver => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_PD
    generic, public :: &
      Initialize => InitializeAllocate_PD
    procedure, public, pass :: &
      Show => Show_FS
    procedure, private, pass :: &
      TimerKernel
    procedure, public, pass :: &
      CloneTimers
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_PD_Form

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


  subroutine InitializeAllocate_PD ( S, RS, SuffixOption )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    character ( * ), intent ( in ), optional :: &
      SuffixOption    

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_PD' 
    
    if ( S % TimerName  ==  '' ) &
      S % TimerName  =  'Slp_DFV_PD_' // trim ( RS % CurrentSet % Name )

    Name  =  'Slp_DFV_PD_' // trim ( RS % CurrentSet % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    S % CurrentSet     =>  RS % CurrentSet
    S % RiemannSolver  =>  RS

    associate ( CS  =>  RS % CurrentSet )

    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

    end associate !-- CS

  end subroutine InitializeAllocate_PD


  subroutine Show_FS ( FS )

    class ( Slope_DFV_PD_Form ), intent ( in ) :: &
      FS

    call FS % Slope_H_Form % Show ( )
    call FS % RiemannSolver % Show ( )

  end subroutine Show_FS


  function TimerKernel ( S, LevelOption ) result ( T )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  S % iTimerKernel )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % TimerName ) // '_Krnl' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function TimerKernel


  subroutine CloneTimers ( S, S_S )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    class ( Slope_H_Form ), intent ( in ) :: &
      S_S  !-- S_Source

    integer ( KDI ) :: &
      iC  !-- iComponent

    call S % Slope_H_Form % CloneTimers ( S_S )

    select type ( S_S )
    class is ( Slope_DFV_PD_Form )

    S % iTimerKernel  =  S_S % iTimerKernel

    end select !-- S_S

  end subroutine CloneTimers


  subroutine Compute ( S, T_Option, iS_Option )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I, &
      V
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      S_4D, &
      F_I
    type ( TimerForm ), pointer :: &
      T_RS, &
      T_K

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 3 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 3 )

    associate &
      ( RS  =>  S % RiemannSolver, &
         G  =>  S % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_RS  =>  RS % Timer       ( LevelOption = T_Option % Level + 1 )
      T_K   =>   S % TimerKernel ( LevelOption = T_Option % Level + 1 )
    else
      T_RS  =>  null ( )
      T_K   =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )
    call S % Clear ( )
    if ( associated ( T_K ) ) call T_K % Stop ( )

    do iC  =  1,  RS % Atlas % nCharts
       
      associate ( C  =>  S % Atlas % Chart ( 1 ) % Element )
      do iD  =  1, C % nDimensions

        if ( associated ( T_RS ) ) then
          call T_RS % Start ( )
          call RS % Compute &
                 ( iC, iD, T_Option = T_RS, iS_Option = iS_Option )
          call T_RS % Stop ( )
        else
          call RS % Compute &
                 ( iC, iD, iS_Option = iS_Option )
        end if

        if ( associated ( T_K ) ) call T_K % Start ( )

        associate &
          (  SV  =>   S % Storage ( iC ) % Value, &
            RSV  =>  RS % Storage ( iC ) % Value, &
             GV  =>   G % Storage ( iC ) % Value )

        select type ( C )
        class is ( Chart_GS_Form )

          call C % SetFieldPointer (  SV ( :, : ), S_4D )
          call C % SetFieldPointer ( RSV ( :, : ), F_I )
          call C % SetFieldPointer (  GV ( :, G % AREA_I_D ( iD ) ), A_I )
          call C % SetFieldPointer (  GV ( :, G % VOLUME ), V )

          call ComputeKernel &
                 ( S_4D, F_I, A_I, V, iD, C % nGhostLayers ( iD ), &
                   UseDeviceOption = S % DeviceMemory )

        class default
          call Show ( 'Chart type not recognized', CONSOLE % ERROR )
          call Show ( 'Slope_DFV_PD__Form', 'module', CONSOLE % ERROR )
          call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- C

        end associate !-- SV, etc.

        if ( associated ( T_K ) ) call T_K % Stop ( )

      end do !-- iD
      end associate !-- C

    end do !-- iC
    end associate !-- RS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S

    nullify ( S % RiemannSolver )
    nullify ( S % CurrentSet )
    
  end subroutine Finalize


end module Slope_DFV_PD__Form
