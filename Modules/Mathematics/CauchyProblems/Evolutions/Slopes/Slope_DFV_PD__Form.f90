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
      iTimer_K = 0
    class ( DivergencePart_CS_Form ), pointer :: &
      DivergencePart => null ( )
    class ( RiemannSolver_HLL_Form ), pointer :: &
      RiemannSolver => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_PD
    generic, public :: &
      Initialize => InitializeAllocate_PD
    procedure, private, pass :: &
      Timer_K
    procedure, public, pass :: &
      CloneTimers
    procedure, public, pass :: &
      ComputeDimension
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_PD_Form

    private :: &
      RecordBoundaryFluence_SCG

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

      module subroutine RecordBoundaryFluence_SCG_Kernel &
               ( BF, F, Factor, nB, oB, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          BF
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          F
        real ( KDR ), intent ( in ) :: &
          Factor
        integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
          nB, &
          oB
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine RecordBoundaryFluence_SCG_Kernel
      
    end interface


contains


  subroutine InitializeAllocate_PD &
               ( S, RS, DP, SuffixOption, IgnorabilityOption )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( in ), target :: &
      DP
    character ( * ), intent ( in ), optional :: &
      SuffixOption    
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_PD'

    associate ( CS  =>  RS % CurrentSet )

    if ( S % TimerName  ==  '' ) &
      S % TimerName  &
        =  'S_DFV_PD_' // trim ( DP % Name ) // '_' // trim ( CS % Name )

    Name  =  'S_DFV_PD_' // trim ( DP % Name ) // '_' // trim ( CS % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    S % DivergencePart  =>  DP
    S % RiemannSolver   =>  RS

    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = IgnorabilityOption )

    end associate !-- CS

  end subroutine InitializeAllocate_PD


  function Timer_K ( S, LevelOption ) result ( T )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  S % iTimer_K )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % TimerName ) // '_K' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_K


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

    S % iTimer_K  =  S_S % iTimer_K

    end select !-- S_S

  end subroutine CloneTimers


  subroutine ComputeDimension ( S, iC, iD, iS, T_Option )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD, &  !-- iDimension
      iS     !-- iStage
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I, &
      V
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      S_4D, &
      F_I
    type ( TimerForm ), pointer :: &
      T_RS, &
      T_K

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate &
      ( RS  =>  S % RiemannSolver, &
        DP  =>  S % DivergencePart, &
         G  =>  S % RiemannSolver % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_RS  =>  RS % Timer_C     ( LevelOption = T_Option % Level + 1 )
      T_K   =>   S % Timer_K ( LevelOption = T_Option % Level + 1 )
    else
      T_RS  =>  null ( )
      T_K   =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )
    if ( iD  ==  1 ) &
      call S % Clear ( )
    if ( associated ( T_K ) ) call T_K % Stop ( )

    associate ( C  =>  S % Atlas % Chart ( iC ) % Element )

    if ( associated ( T_RS ) ) then
      call T_RS % Start ( )
      call RS % ComputeFlux &
             ( DP, iC, iD, T_Option = T_RS, iS_Option = iS )
      call T_RS % Stop ( )
    else
      call RS % ComputeFlux &
             ( DP, iC, iD, iS_Option = iS )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )

    call S % Storage ( iC ) % ReassociateHost &
           ( AssociateVariablesOption = .false. )
      
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
      
    call S % Storage ( iC ) % ReassociateHost &
           ( AssociateVariablesOption = .true. )

    if ( associated ( T_K ) ) call T_K % Stop ( )

    end associate !-- C
    end associate !-- RS, etc.

  end subroutine ComputeDimension


  subroutine Compute ( S, iS, T_Option )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

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

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate &
      ( RS  =>  S % RiemannSolver, &
        DP  =>  S % DivergencePart, &
         G  =>  S % RiemannSolver % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_RS  =>  RS % Timer_C ( LevelOption = T_Option % Level + 1 )
      T_K   =>   S % Timer_K ( LevelOption = T_Option % Level + 1 )
    else
      T_RS  =>  null ( )
      T_K   =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )
    call S % Clear ( )
    if ( associated ( T_K ) ) call T_K % Stop ( )

    do iC  =  1,  S % Atlas % nCharts
       
      associate ( C  =>  S % Atlas % Chart ( iC ) % Element )
      do iD  =  1, C % nDimensions

        if ( associated ( T_RS ) ) then
          call T_RS % Start ( )
          call RS % Prepare &
                 ( iC, iD, T_Option = T_RS, iS_Option = iS )
          call RS % Compute &
                 ( DP, iC, iD, T_Option = T_RS, iS_Option = iS )
          call T_RS % Stop ( )
        else
          call RS % Compute &
                 ( DP, iC, iD, iS_Option = iS )
        end if

        if ( associated ( T_K ) ) call T_K % Start ( )

        call S % Storage ( iC ) % ReassociateHost &
               ( AssociateVariablesOption = .false. )
        
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
        
        call S % Storage ( iC ) % ReassociateHost &
               ( AssociateVariablesOption = .true. )

        if ( associated ( T_K ) ) call T_K % Stop ( )

      end do !-- iD
      end associate !-- C

    end do !-- iC
    end associate !-- RS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S

    nullify ( S % DivergencePart )
    nullify ( S % RiemannSolver )
    
  end subroutine Finalize


  subroutine RecordBoundaryFluence_SCG &
               ( BF, C, F_I, Weight_RK, dT, iD, iC )

    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_I
    real ( KDR ), intent ( in ) :: &
      Weight_RK, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iC     !-- iConserved

    integer ( KDI ) :: &
      jD, kD, &   !-- jDimension, kDimension
      nCells
    integer ( KDI ), dimension ( 3 ) :: &
      oB, & !-- oBoundary
      nB    !-- nBoundary
    logical ( KDL ) :: &
      RecordInner, &
      RecordOuter

    RecordInner  =  ( C % iaBrick ( iD )  ==  1 )
    RecordOuter  =  ( C % iaBrick ( iD )  ==  C % nBricks ( iD ) )
    nCells  =  C % nCellsBrick ( iD )

    jD  =  mod ( iD, 3 ) + 1
    kD  =  mod ( jD, 3 ) + 1
    
    nB ( iD )  =  1
    nB ( jD )  =  C % nCellsBrick ( jD )
    nB ( kD )  =  C % nCellsBrick ( kD )

    if ( RecordInner ) then
      associate ( iCI  =>  C % Connectivity % iaInner ( iD ) )
      associate ( BF_Inner  =>  BF ( iC, iCI ) % Value )
      oB  =  C % nGhostLayers
      call RecordBoundaryFluence_SCG_Kernel &
             ( BF_Inner, F_I, Weight_RK * dT, nB, oB, &
               UseDeviceOption = BF ( iC, iCI ) % AllocatedDevice )
      end associate !-- BF_Inner
      end associate !-- iCI
    end if !-- iaBrick ( iD ) == 1

    if ( RecordOuter ) then
      associate ( iCO  =>  C % Connectivity % iaOuter ( iD ) )
      associate ( BF_Outer  =>  BF ( iC, iCO ) % Value )
      oB         =  C % nGhostLayers
      oB ( iD )  =  oB ( iD )  +  nCells
      call RecordBoundaryFluence_SCG_Kernel &
             ( BF_Outer, F_I, Weight_RK * dT, nB, oB, &
               UseDeviceOption = BF ( iC, iCO ) % AllocatedDevice )
      end associate !-- BF_Outer
      end associate !-- iCO
    end if !-- iaBrick ( iD ) == nBricks ( iD )

  end subroutine RecordBoundaryFluence_SCG


end module Slope_DFV_PD__Form
