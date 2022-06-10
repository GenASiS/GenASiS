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
      iTimer_C = 0, &  !-- Clear
      iTimer_K = 0     !-- Kernel
    class ( DivergencePart_CS_Form ), pointer :: &
      DivergencePart => null ( )
    class ( RiemannSolver_HLL_Form ), pointer :: &
      RiemannSolver => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_PD
    generic, public :: &
      Initialize => InitializeAllocate_PD
    procedure, public, pass :: &
      ComputeDimension
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_PD_Form

    private :: &
      RecordBoundaryFlux_SCG

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

      module subroutine RecordBoundaryFlux_SCG_Kernel &
               ( F, nB, oB, BF, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          F
        integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
          nB, &
          oB
        real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
          BF
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine RecordBoundaryFlux_SCG_Kernel
      
    end interface


contains


  subroutine InitializeAllocate_PD ( S, RS, DP, IgnorabilityOption )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( in ), target :: &
      DP
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_PD'

    associate ( CS  =>  RS % CurrentSet )

    Name  =  trim ( CS % Name ) // '_Slp_DFV_PD_' // trim ( DP % Name )

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


  subroutine ComputeDimension ( S, iC, iD, T_Option )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iF  !-- iFlux
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I, &
      V
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_I_3D
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      S_4D, &
      F_I_4D
    type ( TimerForm ), pointer :: &
      T_RS, &
      T_K

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )
    call Show ( iD, 'Dimension', S % IGNORABILITY + 2 )

    associate &
      ( RS  =>  S % RiemannSolver, &
        DP  =>  S % DivergencePart, &
        CS  =>  S % RiemannSolver % CurrentSet, &
         G  =>  S % RiemannSolver % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_RS  =>  RS % Timer_C ( Level = T_Option % Level + 1 )
      T_K   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_K, &
                    Name = trim ( S % Name ) // '_K', &
                    Level = T_Option % Level + 1 )
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
      call RS % ComputeFlux ( DP, iC, iD, T_Option = T_RS )
      call T_RS % Stop ( )
    else
      call RS % ComputeFlux ( DP, iC, iD )
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
      call C % SetFieldPointer ( RSV ( :, : ), F_I_4D )
      call C % SetFieldPointer (  GV ( :, G % AREA_I_D ( iD ) ), A_I )
      call C % SetFieldPointer (  GV ( :, G % VOLUME ), V )

      call ComputeKernel &
             ( S_4D, F_I_4D, A_I, V, iD, C % nGhostLayers ( iD ), &
               UseDeviceOption = S % DeviceMemory )

      do iF  =  1,  size ( S_4D, dim = 4 ) 
        call C % SetFieldPointer ( RSV ( :, iF ), F_I_3D )
        call RecordBoundaryFlux_SCG &
               ( CS % BoundaryFlux_SCG, C, F_I_3D, iD, iF )
      end do !-- iF

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


  subroutine Compute ( S, T_Option )

    class ( Slope_DFV_PD_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iD, &  !-- iDimension
      iF     !-- iFlux
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I, &
      V
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_I_3D
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      S_4D, &
      F_I_4D
    type ( TimerForm ), pointer :: &
      T_C, &
      T_RSP, &
      T_RSC, &
      T_K

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate &
      ( RS  =>  S % RiemannSolver, &
        DP  =>  S % DivergencePart, &
        CS  =>  S % RiemannSolver % CurrentSet, &
         G  =>  S % RiemannSolver % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_C   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_C, &
                    Name = trim ( S % Name ) // '_Clr', &
                    Level = T_Option % Level + 1 )
    else
      T_C   =>  null ( )
    end if
    if ( associated ( T_C ) ) call T_C % Start ( )
    call S % Clear ( )
    if ( associated ( T_C ) ) call T_C % Stop ( )

    do iC  =  1,  S % Atlas % nCharts
       
      associate ( C  =>  S % Atlas % Chart ( iC ) % Element )
      do iD  =  1, C % nDimensions

        if ( present ( T_Option ) ) then
          T_RSP  =>  RS % Timer_P ( Level = T_Option % Level + 1 )
        else
          T_RSP  =>  null ( )
        end if
        if ( associated ( T_RSP ) ) then
          call T_RSP % Start ( )
          call RS % Prepare ( iC, iD, T_Option = T_RSP )
          call T_RSP % Stop ( )
        else
          call RS % Prepare ( iC, iD )
        end if

        if ( present ( T_Option ) ) then
          T_RSC  =>  RS % Timer_C ( Level = T_Option % Level + 1 )
        else
          T_RSC  =>  null ( )
        end if
        if ( associated ( T_RSC ) ) then
          call T_RSC % Start ( )
          call RS % Compute ( DP, iC, iD, T_Option = T_RSC )
          call T_RSC % Stop ( )
        else
          call RS % Compute ( DP, iC, iD )
        end if

        if ( present ( T_Option ) ) then
          T_K   =>  PROGRAM_HEADER % Timer &
                      ( Handle = S % iTimer_K, &
                        Name = trim ( S % Name ) // '_Krnl', &
                        Level = T_Option % Level + 1 )
        else
          T_K   =>  null ( )
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
          call C % SetFieldPointer ( RSV ( :, : ), F_I_4D )
          call C % SetFieldPointer (  GV ( :, G % AREA_I_D ( iD ) ), A_I )
          call C % SetFieldPointer (  GV ( :, G % VOLUME ), V )
          
          call ComputeKernel &
                 ( S_4D, F_I_4D, A_I, V, iD, C % nGhostLayers ( iD ), &
                   UseDeviceOption = S % DeviceMemory )

          do iF  =  1,  size ( S_4D, dim = 4 ) 
            call C % SetFieldPointer ( RSV ( :, iF ), F_I_3D )
            call RecordBoundaryFlux_SCG &
                   ( CS % BoundaryFlux_SCG, C, F_I_3D, iD, iF )
          end do !-- iF

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

    nullify ( S % RiemannSolver )
    nullify ( S % DivergencePart )
  end subroutine Finalize


  subroutine RecordBoundaryFlux_SCG ( BF, C, F_I, iD, iF )

    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      BF
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_I
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iF     !-- iFlux

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
      associate ( BF_Inner  =>  BF ( iF, iCI ) % Value )
      oB  =  C % nGhostLayers
      call RecordBoundaryFlux_SCG_Kernel &
             ( F_I, nB, oB, BF_Inner, &
               UseDeviceOption = BF ( iF, iCI ) % AllocatedDevice )
      end associate !-- BF_Inner
      end associate !-- iCI
    end if !-- iaBrick ( iD ) == 1

    if ( RecordOuter ) then
      associate ( iCO  =>  C % Connectivity % iaOuter ( iD ) )
      associate ( BF_Outer  =>  BF ( iF, iCO ) % Value )
      oB         =  C % nGhostLayers
      oB ( iD )  =  oB ( iD )  +  nCells
      call RecordBoundaryFlux_SCG_Kernel &
             ( F_I, nB, oB, BF_Outer, &
               UseDeviceOption = BF ( iF, iCO ) % AllocatedDevice )
      end associate !-- BF_Outer
      end associate !-- iCO
    end if !-- iaBrick ( iD ) == nBricks ( iD )

  end subroutine RecordBoundaryFlux_SCG


end module Slope_DFV_PD__Form
