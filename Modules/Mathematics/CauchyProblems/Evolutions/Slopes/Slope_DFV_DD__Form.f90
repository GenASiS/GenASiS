module Slope_DFV_DD__Form

  !-- Slope_DivergenceFiniteVolume_DivergenceDiffusion__Form

  use Basics
  use Manifolds
  use Fields
  use RiemannSolver_HLL__Form
  use Slope_H__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_DD_Form
    integer ( KDI ) :: &
      iTimer_K = 0
    class ( RiemannSolver_HLL_Form ), pointer :: &
      RiemannSolver => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_DD
    generic, public :: &
      Initialize => InitializeAllocate_DD
    procedure, public, pass :: &
      ComputeDimension
    final :: &
      Finalize
  end type Slope_DFV_DD_Form

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


  subroutine InitializeAllocate_DD ( S, RS, IgnorabilityOption )

    class ( Slope_DFV_DD_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_DD'

    associate ( CS  =>  RS % CurrentSet )

    Name  =  trim ( CS % Name ) // '_Slp_DFV_DD'

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

  end subroutine InitializeAllocate_DD


  subroutine ComputeDimension ( S, iC, iD, T_Option )

    class ( Slope_DFV_DD_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
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
         G  =>  S % RiemannSolver % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_RS  =>  RS % Timer_C ( Level = T_Option % Level + 1 )
      T_K   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_K, &
                    Name = trim ( S % Name ) // '_Krnl', &
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
      call RS % ComputeDiffusion ( iC, iD, T_Option = T_RS )
      call T_RS % Stop ( )
    else
      call RS % ComputeDiffusion ( iC, iD )
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
      call Show ( 'Slope_DFV_DD__Form', 'module', CONSOLE % ERROR )
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


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_DD_Form ), intent ( inout ) :: &
      S

    nullify ( S % RiemannSolver )
    
  end subroutine Finalize


end module Slope_DFV_DD__Form
