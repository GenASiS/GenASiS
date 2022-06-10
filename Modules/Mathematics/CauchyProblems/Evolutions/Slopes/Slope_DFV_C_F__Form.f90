module Slope_DFV_C_F__Form

  !-- Slope_DivergenceFiniteVolume_Connection_Flat__Form

  use Basics
  use Manifolds
  use Fields
  use Slope_H__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_C_F_Form
    integer ( KDI ) :: &
      iTimer_S = 0, &  !-- Stresses
      iTimer_K = 0     !-- Kernel
    type ( FieldSetForm ), allocatable :: &
      Stress_UD
    class ( DivergencePart_CS_Form ), pointer :: &
      DivergencePart => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_C_F
    generic, public :: &
      Initialize => InitializeAllocate_C_F
    procedure, public, pass :: &
      ComputeChart
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_C_F_Form

    private :: &
      Compute_C_Kernel, &
      Compute_S_Kernel

    interface

      module subroutine Compute_C_Kernel &
               ( S_M_1, S_UD_33, A_I_1, V, oV, UseDeviceOption )
        !-- Compute_Cylindrical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          S_M_1
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          S_UD_33, &
          A_I_1, &
          V
        integer ( KDI ), intent ( in ) :: &
          oV
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_C_Kernel

      module subroutine Compute_S_Kernel &
               ( S_M_1, S_M_2, S_UD_22, S_UD_33, A_I_1, A_I_2, V, oV, &
                 UseDeviceOption )
        !-- Compute_Spherical_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          S_M_1, S_M_2
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          S_UD_22, S_UD_33, &
          A_I_1, A_I_2, &
          V
        integer ( KDI ), intent ( in ) :: &
          oV
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_S_Kernel

    end interface


contains


  subroutine InitializeAllocate_C_F ( S, DP, IgnorabilityOption )

    class ( Slope_DFV_C_F_Form ), intent ( inout ) :: &
      S
    class ( DivergencePart_CS_Form ), intent ( in ), target :: &
      DP
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_C_F' 
    
    associate ( CS  =>  DP % CurrentSet )

    Name  =  trim ( CS % Name ) // '_Slp_DFV_C_F_' // trim ( DP % Name )

    S % DivergencePart  =>  DP

    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = IgnorabilityOption )

    allocate ( S % Stress_UD )
    associate ( S_UD  =>  S % Stress_UD )
    call S_UD % Initialize &
           ( CS % Atlas, &
             FieldOption = [ 'Stress_UD_22', 'Stress_UD_33' ], &
             NameOption = 'Stress_UD', &
             DeviceMemoryOption = CS % DeviceMemory, &
             nFieldsOption = 2, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- S_UD

    end associate !-- CS

  end subroutine InitializeAllocate_C_F


  subroutine ComputeChart ( S, iC, T_Option )

    class ( Slope_DFV_C_F_Form ), intent ( inout ) :: &
      S
    integer ( KDI ) :: &
      iC  !-- iChart
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iMomentum_1, &
      iMomentum_2
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I_1, A_I_2, &
      V, &
      S_UD_22, &
      S_UD_33, &
      S_M_1, &
      S_M_2
    type ( TimerForm ), pointer :: &
      T_C, &
      T_S, &
      T_K

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )
    call Show ( iC, 'iChart', S % IGNORABILITY + 2 )

    associate &
      ( DP  =>  S % DivergencePart ) 
    associate &
      ( S_UD  =>   S % Stress_UD, &
        G     =>  DP % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_S   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_S, &
                    Name = trim ( S % Name ) // '_Strss', &
                    Level = T_Option % Level + 1 )
      T_K   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_K, &
                    Name = trim ( S % Name ) // '_Krnl', &
                    Level = T_Option % Level + 1 )
    else
      T_S  =>  null ( )
      T_K  =>  null ( )
    end if

    associate ( C  =>  S % Atlas % Chart ( iC ) % Element )

    iMomentum_1  =  -1
    iMomentum_2  =  -1

    if ( associated ( T_S ) ) call T_S % Start ( )
    call DP % ComputeStresses ( S_UD, iC, iMomentum_1, iMomentum_2 )
    if ( associated ( T_S ) ) call T_S % Stop ( )

    if ( iMomentum_1  <  0  .or.  iMomentum_2  <  0 ) &
      return

    if ( associated ( T_K ) ) call T_K % Start ( )

    associate &
      ( SV      =>   S    % Storage ( iC ) % Value, &
        S_UD_V  =>   S_UD % Storage ( iC ) % Value, &
        GV      =>   G    % Storage ( iC ) % Value )

    select type ( C )
    class is ( Chart_GS_Form )

      call C % SetFieldPointer ( SV     ( :, iMomentum_1 ),    S_M_1 )
      call C % SetFieldPointer ( SV     ( :, iMomentum_2 ),    S_M_2 )
      call C % SetFieldPointer ( S_UD_V ( :, 1 ),              S_UD_22 )
      call C % SetFieldPointer ( S_UD_V ( :, 2 ),              S_UD_33 )
      call C % SetFieldPointer ( GV     ( :, G % AREA_I_D_1 ), A_I_1 )
      call C % SetFieldPointer ( GV     ( :, G % AREA_I_D_2 ), A_I_2 )
      call C % SetFieldPointer ( GV     ( :, G % VOLUME ),     V )

      select case ( trim ( C % CoordinateSystem ) )
      case ( 'CYLINDRICAL' )
        call Compute_C_Kernel &
               ( S_M_1, S_UD_33, A_I_1, V, C % nGhostLayers ( 1 ), &
                 UseDeviceOption = S % DeviceMemory )
      case ( 'SPHERICAL' )
        call Compute_S_Kernel &
               ( S_M_1, S_M_2, S_UD_22, S_UD_33, A_I_1, A_I_2, V, &
                 C % nGhostLayers ( 1 ), UseDeviceOption = S % DeviceMemory )
      end select !-- CoordinateSystem

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Slope_DFV_C_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- SV, etc.
      
    if ( associated ( T_K ) ) call T_K % Stop ( )

    end associate !-- C

    end associate !-- S_UD, etc.
    end associate !-- DP

  end subroutine ComputeChart


  subroutine Compute ( S, T_Option )

    class ( Slope_DFV_C_F_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iMomentum_1, &
      iMomentum_2
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      A_I_1, A_I_2, &
      V, &
      S_UD_22, &
      S_UD_33, &
      S_M_1, &
      S_M_2
    type ( TimerForm ), pointer :: &
      T_S, &
      T_K

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate &
      ( DT  =>  S % DivergencePart ) 
    associate &
      ( S_UD  =>   S % Stress_UD, &
        G     =>  DT % CurrentSet % Geometry )

    if ( present ( T_Option ) ) then
      T_S   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_S, &
                    Name = trim ( S % Name ) // '_Strss', &
                    Level = T_Option % Level + 1 )
      T_K   =>  PROGRAM_HEADER % Timer &
                  ( Handle = S % iTimer_K, &
                    Name = trim ( S % Name ) // '_Krnl', &
                    Level = T_Option % Level + 1 )
    else
      T_S  =>  null ( )
      T_K  =>  null ( )
    end if

    do iC  =  1,  S % Atlas % nCharts
       
      associate ( C  =>  S % Atlas % Chart ( iC ) % Element )

      iMomentum_1  =  -1
      iMomentum_2  =  -1

      if ( associated ( T_S ) ) call T_S % Start ( )
      call DT % ComputeStresses ( S_UD, iC, iMomentum_1, iMomentum_2 )
      if ( associated ( T_S ) ) call T_S % Stop ( )

      if ( iMomentum_1  <  0  .or.  iMomentum_2  <  0 ) &
        return

      if ( associated ( T_K ) ) call T_K % Start ( )

      associate &
        ( SV      =>   S    % Storage ( iC ) % Value, &
          S_UD_V  =>   S_UD % Storage ( iC ) % Value, &
          GV      =>   G    % Storage ( iC ) % Value )

      select type ( C )
      class is ( Chart_GS_Form )

        call C % SetFieldPointer ( SV     ( :, iMomentum_1 ),    S_M_1 )
        call C % SetFieldPointer ( SV     ( :, iMomentum_2 ),    S_M_2 )
        call C % SetFieldPointer ( S_UD_V ( :, 1 ),              S_UD_22 )
        call C % SetFieldPointer ( S_UD_V ( :, 2 ),              S_UD_33 )
        call C % SetFieldPointer ( GV     ( :, G % AREA_I_D_1 ), A_I_1 )
        call C % SetFieldPointer ( GV     ( :, G % AREA_I_D_2 ), A_I_2 )
        call C % SetFieldPointer ( GV     ( :, G % VOLUME ),     V )

        select case ( trim ( C % CoordinateSystem ) )
        case ( 'CYLINDRICAL' )
          call Compute_C_Kernel &
                 ( S_M_1, S_UD_33, A_I_1, V, C % nGhostLayers ( 1 ), &
                   UseDeviceOption = S % DeviceMemory )
        case ( 'SPHERICAL' )
          call Compute_S_Kernel &
                 ( S_M_1, S_M_2, S_UD_22, S_UD_33, A_I_1, A_I_2, V, &
                   C % nGhostLayers ( 1 ), UseDeviceOption = S % DeviceMemory )
        end select !-- CoordinateSystem

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_DFV_C_F__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C

      end associate !-- SV, etc.
        
      if ( associated ( T_K ) ) call T_K % Stop ( )

      end associate !-- C

    end do !-- iC

    end associate !-- S_UD, etc.
    end associate !-- DT

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_C_F_Form ), intent ( inout ) :: &
      S

    nullify ( S % DivergencePart )

    if ( allocated ( S % Stress_UD ) ) &
      deallocate ( S % Stress_UD )

  end subroutine Finalize


end module Slope_DFV_C_F__Form
