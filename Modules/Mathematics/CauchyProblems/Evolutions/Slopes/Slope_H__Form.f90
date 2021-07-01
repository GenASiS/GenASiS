module Slope_H__Form

  !-- Slope_Header__Form

  use Basics
  use Algebra
  use Manifolds
  use Fields

  implicit none
  private

    integer, private, parameter :: &
      MAX_COMPONENTS = 16

  type, public, extends ( FieldSetForm ) :: Slope_H_Form
    integer ( KDI ) :: &
      nComponents = 0
    integer ( KDI ) :: &
      iTimer            = 0, &
      iTimerMultiplyAdd = 0
    type ( Slope_H_Element ), dimension ( : ), pointer :: &
      Component => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass ( S ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Timer
    procedure, private, pass :: &
      TimerMultiplyAdd
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      Increment
    final :: &
      Finalize
  end type Slope_H_Form

  type, public :: Slope_H_Element
    class ( Slope_H_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Slope_H_Element


contains


  subroutine InitializeAllocate_FS &
               ( FS, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Slope_H_Form ), intent ( inout ), target :: &
      FS
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a Slope_H' 
    
    call FS % FieldSetForm % Initialize &
           ( A, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
             nFieldsOption, IgnorabilityOption )

    allocate ( FS % Component ( MAX_COMPONENTS ) )

  end subroutine InitializeAllocate_FS


  subroutine SetStream ( Sm, S )

    class ( StreamForm ), intent ( inout ) :: &
      Sm
    class ( Slope_H_Form ), intent ( in ) :: &
      S

    integer ( KDI ) :: &
      iC

    call Sm % AddFieldSet ( S )

    do iC  =  1, S % nComponents
      associate ( SC  =>  S % Component ( iC ) % Element )
      call SC % SetStream ( Sm )
      end associate !-- SCA
    end do

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( Slope_H_Form ), intent ( in ) :: &
      FS

    integer ( KDI ) :: &
      iC

    call FS % FieldSetForm % Show ( )

    call Show ( FS % nComponents, 'nComponents', FS % IGNORABILITY ) 
    do iC  =  1, FS % nComponents
      associate ( SC  =>  FS % Component ( iC ) % Element )
      call Show ( SC % Name, 'Component', FS % IGNORABILITY )
      end associate !-- SC
    end do !-- iC

    do iC  =  1, FS % nComponents
      associate ( SC  =>  FS % Component ( iC ) % Element )
      call SC % Show ( )
      end associate !-- SC
    end do !-- iC

  end subroutine Show_FS


  function Timer ( S, LevelOption ) result ( T )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimer )

    if ( iT == 0 ) then
      TimerName  =  S % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  function TimerMultiplyAdd ( S, LevelOption ) result ( T )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  S % iTimerMultiplyAdd )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % Name ) // '_MltplyAdd' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function TimerMultiplyAdd


  subroutine Compute ( S, T_Option, iS_Option )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), pointer, optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    integer ( KDI ) :: &
      iC, &  !-- iComponent
      iCrt   !-- iChart
    type ( TimerForm ), pointer :: &
      T_MA, &
      T_C

    if ( S % nComponents  >  0 ) then

      if ( present ( T_Option ) ) then
        T_MA  =>  S % TimerMultiplyAdd ( LevelOption = T_Option % Level + 1 )
      else
        T_MA  =>  null ( )
      end if
    
      if ( associated ( T_MA ) ) call T_MA % Start ( )
      call S % Clear ( )
      if ( associated ( T_MA ) ) call T_MA % Stop ( )

      do iC  =  1, S % nComponents
        associate &
          ( SC  =>  S % Component ( iC ) % Element )

        if ( present ( T_Option ) ) then
          T_C  =>  SC % Timer ( LevelOption = T_Option % Level + 1 )
          call T_C % Start ( )
          call SC % Compute ( T_Option = T_C, iS_Option = iS_Option )
          call T_C % Stop ( )
        else
          call SC % Compute ( iS_Option = iS_Option)
        end if

        do iCrt  =  1,  S % Atlas % nCharts
          associate &
            (  SV  =>   S % Storage ( iCrt ) % Value, &
              SCV  =>  SC % Storage ( iCrt ) % Value )

          if ( associated ( T_MA ) ) call T_MA % Start ( )
          call MultiplyAdd &
                 ( SV, SCV, 1.0_KDR, &
                   UseDeviceOption = S % DeviceMemory )
          if ( associated ( T_MA ) ) call T_MA % Stop ( )

        end associate !-- SV, etc.
        end do !-- iCrt

        end associate !-- SC
      end do !-- iC

    else
      call Show ( 'Slope has no components', CONSOLE % ERROR )
      call Show ( 'Compute must be overridden', CONSOLE % ERROR )
      call Show ( S % Name, 'Name', CONSOLE % ERROR )
      call Show ( 'Slope_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine Compute


  subroutine Increment ( S, SS, B, iS )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S  !-- Slope
    class ( Slope_H_Form ), intent ( in ) :: &
      SS  !-- SlopeStage
    real ( KDR ) :: &
      B  !-- RungeKutta weight
    integer ( KDI ) :: &
      iS  !-- RungeKutta iStage

    integer ( KDI ) :: &
      iC, &  !-- iComponent
      iCrt   !-- iChart

    !-- This slope

    call S % Clear ( )

    do iCrt  =  1,  S % Atlas % nCharts
      associate &
        (  SV  =>   S % Storage ( iCrt ) % Value, &
          SSV  =>  SS % Storage ( iCrt ) % Value )

      call MultiplyAdd &
             ( SV, SSV, B, &
               UseDeviceOption = S % DeviceMemory )
    
      end associate !-- SV, etc.
    end do !-- iCrt

    !-- Component slopes

    do iC  =  1, S % nComponents
      associate &
        (  SC  =>   S % Component ( iC ) % Element, &
          SSC  =>  SS % Component ( iC ) % Element )

      call SC % Increment ( SSC, B, iS )

      end associate !-- SCA
    end do !-- iC

  end subroutine Increment


  impure elemental subroutine Finalize ( S )

    type ( Slope_H_Form ), intent ( inout ) :: &
      S

    if ( associated ( S % Component ) ) &
      deallocate ( S % Component )

  end subroutine Finalize


  impure elemental subroutine Finalize_E ( SE )
    
    type ( Slope_H_Element ), intent ( inout ) :: &
      SE

    if ( allocated ( SE % Element ) ) &
      deallocate ( SE % Element )

  end subroutine Finalize_E


end module Slope_H__Form
