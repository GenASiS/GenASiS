module Slope_H__Form

  !-- Slope_Header__Form

  use Basics
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
    character ( LDL ) :: &
      TimerName = ''
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
      CloneTimers
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      ClearRecursive
    procedure, public, pass :: &
      MultiplyAddRecursive
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
                 DevicesCommunicateOption, AssociateFieldsOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

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
      DevicesCommunicateOption, &
      AssociateFieldsOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      Ignorability

    Ignorability  =  A % IGNORABILITY + 1
    if ( present ( IgnorabilityOption ) ) &
      Ignorability  =  IgnorabilityOption

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a Slope_H' 
    
    call FS % FieldSetForm % Initialize &
           ( A, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, AssociateFieldsOption, &
             UnitOption, VectorIndicesOption, nFieldsOption, &
             IgnorabilityOption = Ignorability )

    if ( FS % TimerName  ==  '' ) &
      FS % TimerName  =  FS % Name

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

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  S % iTimer )

    if ( iT == 0 ) then
      TimerName  =  S % TimerName
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

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  S % iTimerMultiplyAdd )

    if ( iT == 0 ) then
      TimerName  =  trim ( S % TimerName ) // '_MltplyAdd' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function TimerMultiplyAdd


  subroutine CloneTimers ( S, S_S )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    class ( Slope_H_Form ), intent ( in ) :: &
      S_S  !-- S_Source

    integer ( KDI ) :: &
      iC  !-- iComponent

    S % iTimer             =  S_S % iTimer
    S % iTimerMultiplyAdd  =  S_S % iTimerMultiplyAdd

    call S % CloneGhostTimers ( S_S )

    do iC  =  1, S % nComponents
      associate &
        ( SC    =>  S   % Component ( iC ) % Element, &
          SC_S  =>  S_S % Component ( iC ) % Element )
      call SC % CloneTimers ( SC_S )
      end associate !-- SC, etc.
    end do !-- iC

  end subroutine CloneTimers


  subroutine Compute ( S, T_Option, iS_Option )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    integer ( KDI ) :: &
      iC  !-- iComponent
    type ( TimerForm ), pointer :: &
      T_MA, &
      T_C

    call Show ( 'Computing a Slope_H', S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

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
          call SC % Compute ( iS_Option = iS_Option )
        end if

        if ( associated ( T_MA ) ) call T_MA % Start ( )
        call S % MultiplyAdd ( SC, 1.0_KDR )
        if ( associated ( T_MA ) ) call T_MA % Stop ( )

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


  subroutine ClearRecursive ( S )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC  !-- iComponent

    !-- This slope

    call S % Clear ( )

    !-- Component slopes

    do iC  =  1, S % nComponents
      associate ( SC  =>   S % Component ( iC ) % Element )

      call SC % ClearRecursive ( )

      end associate !-- SCA
    end do !-- iC

  end subroutine ClearRecursive


  subroutine MultiplyAddRecursive ( S, SS, B )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S  !-- Slope
    class ( Slope_H_Form ), intent ( in ) :: &
      SS  !-- SlopeStage
    real ( KDR ) :: &
      B  !-- RungeKutta weight

    integer ( KDI ) :: &
      iC  !-- iComponent

    !-- This slope

    call S % MultiplyAdd ( SS, B )

    !-- Component slopes

    do iC  =  1, S % nComponents
      associate &
        (  SC  =>   S % Component ( iC ) % Element, &
          SSC  =>  SS % Component ( iC ) % Element )

      call SC % MultiplyAddRecursive ( SSC, B )

      end associate !-- SCA
    end do !-- iC

  end subroutine MultiplyAddRecursive


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
