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
      iTimer    = 0, &
      iTimer_MA = 0  !-- MultiplyAdd
    logical ( KDL ) :: &
      StreamComponents = .true.
    type ( Slope_H_Element ), dimension ( : ), pointer :: &
      Component => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      AddComponents
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
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
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

    allocate ( FS % Component ( MAX_COMPONENTS ) )

  end subroutine InitializeAllocate_FS


  subroutine SetStream ( S, Sm )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    class ( StreamForm ), intent ( inout ) :: &
      Sm

    integer ( KDI ) :: &
      iC

    call Sm % AddFieldSet ( S )

    if ( .not. S % StreamComponents ) &
      return

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


  function Timer ( S, Level ) result ( T )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimer, &
               Name = S % Name, &
               Level = Level )

  end function Timer


  subroutine Compute ( S, T_Option )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iComponent
    type ( TimerForm ), pointer :: &
      T_C

    call Show ( 'Computing a Slope_H', S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    if ( S % nComponents  >  0 ) then

      do iC  =  1, S % nComponents
        associate &
          ( SC  =>  S % Component ( iC ) % Element )

        if ( present ( T_Option ) ) then
          T_C  =>  SC % Timer ( Level = T_Option % Level + 1 )
          call T_C % Start ( )
          call SC % Compute ( T_Option = T_C )
          call T_C % Stop ( )
        else
          call SC % Compute ( )
        end if

        end associate !-- SC
      end do !-- iC

      call S % AddComponents ( T_Option )

    else
      call Show ( 'Slope has no components', CONSOLE % ERROR )
      call Show ( 'Compute must be overridden', CONSOLE % ERROR )
      call Show ( S % Name, 'Name', CONSOLE % ERROR )
      call Show ( 'Slope_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine Compute


  subroutine AddComponents ( S, T_Option )

    class ( Slope_H_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iComponent
    type ( TimerForm ), pointer :: &
      T_MA

    call Show ( 'Adding components of a Slope_H', S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    if ( S % nComponents  >  0 ) then

      if ( present ( T_Option ) ) then
        T_MA  =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_MA, &
                      Name = trim ( S % Name ) // '_MA', &
                      Level = T_Option % Level + 1 )
      else
        T_MA  =>  null ( )
      end if
    
      if ( associated ( T_MA ) ) call T_MA % Start ( )
      call S % Clear ( )
      if ( associated ( T_MA ) ) call T_MA % Stop ( )

      do iC  =  1, S % nComponents
        associate &
          ( SC  =>  S % Component ( iC ) % Element )

        if ( associated ( T_MA ) ) call T_MA % Start ( )
        call S % MultiplyAdd ( SC, 1.0_KDR )
        if ( associated ( T_MA ) ) call T_MA % Stop ( )

        end associate !-- SC
      end do !-- iC

    end if

  end subroutine AddComponents


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
