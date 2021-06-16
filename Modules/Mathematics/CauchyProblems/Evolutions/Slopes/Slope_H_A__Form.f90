module Slope_H_A__Form

  !-- Slope_Header_Atlas_Form

  use Basics
  use Manifolds
  use Fields
  use Slope_H_C__Form

  implicit none
  private

    integer, private, parameter :: &
      MAX_COMPONENTS = 16

  type, public, extends ( FieldSet_A_Form ) :: Slope_H_A_Form
    integer ( KDI ) :: &
      nComponents = 0
    class ( Slope_H_A_Element ), dimension ( : ), pointer :: &
      Component_A => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_H_A_Form

  type, public :: Slope_H_A_Element
    class ( Slope_H_A_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Slope_H_A_Element


contains


  subroutine InitializeAllocate_FS &
               ( FSA, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Slope_H_A_Form ), intent ( inout ), target :: &
      FSA
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
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a Slope_H_A' 
    
    call FSA % FieldSet_A_Form % Initialize &
           ( A, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
             nFieldsOption, IgnorabilityOption )

    allocate ( FSA % Component_A ( MAX_COMPONENTS ) )

  end subroutine InitializeAllocate_FS


  subroutine Show_FS ( FSA )

    class ( Slope_H_A_Form ), intent ( in ) :: &
      FSA

    integer ( KDI ) :: &
      iC

    call FSA % Show ( )
 
    do iC  =  1, FSA % nComponents
      associate ( SCA  =>  FSA % Component_A ( iC ) % Element )
      call SCA % Show ( )
      end associate !-- SCA
    end do !-- iC

  end subroutine Show_FS


  subroutine Compute ( SA, TimerLevelOption, iS_Option )

    class ( Slope_H_A_Form ), intent ( inout ) :: &
      SA
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption, &
      iS_Option

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( SA % FieldSet_C )
      select type ( SC  =>  SA % FieldSet_C ( iC ) % Element )
        class is ( Slope_H_C_Form )
      call SC % Compute ( TimerLevelOption, iS_Option )
      end select !-- SC
    end do !-- iC

  end subroutine Compute


  subroutine Finalize ( SA )

    type ( Slope_H_A_Form ), intent ( inout ) :: &
      SA

    if ( associated ( SA % Component_A ) ) &
      deallocate ( SA % Component_A )

  end subroutine Finalize


  impure elemental subroutine Finalize_E ( SE )
    
    type ( Slope_H_A_Element ), intent ( inout ) :: &
      SE

    if ( allocated ( SE % Element ) ) &
      deallocate ( SE % Element )

  end subroutine Finalize_E


end module Slope_H_A__Form
