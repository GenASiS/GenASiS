module Slope_H_A__Form

  !-- Slope_Header_Atlas_Form

  use Basics
  use Manifolds
  use Fields
  use Slope_H_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: Slope_H_A_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_H_A_Form


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

  end subroutine InitializeAllocate_FS


  subroutine Compute ( SA, TimerLevelOption )

    class ( Slope_H_A_Form ), intent ( inout ) :: &
      SA
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( SA % FieldSet_C )
      select type ( SC  =>  SA % FieldSet_C ( iC ) % Element )
      class is ( Slope_H_C_Form )
      call SC % Compute ( TimerLevelOption )
      end select !-- SC
    end do !-- iC

  end subroutine Compute


  subroutine Finalize ( SA )

    type ( Slope_H_A_Form ), intent ( inout ) :: &
      SA

  end subroutine Finalize


end module Slope_H_A__Form
