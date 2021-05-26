module Slope_H_C__Form

  !-- Slope_Header_Chart_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

  type, public, extends ( FieldSet_C_Form ) :: Slope_H_C_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_H_C_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSC, C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Slope_H_C_Form ), intent ( inout ) :: &
      FSC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
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

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a Slope_H_C' 
    
    call FSC % FieldSet_C_Form % Initialize &
           ( C, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
             nFieldsOption, IgnorabilityOption )

  end subroutine InitializeAllocate_FS


  subroutine Compute ( SC, TimerLevelOption, iS_Option )

    class ( Slope_H_C_Form ), intent ( inout ) :: &
      SC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption, &
      iS_Option

    call Show ( 'Compute must be overridden', CONSOLE % ERROR )
    call Show ( 'Slope_H_C_Form', 'module', CONSOLE % ERROR )
    call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine Compute


  subroutine Finalize ( SC )

    type ( Slope_H_C_Form ), intent ( inout ) :: &
      SC

  end subroutine Finalize


end module Slope_H_C__Form
