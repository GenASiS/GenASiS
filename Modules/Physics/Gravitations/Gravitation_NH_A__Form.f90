module Gravitation_NH_A__Form

  !-- Gravitation_Galileo_Atlas_Form

  use Basics
  use Mathematics
  use Gravitation_NH_C__Form

  implicit none
  private

  type, public, extends ( Geometry_F_A_Form ) :: Gravitation_NH_A_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    final :: &
      Finalize
  end type Gravitation_NH_A_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSA, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Gravitation_NH_A_Form ), intent ( inout ), target :: &
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

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a Gravitation_NH_A'

    Name  =  'Gravitation'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    associate ( nC  =>  A % nCharts )

    if ( allocated ( FSA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( FSA % FieldSet_C ( nC ) )
    end if

    call FSA % Geometry_F_A_Form % Initialize &
           ( A, &
             NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Gravitation_NH_C_Form :: FSA % FieldSet_C ( iC ) % Element ) 
        select type ( GC  =>  FSA % FieldSet_C ( iC ) % Element )
          class is ( Gravitation_NH_C_Form )
        associate ( C  =>  A % Chart ( iC ) % Element )

        call GC % Initialize &
               ( C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

        end associate !-- C
        end select !-- GC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_FS


  impure elemental subroutine Finalize ( GA )

    type ( Gravitation_NH_A_Form ), intent ( inout ) :: &
      GA

  end subroutine Finalize


end module Gravitation_NH_A__Form
