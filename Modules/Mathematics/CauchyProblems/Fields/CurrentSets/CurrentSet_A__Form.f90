module CurrentSet_A__Form

  !-- CurrentSet_Atlas_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams
  use CurrentSet_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: CurrentSet_A_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
    procedure, public, pass ( CSA ) :: &
      SetStream
    final :: &
      Finalize
  end type CurrentSet_A_Form


contains


  subroutine InitializeAllocate_CS &
               ( CSA, A, Velocity_U_Unit, FieldOption, VectorOption, &
                 NameOption, DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( CurrentSet_A_Form ), intent ( inout ), target :: &
      CSA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit
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
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( CSA % Type  ==  '' ) &
      CSA % Type  =  'a CurrentSet_A'

    Name  =  'CurrentSet'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    associate ( nC  =>  A % nCharts )

    if ( allocated ( CSA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( CSA % FieldSet_C ( nC ) )
    end if

    call CSA % FieldSet_A_Form % Initialize &
           ( A, &
             NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( CurrentSet_C_Form :: CSA % FieldSet_C ( iC ) % Element ) 
        select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )

        associate ( C  =>  A % Chart ( iC ) % Element )
        call CSC % Initialize &
               ( C, Velocity_U_Unit, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )
        end associate !-- C

        end select !-- CSC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_CS


  subroutine SetStream ( SA, CSA )

    class ( Stream_A_Form ), intent ( inout ) :: &
      SA
    class ( CurrentSet_A_Form ), intent ( in ) :: &
      CSA

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, CSA % Atlas % nCharts
      associate ( SC  =>  SA % Stream_C ( iC ) % Element )
      select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
      class is ( CurrentSet_C_Form )
      call CSC % SetStream ( SC )
      end select !-- CSC
      end associate !-- SC
    end do !-- iC

  end subroutine SetStream


  impure elemental subroutine Finalize ( CSA )

    type ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA

  end subroutine Finalize


end module CurrentSet_A__Form
