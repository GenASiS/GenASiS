module Geometry_F_A__Form

  !-- Geometry_Flat_Atlas_Form

  use Basics
  use Manifolds
  use FieldSets
  use Streams
  use Geometry_F_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: Geometry_F_A_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass ( GA ) :: &
      SetStream
    final :: &
      Finalize
  end type Geometry_F_A_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSA, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Geometry_F_A_Form ), intent ( inout ), target :: &
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
      FSA % Type  =  'a Geometry_F_A'

    Name  =  'Geometry'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    associate ( nC  =>  A % nCharts )

    if ( allocated ( FSA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( FSA % FieldSet_C ( nC ) )
    end if

    call FSA % FieldSet_A_Form % Initialize &
           ( A, &
             NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Geometry_F_C_Form :: FSA % FieldSet_C ( iC ) % Element ) 
        select type ( GC  =>  FSA % FieldSet_C ( iC ) % Element )
        class is ( Geometry_F_C_Form )

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


  subroutine SetStream ( SA, GA )

    class ( Stream_A_Form ), intent ( inout ) :: &
      SA
    class ( Geometry_F_A_Form ), intent ( in ) :: &
      GA

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, GA % Atlas % nCharts
      associate ( SC  =>  SA % Stream_C ( iC ) % Element )
      select type ( GC  =>  GA % FieldSet_C ( iC ) % Element )
      class is ( Geometry_F_C_Form )
      call GC % SetStream ( SC )
      end select !-- GC
      end associate !-- SC
    end do !-- iC

  end subroutine SetStream


  impure elemental subroutine Finalize ( GA )

    type ( Geometry_F_A_Form ), intent ( inout ) :: &
      GA

  end subroutine Finalize


end module Geometry_F_A__Form
