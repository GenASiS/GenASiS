module CurrentSet_A__Form

  !-- CurrentSet_Atlas_Form

  use Basics
  use FieldSets
  use Streams
  use Geometries
  use CurrentSet_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: CurrentSet_A_Form
    class ( Geometry_F_A_Form ), pointer :: &
      Geometry_A => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_CS
    generic, public :: &
      Initialize => InitializeAllocate_CS
    procedure, public, pass ( CSA ) :: &
      SetStream
    procedure, public, pass :: &
      ComputeFromInitial
    final :: &
      Finalize
  end type CurrentSet_A_Form


contains


  subroutine InitializeAllocate_CS &
               ( CSA, GA, FieldOption, VectorOption, NameOption, UnitOption, &
                 DensityUnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( CurrentSet_A_Form ), intent ( inout ), target :: &
      CSA
    class ( Geometry_F_A_Form ), intent ( in ), target :: &
      GA
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      DensityUnitOption
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

    CSA % Geometry_A  =>  GA

    associate ( nC  =>  GA % Atlas % nCharts )

    if ( allocated ( CSA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( CSA % FieldSet_C ( nC ) )
    end if

    call CSA % FieldSet_A_Form % Initialize &
           ( GA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( CurrentSet_C_Form :: CSA % FieldSet_C ( iC ) % Element ) 
        select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
          class is ( CurrentSet_C_Form )
        select type ( GC  =>  GA % FieldSet_C ( iC ) % Element )
          class is ( Geometry_F_C_Form )

        call CSC % Initialize &
               ( GC, FieldOption, VectorOption, NameOption, &
                 UnitOption, DensityUnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )

        end select !-- GC
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


  subroutine ComputeFromInitial ( CSA )

    class ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, CSA % Atlas % nCharts
      select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )
      call CSC % ComputeFromInitial ( )
      end select !-- CSC
    end do !-- iC

  end subroutine ComputeFromInitial


  impure elemental subroutine Finalize ( CSA )

    type ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA

    nullify ( CSA % Geometry_A )

  end subroutine Finalize


end module CurrentSet_A__Form
