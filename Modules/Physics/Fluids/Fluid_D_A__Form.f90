module Fluid_D_A__Form

  !-- Fluid_Dust_Atlas_Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_D_C__Form

  implicit none
  private

  type, public, extends ( CurrentSet_A_Form ) :: Fluid_D_A_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    final :: &
      Finalize
  end type Fluid_D_A_Form


contains


  subroutine InitializeAllocate_F &
               ( FA, GA, Units_F, NameOption, IgnorabilityOption )

    class ( Fluid_D_A_Form ), intent ( inout ) :: &
      FA
    class ( Geometry_F_A_Form ), intent ( in ) :: &
      GA
    class ( Units_F_Form ), intent ( in ) :: &
      Units_F
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( FA % Type  ==  '' ) &
      FA % Type  =  'a Fluid_D_A'

    Name  =  'Fluid'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    associate ( nC  =>  GA % Atlas % nCharts )

    if ( allocated ( FA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( FA % FieldSet_C ( nC ) )
    end if

    call FA % CurrentSet_A_Form % Initialize &
           ( GA, &
             NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Fluid_D_C_Form :: FA % FieldSet_C ( iC ) % Element ) 
        select type ( FC  =>  FA % FieldSet_C ( iC ) % Element )
          class is ( Fluid_D_C_Form )
        select type ( GC  =>  GA % FieldSet_C ( iC ) % Element )
          class is ( Geometry_F_C_Form )

        call FC % Initialize &
               ( GC, &
                 NameOption = NameOption, &
                 IgnorabilityOption = IgnorabilityOption )

        end select !-- GC
        end select !-- FC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_F


  impure elemental subroutine Finalize ( FA )

    type ( Fluid_D_A_Form ), intent ( inout ) :: &
      FA

  end subroutine Finalize


end module Fluid_D_A__Form
