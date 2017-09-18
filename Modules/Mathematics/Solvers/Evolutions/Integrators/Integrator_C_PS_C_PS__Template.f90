!-- Integrator_C_PS_C_PS is a template for operator-split time evolution of two
!   conserved currents on position space.

module Integrator_C_PS_C_PS__Template

  !-- Integrator_Current_PositionSpace_Current_PositionSpace__Template

  use Basics
  use Fields
  use Integrator_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    Integrator_C_PS_C_PS_Template
      type ( Current_ASC_ElementForm ), allocatable :: &
        Current_ASC_1, &
        Current_ASC_2
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C_PS_C_PS
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C_PS_C_PS
  end type Integrator_C_PS_C_PS_Template

contains


  subroutine InitializeTemplate_C_PS_C_PS &
               ( I, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_C_PS_C_PS'

    if ( .not. allocated ( I % Current_ASC_1 ) ) then
      call Show ( 'Current_ASC_1 not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    if ( .not. allocated ( I % Current_ASC_2 ) ) then
      call Show ( 'Current_ASC_2 not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

  end subroutine InitializeTemplate_C_PS_C_PS



  subroutine FinalizeTemplate_C_PS_C_PS ( I )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % TimeSeries ) ) &
     deallocate ( I % TimeSeries )
   if ( allocated ( I % Current_ASC_2 ) ) &
     deallocate ( I % Current_ASC_2 )
   if ( allocated ( I % Current_ASC_1 ) ) &
     deallocate ( I % Current_ASC_1 )

    call I % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C_PS_C_PS


end module Integrator_C_PS_C_PS__Template
