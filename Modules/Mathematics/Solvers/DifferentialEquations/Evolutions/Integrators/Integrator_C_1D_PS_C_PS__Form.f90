!-- Integrator_C_1D_PS_C_PS is a template for time evolution of multiple
!   similar conserved currents on position space, and an additional conserved
!   current on position space.

module Integrator_C_1D_PS_C_PS__Form

  !-- Integrator_Current_1D_PositionSpace_Current_PositionSpace__Form

  use Basics
  use Fields
  use EvolutionBasics
  use TimeSeries_C_1D_C__Form
  use Integrator_C_1D_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_1D_C_PS_Template ) :: &
    Integrator_C_1D_PS_C_PS_Form
      type ( Current_ASC_ElementForm ), dimension ( : ), allocatable :: &
        Current_ASC_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_I
    final :: &  !-- 1
      Finalize
    procedure, private, pass :: &  !-- 3
      ComputeTally_1D
    procedure, private, pass :: &
      Current_ASC_Pointer   
  end type Integrator_C_1D_PS_C_PS_Form

contains


  subroutine Initialize_I &
               ( I, U, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_1D_PS_C_PS_Form ), intent ( inout ) :: &
      I
    class ( UniverseHeaderForm ), intent ( in ) :: &
      U
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
      I % Type = 'an Integrator_C_1D_PS_C_PS'

    if ( .not. allocated ( I % Current_ASC_1D ) ) then
      call Show ( 'Current_ASC_1D not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_1D_PS_C_PS__Form', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'Initialize', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_1D_C_PS &
           ( U, Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

    select type ( TS  =>  I % TimeSeries )
    type is ( TimeSeries_C_1D_C_Form )
      associate ( CA_1D  =>  I % Current_ASC_1D )
      call TS % Initialize ( I, CA_1D )
      end associate !-- CA_1D
    end select !-- TS

  end subroutine Initialize_I


  subroutine Finalize ( I )

    type ( Integrator_C_1D_PS_C_PS_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % Current_ASC_1D ) ) &
      deallocate ( I % Current_ASC_1D )

    call I % FinalizeTemplate_C_1D_C_PS ( )

  end subroutine Finalize


  subroutine ComputeTally_1D ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_1D_PS_C_PS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC  !-- iCurrent

    if ( .not. allocated ( I % Current_ASC_1D ) ) &
      return

    do iC = 1, I % N_CURRENTS_1D
      associate ( CA => I % Current_ASC_1D ( iC ) % Element )
      call CA % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption = IgnorabilityOption )
      end associate !-- CA
    end do !-- iC

  end subroutine ComputeTally_1D


  function Current_ASC_Pointer ( I, iC ) result ( CA )

    class ( Integrator_C_1D_PS_C_PS_Form ), intent ( inout ), target :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iC
    class ( Current_ASC_Template ), pointer :: &
      CA

    CA => I % Current_ASC_1D ( iC ) % Element

  end function Current_ASC_Pointer

  
end module Integrator_C_1D_PS_C_PS__Form
