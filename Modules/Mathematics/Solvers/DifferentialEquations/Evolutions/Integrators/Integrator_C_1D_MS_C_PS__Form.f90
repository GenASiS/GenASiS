!-- Integrator_C_1D_MS_C_PS is a template for time evolution of multiple
!   conserved currents on momentum space (and position space), and a conserved
!   current on position space.

module Integrator_C_1D_MS_C_PS__Form

  !-- Integrator_Current_1D_MomentumSpace_Current_PositionSpace__Form

  use Basics
  use Fields
  use Integrator_C_1D_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_1D_C_PS_Template ) :: &
    Integrator_C_1D_MS_C_PS_Form
      type ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
        allocatable :: &
          Current_BSLL_ASC_CSLD_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_1D
    final :: &  !-- 1
      Finalize
    procedure, private, pass :: &  !-- 3
      ComputeTally_1D
    procedure, private, pass :: &
      Current_ASC_Pointer   
  end type Integrator_C_1D_MS_C_PS_Form

contains


  subroutine Initialize_1D &
               ( I, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_1D_MS_C_PS_Form ), intent ( inout ) :: &
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
      I % Type = 'an Integrator_C_1D_MS_C_PS'

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) then
      call Show ( 'Current_BSLL_ASC_CSLD_1D not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_1D_MS_C_PS__Form', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'Initialize', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_1D_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

  end subroutine Initialize_1D


  subroutine Finalize ( I )

    type ( Integrator_C_1D_MS_C_PS_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) &
      deallocate ( I % Current_BSLL_ASC_CSLD_1D )

    call I % FinalizeTemplate_C_1D_C_PS ( )

  end subroutine Finalize


  subroutine ComputeTally_1D ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_1D_MS_C_PS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC  !-- iCurrent

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) &
      return

    do iC = 1, I % N_CURRENTS_1D
      associate ( CB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
      call CB % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption  = IgnorabilityOption )
      end associate !-- CB
    end do !-- iC

  end subroutine ComputeTally_1D


  function Current_ASC_Pointer ( I, iC ) result ( CA )

    class ( Integrator_C_1D_MS_C_PS_Form ), intent ( inout ), target :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iC
    class ( Current_ASC_Template ), pointer :: &
      CA
 
    associate ( CB  =>  I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
    select type ( CAE => CB % Section % Atlas ( 1 ) % Element )
    class is ( Current_ASC_Template )
      CA => CAE
    end select !-- CAE
    end associate !-- CB

  end function Current_ASC_Pointer

  
end module Integrator_C_1D_MS_C_PS__Form
