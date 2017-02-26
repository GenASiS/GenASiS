!-- Current_BSLL_ASC_CSLD is a template for a conserved current on a 
!   Bundle_SLL_ASC_CSLD.

module Current_BSLL_ASC_CSLD__Template

  !-- Current_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Template

  use Basics
  use Manifolds
  use Current_Template
  use Current_ASC__Template

  implicit none
  private

  type, public, extends ( Field_BSLL_ASC_CSLD_Template ), abstract :: &
    Current_BSLL_ASC_CSLD_Template
  contains
    procedure, public, pass :: &
      CurrentFiber
    procedure, public, pass :: &
      CurrentSection
  end type Current_BSLL_ASC_CSLD_Template

contains


  function CurrentFiber ( CB, iFiber ) result ( CF )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      CB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( CurrentTemplate ), pointer :: &
      CF

    class ( VariableGroupForm ), pointer :: &
      F

    CF => null ( )

    F => CB % FieldFiber ( iFiber )
    select type ( F )
    class is ( CurrentTemplate )
      CF => F 
    end select !-- F

    nullify ( F )

  end function CurrentFiber


  function CurrentSection ( CB, iSection ) result ( CS )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      CB
    integer ( KDI ), intent ( in ) :: &
      iSection
    class ( CurrentTemplate ), pointer :: &
      CS

    class ( VariableGroupForm ), pointer :: &
      F

    CS => null ( )

    F => CB % FieldSection ( iSection )
    select type ( F )
    class is ( CurrentTemplate )
      CS => F 
    end select !-- F

    nullify ( F )

  end function CurrentSection


end module Current_BSLL_ASC_CSLD__Template
