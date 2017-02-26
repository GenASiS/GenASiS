!-- Field_BSLL_ASC_CSLD is a template for a set of fields on a 
!   Bundle_SLL_ASC_CSLD.

module Field_BSLL_ASC_CSLD__Template

  use Basics
  use Atlases
  use BundleHeader_Form
  use FieldBundle_Template

  !-- Field_BundleSingleLevelLocal_AtlasSingleChart_ChartSingleLevelDistributed
  !   _Template

  implicit none
  private

  type, public, extends ( FieldBundleTemplate ), abstract :: &
    Field_BSLL_ASC_CSLD_Template
      class ( Field_ASC_1D_Form ), allocatable :: &
        Fiber
  contains
    procedure, public, pass :: &
      InitializeTemplate_BSLL_ASC_CSLD
    procedure, public, pass :: &
      FieldFiber
    procedure, public, pass :: &
      FinalizeTemplate_BSLL_ASC_CSLD
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type Field_BSLL_ASC_CSLD_Template

    abstract interface 
      subroutine SF ( FB, NameOutputOption )
        import Field_BSLL_ASC_CSLD_Template
        class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
          FB
        character ( * ), intent ( in ), optional :: &
          NameOutputOption
      end subroutine
    end interface

  type, public :: Field_BSLL_ASC_CSLD_Pointer
    class ( Field_BSLL_ASC_CSLD_Template ), pointer :: &
      Pointer => null ( )
  end type Field_BSLL_ASC_CSLD_Pointer

contains


  subroutine InitializeTemplate_BSLL_ASC_CSLD ( FB, B, NameOutputOption )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB
    class ( BundleHeaderForm ), intent ( in ), target :: &
      B
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    call FB % InitializeTemplate ( B, NameOutputOption = NameOutputOption )

    if ( FB % NameOutput == '' ) then
      call FB % SetField ( )
    else
      call FB % SetField ( NameOutputOption = FB % NameOutput )
    end if

  end subroutine InitializeTemplate_BSLL_ASC_CSLD


  function FieldFiber ( FB, iFiber ) result ( FF )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( VariableGroupForm ), pointer :: &
      FF
    
    class ( FieldChartTemplate ), pointer :: &
      FC 
    
    associate ( FA => FB % Fiber % Atlas ( iFiber ) % Element )
    
    FC => FA % Chart
    select type ( FC )
    class is ( Field_CSL_Template )   
      FF => FC % Field 
    end select !-- FC
    end associate !-- FA

  end function FieldFiber


  impure elemental subroutine FinalizeTemplate_BSLL_ASC_CSLD ( FB )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB

    if ( allocated ( FB % Fiber ) ) &
      deallocate ( FB % Fiber )

    call FB % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_BSLL_ASC_CSLD


end module Field_BSLL_ASC_CSLD__Template
