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
      integer ( KDI ) :: &
        nSections
      class ( Bundle_SLL_ASC_CSLD_Form ), pointer :: &
        Bundle_SLL_ASC_CSLD => null ( )
      type ( Current_ASC_ElementForm ), dimension ( : ), allocatable :: &
        Section_ASC
  contains
    procedure, public, pass :: &
      InitializeTemplate_BSLL_ASC_CSLD_C
    procedure, public, pass :: &
      LoadSections
    procedure, public, pass :: &
      StoreSections
    procedure, public, pass :: &
      FinalizeTemplate_BSLL_ASC_CSLD_C
  end type Current_BSLL_ASC_CSLD_Template

contains


  subroutine InitializeTemplate_BSLL_ASC_CSLD_C ( CB, B, NameOutputOption )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      CB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    call CB % InitializeTemplate_BSLL_ASC_CSLD &
           ( B, NameOutputOption = NameOutputOption )

    CB % nSections = B % nFibers

    CB % Bundle_SLL_ASC_CSLD => B

  end subroutine InitializeTemplate_BSLL_ASC_CSLD_C


  subroutine LoadSections ( CB )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      CB

    integer ( KDI ) :: &
      iS
    class ( CurrentTemplate ), pointer :: &
      C

    associate ( B => CB % Bundle_SLL_ASC_CSLD )

    do iS = 1, CB % nSections
      associate ( CA => CB % Section_ASC ( iS ) % Element )
      C => CA % Current ( )
      call B % LoadSection ( C, CB, iS, GhostExchangeOption = .true. )
      end associate !-- CA
    end do !-- iS

    end associate !-- B
    nullify ( C )

  end subroutine LoadSections


  subroutine StoreSections ( CB )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      CB

    integer ( KDI ) :: &
      iS
    class ( CurrentTemplate ), pointer :: &
      C

    associate ( B => CB % Bundle_SLL_ASC_CSLD )

    do iS = 1, CB % nSections
      associate ( CA => CB % Section_ASC ( iS ) % Element )
      C => CA % Current ( )
      call B % StoreSection ( CB, C, iS )
      end associate !-- CA
    end do !-- iS

    end associate !-- B
    nullify ( C )

  end subroutine StoreSections


  impure elemental subroutine FinalizeTemplate_BSLL_ASC_CSLD_C ( CB )

    class ( Current_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      CB

    if ( allocated ( CB % Section_ASC ) ) &
      deallocate ( CB % Section_ASC )

    call CB % FinalizeTemplate_BSLL_ASC_CSLD ( )

  end subroutine FinalizeTemplate_BSLL_ASC_CSLD_C


end module Current_BSLL_ASC_CSLD__Template
