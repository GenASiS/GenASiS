!-- Field_BSLL_ASC_CSLD is a template for a set of fields on a 
!   Bundle_SLL_ASC_CSLD.

module Field_BSLL_ASC_CSLD__Template

  use Basics
  use Atlases
  use BundleHeader_Form
  use FieldBundle_Template
  use Bundle_SLL_ASC_CSLD__Form

  !-- Field_BundleSingleLevelLocal_AtlasSingleChart_ChartSingleLevelDistributed
  !   _Template

  implicit none
  private

  type, public, extends ( FieldBundleTemplate ), abstract :: &
    Field_BSLL_ASC_CSLD_Template
      integer ( KDI ) :: &
        nFibers   = 0, &
        nSections = 0
      class ( Bundle_SLL_ASC_CSLD_Form ), pointer :: &
        Bundle_SLL_ASC_CSLD => null ( )
      class ( FieldAtlas_1D_Form ), allocatable :: &
        Fiber
      class ( FieldAtlas_1D_Form ), allocatable :: &
        Section
  contains
    procedure, public, pass :: &
      InitializeTemplate_BSLL_ASC_CSLD
    procedure, public, pass :: &
      LoadSections
    procedure, public, pass :: &
      StoreSections
    procedure, public, pass :: &
      LoadSection
    procedure, public, pass :: &
      StoreSection
    procedure, public, pass :: &
      FieldFiber
    procedure, public, pass :: &
      FieldFiber_CSL
    procedure, public, pass :: &
      FieldSection
    procedure, public, pass :: &
      FieldSection_CSL
    procedure, public, pass :: &
      FinalizeTemplate_BSLL_ASC_CSLD
  end type Field_BSLL_ASC_CSLD_Template

  type, public :: Field_BSLL_ASC_CSLD_Pointer
    class ( Field_BSLL_ASC_CSLD_Template ), pointer :: &
      Pointer => null ( )
  end type Field_BSLL_ASC_CSLD_Pointer

contains


  subroutine InitializeTemplate_BSLL_ASC_CSLD &
               ( FB, B, NameShort, UsePinnedMemory_S_Option, &
                 UsePinnedMemory_F_Option, IgnorabilityOption )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( in ), target :: &
      B
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemory_S_Option, &
      UsePinnedMemory_F_Option
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FB % Bundle_SLL_ASC_CSLD => B
    FB % nFibers   = B % nFibers
    FB % nSections = B % nSections

    call FB % InitializeTemplate &
           ( B, NameShort, UsePinnedMemory_S_Option, &
             UsePinnedMemory_F_Option, IgnorabilityOption )

  end subroutine InitializeTemplate_BSLL_ASC_CSLD


  subroutine LoadSections ( FB, GhostExchangeOption )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB
    logical ( KDL ), intent ( in ), optional :: &
      GhostExchangeOption

    integer ( KDI ) :: &
      iS

    do iS = 1, FB % nSections
      call FB % LoadSection ( iS, GhostExchangeOption )
    end do !-- iS

  end subroutine LoadSections


  subroutine StoreSections ( FB )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB

    integer ( KDI ) :: &
      iS

    do iS = 1, FB % nSections
      call FB % StoreSection ( iS )
    end do !-- iS

  end subroutine StoreSections


  subroutine LoadSection ( FB, iFC, GhostExchangeOption )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iFC  !-- iFiberCell
    logical ( KDL ), intent ( in ), optional :: &
      GhostExchangeOption

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iV     !-- iVariable
    logical ( KDL ) :: &
      GhostExchange
    class ( StorageForm ), pointer :: &
      FF, &
      FS

    GhostExchange = .true.
    if ( present ( GhostExchangeOption ) ) &
      GhostExchange = GhostExchangeOption

    associate ( B => FB % Bundle_SLL_ASC_CSLD )

    do iF = 1, FB % nFibers

      FF => FB % FieldFiber ( iF )
      FS => FB % FieldSection ( iFC )

      associate ( iBC => B % iaBaseCell ( iF ) )
      do iV = 1, FF % nVariables
        FS % Value ( iBC, iV ) = FF % Value ( iFC, iV )
      end do 
      end associate !-- iBC

    end do !-- iF
    
    if ( GhostExchange ) then

      select type ( FA => FB % Section % Atlas ( iFC ) % Element )
      class is ( Field_ASC_Template )
        
      select type ( FC => FA % Chart )
      class is ( Field_CSL_Template )   
        
      associate ( CB => B % Base_CSLD )
      call CB % ExchangeGhostData ( FC )
      end associate
        
      end select !-- FC
      end select !-- FA
      
    end if

    end associate !-- B
    nullify ( FF, FS )

  end subroutine LoadSection


  subroutine StoreSection ( FB, iFC )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iFC  !-- iFiberCell

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iV     !-- iVariable
    class ( StorageForm ), pointer :: &
      FF, &
      FS

    associate ( B => FB % Bundle_SLL_ASC_CSLD )

    do iF = 1, FB % nFibers

      FF => FB % FieldFiber ( iF )
      FS => FB % FieldSection ( iFC )

      associate ( iBC => B % iaBaseCell ( iF ) )
      do iV = 1, FS % nVariables
        FF % Value ( iFC, iV ) = FS % Value ( iBC, iV )
      end do 
      end associate !-- iBC

    end do !-- iF

    end associate !-- B
    nullify ( FF, FS )

  end subroutine StoreSection


  function FieldFiber ( FB, iFiber ) result ( FF )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( StorageForm ), pointer :: &
      FF
    
    class ( FieldChartTemplate ), pointer :: &
      FC 
    
    select type ( FA => FB % Fiber % Atlas ( iFiber ) % Element )
    class is ( Field_ASC_Template ) 
      FC => FA % Chart
      select type ( FC )
      class is ( Field_CSL_Template )   
        FF => FC % Field 
      end select !-- FC
    end select !-- FA

  end function FieldFiber


  function FieldFiber_CSL ( FB, iFiber ) result ( FF )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iFiber
    class ( Field_CSL_Template ), pointer :: &
      FF
    
    class ( FieldChartTemplate ), pointer :: &
      FC 
    
    select type ( FA => FB % Fiber % Atlas ( iFiber ) % Element )
    class is ( Field_ASC_Template ) 
      FC => FA % Chart
      select type ( FC )
      class is ( Field_CSL_Template )   
        FF => FC
      end select !-- FC
    end select !-- FA

  end function FieldFiber_CSL


  function FieldSection ( FB, iSection ) result ( FS )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iSection
    class ( StorageForm ), pointer :: &
      FS
    
    class ( FieldChartTemplate ), pointer :: &
      FC 

    select type ( FA => FB % Section % Atlas ( iSection ) % Element )
    class is ( Field_ASC_Template )
      FC => FA % Chart
      select type ( FC )
      class is ( Field_CSL_Template )   
        FS => FC % Field 
      end select !-- FC
    end select !-- FA

  end function FieldSection


  function FieldSection_CSL ( FB, iSection ) result ( FS )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( in ), target :: &
      FB
    integer ( KDI ), intent ( in ) :: &
      iSection
    class ( Field_CSL_Template ), pointer :: &
      FS
    
    class ( FieldChartTemplate ), pointer :: &
      FC 

    select type ( FA => FB % Section % Atlas ( iSection ) % Element )
    class is ( Field_ASC_Template )
      FC => FA % Chart
      select type ( FC )
      class is ( Field_CSL_Template )   
        FS => FC
      end select !-- FC
    end select !-- FA

  end function FieldSection_CSL


  impure elemental subroutine FinalizeTemplate_BSLL_ASC_CSLD ( FB )

    class ( Field_BSLL_ASC_CSLD_Template ), intent ( inout ) :: &
      FB

    nullify ( FB % Bundle_SLL_ASC_CSLD )

    if ( allocated ( FB % Section ) ) &
      deallocate ( FB % Section )
    if ( allocated ( FB % Fiber ) ) &
      deallocate ( FB % Fiber )

    call FB % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_BSLL_ASC_CSLD


end module Field_BSLL_ASC_CSLD__Template
