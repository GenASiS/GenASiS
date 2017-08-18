module ProtoCurrent_ASC__Form

  use Basics
  use Manifolds
  use Current_ASC__Template
  use ProtoCurrent_Form
  use ProtoCurrentSources_CSL__Form
  use ProtoCurrentSources_ASC__Form
  use ProtoCurrent_CSL__Form

  implicit none
  private
  
  type, public, extends ( Current_ASC_Template ) :: ProtoCurrent_ASC_Form
    type ( ProtoCurrentSources_ASC_Form ), allocatable :: &
      Sources_ASC
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      ProtoCurrent_CSL
    generic, public :: &
      ProtoCurrent => ProtoCurrent_CSL
    final :: &
      Finalize
    procedure, public, pass :: &
      SetField
  end type ProtoCurrent_ASC_Form

!     integer ( KDI ), private, parameter :: &
!       iExtent   = 1!, &
!     !   iChartExcision = 2, &
!     !   nBoundaries    = 2
!     ! integer ( KDI ), dimension ( nBoundaries ), private, parameter :: &
!     !   iaB = [ BOUNDARY % CHART_EXTENT, BOUNDARY % CHART_EXCISION ]
!     ! character ( LDL ), dimension ( nBoundaries ), private, parameter :: &
!     !   BoundaryLabel = [ BOUNDARY % LABEL ( BOUNDARY % CHART_EXTENT ), &
!     !                     BOUNDARY % LABEL ( BOUNDARY % CHART_EXCISION ) ]

contains


  subroutine Initialize ( PCA, A, NameShortOption, IgnorabilityOption )

    class ( ProtoCurrent_ASC_Form ), intent ( inout ) :: &
      PCA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameShortOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

!     integer ( KDI ) :: &
!       iB  !-- iBoundary
    character ( LDL ) :: &
      NameShort

    PCA % Type = 'a ProtoCurrent_ASC'

    NameShort = 'Fluid'
    if ( present ( NameShortOption ) ) &
      NameShort = NameShortOption

!     select case ( trim ( FluidType ) )
!     case ( 'DUST' )
!       allocate ( DustTallyForm :: FA % TallyInterior )
!       allocate ( DustTallyForm :: FA % TallyTotal )
!       allocate ( DustTallyForm :: FA % TallyChange )
!       allocate ( FA % TallyBoundary ( A % nBoundaries ) )
!       do iB = 1, A % nBoundaries 
!         allocate ( DustTallyForm :: FA % TallyBoundary ( iB ) % Element )
!       end do !-- iB
!     case ( 'POLYTROPIC', 'MEAN_HEAVY_NUCLEUS' )
!       allocate ( PerfectFluidPolytropicTallyForm :: FA % TallyInterior )
!       allocate ( PerfectFluidPolytropicTallyForm :: FA % TallyTotal )
!       allocate ( PerfectFluidPolytropicTallyForm :: FA % TallyChange )
!       allocate ( FA % TallyBoundary ( A % nBoundaries ) )
!       do iB = 1, A % nBoundaries
!         allocate ( PerfectFluidPolytropicTallyForm &
!                      :: FA % TallyBoundary ( iB ) % Element )
!       end do !-- iB
!     end select !-- FluidType   

!     call FA % TallyInterior % Initialize ( A )
!     call FA % TallyTotal % Initialize ( A )
!     call FA % TallyChange % Initialize ( A )
!     do iB = 1, A % nBoundaries
!       call FA % TallyBoundary ( iB ) % Element % Initialize ( A )
!     end do !-- iB

    call PCA % InitializeTemplate_ASC_C ( A, NameShort, IgnorabilityOption )

    allocate ( PCA % Sources_ASC )
    associate ( PCSA => PCA % Sources_ASC )
    call PCSA % Initialize &
           ( PCA, NameShortOption = trim ( NameShort ) // '_Sources', &
             IgnorabilityOption = IgnorabilityOption )
    select type ( PCSC => PCSA % Chart )
    class is ( ProtoCurrentSources_CSL_Form )
      select type ( PC => PCA % Chart )
      class is ( ProtoCurrent_CSL_Form )
        call PC % SetSources ( PCSC )
      end select !-- PC
    end select !-- PCSC
    end associate !-- PCSA

  end subroutine Initialize


  function ProtoCurrent_CSL ( PCA ) result ( PC )

    class ( ProtoCurrent_ASC_Form ), intent ( in ) :: &
      PCA
    class ( ProtoCurrentForm ), pointer :: &
      PC

    select type ( PCC => PCA % Chart )
    class is ( ProtoCurrent_CSL_Form )
      PC => PCC % ProtoCurrent ( )
    class default
      call Show ( 'ProtoCurrent type not recognized', CONSOLE % ERROR )
      call Show ( 'ProtoCurrent_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ProtoCurrent_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- PCC

  end function ProtoCurrent_CSL


  impure elemental subroutine Finalize ( PCA )

    type ( ProtoCurrent_ASC_Form ), intent ( inout ) :: &
      PCA

    if ( allocated ( PCA % Sources_ASC ) ) &
      deallocate ( PCA % Sources_ASC )

    call PCA % FinalizeTemplate_ASC_C ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( ProtoCurrent_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( ProtoCurrent_CSL_Form :: FA % Chart )

    select type ( PCC => FA % Chart )
    class is ( ProtoCurrent_CSL_Form )
      call PCC % Initialize &
             ( C, FA % NameShort, nValues, &
               IgnorabilityOption = FA % Ignorability )
    end select !-- F

    call A % AddField ( FA )

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


end module ProtoCurrent_ASC__Form
