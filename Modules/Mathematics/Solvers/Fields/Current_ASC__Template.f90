!-- Current_ASC is a template for a conserved current on an Atlas_SC.

module Current_ASC__Template

  !-- Current_AtlasSingleChart_Template

  use Basics
  use Manifolds
  use Current_Template
  use Tally_C__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ), abstract :: &
    Current_ASC_Template
      class ( Tally_C_Form ), allocatable :: &
        TallyInterior, &
        TallyTotal, &
        TallyChange
      class ( Tally_C_Element_Form ), dimension ( : ), allocatable :: &
        TallyBoundaryLocal, &
        TallyBoundaryGlobal
!         TallySources
  contains
    procedure, public, pass :: &
      InitializeTemplate_ASC_C
    procedure, private, pass :: &
      Current_CSL
    generic, public :: &
      Current => Current_CSL
    procedure, private, pass :: &
      AccumulateBoundaryTally_CSL
    generic, public :: &
      AccumulateBoundaryTally &
        => AccumulateBoundaryTally_CSL
    procedure, public, pass :: &
      ComputeTally
    procedure, public, pass :: &
      FinalizeTemplate_ASC_C
  end type Current_ASC_Template

contains


  subroutine InitializeTemplate_ASC_C ( CA, A, NameOutputOption )

    class ( Current_ASC_Template ), intent ( inout ) :: &
      CA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    integer ( KDI ) :: &
      iB  !-- iBoundary
    class ( CurrentTemplate ), pointer :: &
      C

    call CA % InitializeTemplate_ASC &
           ( A, NameOutputOption = NameOutputOption )

    if ( .not. allocated ( CA % TallyInterior ) ) then

      allocate ( CA % TallyInterior )
      allocate ( CA % TallyTotal )
      allocate ( CA % TallyChange )
      allocate ( CA % TallyBoundaryLocal  ( A % nBoundaries ) )
      allocate ( CA % TallyBoundaryGlobal ( A % nBoundaries ) )
      do iB = 1, A % nBoundaries 
        allocate ( CA % TallyBoundaryLocal  ( iB ) % Element )
        allocate ( CA % TallyBoundaryGlobal ( iB ) % Element )
      end do !-- iB

      C => CA % Current ( )
      call CA % TallyInterior % Initialize ( C, A )
      call CA % TallyTotal % Initialize ( C, A )
      call CA % TallyChange % Initialize ( C, A )
      do iB = 1, A % nBoundaries
        call CA % TallyBoundaryLocal  ( iB ) % Element % Initialize ( C, A )
        call CA % TallyBoundaryGlobal ( iB ) % Element % Initialize ( C, A )
      end do !-- iB

    end if

    nullify ( C )

  end subroutine InitializeTemplate_ASC_C


  function Current_CSL ( CA ) result ( C )

    class ( Current_ASC_Template ), intent ( in ), target :: &
      CA
    class ( CurrentTemplate ), pointer :: &
      C
    
    class ( VariableGroupForm ), pointer :: &
      F
    class ( FieldChartTemplate ), pointer :: &
      CC 
    
    CC => CA % Chart
    select type ( CC )
    class is ( Field_CSL_Template )
      F => CC % Field
      select type ( F )
      class is ( CurrentTemplate )
        C => F
      end select
    class default
      call Show ( 'Current type not recognized', CONSOLE % ERROR )
      call Show ( 'Current_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Current_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CC

  end function Current_CSL


  subroutine AccumulateBoundaryTally_CSL ( CA, BoundaryFluence_CSL )

    class ( Current_ASC_Template ), intent ( inout ) :: &
      CA
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence_CSL  !-- boundary slab

    select type ( FC => CA % Chart )
    class is ( Field_CSL_Template )

    associate ( iExtent => 1 )  !-- only boundary for CSL
    call CA % TallyBoundaryLocal ( iExtent ) % Element &
         % ComputeBoundary ( FC, BoundaryFluence_CSL )
    end associate !-- iExtent
      
    end select !-- FC
    
  end subroutine AccumulateBoundaryTally_CSL


  subroutine ComputeTally ( CA, ComputeChangeOption )
    
    class ( Current_ASC_Template ), intent ( inout ) :: &
      CA
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    
    integer ( KDI ) :: &
      iB  !-- iBoundary
    real ( KDR ), dimension ( CA % TallyInterior % N_INTEGRALS ) :: &
      OldTotal
    logical ( KDL ) :: &
      ComputeChange
    type ( CollectiveOperation_R_Form ) :: &
      CO

    ComputeChange = .true.
    if ( present ( ComputeChangeOption ) ) ComputeChange = ComputeChangeOption

    OldTotal = CA % TallyTotal % Value

    select type ( A => CA % Atlas )
    class is ( Atlas_SC_Form )

    associate ( nI => CA % TallyInterior % N_INTEGRALS )

    !-- Interior

    select type ( FC => CA % Chart )
    class is ( Field_CSL_Template )
      call CA % TallyInterior % ComputeInterior ( FC )
    end select !-- FC

    !-- Boundaries: sum local accumulations

    call CO % Initialize &
           ( A % Communicator, &
             nOutgoing = [ A % nBoundaries * nI ], &
             nIncoming = [ A % nBoundaries * nI ] )
    do iB = 1, A % nBoundaries
      CO % Outgoing % Value ( ( iB - 1 ) * nI + 1  :  iB * nI ) &
        = CA % TallyBoundaryLocal ( iB ) % Element % Value
    end do !-- iB
    call CO % Reduce ( REDUCTION % SUM )

    !-- Total

    CA % TallyTotal % Value = CA % TallyInterior % Value

    do iB = 1, A % nBoundaries
      CA % TallyBoundaryGlobal ( iB ) % Element % Value &
        = CO % Incoming % Value ( ( iB - 1 ) * nI + 1  :  iB * nI )
      CA % TallyTotal % Value &
        = CA % TallyTotal % Value  &
          + CA % TallyBoundaryGlobal ( iB ) % Element % Value
    end do !-- iB

    !-- Change

    if ( ComputeChange ) then
      CA % TallyChange % Value &
        = CA % TallyChange % Value &
           + ( CA % TallyTotal % Value - OldTotal )
    end if
  
    !-- Display

    call CA % TallyInterior % Show ( 'Interior Fluid Tally' )

    BoundaryLoop: do iB = 1, A % nBoundaries
      call CA % TallyBoundaryGlobal ( iB ) % Element % Show &
             ( 'Boundary ' // trim ( A % BoundaryName ( iB ) ) &
             // ' Fluid Tally', CONSOLE % INFO_2 )
    end do BoundaryLoop

    call CA % TallyTotal % Show ( 'Total Fluid Tally' )

    if ( ComputeChange ) then
      call CA % TallyChange % Show &
             ( 'Change in Total Fluid Tally' )
    end if

    end associate !-- nI

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Current_ASC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeTally', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine ComputeTally
  
  
  subroutine FinalizeTemplate_ASC_C ( CA )

    class ( Current_ASC_Template ), intent ( inout ) :: &
      CA

    call CA % FinalizeTemplate_ASC ( )

  end subroutine FinalizeTemplate_ASC_C


end module Current_ASC__Template
