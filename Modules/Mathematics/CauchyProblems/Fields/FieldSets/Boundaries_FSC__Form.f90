module Boundaries_FSC__Form

  !-- Boundaries_FieldSetChart_Form

  use Basics
  use Manifolds
  use Storage_FSC__Form

  implicit none
  private

  type, public :: Boundaries_FSC_Form
    integer ( KDI ) :: &
      nBoundaries
    character ( LDL ), dimension ( : ), allocatable :: &
      Boundary
    character ( LDL ), dimension ( :, : ), allocatable :: &
      BoundaryCondition
  contains
    procedure, private, pass :: &
      InitializeAllocate
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeAllocate, InitializeClone
    procedure, public, pass :: &
      SetFace
    procedure, public, pass :: &
      SetEdge
    procedure, public, pass :: &
      Show => Show_BFSC
    final :: &
      Finalize
  end type Boundaries_FSC_Form


contains


  subroutine InitializeAllocate ( BFSC, C, nExcisionsOption )

    class ( Boundaries_FSC_Form ), intent ( inout ) :: &
      BFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption

    BFSC % nBoundaries = 1
    if ( present ( nExcisionsOption ) ) &
      BFSC % nBoundaries = 1 + nExcisionsOption 

    associate &
      ( Cy => C % Connectivity )
    allocate &
      ( BFSC % BoundaryCondition ( Cy % nConnections, BFSC % nBoundaries ) )
    allocate &
      ( BFSC % Boundary ( BFSC % nBoundaries ) )

    BFSC % Boundary = ''
    BFSC % Boundary ( 1 ) = 'Extent' 

    BFSC % BoundaryCondition = ''
    BFSC % BoundaryCondition ( :, 1 ) = 'PERIODIC'

    end associate !-- Cy

  end subroutine InitializeAllocate


  subroutine InitializeClone ( BFSC, BFSC_S )

    class ( Boundaries_FSC_Form ), intent ( inout ) :: &
      BFSC
    class ( Boundaries_FSC_Form ), intent ( in ) :: &
      BFSC_S

      BFSC % nBoundaries  =  BFSC_S % nBoundaries

      allocate &
        ( BFSC % Boundary, source = BFSC_S % Boundary )
      allocate &
        ( BFSC % BoundaryCondition, source = BFSC_S % BoundaryCondition )

  end subroutine InitializeClone


  subroutine SetFace &
               ( BFSC, C, BoundaryCondition, iDimension, BoundaryOption, &
                 iBoundaryOption )

    class ( Boundaries_FSC_Form ), intent ( inout ) :: &
      BFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    character ( * ), dimension ( 2 ), intent ( in ) :: &
      BoundaryCondition  !-- [ Inner, Outer ]
    integer ( KDI ), intent ( in ) :: &
      iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    integer ( KDI ) :: &
      iB  !-- iBoundary

    if ( C % Connectivity % nFaces == 0 ) then
      call Show ( 'Faces not included in Connectivity', CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsFace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iDimension > C % nDimensions ) then
      call Show ( 'Selected iDimension > nDimensions', CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsFace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    iB = 1
    if ( present ( iBoundaryOption ) ) then
      if ( iBoundaryOption > BFSC % nBoundaries ) then
        call Show ( 'Selected iBoundary > nBoundaries', CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iBoundaryOption == 1 ) then
        if ( present ( BoundaryOption ) ) then
          call Show ( 'Boundary name not allowed for iBoundary == 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      else
        if ( .not. present ( BoundaryOption ) ) then
          call Show ( 'Boundary name required for iBoundary > 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      end if
      iB = iBoundaryOption
    end if

    if ( present ( BoundaryOption ) ) then
      if ( .not.present ( iBoundaryOption ) ) then
        call Show ( 'Argument iBoundary required when Boundary name present', &
                    CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )         
      else
        BFSC % Boundary ( iBoundaryOption ) = BoundaryOption
      end if
    end if

    associate &
      ( Cy  => C % Connectivity, &
        iD => iDimension )
    BFSC % BoundaryCondition ( Cy % iaInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    BFSC % BoundaryCondition ( Cy % iaOuter ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    end associate !-- Cy, etc.

  end subroutine SetFace


  subroutine SetEdge &
               ( BFSC, C, BoundaryCondition, iDimension, BoundaryOption, &
                 iBoundaryOption )

    class ( Boundaries_FSC_Form ), intent ( inout ) :: &
      BFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    character ( * ), dimension ( 4 ), intent ( in ) :: &
      BoundaryCondition  !-- [ InnerInner, OuterInner, InnerOuter, OuterOuter ]
    integer ( KDI ), intent ( in ) :: &
      iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- jDimension, etc.
      iB  !-- iBoundary

    iD = iDimension
    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1
 
    if ( C % Connectivity % nEdges == 0 ) then
      call Show ( 'Edges not included in Connectivity', CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsEdge', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( jD > C % nDimensions .or. kD > C % nDimensions ) then
      call Show ( 'Selected jDimension or kDimension > nDimensions', &
                  CONSOLE % ERROR )
      call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsEdge', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    iB = 1
    if ( present ( iBoundaryOption ) ) then
      if ( iBoundaryOption > BFSC % nBoundaries ) then
        call Show ( 'Selected iBoundary > nBoundaries', CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iBoundaryOption == 1 ) then
        if ( present ( BoundaryOption ) ) then
          call Show ( 'Boundary name not allowed for iBoundary == 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      else
        if ( .not. present ( BoundaryOption ) ) then
          call Show ( 'Boundary name required for iBoundary > 1', &
                      CONSOLE % ERROR )
          call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      end if
      iB = iBoundaryOption
    end if

    if ( present ( BoundaryOption ) ) then
      if ( .not.present ( iBoundaryOption ) ) then
        call Show ( 'Argument iBoundary required when BoundaryName present', &
                    CONSOLE % ERROR )
        call Show ( 'FieldSet_C_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )         
      else
        BFSC % Boundary ( iBoundaryOption ) = BoundaryOption
      end if
    end if

    associate ( Cy => C % Connectivity )
    BFSC % BoundaryCondition ( Cy % iaInnerInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    BFSC % BoundaryCondition ( Cy % iaOuterInner ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    BFSC % BoundaryCondition ( Cy % iaInnerOuter ( iD ), iB ) &
      = BoundaryCondition ( 3 )
    BFSC % BoundaryCondition ( Cy % iaOuterOuter ( iD ), iB ) &
      = BoundaryCondition ( 4 )
    end associate !-- Cy

  end subroutine SetEdge


  subroutine Show_BFSC ( BFSC, C, Ignorability )

    class ( Boundaries_FSC_Form ), intent ( in ) :: &
      BFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      Ignorability

    integer ( KDI ) :: &
      iB, &  !-- iBoundary
      iD, jD, kD  !-- iDimension, etc.

    associate &
      ( Cy  =>  C % Connectivity, &
        nD  =>  C % nDimensions, &
        BC  =>  BFSC % BoundaryCondition ( :, : ), &
        BN  =>  BFSC % Boundary ( : ) )

!    call Show ( 'Boundary conditions', Ignorability )
    call Show ( BFSC % nBoundaries, 'nBoundaries', Ignorability )

    do iB = 1, BFSC % nBoundaries
      call Show ( BN ( iB ), 'Boundary', Ignorability )
      call Show ( iB, 'iBoundary', Ignorability )
  
      if ( Cy % nFaces > 0 ) then
          do iD = 1, nD
            call Show ( iD, 'Faces, iDimension', Ignorability )
            associate &
              ( iaI => Cy % iaInner ( iD ), &
                iaO => Cy % iaOuter ( iD ) )
            call Show ( [ BC ( iaI, iB ), BC ( iaO, iB ) ], &
                        '[ Inner, Outer ]', Ignorability )
            end associate !-- iaI, etc.
          end do !-- iD
      end if

      if ( Cy % nEdges > 0 ) then
          do iD = 1, nD
            jD = mod ( iD, 3 ) + 1
            kD = mod ( jD, 3 ) + 1
            if ( jD > nD .or. kD > nD ) &
              cycle
            call Show ( iD, 'Edges parallel to iDimension', Ignorability )
            associate &
              ( iaII => Cy % iaInnerInner ( iD ), &
                iaOI => Cy % iaOuterInner ( iD ), &
                iaIO => Cy % iaInnerOuter ( iD ), &
                iaOO => Cy % iaOuterOuter ( iD ) )
            call Show ( [ BC ( iaII, iB ), BC ( iaOI, iB ), &
                          BC ( iaIO, iB ), BC ( iaOO, iB ) ], &
                        '[ InnerInner, OuterInner, InnerOuter, OuterOuter ]', &
                        Ignorability )
            end associate !-- iaII, etc.
          end do !-- iD
      end if

    end do !-- iB

    end associate !-- Cy, etc.

  end subroutine Show_BFSC


  impure elemental subroutine Finalize ( BFSC )

    type ( Boundaries_FSC_Form ), intent ( inout ) :: &
      BFSC

    if ( allocated ( BFSC % BoundaryCondition ) ) &
      deallocate ( BFSC % BoundaryCondition )
    if ( allocated ( BFSC % Boundary ) ) &
      deallocate ( BFSC % Boundary )

  end subroutine Finalize


end module Boundaries_FSC__Form
