module Boundaries_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: BoundariesForm
    integer ( KDI ) :: &
      nBoundaries = 0
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
      Show => Show_B
    procedure, private, pass ( B ):: &
      ApplyAll
    procedure, private, pass ( B ):: &
      ApplyDimensionConnectivity
    generic, public :: &
      Apply => ApplyAll, ApplyDimensionConnectivity
    final :: &
      Finalize
  end type BoundariesForm

    private :: &
      CopyField, &
      ReverseField

      private :: &
        SetLimits, &
        CopyFieldKernel, &
        ReverseFieldKernel
      
    interface

      module subroutine CopyFieldKernel &
               ( V, nB, dBE, dBI, oBE, oBI, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          V
        integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
          nB,  & 
          dBE, &
          dBI, &
          oBE, &
          oBI
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine CopyFieldKernel
      
      module subroutine ReverseFieldKernel &
               ( V, nB, dBE, oBE, UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          V
        integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
          nB,  & 
          dBE, &
          oBE
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ReverseFieldKernel
    
    end interface


contains


  subroutine InitializeAllocate ( B, C, nExcisionsOption )

    class ( BoundariesForm ), intent ( inout ) :: &
      B
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption

    B % nBoundaries = 1
    if ( present ( nExcisionsOption ) ) &
      B % nBoundaries = 1 + nExcisionsOption 

    associate &
      ( Cy => C % Connectivity )
    allocate &
      ( B % BoundaryCondition ( Cy % nConnections, B % nBoundaries ) )
    allocate &
      ( B % Boundary ( B % nBoundaries ) )

    B % Boundary = ''
    B % Boundary ( 1 ) = 'Extent' 

    B % BoundaryCondition = ''

    end associate !-- Cy

  end subroutine InitializeAllocate


  subroutine InitializeClone ( B, B_S )

    class ( BoundariesForm ), intent ( inout ) :: &
      B
    class ( BoundariesForm ), intent ( in ) :: &
      B_S

      B % nBoundaries  =  B_S % nBoundaries

      allocate &
        ( B % Boundary, source = B_S % Boundary )
      allocate &
        ( B % BoundaryCondition, source = B_S % BoundaryCondition )

  end subroutine InitializeClone


  subroutine SetFace &
               ( B, C, BoundaryCondition, iDimension, BoundaryOption, &
                 iBoundaryOption )

    class ( BoundariesForm ), intent ( inout ) :: &
      B
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
      if ( iBoundaryOption > B % nBoundaries ) then
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
        B % Boundary ( iBoundaryOption ) = BoundaryOption
      end if
    end if

    associate &
      ( Cy  => C % Connectivity, &
        iD => iDimension )
    B % BoundaryCondition ( Cy % iaInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    B % BoundaryCondition ( Cy % iaOuter ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    end associate !-- Cy, etc.

  end subroutine SetFace


  subroutine SetEdge &
               ( B, C, BoundaryCondition, iDimension, BoundaryOption, &
                 iBoundaryOption )

    class ( BoundariesForm ), intent ( inout ) :: &
      B
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
      if ( iBoundaryOption > B % nBoundaries ) then
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
        B % Boundary ( iBoundaryOption ) = BoundaryOption
      end if
    end if

    associate ( Cy => C % Connectivity )
    B % BoundaryCondition ( Cy % iaInnerInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    B % BoundaryCondition ( Cy % iaOuterInner ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    B % BoundaryCondition ( Cy % iaInnerOuter ( iD ), iB ) &
      = BoundaryCondition ( 3 )
    B % BoundaryCondition ( Cy % iaOuterOuter ( iD ), iB ) &
      = BoundaryCondition ( 4 )
    end associate !-- Cy

  end subroutine SetEdge


  subroutine Show_B ( B, C, Ignorability )

    class ( BoundariesForm ), intent ( in ) :: &
      B
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
        BC  =>  B % BoundaryCondition ( :, : ), &
        BN  =>  B % Boundary ( : ) )

    if ( trim ( BC ( 1, 1 ) )  ==  '' ) &
      return

!    call Show ( 'Boundary conditions', Ignorability )
    call Show ( B % nBoundaries, 'nBoundaries', Ignorability )

    do iB = 1, B % nBoundaries
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

  end subroutine Show_B


  subroutine ApplyAll ( S, B, C )

    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( BoundariesForm ), intent ( in ) :: &
      B
    class ( Chart_H_Form ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iD  !-- iDimension

    do iD  =  1,  C % nDimensions
      associate &
        ( iC_I  =>  C % Connectivity % iaInner ( iD ), &
          iC_O  =>  C % Connectivity % iaOuter ( iD ) )
      call B % Apply ( S, C, iD, iC_I )
      call B % Apply ( S, C, iD, iC_O )
      end associate !-- iC_I, etc.
    end do

  end subroutine ApplyAll


  subroutine ApplyDimensionConnectivity ( S, B, C, iD, iC )

    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( BoundariesForm ), intent ( in ) :: &
      B
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iC     !-- iConnection

    integer ( KDI ) :: &
      iB, &  !-- iBoundary
      iS, &  !-- iSelected
      iF, &  !-- iField
      iV     !-- iVector

    do iB  =  1,  B % nBoundaries

      select type ( C )
      class is ( Chart_GS_Form )

        associate ( BC  =>  B % BoundaryCondition ( iC, iB ) )
        select case ( trim ( BC ) )
        case ( 'PERIODIC', 'INFLOW' )

          cycle

        case ( 'OUTFLOW' )

          do iS  =  1,  S % nVariables
            iF  =  S % iaSelected ( iS )
            call CopyField ( S, C, iF, iD, iC )
          end do !-- iS

        case ( 'REFLECTING' )

          do iS  =  1,  S % nVariables
            iF  =  S % iaSelected ( iS )
            call CopyField ( S, C, iF, iD, iC )
            do iV = 1, S % nVectors
              if ( iF == S % VectorIndices ( iV ) % Value ( iD ) ) &
                call ReverseField ( S, C, iF, iD, iC )
            end do !-- iV
          end do !-- iS

        case default
          call Show ( 'BoundaryCondition not recognized', CONSOLE % ERROR )
          call Show ( BC, 'BoundaryCondition', CONSOLE % ERROR )
          call Show ( 'Boundaries_Form', 'module', CONSOLE % ERROR )
          call Show ( 'Apply', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- BC
        end associate !-- BC

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Boundaries_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Apply', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C

    end do !-- iB

  end subroutine ApplyDimensionConnectivity


  impure elemental subroutine Finalize ( B )

    type ( BoundariesForm ), intent ( inout ) :: &
      B

    if ( allocated ( B % BoundaryCondition ) ) &
      deallocate ( B % BoundaryCondition )
    if ( allocated ( B % Boundary ) ) &
      deallocate ( B % Boundary )

  end subroutine Finalize


  subroutine CopyField ( S, C, iF, iD, iC )

    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iF, &  !-- iField
      iD, &  !-- iDimension
      iC     !-- iConnection

    integer ( KDI ), dimension ( 3 ) :: &
      oBI, &  !-- oBoundaryInterior
      oBE, &  !-- oBoundaryExterior
      dBI, &  !-- dBoundaryInterior, i.e. direction
      dBE, &  !-- dBoundaryExterior, i.e. direction
      nB      !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F

    associate ( Cy => C % Connectivity )
    if ( iC == Cy % iaInner ( iD ) &
         .and. C % iaBrick ( iD ) /= 1 ) return
    if ( iC == Cy % iaOuter ( iD ) &
         .and. C % iaBrick ( iD ) /= C % nBricks ( iD ) ) return
    end associate !-- Cy

    associate ( nCB  =>  C % nCellsBrick )
    call SetLimits &
           ( C, nCB, iD, iC, nB, dBE, dBI, oBE, oBI )
    end associate !-- nCB

    call C % SetFieldPointer ( S % Value ( :, iF ), F )

    call CopyFieldKernel &
           ( F, nB, dBE, dBI, oBE, oBI, &
             UseDeviceOption = S % AllocatedDevice )

    nullify ( F )

  end subroutine CopyField


  subroutine ReverseField ( S, C, iF, iD, iC )

    class ( StorageForm ), intent ( inout ) :: &
      S
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iF, &  !-- iField
      iD, &  !-- iDimension
      iC     !-- iConnection

    integer ( KDI ), dimension ( 3 ) :: &
      oBI, &  !-- oBoundaryInterior
      oBE, &  !-- oBoundaryExterior
      dBI, &  !-- dBoundaryInterior, i.e. direction
      dBE, &  !-- dBoundaryExterior, i.e. direction
      nB      !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F

    associate ( Cy => C % Connectivity )
    if ( iC == Cy % iaInner ( iD ) &
         .and. C % iaBrick ( iD ) /= 1 ) return
    if ( iC == Cy % iaOuter ( iD ) &
         .and. C % iaBrick ( iD ) /= C % nBricks ( iD ) ) return
    end associate !-- Cy

    associate ( nCB  =>  C % nCellsBrick )
    call SetLimits &
           ( C, nCB, iD, iC, nB, dBE, dBI, oBE, oBI )
    end associate !-- nCB

    call C % SetFieldPointer ( S % Value ( :, iF ), F )

    call ReverseFieldKernel &
           ( F, nB, dBE, oBE, &
             UseDeviceOption = S % AllocatedDevice )

    nullify ( F )

  end subroutine ReverseField


  subroutine SetLimits &
               ( C, nCells, iDimension, iConnection, nB, dBE, dBI, oBE, oBI )

    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection
    integer ( KDI ), dimension ( 3 ), intent ( out ) :: &
      nB, &   !-- nBoundary
      dBI, &  !-- dBoundaryInterior, i.e. direction
      dBE, &  !-- dBoundaryExterior, i.e. direction
      oBI, &  !-- oBoundaryInterior
      oBE     !-- oBoundaryExterior

    integer ( KDI ) :: &
      jD, kD   !-- jDimension, kDimension

    associate &
      ( Cy => C % Connectivity, &
        iD => iDimension, &
        iC => iConnection )

    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1

    !-- In setting oBI and oBE, note kernel routine does not inherit lbound

    oBI = C % nGhostLayers
    dBI = +1
    if ( iC == Cy % iaOuter ( iD ) ) then
      oBI ( iD ) = oBI ( iD ) + nCells ( iD ) + 1
      dBI ( iD ) = -1
    end if !-- iC

    oBE = oBI
    dBE = dBI
    dBE ( iD ) = -dBI ( iD )
    if ( iC == Cy % iaInner ( iD ) ) then
      oBE ( iD ) = oBE ( iD ) + 1
    else if ( iC == Cy % iaOuter ( iD ) ) then
      oBE ( iD ) = oBE ( iD ) - 1
    end if !-- iC

    nB ( iD ) = C % nGhostLayers ( iD )
    nB ( jD ) = nCells ( jD )
    nB ( kD ) = nCells ( kD )

    end associate !-- Connectivity, etc.
    
  end subroutine SetLimits

  
end module Boundaries_Form
