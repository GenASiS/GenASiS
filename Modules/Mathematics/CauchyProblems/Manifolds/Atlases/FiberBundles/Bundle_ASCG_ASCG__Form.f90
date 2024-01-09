module Bundle_ASCG_ASCG__Form

  !-- Bundle_AtlasSingleChartGrid_AtlasSingleChartGrid__Form

  use Basics
  use Charts
  use BaseManifolds
  use Bundle_H__Form

  implicit none
  private

  type, public, extends ( Bundle_H_Form ) :: Bundle_ASCG_ASCG_Form
    !-- Base manifold ( Distributed on each of several copies )
    class ( Atlas_SCG_Form ), pointer :: &
      Atlas_SCG_Base  => null ( )
    class ( Chart_GS_Form ), pointer :: &
      Chart_GS_Base  => null ( )
    !-- Typical fiber ( Local )
    class ( Atlas_SCG_Form ), pointer :: &
      Atlas_SCG_Fiber => null ( )
    class ( Chart_GS_Form ), pointer :: &
      Chart_GS_Fiber => null ( )
    !-- Members for transposes between:
    !     * Fiber-centric storage ( for operations on fibers ), and 
    !     * Section-centric storage ( for operations on a copy of the 
    !                                 base manifold )
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( PortalHeaderForm ), allocatable :: &
      Portal_F_S, &  !-- Fiber-centric to Section-centric
      Portal_S_F     !-- Section-centric to Fiber-centric
    !-- Metadata for Fiber to Section decomposition
    integer ( KDI ) :: &
      nProcesses, &
      nCopiesBase, &
      nProcessesBase
    !-- Metadata for distribution of fiber operations for different 
    !   base manifold cells
    !-- ( An operation on a fiber for a given base manifold cell is itself 
    !     local. )
    integer ( KDI ) :: &
      nFibers, &
      nMyFibers
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaBrickFirst, iaBrickLast, &
      iaCellFirst, iaCellLast, &
      nFibersGlobal
    integer ( KDI ), dimension ( : ), allocatable :: &
      iFiberExchangeFirst, iFiberExchangeLast
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      iaBinExchangeFirst, iaBinExchangeLast
    !-- Metadata for distribution of position space operations for different 
    !   fiber bins
    !-- ( A position space operation for a given fiber bin is itself 
    !     distributed on a particular copy of the base manifold. )
    integer ( KDI ) :: &
      nSections, &
      nMySections, &
      iCopyBase, &
      iFiberFirst, iFiberLast
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaBinFirst, iaBinLast, &
      nSectionsGlobal
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      iaCellExchangeFirst, iaCellExchangeLast
  contains
    procedure, private, pass :: &
      Initialize_ASCG_ASCG
    generic, public :: &
      Initialize => Initialize_ASCG_ASCG
    procedure, private, pass :: &
      Show_B
    final :: &
      Finalize
  end type Bundle_ASCG_ASCG_Form

    private :: &
      SetDecomposition

      private :: &
        SetPortals

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = 3


contains


  subroutine Initialize_ASCG_ASCG &
               ( B, Base, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, nBricksOption, &
                 nDimensionsOption )

    class ( Bundle_ASCG_ASCG_Form ), intent ( inout ), target :: &
      B
    class ( Atlas_SCG_Form ), intent ( in ), target :: &
      Base
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption, &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption, &
      nBricksOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption

    if ( B % Type  ==  '' ) &
      B % Type  =  'a Bundle_ASCG_ASCG'

    call B % Initialize_H ( Base, NameOption )

    !-- Base

    B % Base            =>  Base
    B % Atlas_SCG_Base  =>  Base
    B % Chart_GS_Base   =>  Base % Chart_GS

    associate ( CB  =>  B % Chart_GS_Base )
    B % nFibers  =  product ( CB % nCells )
    end associate !-- CB

    !-- Fibers

    allocate ( Atlas_SCG_Form :: B % Fiber )
    select type ( AF  =>  B % Fiber )
      class is ( Atlas_SCG_Form )

    call AF % Initialize &
           ( SpacingOption = SpacingOption, &
             CoordinateLabelOption = CoordinateLabelOption, &
             CoordinateSystemOption = CoordinateSystemOption, &
             NameOption = NameOption, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             RatioOption = RatioOption, &
             ScaleOption = ScaleOption, &
             nCellsOption = nCellsOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nBricksOption = nBricksOption, &
             nDimensionsOption = nDimensionsOption, &
             iDimensionalityOption = 2 )

    B % Atlas_SCG_Fiber  =>  AF
    B % Chart_GS_Fiber   =>  AF % Chart_GS

    associate ( CF  =>  B % Chart_GS_Fiber )
    B % nSections  =  product ( CF % nCells )
    end associate !-- CF

    end select !-- AF

    !-- Decomposition

    call SetDecomposition ( B )

  end subroutine Initialize_ASCG_ASCG


  subroutine Show_B ( B )

    class ( Bundle_ASCG_ASCG_Form ), intent ( in ) :: &
      B

    integer ( KDI ) :: &
      iPE  !-- iProcessExchange

    call B % Bundle_H_Form % Show ( )

    call Show ( 'Bundle_ASCG_ASCG Proper Parameters' )

    call Show ( B % nProcesses,     'nProcesses',      B % IGNORABILITY )
    call Show ( B % nCopiesBase,    'nCopiesBase',     B % IGNORABILITY )
    call Show ( B % nProcessesBase, 'nProcessesBase',  B % IGNORABILITY )

    call Show ( 'Fiber-centric storage parameters' )
    call Show ( B % nFibers,       'nFibers',       B % IGNORABILITY )
    call Show ( B % nMyFibers,     'nMyFibers',     B % IGNORABILITY )
    call Show ( B % iaBrickFirst,  'iaBrickFirst',  B % IGNORABILITY )
    call Show ( B % iaBrickLast,   'iaBrickLast',   B % IGNORABILITY )
    call Show ( B % iaCellFirst,   'iaCellFirst',   B % IGNORABILITY )
    call Show ( B % iaCellLast,    'iaCellLast',    B % IGNORABILITY )
    call Show ( B % nFibersGlobal, 'nFibersGlobal', B % IGNORABILITY + 2 )
    call Show ( B % iFiberExchangeFirst, 'iFiberExchangeFirst', &
                B % IGNORABILITY + 2 )
    call Show ( B % iFiberExchangeLast, 'iFiberExchangeLast', &
                B % IGNORABILITY + 2 )
    do iPE  =  1, size ( B % iaBinExchangeFirst )
      call Show ( iPE, 'iProcessExchange', B % IGNORABILITY + 2 )
      call Show ( B % iaBinExchangeFirst ( iPE ) % Value, &
                  'iaBinExchangeFirst', B % IGNORABILITY + 2 )
      call Show ( B % iaBinExchangeLast ( iPE ) % Value, &
                  'iaBinExchangeLast', B % IGNORABILITY + 2 )
    end do

    call Show ( 'Section-centric storage parameters' )
    call Show ( B % nSections,       'nSections',       B % IGNORABILITY )
    call Show ( B % nMySections,     'nMySections',     B % IGNORABILITY )
    call Show ( B % iCopyBase,       'iCopyBase',       B % IGNORABILITY )
    call Show ( B % iFiberFirst,     'iFiberFirst',     B % IGNORABILITY )
    call Show ( B % iFiberLast,      'iFiberLast',      B % IGNORABILITY )
    call Show ( B % iaBinFirst,      'iaBinFirst',      B % IGNORABILITY )
    call Show ( B % iaBinLast,       'iaBinLast',       B % IGNORABILITY )
    call Show ( B % nSectionsGlobal, 'nSectionsGlobal', B % IGNORABILITY + 2 )
    do iPE  =  1, size ( B % iaCellExchangeFirst )
      call Show ( iPE, 'iProcessExchange', B % IGNORABILITY + 2 )
      call Show ( B % iaCellExchangeFirst ( iPE ) % Value, &
                  'iaCellExchangeFirst', B % IGNORABILITY + 2 )
      call Show ( B % iaCellExchangeLast ( iPE ) % Value, &
                  'iaCellExchangeLast', B % IGNORABILITY + 2 )
    end do

    call B % Portal_F_S % Show ( 'Portal_F_S', B % IGNORABILITY + 2 )
    call B % Portal_S_F % Show ( 'Portal_S_F', B % IGNORABILITY + 2 )

  end subroutine Show_B


  subroutine Finalize ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    if ( allocated ( B % iaCellExchangeLast ) ) &
      deallocate ( B % iaCellExchangeLast )
    if ( allocated ( B % iaCellExchangeFirst ) ) &
      deallocate ( B % iaCellExchangeFirst )
    if ( allocated ( B % nSectionsGlobal ) ) &
      deallocate ( B % nSectionsGlobal )
    if ( allocated ( B % iaBinLast ) ) &
      deallocate ( B % iaBinLast )
    if ( allocated ( B % iaBinFirst ) ) &
      deallocate ( B % iaBinFirst )

    if ( allocated ( B % iaBinExchangeLast ) ) &
      deallocate ( B % iaBinExchangeLast )
    if ( allocated ( B % iaBinExchangeFirst ) ) &
      deallocate ( B % iaBinExchangeFirst )
    if ( allocated ( B % iFiberExchangeLast ) ) &
      deallocate ( B % iFiberExchangeLast )
    if ( allocated ( B % iFiberExchangeFirst ) ) &
      deallocate ( B % iFiberExchangeFirst )
    if ( allocated ( B % nFibersGlobal ) ) &
      deallocate ( B % nFibersGlobal )
    if ( allocated ( B % iaCellLast ) ) &
      deallocate ( B % iaCellLast )
    if ( allocated ( B % iaCellFirst ) ) &
      deallocate ( B % iaCellFirst )
    if ( allocated ( B % iaBrickLast ) ) &
      deallocate ( B % iaBrickLast )
    if ( allocated ( B % iaBrickFirst ) ) &
      deallocate ( B % iaBrickFirst )

    if ( allocated ( B % Portal_S_F ) ) &
      deallocate ( B % Portal_S_F )
    if ( allocated ( B % Portal_F_S ) ) &
      deallocate ( B % Portal_F_S )
    nullify ( B % Communicator )

    nullify ( B % Chart_GS_Fiber )
    nullify ( B % Atlas_SCG_Fiber )
    nullify ( B % Chart_GS_Base )
    nullify ( B % Atlas_SCG_Base )

  end subroutine Finalize


  subroutine SetDecomposition ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    integer ( KDI ) :: &
      iB, jB, kB, &     !-- iBrick, etc.
      iC, jC, kC, &     !-- iCell, etc.
      iCB, &            !-- iCopyBase
      iP, &             !-- iProcess
      iF, &             !-- iFiber
      iS, &             !-- iSection
      nC_1, nC_2, nC_3  !-- nCells_1, etc.
    logical ( KDL ) :: &
      NextProcess, &
      NextCopyBase

    !--  The base manifold is distributed, and there are nCopiesBase copies 
    !      of this distributed base manifold, in order to perform base manifold 
    !      operations in parallel for separate groups of fiber cells.
    !--  The fiber is not distributed; a fiber operation for a given
    !      base manifold cell is performed locally.

    if ( associated ( B % Chart_GS_Base % Communicator % Parent ) ) then
      B % Communicator  =>  B % Chart_GS_Base % Communicator % Parent
    else
      B % Communicator  =>  B % Chart_GS_Base % Communicator
    end if

    associate &
      ( C    =>  B % Communicator, &
        CB   =>  B % Chart_GS_Base % Communicator, &
        nP   =>  B % nProcesses, &
        nPB  =>  B % nProcessesBase, &
        nCB  =>  B % nCopiesBase, &
        nF   =>  B % nFibers, &
        nS   =>  B % nSections )

    nP   =  C  % Size
    nPB  =  CB % Size
    nCB  =  nP  /  nPB
    if ( nPB * nCB  /=  nP ) then
      call Show ( 'Size of the base manifold communicator must evenly ' // &
                  'divide the size of its parent communicator', &
                  CONSOLE % ERROR )
      call Show ( CB % Size, 'Base communicator size', &
                  CONSOLE % ERROR )
      call Show ( C % Size, 'Parent communicator size', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    !-- My Fibers

    allocate ( B % nFibersGlobal ( 0 : nP - 1 ) )

    associate &
      ( MyRank  =>  C % Rank, &
           nFG  =>  B % nFibersGlobal, &
            nD  =>  B % Chart_GS_Base % nDimensions, &
            nB  =>  B % Chart_GS_Base % nBricks, &
          nCBG  =>  B % Chart_GS_Base % nCellsBrickGlobal )

    nFG  =  nF  /  nP
    do iP  =  0,  mod ( nF, nP )  -  1
      nFG ( iP )  =  nFG ( iP )  +  1 
    end do !-- iP
    B % nMyFibers  =  nFG ( MyRank )

    allocate ( B % iaBrickFirst ( MAX_DIMENSIONS ) )
    allocate ( B % iaBrickLast  ( MAX_DIMENSIONS ) )
    allocate ( B % iaCellFirst  ( MAX_DIMENSIONS ) )
    allocate ( B % iaCellLast   ( MAX_DIMENSIONS ) )
    B % iaBrickFirst  =  0
    B % iaBrickLast   =  0
    B % iaCellFirst   =  0
    B % iaCellLast    =  0
    iP  =  0
    iF  =  0
    if ( MyRank  == 0 ) then
      NextProcess  =  .true.
    else
      NextProcess  =  .false.
    end if
    do kB  =  1,  nB ( 3 )
      do jB  =  1,  nB ( 2 )
        do iB  =  1,  nB ( 1 )
          nC_1  =  nCBG ( 1 ) % Value ( iB )
          nC_2  =  1
          nC_3  =  1
          if ( nD  >  1 ) &
            nC_2  =  nCBG ( 2 ) % Value ( jB )   
          if ( nD  >  2 ) &
            nC_3  =  nCBG ( 3 ) % Value ( kB )   
          do kC  =  1,  nC_3
            do jC  =  1,  nC_2
              do iC  =  1,  nC_1
                iF  =  iF  +  1
                if ( NextProcess ) then
                  if ( iP  ==  MyRank ) then
                     B % iaBrickFirst  =  [ iB, jB, kB ]
                     B % iaCellFirst   =  [ iC, jC, kC ]
                  end if
                  NextProcess  =  .false.
                end if
                if ( iF  ==  nFG ( iP ) ) then
                  if ( iP  ==  MyRank ) then
                    B % iaBrickLast  =  [ iB, jB, kB ]
                    B % iaCellLast   =  [ iC, jC, kC ]
                  end if
                  iP  =  iP  +  1
                  iF  =  0
                  NextProcess  =  .true.
                end if
              end do !-- iC
            end do !-- jC
          end do !-- kC
        end do !-- iB
      end do !-- jB
    end do !-- kB

    end associate !-- MyRank, etc.

    !-- My Sections

    allocate ( B % nSectionsGlobal ( nCB ) )

    !-- WARNING: MyCopyBase assumes Base ranks are contiguous in the parent
    B % iCopyBase  =  C % Rank  /  nPB  +  1
    associate &
      ( MyCopyBase  =>  B % iCopyBase, &
               nSG  =>  B % nSectionsGlobal, &
                nC  =>  B % Chart_GS_Fiber % nCells )

    nSG  =  nS  /  nCB
    !-- WARNING: Assumes Base ranks are contiguous in the parent
    do iCB  =  1,  mod ( nS, nCB )
      nSG ( iCB )  =  nSG ( iCB )  +  1 
    end do !-- iP
    B % nMySections  =  nSG ( MyCopyBase )

    allocate ( B % iaBinFirst ( MAX_DIMENSIONS ) )
    allocate ( B % iaBinLast  ( MAX_DIMENSIONS ) )
    B % iaBinFirst  =  0
    B % iaBinLast   =  0
    iCB  =  1
    iS   =  0
    if ( MyCopyBase  == 1 ) then
      NextCopyBase  =  .true.
    else
      NextCopyBase  =  .false.
    end if
    do kC  =  1,  nC ( 3 )
      do jC  =  1,  nC ( 2 )
        do iC  =  1,  nC ( 1 )
          iS  =  iS  +  1
          if ( NextCopyBase ) then
            if ( iCB  ==  MyCopyBase ) then
              B % iaBinFirst  =  [ iC, jC, kC ]
            end if
            NextCopyBase  =  .false.
          end if
          if ( iS  ==  nSG ( iCB ) ) then
            if ( iCB  ==  MyCopyBase ) then
              B % iaBinLast  =  [ iC, jC, kC ]
            end if
            iCB  =  iCB  +  1
            iS  =  0
            NextCopyBase  =  .true.
          end if
        end do !-- iC
      end do !-- jC
    end do !-- kC

    end associate !-- MyCopyBase, etc.

    associate &
      (     nD  =>  B % Chart_GS_Base % nDimensions, &
            nB  =>  B % Chart_GS_Base % nBricks, &
           iaB  =>  B % Chart_GS_Base % iaBrick, &
          nCBG  =>  B % Chart_GS_Base % nCellsBrickGlobal )
    iF  =  0
    do kB  =  1,  nB ( 3 )
      do jB  =  1,  nB ( 2 )
        do iB  =  1,  nB ( 1 )
          nC_1  =  nCBG ( 1 ) % Value ( iB )
          nC_2  =  1
          nC_3  =  1
          if ( nD  >  1 ) &
            nC_2  =  nCBG ( 2 ) % Value ( jB )   
          if ( nD  >  2 ) &
            nC_3  =  nCBG ( 3 ) % Value ( kB )
          if ( all ( [ iB, jB, kB ]  ==  iaB ) ) then
            B % iFiberFirst  =  iF  +  1
            B % iFiberLast   =  iF  +  nC_1 * nC_2 * nC_3 
          end if
          iF  =  iF  +  nC_1 * nC_2 * nC_3
        end do !-- iB
      end do !-- jB
    end do !-- kB
    end associate !-- nD, etc.

    end associate !-- C, etc.

    call SetPortals ( B )

  end subroutine SetDecomposition


  subroutine SetPortals ( B )

    !-- WARNING: Assumes Base ranks are contiguous in the parent

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    integer ( KDI ) :: &
      iPS, iPF, &    !-- iProcessSection, iProcessFiber
      iPFF, iPFL, &  !-- iProcessFiberFirst, iProcessFiberLast
      iB, jB, kB, &  !-- iBrick, etc.
      iC, jC, kC, &  !-- iCell, etc.
      iCB, &         !-- iCopyBase
      iF, &          !-- iFiber
      iS, &          !-- iSection
      oF, &          !-- oFiber
      oPS, &         !-- iProcessSection
      nC_1, nC_2, nC_3, &  !-- nCells_1, etc.
      nProcessesExchange_F, &
      nProcessesExchange_S
    integer ( KDI ), dimension ( : ), allocatable :: &
      Source_F_S, Target_F_S, &
      nChunksFrom_F_S, nChunksTo_F_S, &
      Source_S_F, Target_S_F, &
      nChunksFrom_S_F, nChunksTo_S_F
    integer ( KDI ), dimension ( :, :, :, : ), allocatable :: &
      Process_S
    logical ( KDL ) :: &
      MyFibers

    !-- Fiber-centric perspective

    associate &
      (   nP  =>  B % nProcesses, &
         nCB  =>  B % nCopiesBase, &
          nB  =>  B % Chart_GS_Base % nBricks, &
        iaBF  =>  B % iaBrickFirst, &
        iaBL  =>  B % iaBrickLast, &
         nMF  =>  B % nMyFibers, &
         nSG  =>  B % nSectionsGlobal, &
        nPEF  =>  nProcessesExchange_F )

    nPEF  =  0
    do kB  =  iaBF ( 3 ),  iaBL ( 3 )
      do jB  =  iaBF ( 2 ),  iaBL ( 2 )
        do iB  =  iaBF ( 1 ),  iaBL ( 1 )
          nPEF  =  nPEF  +  nCB
        end do !-- iB
      end do !-- jB
    end do !-- kB

    allocate ( Process_S ( nB ( 1 ), nB ( 2 ), nB ( 3 ), nCB ) )
    iPS  =  0
    do iCB  =  1,  nCB
      do kB  =  1,  nB ( 3 )
        do jB  =  1,  nB ( 2 )
          do iB  =  1,  nB ( 1 )
            Process_S ( iB, jB, kB, iCB )  =  iPS
            iPS  =  iPS  +  1
          end do !-- iB
        end do !-- jB
      end do !-- kB
    end do !-- iCB

    allocate ( Source_S_F ( nPEF ) )
    allocate ( Target_F_S ( nPEF ) )
    associate &
      ( iaBF  =>  B % iaBrickFirst, &
        iaBL  =>  B % iaBrickLast )
    iPS  =  1
    do iCB  =  1,  nCB
      do kB  =  iaBF ( 3 ),  iaBL ( 3 )
        do jB  =  iaBF ( 2 ),  iaBL ( 2 )
          do iB  =  iaBF ( 1 ),  iaBL ( 1 )
            Source_S_F ( iPS )  =  Process_S ( iB, jB, kB, iCB )
            Target_F_S ( iPS )  =  Process_S ( iB, jB, kB, iCB )
            iPS  =  iPS  +  1
          end do !-- iB
        end do !-- jB
      end do !-- kB
    end do !-- iCB
    end associate !-- iaBF, etc.

    allocate ( nChunksFrom_S_F ( nPEF ) )
    allocate ( nChunksTo_F_S ( nPEF ) )
    allocate ( B % iFiberExchangeFirst ( nPEF ) )
    allocate ( B % iFiberExchangeLast ( nPEF ) )

    associate &
      ( MyRank  =>  B % Communicator % Rank, &
            nD  =>  B % Chart_GS_Base % nDimensions, &
            nB  =>  B % Chart_GS_Base % nBricks, &
          nCBG  =>  B % Chart_GS_Base % nCellsBrickGlobal, &
           nFG  =>  B % nFibersGlobal, &
          iFEF  =>  B % iFiberExchangeFirst, &
          iFEL  =>  B % iFiberExchangeLast, &
           nCF  =>  nChunksFrom_S_F, &
           nCT  =>  nChunksTo_F_S )
    nCF  =  0
    nCT  =  0
    iPS  =  1
    iFEF  =  huge ( 1_KDI )
    iFEL  =  0
    do iCB  =  1,  nCB
      iPF  =  0
      iF   =  0
      MyFibers  =  .false.
      do kB  =  1,  nB ( 3 )
        do jB  =  1,  nB ( 2 )
          do iB  =  1,  nB ( 1 )
            nC_1  =  nCBG ( 1 ) % Value ( iB )
            nC_2  =  1
            nC_3  =  1
            if ( nD  >  1 ) &
              nC_2  =  nCBG ( 2 ) % Value ( jB )   
            if ( nD  >  2 ) &
              nC_3  =  nCBG ( 3 ) % Value ( kB )
            do kC  =  1,  nC_3
              do jC  =  1,  nC_2
                do iC  =  1,  nC_1
                  iF  =  iF  +  1
                  if ( iPF  ==  MyRank ) then
                    MyFibers  =  .true.
                    iFEF ( iPS )  =  min ( iFEF ( iPS ), iF )
                    nCF  ( iPS )  =  nCF ( iPS )  +  1
                    nCT  ( iPS )  =  nCT ( iPS )  +  1
                  end if
                  if ( iF  ==  nFG ( iPF ) ) then
                    if ( MyFibers ) then
                      MyFibers  =  .false.
                      iFEL ( iPS )  =  iF
                      iPS  =  iPS  +  1  !-- increment for next Base Copy
                    end if
                    iPF  =  iPF  +  1
                    iF   =  0
                  end if
                end do !-- iC
              end do !-- jC
            end do !-- kC
            if ( MyFibers ) then
              iFEL ( iPS )  =  iF
              iPS  =  iPS  +  1
            end if
          end do !-- iB
        end do !-- jB
      end do !-- kB
    end do !-- iCB
    end associate !-- nD, etc.

    associate &
      ( nPB  =>  B % nProcessesBase, &
          S  =>  Source_S_F, &
          T  =>  Target_F_S, &
        nCF  =>  nChunksFrom_S_F, &
        nCT  =>  nChunksTo_F_S )
    do iPS  =  1,  nPEF

      iCB  =  S ( iPS )  /  nPB  +  1
      nCF ( iPS )  =  nCF ( iPS )  *  nSG ( iCB )

      iCB  =  T ( iPS )  /  nPB  +  1 
      nCT ( iPS )  =  nCT ( iPS )  *  nSG ( iCB )

    end do !-- iPS

    end associate !-- nPB, etc.

    allocate ( B % iaBinExchangeFirst ( nPEF ) )
    allocate ( B % iaBinExchangeLast  ( nPEF ) )
    associate &
      ( iaBEF  =>  B % iaBinExchangeFirst, &
        iaBEL  =>  B % iaBinExchangeLast, &
        nC     =>  B % Chart_GS_Fiber % nCells )
    iS   =  0
    iCB  =  1
    oPS  =  ( iCB - 1 )  *  nPEF / nCB 
    do kC  =  1,  nC ( 3 )
      do jC  =  1,  nC ( 2 )
        do iC  =  1,  nC ( 1 )
          iS  =  iS  +  1
          if ( iS  ==  1 ) then
            do iPS  =  oPS + 1,  oPS + nPEF / nCB
              call iaBEF ( iPS ) % Initialize ( [ iC, jC, kC ] )
            end do !-- iPS
          end if
          if ( iS  ==  nSG ( iCB ) ) then
            iS  =  0
            do iPS  =  oPS  +  1,  oPS  +  nPEF / nCB
              call iaBEL ( iPS ) % Initialize ( [ iC, jC, kC ] )
            end do !-- iPS
            iCB  =  iCB  +  1
            oPS  =  ( iCB - 1 )  *  nPEF / nCB 
          end if
        end do !-- iC
      end do !-- jC
    end do !-- kC
    end associate !-- iaBEF, etc.

    end associate !-- nP, etc.

    !-- Section-centric perspective

    associate &
      ( MyRank  =>  B % Communicator % Rank, &
           nCB  =>  B % nCopiesBase, &
            nP  =>  B % nProcesses, &
            nD  =>  B % Chart_GS_Base % nDimensions, &
            nB  =>  B % Chart_GS_Base % nBricks, &
          nCBG  =>  B % Chart_GS_Base % nCellsBrickGlobal, &
           iFF  =>  B % iFiberFirst, &
           iFL  =>  B % iFiberLast, &
           nMS  =>  B % nMySections, &
           nFG  =>  B % nFibersGlobal, &
          nPES  =>  nProcessesExchange_S )

    oF    =  0
    do iPF  =  0,  nP - 1
      if ( iFF  >  oF  .and.  iFF  <=  oF  +  nFG ( iPF ) ) &
        iPFF  =  iPF
      if ( iFL  >  oF  .and.  iFL <=  oF  +  nFG ( iPF ) ) &
        iPFL  =  iPF
      oF  =  oF  +  nFG ( iPF )
    end do !-- iP
    nPES  =  iPFL - iPFF + 1

    allocate ( Source_F_S ( nPES ) )
    allocate ( Target_S_F ( nPES ) )
    Source_F_S  =  [ ( iPF, iPF = iPFF, iPFL ) ]
    Target_S_F  =  [ ( iPF, iPF = iPFF, iPFL ) ]

    allocate ( B % iaCellExchangeFirst ( nPES ) )
    allocate ( B % iaCellExchangeLast ( nPES ) )
    allocate ( nChunksFrom_F_S ( nPES ) )
    allocate ( nChunksTo_S_F ( nPES ) )
    associate &
       ( iaCEF  =>  B % iaCellExchangeFirst, &
         iaCEL  =>  B % iaCellExchangeLast, &
          nCF   =>  nChunksFrom_F_S, &
          nCT   =>  nChunksTo_S_F )
    nCF  =  0
    nCT  =  0
    iPS  =  1
    do iCB  =  1,  nCB    
      iPF  =  0
      iF   =  0
      do kB  =  1,  nB ( 3 )
        do jB  =  1,  nB ( 2 )
          do iB  =  1,  nB ( 1 )
            nC_1  =  nCBG ( 1 ) % Value ( iB )
            nC_2  =  1
            nC_3  =  1
            if ( nD  >  1 ) &
              nC_2  =  nCBG ( 2 ) % Value ( jB )   
            if ( nD  >  2 ) &
              nC_3  =  nCBG ( 3 ) % Value ( kB )   
            do kC  =  1,  nC_3
              do jC  =  1,  nC_2
                do iC  =  1,  nC_1
                  iF  =  iF  +  1
                  if ( Process_S ( iB, jB, kB, iCB )  ==  MyRank ) then
                    if ( iF  ==  1  &
                         .or.  all ( [ iC, jC, kC ]  ==  [ 1, 1, 1 ] ) ) &
                    then
                      call iaCEF ( iPS ) % Initialize ( [ iC, jC, kC ] )
                    end if
                    !-- Workaround for GCC 13.2.0
                    !nCF ( iPS )  =  nCF ( iPS )  +  1
                    !nCT ( iPS )  =  nCT ( iPS )  +  1
                    nChunksFrom_F_S ( iPS )  =  nChunksFrom_F_S ( iPS )  +  1
                    nChunksTo_S_F   ( iPS )  =  nChunksTo_S_F   ( iPS )  +  1
                  end if
                  if ( iF  ==  nFG ( iPF )  &
                       .or.  all ( [ iC, jC, kC ] == [ nC_1, nC_2, nC_3 ] ) ) &
                  then
                    if ( iF  ==  nFG ( iPF ) ) then
                      iPF  =  iPF  +  1
                      iF   =  0
                    end if
                    if ( Process_S ( iB, jB, kB, iCB )  ==  MyRank ) then
                      call iaCEL ( iPS ) % Initialize ( [ iC, jC, kC ] )
                      iPS  =  iPS  +  1
                    end if
                  end if
                end do !-- iC
              end do !-- jC
            end do !-- kC
          end do !-- iB
        end do !-- jB
      end do !-- kB
    end do !-- iCB
    nCF  =  nCF  *  nMS
    nCT  =  nCT  *  nMS
    end associate !-- nCF, etc.

    end associate !-- MyRank, etc.

    !-- Initialize portals

! call Show ( Source_F_S, '>>> Source_F_S' )
! call Show ( nChunksFrom_F_S, '>>> nChunksFrom_F_S' )
! call Show ( Target_F_S, '>>> Target_F_S' )
! call Show ( nChunksTo_F_S, '>>> nChunksTo_F_S' )
! call Show ( Source_S_F, '>>> Source_S_F' )
! call Show ( nChunksFrom_S_F, '>>> nChunksFrom_S_F' )
! call Show ( Target_S_F, '>>> Target_S_F' )
! call Show ( nChunksTo_S_F, '>>> nChunksTo_S_F' )

    allocate ( B % Portal_F_S )
    associate ( P => B % Portal_F_S )
    call P % Initialize &
           ( Source_F_S, Target_F_S, nChunksFrom_F_S, nChunksTo_F_S )
    end associate !-- P

    allocate ( B % Portal_S_F )
    associate ( P => B % Portal_S_F )
    call P % Initialize &
           ( Source_S_F, Target_S_F, nChunksFrom_S_F, nChunksTo_S_F )
    end associate !-- P

  end subroutine SetPortals


end module Bundle_ASCG_ASCG__Form
