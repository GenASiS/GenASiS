module Bundle_ASCG_ASCG__Form

  !-- Bundle_AtlasSingleChartGrid_AtlasSingleChartGrid__Form

  use Basics
  use Charts
  use BaseManifolds
  use Bundle_H__Form

  implicit none
  private

  type, public, extends ( Bundle_H_Form ) :: Bundle_ASCG_ASCG_Form
    integer ( KDI ) :: &
      nProcesses, &
      nCopiesBase, &
      nProcessesBase
    integer ( KDI ) :: &
      nFibers, &
      nMyFibers
    integer ( KDI ) :: &
      nSections, &
      nMySections
    integer ( KDI ), dimension ( : ), allocatable :: &
      nFibersGlobal, &
      iaBrickFirst, iaBrickLast, &
      iaCellFirst, iaCellLast
    integer ( KDI ), dimension ( : ), allocatable :: &
      nSectionsGlobal, &
      iaBinFirst, iaBinLast
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( PortalHeaderForm ), allocatable :: &
      Portal_F_S, &
      Portal_S_F
    class ( Chart_GS_Form ), pointer :: &
      Chart_GS_Base  => null ( ), &
      Chart_GS_Fiber => null ( )
    class ( Atlas_SCG_Form ), pointer :: &
      Atlas_SCG_Base  => null ( ), &
      Atlas_SCG_Fiber => null ( )
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
    B % nSections  =  CF % nCellsProper
    end associate !-- CF

    end select !-- AF

    !-- Decomposition

    call SetDecomposition ( B )

  end subroutine Initialize_ASCG_ASCG


  subroutine Show_B ( B )

    class ( Bundle_ASCG_ASCG_Form ), intent ( in ) :: &
      B

    call B % Bundle_H_Form % Show ( )

    call Show ( 'Bundle_ASCG_ASCG Proper Parameters' )

    call Show ( B % nProcesses,     'nProcesses',      B % IGNORABILITY )
    call Show ( B % nCopiesBase,    'nCopiesBase',     B % IGNORABILITY )
    call Show ( B % nProcessesBase, 'nProcessesBase',  B % IGNORABILITY )

    call Show ( B % nFibers,       'nFibers',       B % IGNORABILITY )
    call Show ( B % nMyFibers,     'nMyFibers',     B % IGNORABILITY )
    call Show ( B % nFibersGlobal, 'nFibersGlobal', B % IGNORABILITY + 2 )
    call Show ( B % iaBrickFirst,  'iaBrickFirst',  B % IGNORABILITY + 2 )
    call Show ( B % iaBrickLast,   'iaBrickLast',   B % IGNORABILITY + 2 )
    call Show ( B % iaCellFirst,   'iaCellFirst',   B % IGNORABILITY + 2 )
    call Show ( B % iaCellLast,    'iaCellLast',    B % IGNORABILITY + 2 )

    call Show ( B % nSections,       'nSections',       B % IGNORABILITY )
    call Show ( B % nMySections,     'nMySections',     B % IGNORABILITY )
    call Show ( B % nSectionsGlobal, 'nSectionsGlobal', B % IGNORABILITY + 2 )
    call Show ( B % iaBinFirst,      'iaBinFirst',      B % IGNORABILITY + 2 )
    call Show ( B % iaBinLast,       'iaBinLast',       B % IGNORABILITY + 2 )

  end subroutine Show_B


  subroutine Finalize ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    nullify ( B % Atlas_SCG_Fiber )
    nullify ( B % Atlas_SCG_Base )
    nullify ( B % Chart_GS_Fiber )
    nullify ( B % Chart_GS_Base )
    nullify ( B % Communicator )

    if ( allocated ( B % Portal_S_F ) ) &
      deallocate ( B % Portal_S_F )
    if ( allocated ( B % Portal_F_S ) ) &
      deallocate ( B % Portal_F_S )

    if ( allocated ( B % iaBinLast ) ) &
      deallocate ( B % iaBinLast )
    if ( allocated ( B % iaBinFirst ) ) &
      deallocate ( B % iaBinFirst )
    if ( allocated ( B % nSectionsGlobal ) ) &
      deallocate ( B % nSectionsGlobal )

    if ( allocated ( B % iaCellLast ) ) &
      deallocate ( B % iaCellLast )
    if ( allocated ( B % iaCellFirst ) ) &
      deallocate ( B % iaCellFirst )
    if ( allocated ( B % iaBrickLast ) ) &
      deallocate ( B % iaBrickLast )
    if ( allocated ( B % iaBrickFirst ) ) &
      deallocate ( B % iaBrickFirst )
    if ( allocated ( B % nFibersGlobal ) ) &
      deallocate ( B % nFibersGlobal )

  end subroutine Finalize


  subroutine SetDecomposition ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    integer ( KDI ) :: &
      iB, jB, kB, &  !-- iBrick, etc.
      iC, jC, kC, &  !-- iBrick, etc.
      iCB, &         !-- iCopyBase
      iP, &          !-- iProcess
      iF, &          !-- iFiber
      iS, &          !-- iSection
      oP, &          !-- oProcess
      nC_1, nC_2, nC_3
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
    associate ( nFG  => B % nFibersGlobal ) 
    nFG  =  nF  /  nP
    do iP  =  0,  mod ( nF, nP )  -  1
      nFG ( iP )  =  nFG ( iP )  +  1 
    end do !-- iP
    B % nMyFibers  =  nFG ( C % Rank )
    end associate !-- nSG

    allocate ( B % iaBrickFirst ( MAX_DIMENSIONS ) )
    allocate ( B % iaBrickLast  ( MAX_DIMENSIONS ) )
    allocate ( B % iaCellFirst  ( MAX_DIMENSIONS ) )
    allocate ( B % iaCellLast   ( MAX_DIMENSIONS ) )
    B % iaBrickFirst  =  0
    B % iaBrickLast   =  0
    B % iaCellFirst   =  0
    B % iaCellLast    =  0
    associate &
      ( MyRank  =>  C % Rank, &
           nFG  =>  B % nFibersGlobal, &
            nD  =>  B % Chart_GS_Base % nDimensions, &
            nB  =>  B % Chart_GS_Base % nBricks, &
          nCBG  =>  B % Chart_GS_Base % nCellsBrickGlobal )
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
              end do !-- iB
            end do !-- jB
          end do !-- kB
        end do !-- iB
      end do !-- jB
    end do !-- kB
    end associate !-- MyRank, etc.

    !-- My Sections

    allocate ( B % nSectionsGlobal ( 0 : nP - 1 ) )
    associate ( nSG  => B % nSectionsGlobal ) 
    nSG  =  nS  /  nCB
    !-- WARNING: Assumes Base ranks are contiguous in the parent
    oP  =  0
    do iCB  =  0,  mod ( nS, nCB )  -  1
      oP  =  oP  +  iCB * nPB
      nSG ( oP : oP + nPB - 1 )  =  nSG ( oP : oP + nPB - 1 )  +  1 
    end do !-- iP
    B % nMySections  =  nSG ( C % Rank )
    end associate !-- nSG

    allocate ( B % iaBinFirst ( MAX_DIMENSIONS ) )
    allocate ( B % iaBinLast  ( MAX_DIMENSIONS ) )
    B % iaBinFirst  =  0
    B % iaBinLast   =  0
    associate &
      ( MyCopyBase  =>  C % Rank  /  nPB  +  1, &
               nSG  =>  B % nSectionsGlobal, &
                nC  =>  B % Chart_GS_Fiber % nCells )
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
              B % iaBinFirst   =  [ iC, jC, kC ]
            end if
            NextCopyBase  =  .false.
          end if
          if ( iS  ==  nSG ( iCB ) ) then
            if ( iCB  ==  MyCopyBase ) then
              B % iaBinLast   =  [ iC, jC, kC ]
            end if
            iCB  =  iCB  +  1
            iS  =  0
            NextCopyBase  =  .true.
          end if
        end do !-- iB
      end do !-- jB
    end do !-- kB
    end associate !-- MyCopyBrick, etc.

    end associate !-- C, etc.

  end subroutine SetDecomposition


  subroutine SetPortals ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    integer ( KDI ) :: &
      iP, &          !-- iProcess
      iB, jB, kB, &  !-- iBrick, etc.
      iCB            !-- iCopyBase
    integer ( KDI ), dimension ( : ), allocatable :: &
      Source_F_S, &
      Source_S_F, &
      Target_F_S, &
      Target_S_F
    integer ( KDI ), dimension ( :, :, :, : ), allocatable :: &
      Process

    associate &
      ( nB   =>  B % Chart_GS_Base % nBricks, &
        nCB  =>  B % nCopiesBase, &
        nP   =>  B % Communicator % Size ) 

    allocate ( Process ( nB ( 1 ), nB ( 2 ), nB ( 3 ), nCB ) )

    iP  =  0
    do iCB  =  1,  nCB
      do kB  =  1,  nB ( 3 )
        do jB  =  1,  nB ( 2 )
          do iB  =  1,  nB ( 1 )
            Process ( iB, jB, kB, iCB )  =  iP
            iP  =  iP  +  1
          end do !-- iB
        end do !-- jB
      end do !-- kB
    end do !-- iCB


    allocate ( Source_F_S ( nP ) )
    allocate ( Source_S_F ( nP ) )
    allocate ( Target_F_S ( nP ) )
    allocate ( Target_S_F ( nP ) )
    Source_F_S  =  -1
    Source_S_F  =  -1
    Target_F_S  =  -1
    Target_S_F  =  -1

    end associate !-- nB, etc.

  end subroutine SetPortals


end module Bundle_ASCG_ASCG__Form
