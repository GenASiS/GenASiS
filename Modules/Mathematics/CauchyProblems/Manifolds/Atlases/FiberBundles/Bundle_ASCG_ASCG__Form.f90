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
      nFibers         = 0, &
      nSections       = 0, &
      nProcesses      = 0, &
      nProcessesBase  = 0, &
      nProcessesFiber = 0, &
      nMyFibers       = 0, &
      nMySections     = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      nFibersGlobal, &
      nSectionsGlobal
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
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

    integer ( KDI ) :: &
      iPF, &  !-- iProcessFiber
      iP, &   !-- iProcess
      oP      !-- oProcess

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

    !-- Distribution of Sections and Fibers
    !   ( Solution on the base manifold is distributed, and there are
    !     nProcessesFiber copies of this distributed base manifold. )
    !   ( Solution on a fiber is not distributed, however. )

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
        nPF  =>  B % nProcessesFiber, &
        nF   =>  B % nFibers, &
        nS   =>  B % nSections )

    nP   =  C  % Size
    nPB  =  CB % Size
    nPF  =  nP  /  nPB
    if ( nPB * nPF  /=  C % Size ) then
      call Show ( 'Size of the base manifold communicator must evenly ' // &
                  'divide the size of its parent communicator', &
                  CONSOLE % ERROR )
      call Show ( CB % Size,          'Base communicator size', &
                  CONSOLE % ERROR )
      call Show ( C % Size, 'Parent communicator size', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    allocate ( B % nFibersGlobal ( 0 : nP - 1 ) )
    associate ( nFG  => B % nFibersGlobal ) 
    nFG  =  nF  /  nP
    do iP  =  0,  mod ( nF, nP )  -  1
      nFG ( iP )  =  nFG ( iP )  +  1 
    end do !-- iP
    B % nMyFibers  =  nFG ( C % Rank )
    end associate !-- nSG

    allocate ( B % nSectionsGlobal ( 0 : nP - 1 ) )
    associate ( nSG  => B % nSectionsGlobal ) 
    nSG  =  nS  /  nPF
    !-- WARNING: Assumes Base ranks are contiguous in the parent
    oP  =  0
    do iPF  =  0,  mod ( nS, nPF )  -  1
      oP  =  oP  +  iPF * nPB
      nSG ( oP : oP + nPB - 1 )  =  nSG ( oP : oP + nPB - 1 )  +  1 
    end do !-- iP
    B % nMySections  =  nSG ( C % Rank )
    end associate !-- nSG

    end associate !-- C, etc.

  end subroutine Initialize_ASCG_ASCG


  subroutine Show_B ( B )

    class ( Bundle_ASCG_ASCG_Form ), intent ( in ) :: &
      B

    call B % Bundle_H_Form % Show ( )

    call Show ( B % nFibers,         'nFibers',         B % IGNORABILITY )
    call Show ( B % nSections,       'nSections',       B % IGNORABILITY )
    call Show ( B % nProcesses,      'nProcesses',      B % IGNORABILITY )
    call Show ( B % nProcessesBase,  'nProcessesBase',  B % IGNORABILITY )
    call Show ( B % nProcessesFiber, 'nProcessesFiber', B % IGNORABILITY )
    call Show ( B % nMyFibers,       'nMyFibers',       B % IGNORABILITY )
    call Show ( B % nMySections,     'nMySections',     B % IGNORABILITY )
    call Show ( B % nFibersGlobal,   'nFibersGlobal',   B % IGNORABILITY + 2 )
    call Show ( B % nSectionsGlobal, 'nSectionsGlobal', B % IGNORABILITY + 2 )

  end subroutine Show_B


  subroutine Finalize ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    nullify ( B % Atlas_SCG_Fiber )
    nullify ( B % Atlas_SCG_Base )
    nullify ( B % Chart_GS_Fiber )
    nullify ( B % Chart_GS_Base )
    nullify ( B % Communicator )

    if ( allocated ( B % nSectionsGlobal ) ) &
      deallocate ( B % nSectionsGlobal )

  end subroutine Finalize


end module Bundle_ASCG_ASCG__Form
