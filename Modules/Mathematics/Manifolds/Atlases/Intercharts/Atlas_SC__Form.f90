!-- Atlas_SC represents an Atlas with a single chart.

module Atlas_SC__Form

  !-- Atlas_SingleChart_Form

  use Basics
  use Charts
  use Atlas_SC__Template
  use GeometryFlat_ASC__Form

  implicit none
  private

  type, public, extends ( Atlas_SC_Template ) :: Atlas_SC_Form
    logical ( KDL ) :: &
      AllocatedGeometry = .false.
    class ( GeometryFlat_ASC_Form ), pointer :: &
      Geometry_ASC => null ( )
  contains
    procedure, public, pass :: &
      SetGeometry
    procedure, public, pass :: &
      SetCoarsening   
    procedure, private, pass :: &
      Geometry_CSL
    generic, public :: &
      Geometry => Geometry_CSL
    procedure, public, pass :: &
      OpenStream
    procedure, public, pass :: &
      Write
  !   procedure, public, pass :: &
  !     Read
    procedure, public, pass :: &
      CloseStreams
    final :: &
      Finalize
  end type Atlas_SC_Form

  type, public :: Atlas_SC_ElementForm
    class ( Atlas_SC_Form ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type Atlas_SC_ElementForm

  type, public :: Atlas_SC_1D_Form
    class ( Atlas_SC_ElementForm ), dimension ( : ), allocatable :: &
      Atlas
  contains
    procedure, public, pass :: &
      Initialize => Initialize_1D
    final :: &
      Finalize_1D
  end type Atlas_SC_1D_Form

contains


  subroutine SetGeometry &
               ( A, GeometryOption, UsePinnedMemoryOption, EdgeOption )

    class ( Atlas_SC_Form ), intent ( inout ) :: &
      A
    class ( GeometryFlat_ASC_Form ), intent ( in ), target, optional :: &
      GeometryOption
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    type ( Real_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      EdgeOption

    if ( associated ( A % Geometry_ASC ) ) then
      A % AllocatedGeometry = .true.
    else
      if ( present ( GeometryOption ) ) then
        A % AllocatedGeometry = .false.
        A % Geometry_ASC => GeometryOption
      else
        A % AllocatedGeometry = .true.
        allocate ( A % Geometry_ASC )
        associate ( GA => A % Geometry_ASC )
        call GA % InitializeFlat &
               ( A, UsePinnedMemoryOption = UsePinnedMemoryOption ) 
        end associate !-- GA
      end if
    end if

    call A % AddField ( A % Geometry_ASC )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )

      C % iFieldGeometry = C % nFields

      select type ( GC => A % Geometry_ASC % Chart )
      class is ( GeometryFlat_CSL_Form ) 
      call C % SetGeometry ( GC, EdgeOption )
      end select !-- G_CSL

      call A % Show ( )
      call C % Show ( )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetGeometry', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine SetGeometry


  subroutine SetCoarsening ( A )

    class ( Atlas_SC_Form ), intent ( inout ) :: &
      A

    select type ( C => A % Chart )
    class is ( Chart_SLD_CE_Form )
      call C % SetCoarsening ( )
    class is ( Chart_SLD_CC_Form )
      call C % SetCoarsening ( )
    class default
      call Show ( 'Chart type does not support coarsening', CONSOLE % WARNING )
      call Show ( 'Name', C % Name, CONSOLE % WARNING )
      call Show ( 'Type', C % Type, CONSOLE % WARNING )
      call Show ( 'Atlas_SC__Form', 'module', CONSOLE % WARNING )
      call Show ( 'SetCoarsening', 'subroutine', CONSOLE % WARNING )
    end select !-- C

  end subroutine SetCoarsening


  function Geometry_CSL ( A ) result ( G )

    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
      G => C % Geometry ( )
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Geometry_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end function Geometry_CSL


  subroutine OpenStream ( A, GIS, Label, iStream, VerboseOption )

    class ( Atlas_SC_Form ), intent ( inout ) :: &
      A
    type ( GridImageStreamForm ), intent ( in )  :: &
      GIS
    character ( * ), intent ( in ) :: &
      Label
    integer ( KDI ), intent ( in ) :: &
      iStream
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
      call C % OpenStream ( GIS, Label, iStream, VerboseOption )
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'OpenStream', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine OpenStream


  subroutine Write ( A, iStream, DirectoryOption, TimeOption, &
                     CycleNumberOption )

    class ( Atlas_SC_Form ), intent ( inout ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iStream
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
  
    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
      call C % Write &
             ( iStream, DirectoryOption = DirectoryOption, &
               TimeOption = TimeOption, &
               CycleNumberOption = CycleNumberOption )
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Write', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine Write


!   subroutine Read ( A, iStream, TimeOption, CycleNumberOption )

!     class ( AtlasForm ), intent ( inout ) :: &
!       A
!     integer ( KDI ), intent ( in ) :: &
!       iStream
!     type ( MeasuredValueForm ), intent ( out ), optional :: &
!       TimeOption
!     integer ( KDI ), intent ( out ), optional :: &
!       CycleNumberOption

!     integer ( KDI ) :: &
!       iC  !-- iChart

!     associate &
!       ( Timer => PROGRAM_HEADER % Timer ( A % InputOutputTimer ( iStream ) ) )
!     call Timer % Start ( )

!     do iC = 1, A % nCharts
!       select type ( C => A % Chart ( iC ) )
!       class is ( ChartForm )
!         call C % Read ( iStream, TimeOption, CycleNumberOption )
!       end select
!     end do

!     call Timer % Stop ( )
!     end associate !-- Timer

!   end subroutine Read


  subroutine CloseStreams ( A )

    class ( Atlas_SC_Form ), intent ( inout ) :: &
      A

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
      call C % CloseStreams ( )
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'CloseStreams', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine CloseStreams


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SC_Form ), intent ( inout ) :: &
      A

    if ( A % AllocatedGeometry ) then
      deallocate ( A % Geometry_ASC )
    else
      nullify ( A % Geometry_ASC )
    end if

    call A % FinalizeTemplate ( )

  end subroutine Finalize


  impure elemental subroutine FinalizeElement ( AE )
    
    type ( Atlas_SC_ElementForm ), intent ( inout ) :: &
      AE

    if ( allocated ( AE % Element ) ) deallocate ( AE % Element )

  end subroutine FinalizeElement


  subroutine Initialize_1D ( A_1D, nElements )

    class ( Atlas_SC_1D_Form ), intent ( inout ) :: &
      A_1D
    integer ( KDI ), intent ( in ) :: &
      nElements

    allocate ( A_1D % Atlas ( nElements ) )

  end subroutine Initialize_1D


  impure elemental subroutine Finalize_1D ( A_1D )

    type ( Atlas_SC_1D_Form ), intent ( inout ) :: &
      A_1D
    
    if ( allocated ( A_1D % Atlas ) ) &
      deallocate ( A_1D % Atlas )

  end subroutine Finalize_1D


end module Atlas_SC__Form
