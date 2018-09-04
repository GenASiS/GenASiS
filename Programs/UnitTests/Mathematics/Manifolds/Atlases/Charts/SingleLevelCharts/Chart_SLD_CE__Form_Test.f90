module Chart_SLD_CE__Form_Test__Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLD_CE__Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'Chart_SLD_CE__Form_Test'

  type, public :: Chart_SLD_CE__Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GIS
    type ( AtlasHeaderForm ), allocatable :: &
      AtlasHeader
    type ( GeometryFlat_CSL_Form ), allocatable :: &
      Geometry
!    type ( Chart_SLD_CE_Form ), allocatable :: &
!      Chart
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Chart_SLD_CE__Form_Test_Form

contains


  subroutine Initialize ( CFT )

    class ( Chart_SLD_CE__Form_Test_Form ), intent ( inout ), target :: &
      CFT

    ! allocate &
    !   ( CFT % AtlasHeader, &
    !     CFT % Geometry, &
    !     CFT % Chart, &
    !     CFT % GIS )
    ! associate &
    !   ( A => CFT % AtlasHeader, &
    !     G => CFT % Geometry, &
    !     C => CFT % Chart, &
    !     GIS => CFT % GIS )

    ! call A % Initialize &
    !        ( 'Atlas', CommunicatorOption = PROGRAM_HEADER % Communicator )

    ! call C % Initialize ( A, iChart = 1 )

    ! associate ( nValues => C % nProperCells + C % nGhostCells )
    ! call G % InitializeFlat ( C, 'Geometry', nValues )
    ! call C % AddField ( G )
    ! C % iFieldGeometry = C % nFields
    ! call C % SetGeometry ( G )
    ! end associate !-- nValues

    ! call A % Show ( )
    ! call C % Show ( )

    ! call GIS % Initialize &
    !        ( PROGRAM_HEADER % Name, CommunicatorOption = A % Communicator )
    ! call C % OpenStream ( GIS, '1', iStream = 1 )

    ! call GIS % Open ( GIS % ACCESS_CREATE )
    ! call C % Write ( iStream = 1 )
    ! call GIS % Close ( )

    ! end associate !-- A, etc.

  end subroutine Initialize


  subroutine Finalize ( CFT )

    type ( Chart_SLD_CE__Form_Test_Form ), intent ( inout ) :: &
      CFT

    ! deallocate ( CFT % GIS )
    ! deallocate ( CFT % Chart )
    ! deallocate ( CFT % Geometry )
    ! deallocate ( CFT % AtlasHeader )

  end subroutine Finalize


end module Chart_SLD_CE__Form_Test__Form



program Chart_SLD_CE__Form_Test

  use Basics
  use Chart_SLD_CE__Form_Test__Form
  
  implicit none

  type ( Chart_SLD_CE__Form_Test_Form ), allocatable :: &
    CFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, DimensionalityOption = '2D' )
    
  allocate ( CFT )
  call CFT % Initialize ( )
  deallocate ( CFT )

  deallocate ( PROGRAM_HEADER )

end program Chart_SLD_CE__Form_Test
