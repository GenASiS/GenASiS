program Reconstruction_C__Form_Test

  !-- Reconstruction_Chart__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Slopes

  implicit none

  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    FSC
  type ( Geometry_F_C_Form ), allocatable :: &
    GC
  type ( Reconstruction_C_Form ), allocatable :: &
    RC_0, &
    RC_1, &
    RC_2

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Reconstruction_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( C )
  call C % Initialize &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSC )
  call FSC % Initialize ( C )

  allocate ( GC )
  call GC % Initialize ( C )

  allocate ( RC_0 )
  call RC_0 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_0', &
           StreamedOption = .true., &
           OrderOption = 0 )

  allocate ( RC_1 )
  call RC_1 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_1', &
           StreamedOption = .true., &
           OrderOption = 1 )

  allocate ( RC_2 )
  call RC_2 % Initialize &
         ( GC, FSC, &
           NameOption = 'Reconstruction_2', &
           StreamedOption = .true., &
           OrderOption = 2 )

  call   C   % Show ( )
  call FSC   % Show ( )
  call  GC   % Show ( )
  call  RC_0 % Show ( )
  call  RC_1 % Show ( )
  call  RC_2 % Show ( )

  deallocate ( RC_2 )
  deallocate ( RC_1 )
  deallocate ( RC_0 )
  deallocate ( GC )
  deallocate ( FSC )
  deallocate ( C )
  deallocate ( PROGRAM_HEADER )

end program Reconstruction_C__Form_Test
