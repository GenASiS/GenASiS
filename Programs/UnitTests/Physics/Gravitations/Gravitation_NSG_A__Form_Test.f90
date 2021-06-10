program Gravitation_NSG_A__Form_Test

  !-- Gravitation_NewtonSelfGravity_Atlas__Form_Test

  use Basics
  use Mathematics
  use Gravitations

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Gravitation_NSG_A_Form ), allocatable :: &
    GA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Gravitation_NSG_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )
  call GA % SetStream ( SA )
  call SA % AddFieldSet ( GA % Source_A )

  call  A % Show ( )
  call GA % Show ( )
  call SA % Show ( )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call SA % Write ( )
  call GIS % Close ( )

  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Gravitation_NSG_A__Form_Test
