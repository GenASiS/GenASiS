program Step_RK_H__Form_Test

  !-- Step_RK_Header__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Steps

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( TimerForm ), pointer :: &
    T
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    Sm
  type ( Step_RK_H_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Step_RK_H__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( Sm )
  call Sm % Initialize ( A, GIS )

  allocate ( S )
  call S % Initialize_H ( A )
  call S % SetStream ( Sm )

  call S  % Show ( )
  call Sm % Show ( )

  T  =>  S % Timer ( Level = 1 )
  call T % Start ( )
  call S % Compute ( T = 0.0_KDR, dT = 1.0e-2_KDR, T_Option = T )
  call T % Stop ( )

  deallocate ( S )
  deallocate ( Sm )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Step_RK_H__Form_Test
