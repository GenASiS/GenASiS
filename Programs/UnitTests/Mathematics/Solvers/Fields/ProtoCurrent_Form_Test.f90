program ProtoCurrent_Form_Test

  use Basics
  use ProtoCurrent_Form

  implicit none

  character ( LDF ), parameter :: &
    ProgramName = 'ProtoCurrent_Form_Test'

  integer ( KDI ) :: &
    iF, &  !-- iField, &
    iV     !-- iVector
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    VelocityUnit
  type ( ProtoCurrentForm ), allocatable :: &
    PC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( PC )
  call PC % Initialize ( VelocityUnit, nValues = 8 )
  call PC % ShowPrimitiveConserved ( CONSOLE % INFO_1 )

  call Show ( 'ProtoCurrent variables' )
  call Show ( PC % Name, 'Name' )
  do iF = 1, PC % nVariables
    call Show ( PC % Value ( :, iF ), PC % Unit ( iF ), PC % Variable ( iF ) )
  end do

  call Show ( 'ProtoCurrent vectors' )
  call Show ( PC % Name, 'Name' )
  do iV = 1, PC % nVectors
    call Show ( PC % VectorIndices ( iV ) % Value, PC % Vector ( iV ) )
  end do

  deallocate ( PC )
  deallocate ( PROGRAM_HEADER )

end program ProtoCurrent_Form_Test
