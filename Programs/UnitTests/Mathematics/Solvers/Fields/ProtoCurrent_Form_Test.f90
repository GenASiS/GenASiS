program ProtoCurrent_Form_Test

  use Basics
  use CurrentSources_Form
  use ProtoCurrent_Form

  implicit none

  character ( LDF ), parameter :: &
    ProgramName = 'ProtoCurrent_Form_Test'

  integer ( KDI ) :: &
    iF, &  !-- iField, &
    iV     !-- iVector
  type ( MeasuredValueForm ) :: &
    LengthUnit
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    VelocityUnit
  type ( CurrentSourcesForm ), allocatable :: &
    PCS
  type ( ProtoCurrentForm ), allocatable :: &
    PC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( PC )
  call PC % Initialize ( VelocityUnit, nValues = 8 )
  call PC % SetPrimitiveConserved ( )
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

  allocate ( PCS )
  call PCS % Initialize ( PC, LengthUnit, PC % iaConserved )
  call PC % SetSources ( PCS )

  call Show ( 'ProtoCurrent % Sources variables' )
  call Show ( PC % Sources % Name, 'Name' )
  do iF = 1, PC % Sources % nVariables
    call Show ( PC % Sources % Value ( :, iF ), PC % Sources % Unit ( iF ), &
                PC % Sources % Variable ( iF ) )
  end do

  deallocate ( PCS )
  deallocate ( PC )
  deallocate ( PROGRAM_HEADER )

end program ProtoCurrent_Form_Test
