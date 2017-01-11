program IncrementDivergence_FV__Form_Test

  use Basics
  use DensityWaveIncrement_Form

  implicit none

  type ( DensityWaveIncrementForm ), allocatable :: &
    DW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'IncrementDivergence_FV__Form_Test' )
    
  allocate ( DW )
  call DW % Initialize ( PROGRAM_HEADER % Name )
  call DW % ProtoCurrent % ComputeTally ( ComputeChangeOption = .false. )
  call DW % Write ( )
  deallocate ( DW )

  deallocate ( PROGRAM_HEADER )

end program IncrementDivergence_FV__Form_Test
