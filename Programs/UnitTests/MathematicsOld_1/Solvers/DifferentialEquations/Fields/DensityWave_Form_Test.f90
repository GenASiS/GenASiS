program DensityWave_Form_Test

  use Basics
  use DensityWave_Form

  implicit none

  type ( DensityWaveForm ), allocatable :: &
    DW

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'DensityWave' )
    
  allocate ( DW )
  call DW % Initialize ( PROGRAM_HEADER % Name )
  call DW % ProtoCurrent % ComputeTally ( ComputeChangeOption = .false. )
  call DW % Write ( )
  deallocate ( DW )

  deallocate ( PROGRAM_HEADER )

end program DensityWave_Form_Test
