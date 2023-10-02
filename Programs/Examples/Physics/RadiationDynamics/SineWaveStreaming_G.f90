program SineWaveStreaming_G

  use GenASiS
  use SineWaveStreaming_Form

  implicit none

  type ( SineWaveStreamingForm ), allocatable :: &
    SWS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveStreaming_G', DimensionalityOption = '2D' )

  allocate ( SWS )
  call SWS % Initialize ( 'GREY', PROGRAM_HEADER % Name )
call Show ( 'SineWaveStreaming parameters' )
call SWS % ShowParameters ( )
select type ( I => SWS % Integrator )
class is ( Integrator_CS_Form )
  call I % X % Show ( )
  call I % Geometry_X % Show ( )
  call I % CurrentSet_X % Show ( )
end select
  ! call SWS % Evolve ( )
  ! call SWS % ComputeError ( )
  deallocate ( SWS )

  deallocate ( PROGRAM_HEADER )

end program SineWaveStreaming_G
