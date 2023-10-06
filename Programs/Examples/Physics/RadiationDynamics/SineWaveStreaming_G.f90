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
  call I % Show ( )
  associate ( GIS => I % GridImageStream )
  associate ( S_X  =>  I % Checkpoint_X )
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S_X % Write &
         ( TimeOption  =  I % T  /  I % Unit_T, &
           CycleNumberOption  =  I % iCycle )
  call GIS % Close ( )
  end associate !-- S_X
  end associate !-- GIS
end select
  ! call SWS % Evolve ( )
  ! call SWS % ComputeError ( )
  deallocate ( SWS )

  deallocate ( PROGRAM_HEADER )

end program SineWaveStreaming_G
