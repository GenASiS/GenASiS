program GPU_Rank_Binding

  use OMP_LIB
  use Basics
  
  implicit none
  
  integer ( KDI ) :: &
    nVariables
  real ( KDR ) :: &
    StartTime, &
    TotalTime, &
    DataSize_MB
  logical ( KDL ) :: &
    PinnedOption
  type ( StorageForm ) :: &
    S
  type ( TimerForm ) :: &
    T
    
  allocate ( PROGRAM_HEADER )
  
  associate ( PH => PROGRAM_HEADER )
  
  call PH % Initialize ( 'GPU_Rank_Binding' )
  
  call T % Initialize ( 'UpdateDevice', 0 )
  
  nVariables = 32
  call PH % GetParameter ( nVariables, 'nVariables' )
  
  PinnedOption = .false.
  call PH % GetParameter ( PinnedOption, 'PinnedOption' )
  
  call S % Initialize &
         ( ValueShape = [ 256 ** 3, nVariables ], &
           PinnedOption = PinnedOption )
  call random_number ( S % Value )
  
  DataSize_MB = 1.0_KDR * S % nValues * S % nVariables * 8 &
                / 1.0e6_KDR 
  call Show ( DataSize_MB,  'Data size (MB)')

  call Show ( SelectedDevice ( ), 'Selected Device' )
  
  call S % AllocateDevice ( )
  
  call Show ( SelectedDevice ( ), 'Selected Device' )
  
  call T % Start ( )
  
  call S % UpdateDevice ( )
  
  call T % Stop ( )
  call T % ShowTotal ( CONSOLE % INFO_1 )
  
  call PH % ShowStatistics ( CONSOLE % INFO_1 )
  
  end associate !-- PH
  
  deallocate ( PROGRAM_HEADER )

end program GPU_Rank_Binding