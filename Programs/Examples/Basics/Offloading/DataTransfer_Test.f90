program DataManagement_Test

  !-- Replicate most of Storage_Form_Test for data movement with the addition
  !   of parallelism provided by Basics functionalities
  
  use Basics
  
  implicit none
  
  real ( KDR ) :: &
    DataSize_GiB, &
    Bandwidth, &
    BandwidthMax, &
    BandwidthMin, &
    Diff
  type ( StorageForm ) :: &
    S
  type ( TimerForm ) :: &
    T
  type ( CollectiveOperation_R_Form ) :: &
    CO
    
  allocate ( PROGRAM_HEADER )
  associate ( PH => PROGRAM_HEADER )
  
  call PH % Initialize ( 'DataManagement_Test' )
  
  call T % Initialize ( 'Timer', 0 )
  
  call CO % Initialize &
         ( PH % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )  
  
  call S % Initialize ( ValueShape = [ 256 ** 3, 64 ], PinnedOption = .true. )
  call S % AllocateDevice ( )
  
  call random_number ( S % Value )
  
  DataSize_GiB = 1.0_KDR * S % nValues * S % nVariables * 8 &
                   / ( 1024.0_KDR ** 3 )
  
  call Show ( DataSize_GiB, 'Data size (GiB)' )
  
  call T % Start ( )
  call S % UpdateDevice ( )
  call T % Stop ( )
  Bandwidth = DataSize_GiB / T % TimeInterval % Number
  CO % Outgoing % Value ( 1 ) = Bandwidth
  call CO % Reduce ( REDUCTION % MAX )
  BandwidthMax = CO % Incoming % Value ( 1 )
  call CO % Reduce ( REDUCTION % MIN )
  BandwidthMin = CO % Incoming % Value ( 1 )
  Diff = ( BandwidthMax - BandwidthMin ) / BandwidthMax
  
  call Show ( S % ErrorDevice,  'Device Error' )
  call Show ( T % TimeInterval, 'H-to-D Time (s)' )
  call Show ( Bandwidth,        'H-to-D Bandwidth Process (GiB/s)' )
  call Show ( BandwidthMin,     'H-to-D Bandwidth Min     (GiB/s)' )
  call Show ( BandwidthMax,     'H-to-D Bandwidth Max     (GiB/s)' )
  call Show ( Diff,             'H-to-D Bandwidth Variation' )
  

  call T % Start ( )
  call S % UpdateHost ( )
  call T % Stop ( )
  
  Bandwidth = DataSize_GiB / T % TimeInterval % Number
  CO % Outgoing % Value ( 1 ) = Bandwidth
  call CO % Reduce ( REDUCTION % MAX )
  BandwidthMax = CO % Incoming % Value ( 1 )
  call CO % Reduce ( REDUCTION % MIN )
  BandwidthMin = CO % Incoming % Value ( 1 )
  Diff = ( BandwidthMax - BandwidthMin ) / BandwidthMax
  
  call Show ( S % ErrorDevice,  'Device Error' )
  call Show ( T % TimeInterval, 'D-to-H Time (s)' )
  call Show ( Bandwidth,        'D-to-H Bandwidth Process (GiB/s)' )
  call Show ( BandwidthMin,     'D-to-H Bandwidth Min     (GiB/s)' )
  call Show ( BandwidthMax,     'D-to-H Bandwidth Max     (GiB/s)' )
  call Show ( Diff,             'H-to-D Bandwidth Variation' )
  
  end associate !-- PH

end program DataManagement_Test
