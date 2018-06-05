program GetMemoryUsage_Command_Test

  use VariableManagement
  use Display
  use MessagePassing
  use CommandLineOptions_Form
  use GetMemoryUsage_Command

  implicit none
  
  integer ( KDI ) :: &
    iStrg, &
    nValues, &
    DisplayRank
  type ( MeasuredValueForm ) :: &
    HighWaterMark, &
    AcrossProcessesMinHighWaterMark, &
    AcrossProcessesMaxHighWaterMark, &
    AcrossProcessesMeanHighWaterMark, &
    ResidentSetSize, &
    AcrossProcessesMinResidentSetSize, &
    AcrossProcessesMaxResidentSetSize, &
    AcrossProcessesMeanResidentSetSize
  type ( StorageForm ), dimension ( 10 ) :: &
    Dummy
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( CommandLineOptionsForm ) :: &
    CLO

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  
  call CLO % Initialize ( ) 
  
  DisplayRank = 0
  call CLO % Read ( DisplayRank, 'DisplayRank' )
  call CONSOLE % SetDisplayRank ( DisplayRank )
  
  !-- dummy allocation to eat some memory
  
  nValues = 100
  call CLO % Read ( nValues, 'nValues' )
  call Show ( nValues, 'nValues' )
  
  do iStrg = 1, size ( Dummy ) 
    call Dummy ( iStrg ) % Initialize ( [ nValues, C % Rank + 1 ] )
    call random_number ( Dummy ( iStrg ) % Value )
  end do
  
  !-- get memory usage
  
  call Show ( 'Getting memory usage', CONSOLE % INFO_1 )
  call GetMemoryUsage &
         ( HighWaterMark, ResidentSetSize, CONSOLE % INFO_1, C_Option = C, &
           Max_HWM_Option = AcrossProcessesMaxHighWaterMark, &
           Min_HWM_Option = AcrossProcessesMinHighWaterMark, &
           Mean_HWM_Option = AcrossProcessesMeanHighWaterMark, &
           Max_RSS_Option = AcrossProcessesMaxResidentSetSize, &
           Min_RSS_Option = AcrossProcessesMinResidentSetSize, &
           Mean_RSS_Option = AcrossProcessesMeanResidentSetSize )

  call Show ( HighWaterMark, 'High water mark' )
  call Show ( ResidentSetSize, 'Resident set size' )
  
  call Show ( 'Across processes memory usage', CONSOLE % INFO_1 )
  call Show &
         ( AcrossProcessesMaxHighWaterMark, &
           'Max high water mark' )
  call Show &
         ( AcrossProcessesMinHighWaterMark, &
           'Min high water mark' )
  call Show &
         ( AcrossProcessesMeanHighWaterMark, &
           'Mean high water mark' )
           
  
  call Show &
         ( AcrossProcessesMaxResidentSetSize, &
           'Max resident set size' )
  call Show &
         ( AcrossProcessesMinResidentSetSize, &
           'Min resident set size' )
  call Show &
         ( AcrossProcessesMeanResidentSetSize, &
           'Mean resident set size' )

  deallocate ( C )

end program GetMemoryUsage_Command_Test
