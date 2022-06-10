program Interpolation_Form_Test

  use Basics
  use Calculus
  
  implicit none
  
  integer ( KDI ) :: &
    iV
  real ( KDR ) :: &
    Delta
  real ( KDR ), dimension ( : ), allocatable :: &
    KnownInput, &
    InterpolatedInput
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Table
  character ( LDF ) :: &
    DataFile = ''
  type ( StorageForm ) :: &
    S_Known, &
    S_Interpolated
  type ( GridImageStreamForm ) :: &
    GIS
  type ( TableStreamForm ) :: &
    TS
  type ( CurveImageForm ) :: &
    CI_Known, &
    CI_Interpolated
  type ( InterpolationForm ) :: &
    I
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Interpolation_Form_Test', AppendDimensionalityOption = .false. )
  
  !-- read table of data

  call PROGRAM_HEADER % GetParameter ( DataFile, 'DataFile' )
  
  call TS % Initialize ( DataFile, PROGRAM_HEADER % Communicator % Rank )
  call TS % Read ( Table )
  
  allocate ( KnownInput, source = Table ( :, 1 ) )
  
  !-- initialize spline interpolation

  call I % Initialize ( Table ( :, 1 ), Table ( :, 2 ) )
  
  !-- initialize storage for interpolated values, and interpolate

  allocate ( InterpolatedInput ( 5 * size ( KnownInput ) ) )
  call S_Interpolated % Initialize &
         ( [ size ( InterpolatedInput ), 1 ], &
           VariableOption = [ 'Interpolated' ] )
  Delta &
    = ( maxval ( KnownInput ) - minval ( KnownInput ) ) &
      / S_Interpolated % nValues 
  
  do iV = 1, S_Interpolated % nValues
    InterpolatedInput ( iV ) = KnownInput ( 1 ) + ( iV - 1 ) * Delta 
    call I % Evaluate &
           ( InterpolatedInput ( iV ), S_Interpolated % Value ( iV, 1 ) )
  end do 
  
  call Show ( Table, 'Original value', CONSOLE % INFO_5 )
  call Show ( reshape ( [ InterpolatedInput, &
                          S_Interpolated % Value ( :, 1 ) ], &
              [ S_Interpolated % nValues, 2 ] ), 'Interpolated value', &
              CONSOLE % INFO_5 )

  !-- write

  call S_Known % Initialize &
         ( [ size ( KnownInput ), 1 ], VariableOption = [ 'Known' ] )
  S_Known % Value ( :, 1 ) = Table ( :, 2 )
  
  call GIS % Initialize ( PROGRAM_HEADER % Name )
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call CI_Known % Initialize ( GIS )
  call CI_Known % AddStorage ( S_Known )
  call CI_Known % SetGridWrite &
         ( Directory = '', NodeCoordinate = KnownInput, &
           nProperCells = size ( KnownInput ), oValue = 0 )
  call CI_Known % Write ( )
  
  call CI_Interpolated % Initialize ( GIS )
  call CI_Interpolated % AddStorage ( S_Interpolated )
  call CI_Interpolated % SetGridWrite &
         ( Directory = '', NodeCoordinate = InterpolatedInput, &
           nProperCells = size ( InterpolatedInput ), oValue = 0 )
  call CI_Interpolated % Write ( )
  
  call GIS % Close ( )

  deallocate ( PROGRAM_HEADER )

end program Interpolation_Form_Test
