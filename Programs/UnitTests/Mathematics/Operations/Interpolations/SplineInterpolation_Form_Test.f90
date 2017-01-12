program SplineInterpolation_Form_Test

  use Basics
  use SplineInterpolation_Form
  
  implicit none
  
  character ( LDL ), parameter :: &
    ProgramName = 'SplineInterpolation_Form_Test'
  
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
  type ( VariableGroupForm ) :: &
    VG_Known, &
    VG_Interpolated
  type ( GridImageStreamForm ) :: &
    GIS
  type ( TableStreamForm ) :: &
    TS
  type ( CurveImageForm ) :: &
    CI_Known, &
    CI_Interpolated
  type ( SplineInterpolationForm ) :: &
    SI
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )
  
  !-- read table data of value

  call PROGRAM_HEADER % GetParameter ( DataFile, 'DataFile' )
  
  call TS % Initialize ( DataFile, PROGRAM_HEADER % Communicator % Rank )
  call TS % Read ( Table )
  
!-- FIXME: NAG doesn't like sourced allocation but it should be supported
!  allocate ( KnownInput, source = Table ( :, 1 ) )
  allocate ( KnownInput ( size ( Table ( :, 1 ) ) ) )
  KnownInput = Table ( :, 1 )
  
  !-- initialize spline interpolation

  call SI % Initialize ( Table ( :, 1 ), Table ( :, 2 ) )
  
  !-- initialize storage for interpolated values, and interpolate

  allocate ( InterpolatedInput ( 5 * size ( KnownInput ) ) )
  call VG_Interpolated % Initialize &
         ( [ size ( InterpolatedInput ), 1 ], &
           VariableOption = [ 'Interpolated' ] )
  Delta &
    = ( maxval ( KnownInput ) - minval ( KnownInput ) ) &
      / VG_Interpolated % nValues 
  
  do iV = 1, VG_Interpolated % nValues
    InterpolatedInput ( iV ) = KnownInput ( 1 ) + ( iV - 1 ) * Delta 
    call SI % Evaluate &
           ( InterpolatedInput ( iV ), VG_Interpolated % Value ( iV, 1 ) )
  end do 
  
  call Show ( Table, 'Original value', CONSOLE % INFO_5 )
  call Show &
         ( reshape &
             ( [ InterpolatedInput, VG_Interpolated % Value ( :, 1 ) ], &
               [ VG_Interpolated % nValues, 2 ] ), 'Interpolated value', &
           CONSOLE % INFO_5 )

  !-- write
  call VG_Known % Initialize &
         ( [ size ( KnownInput ), 1 ], VariableOption = [ 'Known' ] )
  VG_Known % Value ( :, 1 ) = Table ( :, 2 )
  
  call GIS % Initialize ( ProgramName )
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call CI_Known % Initialize ( GIS )
  call CI_Known % AddVariableGroup ( VG_Known )
  call CI_Known % SetGrid &
         ( Directory = '', NodeCoordinate = KnownInput, &
           nProperCells = size ( KnownInput ), oValue = 0 )
  call CI_Known % Write ( )
  
  call CI_Interpolated % Initialize ( GIS )
  call CI_Interpolated % AddVariableGroup ( VG_Interpolated )
  call CI_Interpolated % SetGrid &
         ( Directory = '', NodeCoordinate = InterpolatedInput, &
           nProperCells = size ( InterpolatedInput ), oValue = 0 )
  call CI_Interpolated % Write ( )
  
  call GIS % Close ( )

  deallocate ( PROGRAM_HEADER )

end program SplineInterpolation_Form_Test
