program CurveImage_Form_Test

  use ISO_FORTRAN_ENV
  use Specifiers
  use DataManagement
  use Display
  use MessagePassing
  use GridImageStream_Form
  use CurveImage_Form
  
  implicit none
  
  integer ( KDI ) :: &
    iV, &
    iC, &
    oC_1, &
    DisplayRank = 0
  real ( KDR ) :: &
    L1_Error  
  real ( KDR ), dimension ( 20 ) :: &
    NodeCoordinate
  character ( 5 ) :: &
    Encoding
  character ( LDF ) :: &
    Name = 'CurveImage_Form_Test'
  character ( LDL ), dimension ( 1 ) :: &
    VariableName
  type ( StorageForm ) :: &
    S, &
    S_SO
  type ( CommunicatorForm ), allocatable :: &
    C
  type ( GridImageStreamForm ) :: &
    GIS
  type ( CurveImageForm ) :: & 
    CI_Write, &
    CI_Read, &
    CI_Read_SO  !-- Read_StorageOnly
  
!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  call UNIT % Initialize ( )

  allocate ( C )
  call C % Initialize ( )
  
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetDisplayRank ( DisplayRank )

  NodeCoordinate =  [ ( iC * 0.5_KDR, iC = 1, 20 ) ]
  
  oC_1 = C % Rank * maxval ( NodeCoordinate )
  NodeCoordinate = NodeCoordinate  + oC_1
  
  VariableName ( 1 ) = 'VariableName'
  call S % Initialize &
         ( [ 20, 1 ], VariableOption = VariableName, &
           NameOption = 'Storage', &
           UnitOption = [ UNIT % METER ] )
  S % Value ( :, 1 ) &
    = [ ( ( iC + 20 * C % Rank ) * 2.0_KDR, iC = 1, 20 ) ] 
    
  call Show ( NodeCoordinate, 'NodeCoordinate' )
  call Show ( S % Value ( :, 1 ), S % Unit ( 1 ), 'Variable' )
  call Show ( S % Unit, 'Unit' )

  call GIS % Initialize ( Name, CommunicatorOption = C )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  
  call CI_Write % Initialize ( GIS )
  
  call CI_Write % AddStorage ( S )

  call CI_Write % SetGridWrite  &
         ( Directory = 'Curves', NodeCoordinate = NodeCoordinate, &
           nProperCells = 20, oValue = 0 )

  call CI_Write % Write ()

  call GIS % Close ( ) 
  
  
  !-- Test Read
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call CI_Read % Initialize ( GIS )
  
  call CI_Read % SetGridRead &
         ( Directory = 'Curves', nProperCells = 20, oValue = 0 )

  call CI_Read % Read ( )
  
  !-- Note unit not read!
  CI_Read % Storage ( 1 ) % Unit = S % Unit
  do iV = 1, CI_Read % Storage ( 1 ) % nVariables
    CI_Read % Storage ( 1 ) % Value ( :, iV ) &
      = CI_Read % Storage ( 1 ) % Value ( :, iV ) &
          * CI_Read % Storage ( 1 ) % Unit ( iV ) % Number
  end do

  call Show ( CI_Read % NodeCoordinate_1, 'NodeCoordinate_1' )
  call Show ( CI_Read % Storage ( 1 ) % Value ( :, 1 ), &
              CI_Read % Storage ( 1 ) % Unit ( 1 ), &
              CI_Read % Storage ( 1 ) % Variable ( 1 ) )
  call Show ( CI_Read % Storage ( 1 ) % Unit, 'Unit' )
  
  call GIS % Close ( )

  associate &
    ( V_0 => CI_Write % Storage ( 1 ) % Value ( :, 1 ), &
      V_1 => CI_Read  % Storage ( 1 ) % Value ( :, 1 ) )
  L1_Error = sum ( abs ( V_1 - V_0 ) ) / sum ( abs ( V_0 ) )
  call Show ( L1_Error, 'L1_Error Read - Write' )
  call Show ( L1_Error <= epsilon ( V_0 ( 1 ) ) * 10.0, 'Read == Write ?' )
  end associate
  
  
  !-- Test Read StorageOnly
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call CI_Read_SO % Initialize ( GIS )

  call S_SO % Initialize &
         ( shape ( S % Value ), VariableOption = S % Variable, &
           NameOption = S % Name, UnitOption = S % Unit )
  call CI_Read_SO % SetGridRead &
         ( Directory = 'Curves', nProperCells = 20, oValue = 0 )
  call CI_Read_SO % AddStorage ( S_SO )
!  call CI_Read_SO % SetGrid &
!         ( Directory = 'Curves', NodeCoordinate = NodeCoordinate, &
!           nProperCells = 20, oValue = 0 )

  call CI_Read_SO % Read ( StorageOnlyOption = .true. )

  call GIS % Close ( )

  associate &
    ( V_0 => CI_Write % Storage ( 1 ) % Value ( :, 1 ), &
      V_1 => CI_Read_SO  % Storage ( 1 ) % Value ( :, 1 ) )
  L1_Error = sum ( abs ( V_1 - V_0 ) ) / sum ( abs ( V_0 ) )
  call Show ( L1_Error, 'L1_Error Read_SO - Write' )
  call Show ( L1_Error <= epsilon ( V_0 ( 1 ) ) * 10.0, 'Read_SO == Write ?' )
  end associate

  deallocate ( C )
  
end program CurveImage_Form_Test
