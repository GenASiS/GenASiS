program ReadLabelValue_Command_Test

  use ISO_FORTRAN_ENV
  use Specifiers
  use Display
  use MessagePassing
  use ReadLabelValue_Command

  implicit none

  integer ( KDI ) :: &
    nValues, &
    ScalarIntegerValue
  integer ( KDI), dimension ( 10 ) :: & 
    ArrayIntegerValue
  real ( KDR ) :: &
    ScalarRealValue, &
    ScalarRealUnitValue
  real ( KDR ), dimension ( 10 ) :: &
    ArrayRealValue, &
    ArrayRealUnitValue
  type ( QuantityForm ) :: &
    ScalarQuantity, &
    InputUnit
  type ( QuantityForm ), dimension ( 10 ) :: &
    ArrayQuantity, &
    ArrayInputUnit
  logical ( KDL ) :: &
    ScalarLogicalValue, &
    Success
  logical ( KDL ), dimension ( 10 ) :: &
    ArrayLogicalValue
  character ( 5 ) :: &
    Encoding
  character ( LDL ) :: &
    ScalarStringValue, &
    Label
  character ( LDL ), dimension ( 10 ) :: &
    ArrayStringValue
  character ( LDB ) :: &
    CommentedString, &
    EmptyString, &
    ScalarIntegerString, &
    ScalarRealString, &
    ScalarRealUnitString, &
    ScalarQuantityString, &
    ScalarLogicalString, &
    ScalarStringString, &
    ArrayIntegerString, &
    ArrayRealString, &
    ArrayRealUnitString, &
    ArrayQuantityString, &
    ArrayLogicalString, &
    ArrayStringString
  type ( CommunicatorForm ), allocatable :: &
    C

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )

!-- Runtime error with CCE
!  if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!    open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!  else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
  if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    Encoding = 'UTF-8'
    open ( OUTPUT_UNIT, encoding = Encoding )
  end if

  call UNIT % Initialize ( )
    
  CommentedString           = '!-- Commented string should not be read'
  ScalarIntegerString       = 'ScalarInteger=10'
  ScalarRealString          = 'ScalarReal=10.1'
  ScalarRealUnitString      = 'ScalarRealUnit=10.0~KILOMETER'
  ScalarQuantityString = 'ScalarQuantity=20.0~KILOMETER'
  ScalarLogicalString       = 'ScalarLogical=True'
  ScalarStringString        = 'ScalarString=Hello World'
  ArrayIntegerString        = 'ArrayInteger=10,20,30'
  ArrayRealString           = 'ArrayReal=10.5,20.5,30.5,40.5,50.5'
  ArrayRealUnitString  &
    = 'ArrayRealUnit=10.5~SECOND,20.5,30.5~MILLISECOND,40.5,50.5~MILLISECOND'
  ArrayQuantityString  &
    = 'ArrayQuantity=10.5~SECOND,20.5,30.5~MILLISECOND,40.5,' &
      // '50.5~MILLISECOND'
  ArrayLogicalString        = 'ArrayLogical=T,F,False,True,True,F'
  ArrayStringString         = 'ArrayString=Lorem, ipsum, dolor, sit, amet'
  
  Label = ''
  EmptyString = ''
  call ReadLabelValue &
         ( Label, EmptyString, CommentedString, SuccessOption = Success )
  call Show ( Success, 'Is comment parsed ?' )
  
  call ReadLabelValue ( Label, ScalarIntegerValue, ScalarIntegerString )
  call Show ( ScalarIntegerValue, Label )
  
  call ReadLabelValue ( Label, ScalarRealValue, ScalarRealString )
  call Show ( ScalarRealValue, Label )
  
  call ReadLabelValue &
         ( Label, ScalarRealUnitValue, ScalarRealUnitString, &
           InputUnitOption = InputUnit )
  call Show ( ScalarRealUnitValue, InputUnit, trim ( Label ) &
              // ' input units' )
  call Show ( ScalarRealUnitValue, trim ( Label ) // ', program units' )
  
  call ReadLabelValue &
         ( Label, ScalarQuantity, ScalarQuantityString, &
           InputUnitOption = InputUnit )
  call Show ( ScalarQuantity, InputUnit, trim ( Label ) // &
              ' input units' )
  call Show ( ScalarQuantity, trim ( Label ) // ', program units' )
  
  call ReadLabelValue ( Label, ScalarLogicalValue, ScalarLogicalString )
  call Show ( ScalarLogicalValue, Label )
  
  call ReadLabelValue ( Label, ScalarStringValue, ScalarStringString )
  call Show ( ScalarStringValue, Label )
  
  call ReadLabelValue &
         ( Label, ArrayIntegerValue, ArrayIntegerString, &
           nValuesOption = nValues )
  call Show ( ArrayIntegerValue ( : nValues ), Label )
  
  call ReadLabelValue &
         ( Label, ArrayRealValue, ArrayRealString, &
           nValuesOption = nValues )
  call Show ( ArrayRealValue ( : nValues ), Label )
  
  call ReadLabelValue &
         ( Label, ArrayRealUnitValue, ArrayRealUnitString, &
           InputUnitOption = ArrayInputUnit, nValuesOption = nValues )
  call Show &
         ( ArrayRealUnitValue ( : nValues ), ArrayInputUnit ( : nValues ), &
           trim ( Label ) // ', input units' )
  call Show &
         ( ArrayRealUnitValue ( : nValues ), &
           trim ( Label ) // ', program units' )
  
  call ReadLabelValue &
         ( Label, ArrayQuantity, ArrayQuantityString, &
           InputUnitOption = ArrayInputUnit, nValuesOption = nValues )
  call Show &
         ( ArrayQuantity ( : nValues ), ArrayInputUnit ( : nValues ), &
           trim ( Label ) // ', input units' )
  call Show &
         ( ArrayQuantity ( : nValues ), &
           trim ( Label ) // ', program units' )
  
  call ReadLabelValue &
         ( Label, ArrayLogicalValue, ArrayLogicalString, &
           nValuesOption = nValues )
  call Show ( ArrayLogicalValue ( : nValues ), Label )
  
  call ReadLabelValue &
         ( Label, ArrayStringValue, ArrayStringString, &
           nValuesOption = nValues )
  call Show ( ArrayStringValue ( : nValues ), Label )
  
  !-- Test for retaining default value
  ScalarIntegerValue = 42
  Label = 'Default label'
  call ReadLabelValue ( Label, ScalarIntegerValue, 'Invalid string' )
  call Show ( ScalarIntegerValue, Label )
  
  ScalarRealValue = 42.0
  call ReadLabelValue ( Label, ScalarRealValue, 'Invalid string' )
  call Show ( ScalarRealValue, Label )
  
  ScalarStringValue = 'Forty-two'
  call ReadLabelValue ( Label, ScalarStringValue, 'Invalid string' )
  call Show ( ScalarStringValue, Label )
  
  ArrayRealValue ( 1 : 3 ) = [ 100.0_KDR, 200.0_KDR, 300.0_KDR ]
  call ReadLabelValue ( Label, ArrayRealValue, 'Invalid string' )
  call Show ( ArrayRealValue ( 1 : 3 ), Label ) 
  
  ArrayStringValue ( 1 ) = 'Hello'
  ArrayStringValue ( 2 ) = 'World'
  call ReadLabelValue ( Label, ArrayStringValue, 'Invalid string' )
  call Show ( ArrayStringValue ( : 2 ), Label ) 
  
  deallocate ( C )

end program ReadLabelValue_Command_Test
