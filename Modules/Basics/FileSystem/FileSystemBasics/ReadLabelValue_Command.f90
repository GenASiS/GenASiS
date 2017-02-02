!-- Given a string in the general form label=value, ReadLabelValue parses the 
!   string and returns the label and its value as separate variables

module ReadLabelValue_Command
  
  use VariableManagement

  implicit none
  private 
  
  
  public :: &
    ReadLabelValue
  
  interface ReadLabelValue
    module procedure ReadLabelValueInteger_0D
    module procedure ReadLabelValueReal_0D
    module procedure ReadLabelValueMeasuredValue_0D
    module procedure ReadLabelValueLogical_0D
    module procedure ReadLabelValueCharacter_0D
    module procedure ReadLabelValueInteger_1D
    module procedure ReadLabelValueReal_1D
    module procedure ReadLabelValueMeasuredValue_1D
    module procedure ReadLabelValueLogical_1D
    module procedure ReadLabelValueCharacter_1D
  end interface ReadLabelValue

    character ( 1 ), private, parameter :: &
      DELIMITER_VALUE     = '=', &
      DELIMITER_ARRAY     = ',', &
      DELIMITER_COMMENT_1 = '!', &
      DELIMITER_COMMENT_2 = '#', &
      DELIMITER_UNIT      = '~'
      
contains
    

  subroutine ReadLabelValueInteger_0D ( Label, Value, Buffer, SuccessOption )
    
    character ( * ), intent ( inout ) :: &
      Label
    integer ( KDI ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      ValueMarker, &
      Status
    logical ( KDL ) :: & 
      IsValid
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a8)' ) '(a', ValueMarker, ',a1,i10)'
    read ( Buffer, fmt = ReadFormat, iostat = Status ) &
      Label, Delimiter, Value 

    if ( Status /= 0 .and. present ( SuccessOption ) ) then
      SuccessOption = .false.
      return
    end if
    
    if ( present ( SuccessOption ) ) SuccessOption = .true.
     
  end subroutine ReadLabelValueInteger_0D

  
  subroutine ReadLabelValueReal_0D &
               ( Label, Value, Buffer, InputUnitOption, SuccessOption )
    
    character ( * ), intent ( inout ) :: &
      Label
    real ( KDR ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      ValueMarker, &
      UnitMarker, &
      Status
    type ( MeasuredValueForm ) :: &
      InputUnit
    logical ( KDL ) :: & 
      IsValid
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      StrippedBuffer

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if

    UnitMarker = index ( Buffer, DELIMITER_UNIT, back = .true. )
    if ( UnitMarker > 0 ) then
      call UNIT % Select ( trim ( Buffer ( UnitMarker + 1 : ) ), InputUnit )
      StrippedBuffer = Buffer ( : UnitMarker - 1 )
    else 
      InputUnit = UNIT % IDENTITY
      StrippedBuffer = Buffer
    end if
    if ( present ( InputUnitOption ) ) InputUnitOption = InputUnit

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a11)' ) '(a', ValueMarker, ',a1,f30.20)'
    read ( StrippedBuffer, fmt = ReadFormat, iostat = Status ) &
      Label, Delimiter, Value 
    Value = Value * InputUnit % Number

    if ( Status /= 0 .and. present ( SuccessOption ) ) then
      SuccessOption = .false.
      return
    end if
    
    if ( present ( SuccessOption ) ) SuccessOption = .true.
     
  end subroutine ReadLabelValueReal_0D

  
  subroutine ReadLabelValueMeasuredValue_0D &
               ( Label, Value, Buffer, InputUnitOption, SuccessOption )
    
    character ( * ), intent ( inout ) :: &
      Label
    type ( MeasuredValueForm ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      ValueMarker, &
      UnitMarker, &
      Status
    type ( MeasuredValueForm ) :: &
      InputUnit
    logical ( KDL ) :: & 
      IsValid
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      StrippedBuffer

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if

    UnitMarker = index ( Buffer, DELIMITER_UNIT, back = .true. )
    if ( UnitMarker > 0 ) then
      call UNIT % Select ( trim ( Buffer ( UnitMarker + 1 : ) ), InputUnit )
      StrippedBuffer = Buffer ( : UnitMarker - 1 )
    else 
      InputUnit = UNIT % IDENTITY
      StrippedBuffer = Buffer
    end if
    if ( present ( InputUnitOption ) ) InputUnitOption = InputUnit

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a11)' ) '(a', ValueMarker, ',a1,f30.20)'
    read ( StrippedBuffer, fmt = ReadFormat, iostat = Status ) &
      Label, Delimiter, Value % Number
    Value = Value % Number * InputUnit

    if ( Status /= 0 .and. present ( SuccessOption ) ) then
      SuccessOption = .false.
      return
    end if
    
    if ( present ( SuccessOption ) ) SuccessOption = .true.
     
  end subroutine ReadLabelValueMeasuredValue_0D

  
  subroutine ReadLabelValueLogical_0D ( Label, Value, Buffer, SuccessOption )
    
    character ( * ), intent ( inout ) :: &
      Label
    logical ( KDL ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      ValueMarker, &
      Status
    logical ( KDL ) :: & 
      IsValid
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a11)' ) '(a', ValueMarker, ',a1,l10)'
    read ( Buffer, fmt = ReadFormat, iostat = Status ) &
      Label, Delimiter, Value 

    if ( Status /= 0 .and. present ( SuccessOption ) ) then
      SuccessOption = .false.
      return
    end if
    
    if ( present ( SuccessOption ) ) SuccessOption = .true.
     
  end subroutine ReadLabelValueLogical_0D

  
  subroutine ReadLabelValueCharacter_0D &
               ( Label, Value, Buffer, SuccessOption )
    
    character ( * ), intent ( inout ) :: &
      Label
    character ( * ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      ValueMarker, &
      Status
    logical ( KDL ) :: & 
      IsValid
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a5,i4,a1)' ) &
      '(a', ValueMarker, ',a1,a,', len ( Value ), ')'
    read ( Buffer, fmt = ReadFormat, iostat = Status ) &
      Label, Delimiter, Value 

    if ( Status /= 0 .and. present ( SuccessOption ) ) then
      SuccessOption = .false.
      return
    end if
    
    if ( present ( SuccessOption ) ) SuccessOption = .true.
     
  end subroutine ReadLabelValueCharacter_0D

  
  subroutine ReadLabelValueInteger_1D &
               ( Label, Value, Buffer, nValuesOption, SuccessOption )

    character ( * ), intent ( inout ) :: &
      Label
    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      iV, &  !-- iValue
      oChar, &
      nValues, &
      Marker, &
      ValueMarker, &
      Status
    logical ( KDL ) :: & 
      IsValid, &
      Success
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      JoinedValue, &
      SingleValue

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if 
    
    Success = .true. 
    
    !-- Parse Buffer

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a5,i4,a1)' ) &
      '(a', ValueMarker, ',a1,a', LDB, ')'
    read ( Buffer, fmt = ReadFormat ) Label, Delimiter, JoinedValue 
    
    !-- Find out how many values there are
    
    nValues = 1
    oChar = 0
    do
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      if ( Marker > 0 )then
        nValues = nValues + 1
      else
        exit
      end if
      oChar = oChar + Marker
    end do

    if ( present ( nValuesOption ) ) nValuesOption = nValues
    
    !-- Read values
    oChar = 0
    do iV = 1, min ( size ( Value ), nValues )
      
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      
      if ( Marker > 0 ) then
        SingleValue = JoinedValue ( oChar + 1 : oChar + Marker - 1 )
      else
        SingleValue = JoinedValue ( oChar + 1 : )
      end if

      !-- Read value

      read ( SingleValue, fmt = '(i10)', iostat = Status ) Value ( iV )
        
      if ( Status /= 0 ) then
        Success = .false.
        exit 
      end if 
      oChar = oChar + Marker
    
    end do
    
    if ( present ( SuccessOption ) ) SuccessOption = Success
     
  end subroutine ReadLabelValueInteger_1D


  subroutine ReadLabelValueReal_1D &
               ( Label, Value, Buffer, InputUnitOption, nValuesOption, &
                 SuccessOption )

    character ( * ), intent ( inout ) :: &
      Label
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      iV, &  !-- iValue
      oChar, &
      nValues, &
      Marker, &
      ValueMarker, &
      UnitMarker, &
      Status
    type ( MeasuredValueForm ) :: &
      InputUnit
    logical ( KDL ) :: & 
      IsValid, &
      Success
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      JoinedValue, &
      SingleValue, &
      StrippedValue 

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if 
    
    Success = .true. 
    
    !-- Parse Buffer

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a5,i4,a1)' ) &
      '(a', ValueMarker, ',a1,a', LDB, ')'
    read ( Buffer, fmt = ReadFormat ) Label, Delimiter, JoinedValue 
    
    !-- Find out how many values there are
    
    nValues = 1
    oChar = 0
    do
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      if ( Marker > 0 )then
        nValues = nValues + 1
      else
        exit
      end if
      oChar = oChar + Marker
    end do

    if ( present ( nValuesOption ) ) nValuesOption = nValues

    !-- Read values
    oChar = 0
    do iV = 1, min ( size ( Value ), nValues )
      
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      
      if ( Marker > 0 ) then
        SingleValue = JoinedValue ( oChar + 1 : oChar + Marker - 1 )
      else
        SingleValue = JoinedValue ( oChar + 1 : )
      end if

      !-- Read unit and strip it from SingleValue

      UnitMarker = index ( SingleValue, DELIMITER_UNIT, back = .true. )
      if ( UnitMarker > 0 ) then
        call UNIT % Select &
               ( trim ( SingleValue ( UnitMarker + 1 : ) ), InputUnit )
        StrippedValue = SingleValue ( : UnitMarker - 1 )
      else 
        InputUnit = UNIT % IDENTITY
        StrippedValue = SingleValue
      end if
      if ( present ( InputUnitOption ) ) InputUnitOption ( iV ) = InputUnit

      !-- Read value

      read ( StrippedValue, fmt = '(f30.20)', iostat = Status ) Value ( iV )
      Value ( iV ) = Value ( iV ) * InputUnit % Number
        
      if ( Status /= 0 ) then
        Success = .false.
        exit 
      end if 
      oChar = oChar + Marker
    
    end do
    
    if ( present ( SuccessOption ) ) SuccessOption = Success
     
  end subroutine ReadLabelValueReal_1D


  subroutine ReadLabelValueMeasuredValue_1D &
               ( Label, Value, Buffer, InputUnitOption, nValuesOption, &
                 SuccessOption )

    character ( * ), intent ( inout ) :: &
      Label
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), &
      optional :: &
        InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      iV, &  !-- iValue
      oChar, &
      nValues, &
      Marker, &
      ValueMarker, &
      UnitMarker, &
      Status
    type ( MeasuredValueForm ) :: &
      InputUnit
    logical ( KDL ) :: & 
      IsValid, &
      Success
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      JoinedValue, &
      SingleValue, &
      StrippedValue 

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if 
    
    Success = .true. 
    
    !-- Parse Buffer

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a5,i4,a1)' ) &
      '(a', ValueMarker, ',a1,a', LDB, ')'
    read ( Buffer, fmt = ReadFormat ) Label, Delimiter, JoinedValue 
    
    !-- Find out how many values there are
    
    nValues = 1
    oChar = 0
    do
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      if ( Marker > 0 )then
        nValues = nValues + 1
      else
        exit
      end if
      oChar = oChar + Marker
    end do

    if ( present ( nValuesOption ) ) nValuesOption = nValues
    
    !-- Read values
    oChar = 0
    do iV = 1, min ( size ( Value ), nValues )
      
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      
      if ( Marker > 0 ) then
        SingleValue = JoinedValue ( oChar + 1 : oChar + Marker - 1 )
      else
        SingleValue = JoinedValue ( oChar + 1 : )
      end if

      !-- Read unit and strip it from SingleValue

      UnitMarker = index ( SingleValue, DELIMITER_UNIT, back = .true. )
      if ( UnitMarker > 0 ) then
        call UNIT % Select &
               ( trim ( SingleValue ( UnitMarker + 1 : ) ), InputUnit )
        StrippedValue = SingleValue ( : UnitMarker - 1 )
      else 
        InputUnit = UNIT % IDENTITY
        StrippedValue = SingleValue
      end if
      if ( present ( InputUnitOption ) ) InputUnitOption ( iV ) = InputUnit

      !-- Read value

      read ( StrippedValue, fmt = '(f30.20)', iostat = Status ) &
        Value ( iV ) % Number
      Value ( iV ) = Value ( iV ) % Number * InputUnit
        
      if ( Status /= 0 ) then
        Success = .false.
        exit 
      end if 
      oChar = oChar + Marker
    
    end do
    
    if ( present ( SuccessOption ) ) SuccessOption = Success
     
  end subroutine ReadLabelValueMeasuredValue_1D


  subroutine ReadLabelValueLogical_1D &
               ( Label, Value, Buffer, nValuesOption, SuccessOption )

    character ( * ), intent ( inout ) :: &
      Label
    logical ( KDL ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      iV, &  !-- iValue
      oChar, &
      nValues, &
      Marker, &
      ValueMarker, &
      Status
    logical ( KDL ) :: & 
      IsValid, &
      Success
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      JoinedValue, &
      SingleValue

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if 
    
    Success = .true. 
    
    !-- Parse Buffer

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a5,i4,a1)' ) &
      '(a', ValueMarker, ',a1,a', LDB, ')'
    read ( Buffer, fmt = ReadFormat ) Label, Delimiter, JoinedValue 
    
    !-- Find out how many values there are
    
    nValues = 1
    oChar = 0
    do
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      if ( Marker > 0 ) then
        nValues = nValues + 1
      else
        exit
      end if
      oChar = oChar + Marker
    end do

    if ( present ( nValuesOption ) ) nValuesOption = nValues
    
    !-- Read values
    oChar = 0
    do iV = 1, min ( size ( Value ), nValues )
      
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      
      if ( Marker > 0 ) then
        SingleValue = JoinedValue ( oChar + 1 : oChar + Marker - 1 )
      else
        SingleValue = JoinedValue ( oChar + 1 : )
      end if

      !-- Read value

      read ( SingleValue, fmt = '(l10)', iostat = Status ) Value ( iV )
        
      if ( Status /= 0 ) then
        Success = .false.
        exit 
      end if 
      oChar = oChar + Marker
    
    end do
    
    if ( present ( SuccessOption ) ) SuccessOption = Success
     
  end subroutine ReadLabelValueLogical_1D


  subroutine ReadLabelValueCharacter_1D &
               ( Label, Value, Buffer, nValuesOption, SuccessOption )

    character ( * ), intent ( inout ) :: &
      Label
    character ( * ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Buffer
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
      
    integer ( KDI ) :: &
      iV, &  !-- iValue
      oChar, &
      nValues, &
      Marker, &
      ValueMarker, &
      Status
    logical ( KDL ) :: & 
      IsValid, &
      Success
    character ( 1 ) :: &
      Delimiter
    character ( LDL ) :: &
      ReadFormat
    character ( LDB ) :: &
      JoinedValue, &
      SingleValue

    call CheckStringValidity ( Buffer, IsValid, ValueMarker )
    if ( .not. IsValid ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false. 
      return
    end if 
    
    Success = .true. 
    
    !-- Parse Buffer

    ReadFormat = ''
    write ( ReadFormat, fmt = '(a2,i4,a5,i4,a1)' ) &
      '(a', ValueMarker, ',a1,a', LDB, ')'
    read ( Buffer, fmt = ReadFormat ) Label, Delimiter, JoinedValue 
    
    !-- Find out how many values there are
    
    nValues = 1
    oChar = 0
    do
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      if ( Marker > 0 ) then
        nValues = nValues + 1
      else
        exit
      end if
      oChar = oChar + Marker
    end do

    if ( present ( nValuesOption ) ) nValuesOption = nValues
    
    !-- Read values
    oChar = 0
    do iV = 1, min ( size ( Value ), nValues )
      
      Marker = index ( JoinedValue ( oChar + 1 : ), DELIMITER_ARRAY )
      
      if ( Marker > 0 ) then
        SingleValue = JoinedValue ( oChar + 1 : oChar + Marker - 1 )
      else
        SingleValue = JoinedValue ( oChar + 1 : )
      end if

      !-- Read value

      read ( SingleValue, fmt = '(a31)', iostat = Status ) Value ( iV )
      Value ( iV ) = trim ( adjustl ( Value ( iV ) ) )
        
      if ( Status /= 0 ) then
        Success = .false.
        exit 
      end if 
      oChar = oChar + Marker
    
    end do
    
    if ( present ( SuccessOption ) ) SuccessOption = Success
     
  end subroutine ReadLabelValueCharacter_1D


  subroutine CheckStringValidity ( Buffer, IsValid, ValueMarker ) 
  
    character ( * ), intent ( in ) :: &
      Buffer
    logical ( KDL ), intent ( out ) :: &
      IsValid
    integer ( KDI ), intent ( out ) :: &
      ValueMarker
    
    integer ( KDI ) :: &
      iC
      
    if ( Buffer ( 1:1 ) == DELIMITER_COMMENT_1 &
         .or. Buffer ( 1:1 ) == DELIMITER_COMMENT_2 ) then
       IsValid = .false.
       return
    end if
    
    iC = index ( Buffer, DELIMITER_VALUE )
    if ( iC == 0 ) then
      IsValid = .false.
      return
    end if
    
    ValueMarker = iC - 1 
    IsValid = .true.
  
  end subroutine CheckStringValidity
  

end module ReadLabelValue_Command
