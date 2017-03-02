!-- Given an array of strings of the form label=value, FindParameter 
!   finds the one with a particular label (if it exists) and reads the value.

module FindParameter_Command

  use VariableManagement
  use Display
  use ReadLabelValue_Command

  implicit none
  private


  public :: &
    FindParameter

  interface FindParameter
    module procedure FindParameterInteger_0D
    module procedure FindParameterReal_0D
    module procedure FindParameterMeasuredValue_0D
    module procedure FindParameterLogical_0D
    module procedure FindParameterCharacter_0D
    module procedure FindParameterInteger_1D
    module procedure FindParameterReal_1D
    module procedure FindParameterMeasuredValue_1D
    module procedure FindParameterLogical_1D
    module procedure FindParameterCharacter_1D
  end interface FindParameter

contains


  subroutine FindParameterInteger_0D &
               ( Value, Buffer, Source, Name, IgnorabilityOption, &
                 SuccessOption )

    integer ( KDI ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    integer ( KDI ) :: &
      Scratch

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue ( Label, Scratch, Buffer ( iB ), Success ) 
      if ( Success .and. trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      Value = Scratch
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      call Show ( Value, Name, Ignorability )
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterInteger_0D


  subroutine FindParameterReal_0D &
               ( Value, Buffer, Source, Name, InputUnitOption, &
                 IgnorabilityOption, SuccessOption )

    real ( KDR ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    real ( KDR ) :: &
      Scratch
    type ( MeasuredValueForm ) :: &
      InputUnit

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), InputUnit, Success ) 
      if ( Success .and. trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      Value = Scratch
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      if ( InputUnit % Label == '' ) then
        call Show &
               ( Value, trim ( Name ) // ' ( dimensionless )', &
                 Ignorability )
      else
        call Show &
               ( Value, InputUnit, trim ( Name ) // ' ( input units )', &
                 Ignorability )
        call Show &
               ( Value, trim ( Name ) // ' ( program units )', &
                 Ignorability )
      end if
      if ( present ( InputUnitOption ) ) InputUnitOption = InputUnit
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterReal_0D


  subroutine FindParameterMeasuredValue_0D &
               ( Value, Buffer, Source, Name, InputUnitOption, &
                 IgnorabilityOption, ConvertOption, SuccessOption )

    type ( MeasuredValueForm ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( in ), optional :: &
      ConvertOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      Ignorability
    logical ( KDL ) :: &
      Convert, &
      Success
    character ( LDF ) :: &
      Label
    type ( MeasuredValueForm ) :: &
      Scratch
    type ( MeasuredValueForm ) :: &
      InputUnit

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Convert = .true.
    if ( present ( ConvertOption ) ) Convert = ConvertOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), InputUnit, Success ) 
      if ( Success .and. trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      Value = Scratch
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      if ( InputUnit % Label == '' ) then
        call Show &
               ( Value, trim ( Name ) // ' ( dimensionless )', &
                 Ignorability )
      else if ( Convert ) then
        call Show &
               ( Value, InputUnit, trim ( Name ) // ' ( input units )', &
                 Ignorability )
        call Show &
               ( Value, trim ( Name ) // ' ( program units )', &
                 Ignorability )
      else 
        call Value % Initialize &
               ( InputUnit % Label, Value % Number / InputUnit % Number )
        call Show &
               ( Value, trim ( Name ) // ' ( unconverted )', &
                 Ignorability )        
      end if
      if ( present ( InputUnitOption ) ) InputUnitOption = InputUnit
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterMeasuredValue_0D


  subroutine FindParameterLogical_0D &
               ( Value, Buffer, Source, Name, IgnorabilityOption, &
                 SuccessOption )

    logical ( KDL ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    logical ( KDL ) :: &
      Scratch

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue ( Label, Scratch, Buffer ( iB ), Success ) 
      if ( Success .and. trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      Value = Scratch
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      call Show ( Value, Name, Ignorability )
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterLogical_0D


  subroutine FindParameterCharacter_0D &
               ( Value, Buffer, Source, Name, IgnorabilityOption, &
                 SuccessOption )

    character ( * ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    character ( LDB ) :: &
      Scratch

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue ( Label, Scratch, Buffer ( iB ), Success ) 
      if ( Success .and. trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      Value = Scratch
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      call Show ( Value, Name, Ignorability )
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterCharacter_0D


  subroutine FindParameterInteger_1D &
               ( Value, Buffer, Source, Name, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      nValues, &
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    integer ( KDI ), dimension ( size ( Value ) ) :: &
      Scratch

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), nValues, Success )
      if ( Success .and.  trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      if ( nValues > size ( Value ) ) then
        call Show ( 'Ignoring extra input values', CONSOLE % WARNING )
        call Show ( size ( Value ), 'nValues requested', CONSOLE % WARNING )
        call Show ( nValues, 'nValues input', CONSOLE % WARNING )
        nValues = size ( Value )
      end if
      Value ( : nValues ) = Scratch ( : nValues )
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      call Show ( Value ( : nValues ), Name, Ignorability )
      if ( present ( nValuesOption ) ) nValuesOption = nValues
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterInteger_1D


  subroutine FindParameterReal_1D &
               ( Value, Buffer, Source, Name, InputUnitOption, &
                 nValuesOption, IgnorabilityOption, SuccessOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), &
      optional :: &
        InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      nValues, &
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    real ( KDR ), dimension ( size ( Value ) ) :: &
      Scratch
    type ( MeasuredValueForm ), dimension ( size ( Value ) ) :: &
      InputUnit

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), InputUnit, nValues, Success )
      if ( Success .and.  trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      if ( nValues > size ( Value ) ) then
        call Show ( 'Ignoring extra input values', CONSOLE % WARNING )
        call Show ( size ( Value ), 'nValues requested', CONSOLE % WARNING )
        call Show ( nValues, 'nValues input', CONSOLE % WARNING )
        nValues = size ( Value )
      end if
      Value ( : nValues ) = Scratch ( : nValues )
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      if ( all ( InputUnit ( : nValues ) % Label == '' ) ) then
        call Show &
               ( Value ( : nValues ), &
                 trim ( Name ) // ' ( dimensionless )', Ignorability )
      else
        call Show &
               ( Value ( : nValues ), InputUnit ( : nValues ), &
                 trim ( Name ) // ' ( input units )', Ignorability )
        call Show &
               ( Value ( : nValues ), &
                 trim ( Name ) // ' ( program units )', Ignorability )
      end if
      if ( present ( nValuesOption ) ) nValuesOption = nValues
      if ( present ( InputUnitOption ) ) InputUnitOption = InputUnit
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterReal_1D


  subroutine FindParameterMeasuredValue_1D &
               ( Value, Buffer, Source, Name, InputUnitOption, &
                 nValuesOption, IgnorabilityOption, SuccessOption )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), &
      optional :: &
        InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      nValues, &
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    type ( MeasuredValueForm ), dimension ( size ( Value ) ) :: &
      Scratch
    type ( MeasuredValueForm ), dimension ( size ( Value ) ) :: &
      InputUnit

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), InputUnit, nValues, Success )
      if ( Success .and.  trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      if ( nValues > size ( Value ) ) then
        call Show ( 'Ignoring extra input values', CONSOLE % WARNING )
        call Show ( size ( Value ), 'nValues requested', CONSOLE % WARNING )
        call Show ( nValues, 'nValues input', CONSOLE % WARNING )
        nValues = size ( Value )
      end if
      Value ( : nValues ) = Scratch ( : nValues )
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      if ( all ( InputUnit ( : nValues ) % Label == '' ) ) then
        call Show &
               ( Value ( : nValues ), &
                 trim ( Name ) // ' ( dimensionless )', Ignorability )
      else
        call Show &
               ( Value ( : nValues ), InputUnit ( : nValues ), &
                 trim ( Name ) // ' ( input units )', Ignorability )
        call Show &
               ( Value ( : nValues ), &
                 trim ( Name ) // ' ( program units )', Ignorability )
      end if
      if ( present ( nValuesOption ) ) nValuesOption = nValues
      if ( present ( InputUnitOption ) ) InputUnitOption = InputUnit
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterMeasuredValue_1D


  subroutine FindParameterLogical_1D &
               ( Value, Buffer, Source, Name, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    logical ( KDL ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    integer ( KDI ), intent ( out ), optional :: &
      nValuesOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      nValues, &
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    logical ( KDL ), dimension ( size ( Value ) ) :: &
      Scratch

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), nValues, Success )
      if ( Success .and.  trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      if ( nValues > size ( Value ) ) then
        call Show ( 'Ignoring extra input values', CONSOLE % WARNING )
        call Show ( size ( Value ), 'nValues requested', CONSOLE % WARNING )
        call Show ( nValues, 'nValues input', CONSOLE % WARNING )
        nValues = size ( Value )
      end if
      Value ( : nValues ) = Scratch ( : nValues )
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      call Show ( Value ( : nValues ), Name, Ignorability )
      if ( present ( nValuesOption ) ) nValuesOption = nValues
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterLogical_1D


  subroutine FindParameterCharacter_1D &
               ( Value, Buffer, Source, Name, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    character ( * ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), dimension ( : ), intent ( in ) :: &
      Buffer
    character ( * ), intent ( in ) :: &
      Source, &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    integer ( KDI ), intent ( out ), optional :: &
      nValuesOption

    integer ( KDI ) :: &
      iB, &  !-- iBuffer
      nValues, &
      Ignorability
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      Label
    character ( LDB ), dimension ( size ( Value ) ) :: &
      Scratch

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    Success = .false.
    do iB = 1, size ( Buffer )
      call ReadLabelValue &
             ( Label, Scratch, Buffer ( iB ), nValues, Success )
      if ( Success .and.  trim ( Label ) == trim ( Name ) ) exit
      Success = .false.
    end do

    if ( present ( SuccessOption ) ) SuccessOption = Success

    if ( Success ) then
      if ( nValues > size ( Value ) ) then
        call Show ( 'Ignoring extra input values', CONSOLE % WARNING )
        call Show ( size ( Value ), 'nValues requested', CONSOLE % WARNING )
        call Show ( nValues, 'nValues input', CONSOLE % WARNING )
        nValues = size ( Value )
      end if
      Value ( : nValues ) = Scratch ( : nValues )
      call Show ( 'Parameter ' // trim ( Name ) // ' found', Ignorability )
      call Show ( Source, 'Source', Ignorability )
      call Show ( Value ( : nValues ), Name, Ignorability )
      if ( present ( nValuesOption ) ) nValuesOption = nValues
    else
      call Show ( 'Parameter ' // trim ( Name ) // ' not found', &
                  Ignorability )
      call Show ( Source, 'Source', Ignorability )
    end if

  end subroutine FindParameterCharacter_1D


end module FindParameter_Command
