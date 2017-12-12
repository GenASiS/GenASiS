!-- Show_Command implements an overloaded "Show" subroutine that shows either 
!   a message or a piece of data together with a label

module Show_Command

  use VariableManagement
  use CONSOLE_Singleton
 
  implicit none
  private

  public :: &
    Show, &
    ShowCharacter_KBCH

  interface Show
    module procedure ShowInteger
    module procedure ShowInteger_1D
    module procedure ShowInteger_2D
    module procedure ShowBigInteger
    module procedure ShowBigInteger_1D
    module procedure ShowBigInteger_2D
    module procedure ShowTinyInteger
    module procedure ShowTinyInteger_1D
    module procedure ShowReal
    module procedure ShowReal_1D
    module procedure ShowReal_1D_Description
    module procedure ShowReal_2D
    module procedure ShowReal_3D
    module procedure ShowRealUnitized
    module procedure ShowRealUnitized_1D
    module procedure ShowRealUnitized_1D_1D
    module procedure ShowComplex_1D
    module procedure ShowComplex_2D
    module procedure ShowComplex_3D
    module procedure ShowLogical
    module procedure ShowLogical_1D
    module procedure ShowCharacter
!    module procedure ShowCharacter_KBCH
    module procedure ShowCharacterNoDescription
    module procedure ShowCharacter_1D
    module procedure ShowMeasuredValue
    module procedure ShowMeasuredValueConvertUnit
    module procedure ShowMeasuredValue_1D
    module procedure ShowMeasuredValue_1D_ConvertUnit
  end interface Show

    private :: &
      PrepareShow, &
      EndShow

contains


  subroutine ShowInteger &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KDI ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    logical( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35,a3,i10)', trim ( Description ), '  =', Integer

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowInteger


  subroutine ShowInteger_1D &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do i = 1, size ( Integer )
      write ( IndexLabel, fmt = '( i7 )' ) i
      print &
        '(a38,i10)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', Integer ( i )
    end do

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowInteger_1D


  subroutine ShowInteger_2D &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KDI ), dimension ( :, : ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, j
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel_1, IndexLabel_2
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do j = 1, size ( Integer, dim = 2 )
      do i = 1, size ( Integer, dim = 1 )
        write ( IndexLabel_1, fmt = '( i7 )' ) i
        write ( IndexLabel_2, fmt = '( i7 )' ) j
        print &
          '(a38,i10)', &
          '( ' // trim ( adjustl ( IndexLabel_1 ) ) // ', ' &
            // trim ( adjustl ( IndexLabel_2 ) ) // ' ) =', Integer ( i, j )
      end do
    end do

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowInteger_2D


  subroutine ShowBigInteger &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KBI ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    logical( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35,a3,i10)', trim ( Description ), '  =', Integer

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowBigInteger


  subroutine ShowBigInteger_1D &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KBI ), dimension ( : ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do i = 1, size ( Integer )
      write ( IndexLabel, fmt = '( i7 )' ) i
      print &
        '(a38,i10)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', Integer ( i )
    end do

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowBigInteger_1D


  subroutine ShowBigInteger_2D &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KBI ), dimension ( :, : ), intent ( in ) ::&
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, j
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel_1, IndexLabel_2
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do j = 1, size ( Integer, dim = 2 )
      do i = 1, size ( Integer, dim = 1 )
        write ( IndexLabel_1, fmt = '( i7 )' ) i
        write ( IndexLabel_2, fmt = '( i7 )' ) j
        print &
          '(a38,i10)', &
          '( ' // trim ( adjustl ( IndexLabel_1 ) ) // ', ' &
            // trim ( adjustl ( IndexLabel_2 ) ) // ' ) =', Integer ( i, j )
      end do
    end do

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowBigInteger_2D


  subroutine ShowTinyInteger &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KTI ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    logical( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35,a3,i10)', trim ( Description ), '  =', Integer

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowTinyInteger


  subroutine ShowTinyInteger_1D &
               ( Integer, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Integer being
    !   "Show"n is more important than the Description.

    integer ( KTI ), dimension ( : ), intent ( in ) :: &
      Integer
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do i = 1, size ( Integer )
      write ( IndexLabel, fmt = '( i7 )' ) i
      print &
        '(a38,i10)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', Integer ( i )
    end do

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowTinyInteger_1D


  subroutine ShowReal &
               ( Real, Description, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption, nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Description.

    real ( KDR ), intent ( in ) :: &
      Real
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    logical ( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35,a3,es15.6e3)', trim ( Description ), '  =', Real
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowReal
  
  
  subroutine ShowReal_1D &
               ( Real, Description, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption, nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Description.

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Real
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do i = 1, size ( Real )
      write ( IndexLabel, fmt = '( i7 )' ) i
      print &
        '(a38,es15.6e3)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', Real ( i )
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowReal_1D
  
  
  subroutine ShowReal_1D_Description &
               ( Real, Description, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption, nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Description.

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Real
    character ( * ), dimension ( : ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    do i = 1, size ( Real )
      print '(a35,a3,es15.6e3)', trim ( Description ( i ) ), '  =', Real ( i )
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowReal_1D_Description
  
  
  subroutine ShowReal_2D &
               ( Real, Description, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption, nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Description.

    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Real
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, j
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel_1, IndexLabel_2

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do j = 1, size ( Real, dim = 2 )
      do i = 1, size ( Real, dim = 1 )
        write ( IndexLabel_1, fmt = '( i7 )' ) i
        write ( IndexLabel_2, fmt = '( i7 )' ) j
        print &
          '(a38,es15.6e3)', &
          '( ' // trim ( adjustl ( IndexLabel_1 ) ) // ', ' &
            // trim ( adjustl ( IndexLabel_2 ) ) // ' ) =', Real ( i, j )
      end do
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowReal_2D
  
  
  subroutine ShowReal_3D &
               ( Real, Description, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption, nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Description.

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      Real
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, j, k
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel_1, IndexLabel_2, IndexLabel_3

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do k = 1, size ( Real, dim = 3 )
      do j = 1, size ( Real, dim = 2 )
        do i = 1, size ( Real, dim = 1 )
          write ( IndexLabel_1, fmt = '( i7 )' ) i
          write ( IndexLabel_2, fmt = '( i7 )' ) j
          write ( IndexLabel_3, fmt = '( i7 )' ) k
          print &
            '(a38,es15.6e3)', &
            '( ' // trim ( adjustl ( IndexLabel_1 ) ) // ', ' &
                 // trim ( adjustl ( IndexLabel_2 ) ) // ', ' &
                 // trim ( adjustl ( IndexLabel_3 ) ) // ' ) =', &
            Real ( i, j, k )
        end do
      end do
    end do
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowReal_3D
  
  
  subroutine ShowRealUnitized &
               ( Real, Unit, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Unit or the Description.

    real ( KDR ), intent ( in ) :: &
      Real
    type ( MeasuredValueForm ), intent ( in ) :: &
      Unit
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
    
    type ( MeasuredValueForm ) :: &
      MV
    
    if ( KBCH > KDCH ) then
      call MV % Initialize_UCS &
             ( Unit % Label_UCS, Unit % Label, Real / Unit % Number )
    else
      call MV % Initialize ( Unit % Label, Real / Unit % Number )
    end if
  
    call Show &
           ( MV, Description, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption, nTrailingLinesOption )

  end subroutine ShowRealUnitized
  
  
  subroutine ShowRealUnitized_1D &
               ( Real, Unit, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Unit or the Description.

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Real
    type ( MeasuredValueForm ), intent ( in ) :: &
      Unit
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
    
    integer ( KDI ) :: &
      iV
    !-- FIXME: Made this allocatable instead of automatic array due to 
    !          PGI 11.10 bug
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      MV
      
    allocate ( MV ( size ( Real ) ) )
    
    do iV = 1, size ( Real )
      if ( KBCH > KDCH ) then
        call MV ( iV ) % Initialize_UCS &
               ( Unit % Label_UCS, Unit % Label, Real ( iV ) / Unit % Number )
      else
        call MV ( iV ) % Initialize &
               ( Unit % Label, Real ( iV ) / Unit % Number )
      end if
    end do
    
    call Show &
           ( MV, Description, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption, nTrailingLinesOption )

  end subroutine ShowRealUnitized_1D
  
  
  subroutine ShowRealUnitized_1D_1D &
               ( Real, Unit, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Unit or the Description.

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Real
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      Unit
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
    
    integer ( KDI ) :: &
      iV
    type ( MeasuredValueForm ), dimension ( size ( Real ) ) :: &
      MV
    
    do iV = 1, size ( Real )
      if ( KBCH > KDCH ) then
        call MV ( iV ) % Initialize_UCS &
               ( Unit ( iV ) % Label_UCS, Unit ( iV ) % Label, &
                 Real ( iV ) / Unit ( iV ) % Number )
      else
        call MV ( iV ) % Initialize &
               ( Unit ( iV ) % Label, Real ( iV ) / Unit ( iV ) % Number )
      end if
    end do
    
    call Show &
           ( MV, Description, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption, nTrailingLinesOption )

  end subroutine ShowRealUnitized_1D_1D
  
  
  subroutine ShowComplex_1D &
               ( Complex, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Complex being
    !   "Show"n is more important than the Description.

    complex ( KDC ), dimension ( : ), intent ( in ) :: &
      Complex
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do i = 1, size ( Complex )
      write ( IndexLabel, fmt = '( i7 )' ) i
      print &
        '(a38,2es15.6e3)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', Complex ( i )
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowComplex_1D
  
  
  subroutine ShowComplex_2D &
               ( Complex, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Complex being
    !   "Show"n is more important than the Description.

    complex ( KDC ), dimension ( :, : ), intent ( in ) :: &
      Complex
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, j
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel_1, IndexLabel_2

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do j = 1, size ( Complex, dim = 2 )
      do i = 1, size ( Complex, dim = 1 )
        write ( IndexLabel_1, fmt = '( i7 )' ) i
        write ( IndexLabel_2, fmt = '( i7 )' ) j
        print &
          '(a38,2es15.6e3)', &
          '( ' // trim ( adjustl ( IndexLabel_1 ) ) // ', ' &
            // trim ( adjustl ( IndexLabel_2 ) ) // ' ) =', Complex ( i, j )
      end do
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowComplex_2D
  
  
  subroutine ShowComplex_3D &
               ( Complex, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Complex being
    !   "Show"n is more important than the Description.

    complex ( KDC ), dimension ( :, :, : ), intent ( in ) :: &
      Complex
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, j, k
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel_1, IndexLabel_2, IndexLabel_3

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do k = 1, size ( Complex, dim = 3 )
      do j = 1, size ( Complex, dim = 2 )
        do i = 1, size ( Complex, dim = 1 )
          write ( IndexLabel_1, fmt = '( i7 )' ) i
          write ( IndexLabel_2, fmt = '( i7 )' ) j
          write ( IndexLabel_3, fmt = '( i7 )' ) k
          print &
            '(a38,2es15.6e3)', &
            '( ' // trim ( adjustl ( IndexLabel_1 ) ) // ', ' &
                 // trim ( adjustl ( IndexLabel_2 ) ) // ', ' &
                 // trim ( adjustl ( IndexLabel_3 ) ) // ' ) =', &
            Complex ( i, j, k )
        end do
      end do
    end do
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowComplex_3D
  
  
  subroutine ShowLogical &
               ( Logical, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Logical being
    !   "Show"n is more important than the Description.

    logical ( KDL ), intent ( in ) :: &
      Logical
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
      
    logical ( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return
    
    if ( Logical ) then
      print '(a35,a5,a5)', trim ( Description ), '  =  ', 'TRUE'
    else 
      print '(a35,a5,a6)', trim ( Description ), '  =  ', 'FALSE'
    end if

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowLogical
  
  
  subroutine ShowLogical_1D &
               ( Logical, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Logical being
    !   "Show"n is more important than the Description.

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      Logical
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35)', trim ( Description )
    do i = 1, size ( Logical )
      write ( IndexLabel, fmt = '( i7 )' ) i
      if ( Logical ( i ) ) then
        print &
          '(a40,a4)', &
          '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =  ', 'TRUE'
      else
        print &
          '(a40,a5)', &
          '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =  ', 'FALSE'
      end if
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowLogical_1D


  subroutine ShowCharacter &
               ( Character, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    character ( * ), intent ( in ) :: &
      Character, &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    logical( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35,a,a)', trim ( Description ), '  =  ', trim ( Character )

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowCharacter
  
  
  subroutine ShowCharacter_KBCH &
               ( Character, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    character ( *, KBCH ), intent ( in ) :: &
      Character
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    logical( KDL ) :: &
      AbortShow
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    print '(a35,a,a)', trim ( Description ), '  =  ', trim ( Character )

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowCharacter_KBCH
  
  
  subroutine ShowCharacterNoDescription &
               ( Character, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption, nTrailingLinesOption )

    character ( * ), intent ( in ) :: &
      Character
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      iLabel, &
      oChar, &
      lLabel, &
      Ignorability
    integer ( KDI ), dimension ( size ( CONSOLE % LABEL ) ) :: &
      length
    logical ( KDL ) :: &
      AbortShow
    character ( LDL ) :: &
      Label
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    Ignorability = CONSOLE % INFO_1
    if( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption 

    oChar = 2 * ( Ignorability - 1 )
!-- FIXME: This constructor gives ridiculous values, at least on Intel 12.1.2
!    lLabel = maxval ( [ ( len_trim ( CONSOLE % LABEL ( iLabel ) ), &
!                        iLabel = 1, size ( CONSOLE % LABEL ) ) ] ) 
    do iLabel = 1, size ( CONSOLE % LABEL )
      length ( iLabel ) = len_trim ( CONSOLE % LABEL ( iLabel ) )
    end do
    lLabel = maxval ( length )
    Label = ''
    Label ( oChar + 1 : oChar + lLabel + 2 ) &
      = trim ( CONSOLE % LABEL ( Ignorability ) ) // ': '

    print *, trim ( Label ) // ' ' // trim ( Character )

    call EndShow ( nTrailingLinesOption )

  end subroutine ShowCharacterNoDescription
  
  
  subroutine ShowCharacter_1D &
               ( Character, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption, oIndexOption )

    !-- Convention on argument order violated because the Character being
    !   "Show"n is more important than the Description.

    character ( * ), dimension ( : ), intent ( in ) :: &
      Character
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption, &
      oIndexOption

    integer ( KDI ) :: &
      i, &
      oIndex
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    oIndex = 0
    if ( present ( oIndexOption ) ) &
      oIndex = oIndexOption

    print '(a35)', trim ( Description )
    do i = 1, size ( Character )
      write ( IndexLabel, fmt = '( i7 )' ) oIndex + i
      print & 
        '(a40,a)', '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =  ', &
        trim ( Character ( i ) )
    end do
    
    call EndShow ( nTrailingLinesOption )

  end subroutine ShowCharacter_1D


  subroutine ShowMeasuredValue &
              ( MeasuredValue, Description, IgnorabilityOption, &
                DisplayRankOption, nLeadingLinesOption, &
                nTrailingLinesOption )

    !-- Convention on argument order violated because the MeasuredValue being
    !   "Show"n is more important than the Description.

    type ( MeasuredValueForm ), intent ( in ) :: &
      MeasuredValue
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      LenUnit
    logical ( KDL ) :: &
      AbortShow
    character ( LDL ) :: &
      PrintFormat

    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return

    if ( MeasuredValue % Label_UCS == KBCH_'' ) then

      LenUnit = len_trim ( MeasuredValue % Unit_UCS ) + 1
      write ( PrintFormat, fmt = '(a18,i0,a1)' ) &
        '(a35,a3,es15.6e3,a', LenUnit, ')'
    
      print trim ( PrintFormat ), &
        trim ( Description ), KBCH_'  =', MeasuredValue % Number, &
        KBCH_' ' // trim ( MeasuredValue % Unit_UCS )
    
    else

      LenUnit &
        = len_trim ( MeasuredValue % Unit_UCS ) + 1 &
          + len_trim ( MeasuredValue % Label_UCS ) + 5
      write ( PrintFormat, fmt = '(a18,i0,a2)' ) &
        '(a35,a3,es15.6e3,a', LenUnit, ')'
      
      print trim ( PrintFormat ), &
        trim ( Description ), KBCH_'  =', MeasuredValue % Number, &
        KBCH_' ' // trim ( MeasuredValue % Unit_UCS ) &
        // KBCH_' ( ' // trim ( MeasuredValue % Label_UCS ) // KBCH_' )'

    end if

  end subroutine ShowMeasuredValue
  
  
  subroutine ShowMeasuredValueConvertUnit &
               ( MV_Source, Unit, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Unit or the Description.

    type ( MeasuredValueForm ), intent ( in ) :: &
      MV_Source
    type ( MeasuredValueForm ), intent ( in ) :: &
      Unit
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
    
    type ( MeasuredValueForm ) :: &
      MV
    
    if ( KBCH > KDCH ) then
      call MV % Initialize_UCS &
             ( Unit % Label_UCS, Unit % Label, &
               MV_Source % Number / Unit % Number )
    else
      call MV % Initialize &
             ( Unit % Label, MV_Source % Number / Unit % Number )
    end if

    call Show ( MV, Description, IgnorabilityOption, DisplayRankOption, &
                nLeadingLinesOption, nTrailingLinesOption )

  end subroutine ShowMeasuredValueConvertUnit
  
  
  subroutine ShowMeasuredValue_1D &
              ( MeasuredValue, Description, IgnorabilityOption, &
                DisplayRankOption, nLeadingLinesOption, &
                nTrailingLinesOption )

    !-- Convention on argument order violated because the MeasuredValue being
    !   "Show"n is more important than the Description.

    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      MeasuredValue
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption

    integer ( KDI ) :: &
      i, &
      LenUnit
    logical ( KDL ) :: &
      AbortShow
    character ( LDN ) :: &
      IndexLabel
    character ( LDL ) :: &
      PrintFormat
    
    call PrepareShow &
           ( AbortShow, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption )

    if ( AbortShow ) return
    
    print '(a35)', trim ( Description )
    
    do i = 1, size ( MeasuredValue )

      if ( MeasuredValue ( i ) % Label_UCS == KBCH_'' ) then

        LenUnit &
          = len_trim ( MeasuredValue ( i ) % Unit_UCS ) + 1
        write ( PrintFormat, fmt = '(a18,i0,a1)' ) &
          '(a38,es15.6e3,a', LenUnit, ')'
      
        write ( IndexLabel, fmt = '( i7 )' ) i
        print &
          trim ( PrintFormat ), &
          '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', &
          MeasuredValue ( i ) % Number, &
          KBCH_' ' // trim ( MeasuredValue ( i ) % Unit_UCS )

      else

        LenUnit &
          = len_trim ( MeasuredValue ( i ) % Unit_UCS ) + 1 &
            + len_trim ( MeasuredValue ( i ) % Label_UCS ) + 5
        write ( PrintFormat, fmt = '(a18,i0,a2)' ) &
          '(a38,es15.6e3,a', LenUnit, ')'
      
        write ( IndexLabel, fmt = '( i7 )' ) i
        print &
          trim ( PrintFormat ), &
          '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) =', &
          MeasuredValue ( i ) % Number, &
          KBCH_' ' // trim ( MeasuredValue ( i ) % Unit_UCS ) &
          // KBCH_' ( ' // trim ( MeasuredValue ( i ) % Label_UCS ) &
          // KBCH_' )'

      end if
    
    end do
      
  end subroutine ShowMeasuredValue_1D
  
  
  subroutine ShowMeasuredValue_1D_ConvertUnit &
               ( MV_Source, Unit, Description, IgnorabilityOption, &
                 DisplayRankOption, nLeadingLinesOption, &
                 nTrailingLinesOption )

    !-- Convention on argument order violated because the Real being
    !   "Show"n is more important than the Unit or the Description.

    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      MV_Source
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      Unit
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
    
    integer ( KDI ) :: &
      iMV
    type ( MeasuredValueForm ), dimension ( size ( MV_Source ) ) :: &
      MV
    
    do iMV = 1, size ( MV_Source )
      if ( KBCH > KDCH ) then
        call MV ( iMV ) % Initialize_UCS &
               ( Unit ( iMV ) % Label_UCS, Unit ( iMV ) % Label, &
                 MV_Source ( iMV ) % Number / Unit ( iMV ) % Number )
      else
        call MV ( iMV ) % Initialize &
               ( Unit ( iMV ) % Label, &
                 MV_Source ( iMV ) % Number / Unit ( iMV ) % Number )
      end if
    end do
  
    call Show &
           ( MV, Description, IgnorabilityOption, DisplayRankOption, &
             nLeadingLinesOption, nTrailingLinesOption )

  end subroutine ShowMeasuredValue_1D_ConvertUnit
  
  
  subroutine PrepareShow &
               ( AbortShow, IgnorabilityOption, DisplayRankOption, &
                 nLeadingLinesOption )

    logical ( KDL ), intent ( out ) :: &
      AbortShow
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      DisplayRankOption, &
      nLeadingLinesOption

    integer ( KDI ) :: &
      iLine, &
      Ignorability

    Ignorability = CONSOLE % INFO_1
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption 

    AbortShow = .false.
    
    if ( CONSOLE % Muted .or. Ignorability > CONSOLE % Verbosity ) then
      AbortShow = .true.
      return
    end if

    if ( Ignorability > CONSOLE % ERROR ) then
      if ( present ( DisplayRankOption ) ) then
        if ( DisplayRankOption /= CONSOLE % ProcessRank )then
          AbortShow = .true.
          return
        end if
      else 
        if ( CONSOLE % DisplayRank /= CONSOLE % ProcessRank )then
          AbortShow = .true.
          return
        end if
      end if
    end if

    if ( present ( nLeadingLinesOption ) ) then
      do iLine = 1, nLeadingLinesOption
        print *
      end do
    end if
    
  end subroutine PrepareShow


  subroutine EndShow ( nTrailingLinesOption )

    integer ( KDI ), intent ( in ), optional :: &
      nTrailingLinesOption

    integer ( KDI ) :: &
      iLine

    if ( present ( nTrailingLinesOption ) ) then
      do iLine = 1, nTrailingLinesOption
        print *
      end do
    end if

  end subroutine EndShow


end module Show_Command
