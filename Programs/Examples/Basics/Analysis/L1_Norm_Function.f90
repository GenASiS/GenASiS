module L1_Norm_Function

  use Basics
  
  implicit none
  private
  
  public :: &
    L1_Norm
    
  interface L1_Norm
    module procedure L1_Norm_Array_1D
    module procedure L1_Norm_Array_3D
    module procedure L1_Norm_VariableGroup
  end interface

contains


  function L1_Norm_Array_1D &
             ( A_Origin, A_End, CommunicatorOption, RootOption, &
               nValuesOption, oValueOption ) result ( L1 )
               
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      A_Origin, &
      A_End
    type ( CommunicatorForm ), intent ( in ), optional, target :: &
      CommunicatorOption
    integer ( KDI ), intent ( in ), optional :: &
      RootOption, &
      nValuesOption, &
      oValueOption
    real ( KDR ) :: &
      L1
      
    integer ( KDI ) :: &
      oV, &
      nValues
    real ( KDR ) :: &
      MyDistanceSum, &
      MyOriginSum
    type ( CommunicatorForm ), pointer :: &
      C
    type ( CollectiveOperation_R_Form ) :: &
      CO
      
    oV = 0
    if ( present ( oValueOption ) ) oV = oValueOption
    
    nValues = size ( A_Origin )
    if ( present ( nValuesOption ) ) nValues = nValuesOption
    
    C => PROGRAM_HEADER % Communicator
    if ( present ( CommunicatorOption ) ) C => CommunicatorOption
    
    MyDistanceSum &
      = sum ( abs ( A_End ( oV + 1 : oV + nValues ) &
                    - A_Origin ( oV + 1 : oV + nValues ) ) )
    MyOriginSum = sum ( abs ( A_Origin ( oV + 1 : oV + nValues ) ) )
  
    call CO % Initialize ( C, [ 2 ], [ 2 ], RootOption = RootOption )
    
    CO % Outgoing % Value = [ MyDistanceSum, MyOriginSum ]
    call CO % Reduce ( REDUCTION % SUM )

    if ( present ( RootOption ) ) then
      if ( C % Rank == RootOption ) then
        MyDistanceSum = CO % Incoming % Value ( 1 )
        MyOriginSum   = CO % Incoming % Value ( 2 )
        if ( MyOriginSum /= 0.0_KDR ) then
          L1 = MyDistanceSum / MyOriginSum
        else
          L1 = 0.0_KDR
          call Show ( 'Cannot compute L1 error for vanishing variable', &
                      CONSOLE % WARNING )
        end if
      end if
    else
      MyDistanceSum = CO % Incoming % Value ( 1 )
      MyOriginSum   = CO % Incoming % Value ( 2 )
      if ( MyOriginSum /= 0.0_KDR ) then
        L1 = MyDistanceSum / MyOriginSum
      else
        L1 = 0.0_KDR
        call Show ( 'Cannot compute L1 error for vanishing variable', &
                    CONSOLE % WARNING )
      end if
    end if
       
    nullify ( C )

  end function L1_Norm_Array_1D
  
  
  function L1_Norm_Array_3D &
             ( A_Origin, A_End, CommunicatorOption, RootOption, &
               nValuesOption, oValueOption ) result ( L1 )
               
    real ( KDR ), dimension ( :, :, : ), intent ( in ), target :: &
      A_Origin, &
      A_End
    type ( CommunicatorForm ), intent ( in ), optional, target :: &
      CommunicatorOption
    integer ( KDI ), intent ( in ), optional :: &
      RootOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nValuesOption, &
      oValueOption
    real ( KDR ) :: &
      L1
      
    integer ( KDI ), dimension ( 3 ) :: &
      oV, &
      nValues
    real ( KDR ) :: &
      MyDistanceSum, &
      MyOriginSum
    type ( CommunicatorForm ), pointer :: &
      C
    type ( CollectiveOperation_R_Form ) :: &
      CO
      
    oV = 0
    if ( present ( oValueOption ) ) oV = oValueOption
    
    nValues = shape ( A_Origin )
    if ( present ( nValuesOption ) ) nValues = nValuesOption
    
    C => PROGRAM_HEADER % Communicator
    if ( present ( CommunicatorOption ) ) C => CommunicatorOption
    
    call Show ( lbound ( A_Origin ), 'LowerBounds' )
    call Show ( ubound ( A_Origin ), 'UpperBounds' )
    
    MyDistanceSum &
      = sum ( abs ( A_End ( oV ( 1 ) + 1 : oV ( 1 ) + nValues ( 1 ), &
                            oV ( 2 ) + 1 : oV ( 2 ) + nValues ( 2 ), &
                            oV ( 3 ) + 1 : oV ( 3 ) + nValues ( 3 ) ) &
                    - A_Origin ( oV ( 1 ) + 1 : oV ( 1 ) + nValues ( 1 ), &
                                 oV ( 2 ) + 1 : oV ( 2 ) + nValues ( 2 ), &
                                 oV ( 3 ) + 1 : oV ( 3 ) + nValues ( 3 ) ) ) )
    MyOriginSum &
      = sum ( abs ( A_Origin ( oV ( 1 ) + 1 : oV ( 1 ) + nValues ( 1 ), &
                               oV ( 2 ) + 1 : oV ( 2 ) + nValues ( 2 ), &
                               oV ( 3 ) + 1 : oV ( 3 ) + nValues ( 3 ) ) ) )
  
    call CO % Initialize ( C, [ 2 ], [ 2 ], RootOption = RootOption )
    
    CO % Outgoing % Value = [ MyDistanceSum, MyOriginSum ]
    call CO % Reduce ( REDUCTION % SUM )

    if ( present ( RootOption ) ) then
      if ( C % Rank == RootOption ) then
        MyDistanceSum = CO % Incoming % Value ( 1 )
        MyOriginSum   = CO % Incoming % Value ( 2 )
        if ( MyOriginSum /= 0.0_KDR ) then
          L1 = MyDistanceSum / MyOriginSum
        else
          L1 = 0.0_KDR
          call Show ( 'Cannot compute L1 error for vanishing variable', &
                      CONSOLE % WARNING )
        end if
      end if
    else
      MyDistanceSum = CO % Incoming % Value ( 1 )
      MyOriginSum   = CO % Incoming % Value ( 2 )
      if ( MyOriginSum /= 0.0_KDR ) then
        L1 = MyDistanceSum / MyOriginSum
      else
        L1 = 0.0_KDR
        call Show ( 'Cannot compute L1 error for vanishing variable', &
                    CONSOLE % WARNING )
      end if
    end if
       
    nullify ( C )

  end function L1_Norm_Array_3D
  
  
  function L1_Norm_VariableGroup &
             ( VG_Origin, VG_End, CommunicatorOption, iaSelectedOption, &
               RootOption, nValuesOption, oValueOption ) result ( L1 )
             
    class ( VariableGroupForm ), intent ( in ), target :: &
      VG_Origin, &
      VG_End
    type ( CommunicatorForm ), intent ( in ), optional, target :: &
      CommunicatorOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional, target :: &
      iaSelectedOption
    integer ( KDI ), intent ( in ), optional :: &
      RootOption, &
      nValuesOption, &
      oValueOption
    real ( KDR ), dimension ( : ), allocatable :: &
      L1
    
    integer ( KDI ) :: &
      iS, &
      iVrbl, &
      oV, &
      nValues
    integer ( KDI ), dimension ( : ), pointer :: &
      iaSelected
    real ( KDR ), dimension ( : ), allocatable :: &
      MyDistanceSum, &
      MyOriginSum
    type ( CommunicatorForm ), pointer :: &
      C
    type ( CollectiveOperation_R_Form ) :: &
      CO

    oV = 0
    if ( present ( oValueOption ) ) oV = oValueOption
    
    nValues = VG_Origin % nValues
    if ( present ( nValuesOption ) ) nValues = nValuesOption
    
    iaSelected => VG_Origin % iaSelected
    if ( present ( iaSelectedOption ) ) iaSelected => iaSelectedOption
    
    C => PROGRAM_HEADER % Communicator
    if ( present ( CommunicatorOption ) ) C => CommunicatorOption
    
    allocate ( MyDistanceSum ( size ( iaSelected ) ) )
    allocate ( MyOriginSum ( size ( iaSelected ) ) )  
    allocate ( L1 ( size ( iaSelected ) ) )
    
    do iS = 1, size ( iaSelected )
      iVrbl = iaSelected ( iS )
      MyDistanceSum ( iS ) &
        = sum ( abs ( VG_End % Value ( oV + 1 : oV + nValues, iVrbl ) &
                      - VG_Origin % Value ( oV + 1 : oV + nValues, iVrbl ) ) )
      MyOriginSum ( iS ) &
        = sum ( abs ( VG_Origin % Value ( oV + 1 : oV + nValues, iVrbl ) ) )
    end do
    
    call CO % Initialize &
           ( C, [ 2 * size ( iaSelected ) ], [ 2 * size ( iaSelected ) ], &
             RootOption = RootOption )
    
    CO % Outgoing % Value = [ MyDistanceSum, MyOriginSum ]
    call CO % Reduce ( REDUCTION % SUM )

    if ( present ( RootOption ) ) then
      if ( C % Rank == RootOption ) then
        do iS = 1, size ( iaSelected )
          MyDistanceSum ( iS ) = CO % Incoming % Value ( iS )
          MyOriginSum ( iS )  = CO % Incoming % Value ( size ( iaSelected ) + iS )
          if ( MyOriginSum ( iS ) /= 0.0_KDR ) then
            L1 ( iS ) = MyDistanceSum ( iS ) / MyOriginSum ( iS )
          else
            L1 ( iS ) = 0.0_KDR
            call Show ( 'Cannot compute L1 error for vanishing variable', &
                        CONSOLE % WARNING )
          end if
        end do
      end if
    else
      do iS = 1, size ( iaSelected )
        MyDistanceSum ( iS ) = CO % Incoming % Value ( iS )
        MyOriginSum ( iS )  = CO % Incoming % Value ( size ( iaSelected ) + iS )
        if ( MyOriginSum ( iS ) /= 0.0_KDR ) then
          L1 ( iS ) = MyDistanceSum ( iS ) / MyOriginSum ( iS )
        else
          L1 ( iS ) = 0.0_KDR
          call Show ( 'Cannot compute L1 error for vanishing variable', &
                      CONSOLE % WARNING )
        end if
      end do
    end if

    nullify ( C )
    nullify ( iaSelected )

  end function L1_Norm_VariableGroup
  
  
end module L1_Norm_Function
