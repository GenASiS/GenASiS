program Stream_Form_Test

  use Basics
  use Manifolds
  use FieldSets

  implicit none

  integer ( KDI ) :: &
    iD, &  !-- iDimension
    nFields, &
    nGhostExchanges
  type ( Integer_1D_Form ), dimension ( 1 ) :: &
    VectorIndices
  type ( QuantityForm ), dimension ( :, : ), allocatable :: &
    FieldUnit
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( FieldSetForm ), allocatable :: &
    FS,     FS_R, &
    FS_234, FS_234_R, &
    FS_5,   FS_5_R
  type ( StreamForm ), allocatable :: &
    S, &
    S_234, &
    S_5

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'StreamForm_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  nFields  =  5

  allocate ( FieldUnit ( nFields, A % nCharts ) )
  FieldUnit ( 1,     1 )  =  UNIT % MASS_DENSITY_MKS
  FieldUnit ( 2 : 4, 1 )  =  UNIT % SPEED_MKS
  FieldUnit ( 5,     1 )  =  UNIT % JOULE

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3, 4 ] )

  DeviceMemory  =  OffloadEnabled ( )  .and.  NumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory  =  OffloadEnabled ( )  .and.  NumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( PinnedMemory, 'PinnedMemory' )

  DevicesCommunicate  =  OffloadEnabled ( )  .and.  NumberOfDevices ( ) >= 1
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

  allocate ( FS_R )
  allocate ( FS_234_R )
  allocate ( FS_5_R )
  call FS_R % Initialize &
         ( A, &
           NameOption = 'Fields_R', &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FS_234_R % Initialize &
         ( FS_R, iaSelected = [ 2, 3, 4 ], NameOption = 'Fields_234_R' )
  call FS_5_R % Initialize &
         ( FS_R, iaSelected = [ 5 ], NameOption = 'Fields_5_R' )

  allocate ( FS )
  allocate ( FS_234 )
  allocate ( FS_5 )
  call FS % Initialize &
         ( A, &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FS_234 % Initialize &
         ( FS, iaSelected = [ 2, 3, 4 ], NameOption = 'Fields_234' )
  call FS_5 % Initialize &
         ( FS, iaSelected = [ 5 ], NameOption = 'Fields_5' )
  do iD  =  1, 3
    call FS     % SetBoundaryConditionsFace &
                   ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    call FS_234 % SetBoundaryConditionsFace &
                   ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
    call FS_5   % SetBoundaryConditionsFace &
                   ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
  end do !-- iD

  allocate ( S )
  allocate ( S_234 )
  allocate ( S_5 )
  call S     % Initialize ( A, GIS )
  call S_234 % Initialize ( A, GIS, NameOption = 'Stream_234' )
  call S_5   % Initialize ( A, GIS, NameOption = 'Stream_5' )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call S     % AddFieldSet ( FS )
  call S_234 % AddFieldSet ( FS_234 )
  call S_5   % AddFieldSet ( FS_5 )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call  A     % Show ( )
  call FS     % Show ( )
  call FS_234 % Show ( )
  call FS_5   % Show ( )
  call  S     % Show ( )
  call  S_234 % Show ( )
  call  S_5   % Show ( )

  nGhostExchanges  =  1
  call PROGRAM_HEADER % GetParameter ( nGhostExchanges, 'nGhostExchanges' )

  call TestReadWrite ( S,     FS,     FS_R )
  call TestReadWrite ( S_234, FS_234, FS_234_R )
  call TestReadWrite ( S_5,   FS_5,   FS_5_R )

  deallocate ( S_5 )
  deallocate ( S_234 )
  deallocate ( S )
  deallocate ( FS_5 )
  deallocate ( FS_234 )
  deallocate ( FS )
  deallocate ( FS_5_R )
  deallocate ( FS_234_R )
  deallocate ( FS_R )
  deallocate ( A )  
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestReadWrite ( S, FS, FS_R )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( FieldSetForm ), intent ( inout ) :: &
      FS, &
      FS_R

    call SetFieldSet ( FS_R )
    call SetFieldSet ( FS )
    call WriteStream ( S )
    call FS % Clear ( )
    !-- avoid dbopen errors
    call PROGRAM_HEADER % Communicator % Synchronize ( )  
    call ReadStream ( S, FS )
    call CompareFieldSets ( FS, FS_R )

  end subroutine TestReadWrite


  subroutine SetFieldSet ( FS )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS

    integer ( KDI ) :: &
      iS, &          !-- iSelected
      iF, &          !-- iField
      iC, jC, kC, &  !-- iCell, etc.
      iGE            !-- iGhostExchange
    integer ( KDI ), dimension ( 3 ) :: &
      oC
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field
    type ( TimerForm ), pointer :: &
      T_G

    call Show ( 'Ghost exchange' )
    call Show ( FS % Name, 'FieldSet' )
    call Show ( nGhostExchanges, 'nGhostExchanges' )
    call FS % Clear ( )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    associate ( nCB  =>  C % nCellsBrick )

    do iS  =  1, FS % nFields
      iF  =  FS % iaSelected ( iS )
      associate &
        ( F    =>  FS % Storage_GS % Value ( :, iF ), &
          F_U  =>  FS % Storage_GS % Unit ( iF ) )
      call C % SetFieldPointer ( F, F_3D )

      oC  =  ( C % iaBrick  -  1 )  *  nCB
      do kC  =  1,  nCB ( 3 )
        do jC  =  1,  nCB ( 2 )
          do iC  =  1,  nCB ( 1 )
            F_3D ( iC, jC, kC )  &
              =  iS  *  (    1.e0  *  ( oC ( 1 )  +  iC  -  1 )  &
                          +  1.e2  *  ( oC ( 2 )  +  jC  -  1 )  &
                          +  1.e4  *  ( oC ( 3 )  +  kC  -  1 ) )
          end do !-- iC
        end do !-- jC
      end do !-- kC

      call Show ( 'Field before ghost exchange', CONSOLE % INFO_2 )
      call Show ( FS % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )

      end associate !-- F, etc.
    end do !-- iS

    call FS % UpdateDevice ( )

    do iGE  =  1, nGhostExchanges
      T_G  =>  FS % TimerGhost ( Level = 1 )
      call T_G % Start ( )
      call FS % ExchangeGhostData ( )
      call T_G % Stop ( )
    end do !-- iGE

    if ( FS % DevicesCommunicate ) &
      call FS % UpdateHost ( )

    do iS  =  1, FS % nFields
      iF  =  FS % iaSelected ( iS )
      associate &
        ( F    =>  FS % Storage_GS % Value ( :, iF ), &
          F_U  =>  FS % Storage_GS % Unit ( iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchanges', CONSOLE % INFO_2 )
      call Show ( FS % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )
      end associate !-- F, etc.
    end do !-- iF

    end associate !-- nCB
    end associate !-- C
    end select !-- A

  end subroutine SetFieldSet


  subroutine ShowField ( F_3D, F_Unit, nGL, nD )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_3D
    type ( QuantityForm ), intent ( in ) :: &
      F_Unit
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nGL  !-- nGhostLayers
    integer ( KDI ), intent ( in ) :: &
      nD  !-- nDimensions
    
    call Show ( 'Selected X strips', CONSOLE % INFO_2 )
    call Show ( F_3D ( :, nGL ( 2 ) + 1, nGL ( 3 ) + 1 ), &
                F_Unit, 'F_3D ( :, 1, 1 )', CONSOLE % INFO_2 )
    if ( nD > 1 ) &
      call Show ( F_3D ( :, nGL ( 2 ), nGL ( 3 ) + 1 ), &
                  F_Unit, 'F_3D ( :, 0, 1 )', CONSOLE % INFO_2 )
    if ( nD > 2 ) then
      call Show ( F_3D ( :, nGL ( 2 ) + 1, nGL ( 3 ) ), &
                  F_Unit, 'F_3D ( :, 1, 0 )', CONSOLE % INFO_2 )
      call Show ( F_3D ( :, nGL ( 2 ),     nGL ( 3 ) ), &
                  F_Unit, 'F_3D ( :, 0, 0 )', CONSOLE % INFO_2 )
    end if

    if ( nD > 1 ) then
      call Show ( 'Selected Y strips', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ) + 1, :, nGL ( 3 ) + 1 ), &
                  F_Unit, 'F_3D ( 1, :, 1 )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ),     :, nGL ( 3 ) + 1 ), &
                  F_Unit, 'F_3D ( 0, :, 1 )', CONSOLE % INFO_2 )
      if ( nD > 2 ) then
        call Show ( F_3D ( nGL ( 1 ) + 1, :, nGL ( 3 ) ), &
                    F_Unit, 'F_3D ( 1, :, 0 )', CONSOLE % INFO_2 )
        call Show ( F_3D ( nGL ( 1 ) + 1, :, nGL ( 3 ) ), &
                    F_Unit, 'F_3D ( 1, :, 0 )', CONSOLE % INFO_2 )
      end if
    end if

    if ( nD > 2 ) then
      call Show ( 'Selected Z strips', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ) + 1, nGL ( 2 ) + 1, : ), &
                  F_Unit, 'F_3D ( 1, 1, : )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ),     nGL ( 2 ) + 1, : ), &
                  F_Unit, 'F_3D ( 0, 1, : )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ) + 1, nGL ( 2 ),     : ), &
                  F_Unit, 'F_3D ( 1, 0, : )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ),     nGL ( 2 ),     : ), &
                  F_Unit, 'F_3D ( 0, 0, : )', CONSOLE % INFO_2 )
    end if

  end subroutine ShowField


  subroutine WriteStream ( S )

    class ( StreamForm ), intent ( inout ) :: &
      S

    type ( TimerForm ), pointer :: &
      T_W

    T_W  =>  S % TimerWrite ( Level = 1 )
    call T_W % Start ( )
    associate ( GIS  =>  S % GridImageStream ) 
    call GIS % Open ( GIS % ACCESS_CREATE )
    call   S % Write ( )
    call GIS % Close ( )
    end associate !-- GIS
    call T_W % Stop ( )

  end subroutine WriteStream


  subroutine ReadStream ( S, FS )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( FieldSetForm ), intent ( in ) :: &
      FS

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field
    type ( TimerForm ), pointer :: &
      T_R

    T_R  =>  S % TimerRead ( Level = 1 )
    call T_R % Start ( )
    associate ( GIS  =>  S % GridImageStream ) 
    call GIS % Open ( GIS % ACCESS_READ, NumberOption = GIS % Number )
    call  S % Read ( )
    call GIS % Close ( )
    end associate !-- GIS
    call T_R % Stop ( )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    call Show ( 'FieldSet after reading', CONSOLE % INFO_2 )
    call Show ( FS % Name, 'FieldSet', CONSOLE % INFO_2 )
    do iS  =  1, FS % nFields
      iF  =  FS % iaSelected ( iS )
      associate &
        ( F    =>  FS % Storage_GS % Value ( :, iF ), &
          F_U  =>  FS % Storage_GS % Unit ( iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after reading', CONSOLE % INFO_2 )
      call Show ( FS % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )
      end associate !-- F, etc.
    end do !-- iF

    end associate !-- C
    end select !-- A

  end subroutine ReadStream


  subroutine CompareFieldSets ( FS, FS_R )

    class ( FieldSetForm ), intent ( in ) :: &
      FS, &
      FS_R

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Comparing FieldSet with reference' )
    call Show ( FS % Name, 'FieldSet' )

    associate ( nF  =>  FS % nFields )

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    associate ( Cr  =>  C % Communicator )
    call CO % Initialize &
           ( Cr, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cr

    do iS  =  1, nF
      iF  =  FS % iaSelected ( iS )
      associate &
        ( F_R  =>  FS_R % Storage_GS % Value ( :, iF ), &
          F    =>  FS   % Storage_GS % Value ( :, iF ) )

      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F  -  F_R ), mask = C % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = C % ProperCell ) )

!      !-- with ghost cells
!      CO % Outgoing % Value ( iF )  &
!        =  sum ( abs ( F  -  F_R ) )

     end associate !-- F_R, etc.
    end do !-- iS

    end associate !-- C
    end select !-- A

    call CO % Reduce ( REDUCTION % SUM )
    call Show (    CO % Incoming % Value (      1 :      nF )  &
                /  CO % Incoming % Value ( nF + 1 : nF + nF ), 'L1 Error' )

    end associate !-- nF

  end subroutine CompareFieldSets


end program Stream_Form_Test
