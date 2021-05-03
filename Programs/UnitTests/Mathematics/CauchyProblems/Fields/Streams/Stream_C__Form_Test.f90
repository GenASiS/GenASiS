program Stream_C__Form_Test

  !-- Stream_Chart__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none

  integer ( KDI ) :: &
    nFields, &
    nGhostExchanges
  type ( Integer_1D_Form ), dimension ( 1 ) :: &
    VectorIndices
  type ( MeasuredValueForm ), dimension ( 5 ) :: &
    FieldUnit
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    FSC,     FSC_R, &
    FSC_234, FSC_234_R, &
    FSC_5,   FSC_5_R
  type ( Stream_C_Form ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  nFields  =  5

  FieldUnit ( 1 )      =  UNIT % MASS_DENSITY_MKS
  FieldUnit ( 2 : 4 )  =  UNIT % SPEED_MKS
  FieldUnit ( 5 )      =  UNIT % JOULE

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3, 4 ] )

  DeviceMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( PinnedMemory, 'PinnedMemory' )

  DevicesCommunicate  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

  allocate ( FSC_R )
  allocate ( FSC_234_R )
  allocate ( FSC_5_R )
  call FSC_R % Initialize &
         ( C, &
           NameOption = 'Fields_R', &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FSC_234_R % Initialize &
         ( FSC_R, iaSelected = [ 2, 3, 4 ], NameOption = 'Fields_234_R' )
  call FSC_5_R % Initialize &
         ( FSC_R, iaSelected = [ 5 ], NameOption = 'Fields_5_R' )

  allocate ( FSC )
  allocate ( FSC_234 )
  allocate ( FSC_5 )
  call FSC % Initialize &
         ( C, &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FSC_234 % Initialize &
         ( FSC, iaSelected = [ 2, 3, 4 ], NameOption = 'Fields_234' )
  call FSC_5 % Initialize &
         ( FSC, iaSelected = [ 5 ], NameOption = 'Fields_5' )

  allocate ( SC )
  call SC % Initialize ( C, GIS )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SC % AddFieldSet ( FSC )
  call SC % AddFieldSet ( FSC_234 )
  call SC % AddFieldSet ( FSC_5 )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   C     % Show ( )
  call FSC     % Show ( )
  call FSC_234 % Show ( )
  call FSC_5   % Show ( )
  call  SC     % Show ( )

  nGhostExchanges  =  1
  call PROGRAM_HEADER % GetParameter ( nGhostExchanges, 'nGhostExchanges' )

  call TestReadWrite ( SC, FSC,     FSC_R )
  call TestReadWrite ( SC, FSC_234, FSC_234_R )
  call TestReadWrite ( SC, FSC_5,   FSC_5_R )

  deallocate ( SC )
  deallocate ( FSC_5 )
  deallocate ( FSC_234 )
  deallocate ( FSC )
  deallocate ( FSC_5_R )
  deallocate ( FSC_234_R )
  deallocate ( FSC_R )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestReadWrite ( SC, FSC, FSC_R )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC, &
      FSC_R

  call SetFieldSet ( FSC_R )
  call SetFieldSet ( FSC )
  call WriteStream ( SC )
  call ClearFieldSet ( FSC )
  call PROGRAM_HEADER % Communicator % Synchronize ( )  !-- avoid dbopen errors
  call ReadStream ( SC, FSC )
  call CompareFieldSets ( FSC, FSC_R )

  end subroutine TestReadWrite


  subroutine SetFieldSet ( FSC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    integer ( KDI ) :: &
      iS, &          !-- iSelected
      iF, &          !-- iField
      iC, jC, kC, &  !-- iCell, etc.
      iGE            !-- iGhostExchange
    integer ( KDI ), dimension ( 3 ) :: &
      oC
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field

    call Show ( 'Set FieldSet' )
    call Show ( FSC % Name, 'FieldSet' )
    call Clear ( FSC % Storage_FSC % Storage % Value )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    associate ( nCB  =>  C % nCellsBrick )

    do iS  =  1, FSC % nFields
      iF  =  FSC % iaSelected ( iS )
      associate &
        ( F    =>  FSC % Storage_FSC % Storage % Value ( :, iF ), &
          F_U  =>  FSC % Storage_FSC % Storage % Unit ( iF ) )
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
      call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )

      end associate !-- F, etc.
    end do !-- iF

    call FSC % UpdateDevice ( )

    do iGE  =  1, nGhostExchanges
      if ( .not. FSC % GhostExchange_FSC % DevicesCommunicate ) &
        call FSC % UpdateHost ( )
      call FSC % ExchangeGhostData ( )
      if ( .not. FSC % GhostExchange_FSC % DevicesCommunicate ) &
        call FSC % UpdateDevice ( )
    end do !-- iGE

    do iS  =  1, FSC % nFields
      iF  =  FSC % iaSelected ( iS )
      associate &
        ( F    =>  FSC % Storage_FSC % Storage % Value ( :, iF ), &
          F_U  =>  FSC % Storage_FSC % Storage % Unit ( iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchanges', CONSOLE % INFO_2 )
      call Show ( nGhostExchanges, 'nGhostExchanges', CONSOLE % INFO_2 )
      call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )
      end associate !-- F
    end do !-- iF

    if ( FSC % GhostExchange_FSC % DevicesCommunicate ) then
      call FSC % UpdateHost ( )
      do iS  =  1, FSC % nFields
        iF  =  FSC % iaSelected ( iS )
      associate &
        ( F    =>  FSC % Storage_FSC % Storage % Value ( :, iF ), &
          F_U  =>  FSC % Storage_FSC % Storage % Unit ( iF ) )
        call C % SetFieldPointer ( F, F_3D )
        call Show ( 'Field after update host', CONSOLE % INFO_2 )
        call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
        call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )
        end associate !-- F
      end do !-- iF
    end if !-- DevicesCommunicate

    end associate !-- nCB
    end select !-- G
    nullify ( F_3D )

  end subroutine SetFieldSet


  subroutine ShowField ( F_3D, F_Unit, nGL, nD )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_3D
    type ( MeasuredValueForm ), intent ( in ) :: &
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


  subroutine WriteStream ( SC )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC

    associate ( GIS  =>  SC % GridImageStream ) 
    call GIS % Open ( GIS % ACCESS_CREATE )
    call  SC % Write ( )
    call GIS % Close ( )
    end associate !-- GIS

  end subroutine WriteStream


  subroutine ClearFieldSet ( FSC )

    class ( FieldSet_C_Form ), intent ( inout ) :: &
      FSC

    call Clear ( FSC % Storage_FSC % Storage % Value )

  end subroutine ClearFieldSet


  subroutine ReadStream ( SC, FSC )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field

    associate ( GIS  =>  SC % GridImageStream ) 
    call GIS % Open ( GIS % ACCESS_READ, NumberOption = GIS % Number )
    call  SC % Read ( )
    call GIS % Close ( )
    end associate !-- GIS

    call Show ( 'FieldSet after reading', CONSOLE % INFO_2 )
    call Show ( FSC % Name, 'FieldSet', CONSOLE % INFO_2 )
    do iS  =  1, FSC % nFields
      iF  =  FSC % iaSelected ( iS )
      associate &
        ( F    =>  FSC % Storage_FSC % Storage % Value ( :, iF ), &
          F_U  =>  FSC % Storage_FSC % Storage % Unit ( iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after reading', CONSOLE % INFO_2 )
      call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, C % nGhostLayers, C % nDimensions )
      end associate !-- F, etc.
    end do !-- iF

  end subroutine ReadStream


  subroutine CompareFieldSets ( FSC, FSC_R )

    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC, &
      FSC_R

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Comparing FieldSet with reference' )
    call Show ( FSC % Name, 'FieldSet' )

    associate ( nF  =>  FSC % nFields )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    associate ( Cmm  =>  C % Communicator )
    call CO % Initialize &
           ( Cmm, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cmm

    do iS  =  1, nF
      iF  =  FSC % iaSelected ( iS )
      associate &
        ( F_R  =>  FSC_R % Storage_FSC % Storage % Value ( :, iF ), &
          F    =>  FSC   % Storage_FSC % Storage % Value ( :, iF ) )

      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F  -  F_R ), mask = C % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = C % ProperCell ) )

!      !-- with ghost cells
!      CO % Outgoing % Value ( iF )  &
!        =  sum ( abs ( F  -  F_R ) )

      end associate !-- F_R, etc.
    end do !-- iF
    end select !-- G

    call CO % Reduce ( REDUCTION % SUM )
    call Show (    CO % Incoming % Value (      1 :      nF )  &
                /  CO % Incoming % Value ( nF + 1 : nF + nF ), 'L1 Error' )

    end associate !-- nF

  end subroutine CompareFieldSets


end program Stream_C__Form_Test
