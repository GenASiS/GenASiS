program FieldSet_C__Form_Test

  !-- FieldSet_Chart_Form_Test

  use Basics
  use Manifolds
  use FieldSets

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
  type ( Chart_GS_Form ), allocatable :: &
    C
  type ( FieldSet_C_Form ), allocatable :: &
    FSC, &
    FSC_234, &
    FSC_5

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_C__Form_Test', DimensionalityOption = '2D' )

  allocate ( C )
  call C % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

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

  call   C     % Show ( )
  call FSC     % Show ( )
  call FSC_234 % Show ( )
  call FSC_5   % Show ( )

  nGhostExchanges  =  1000
  call PROGRAM_HEADER % GetParameter ( nGhostExchanges, 'nGhostExchanges' )

  call SetFieldSet ( FSC )
  call SetFieldSet ( FSC_234 )
  call SetFieldSet ( FSC_5 )

  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( FSC_5 )
  deallocate ( FSC_234 )
  deallocate ( FSC )
  deallocate ( C )  
  deallocate ( PROGRAM_HEADER )


contains


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

    call Show ( 'Ghost exchange' )
    call Show ( FSC % Name, 'FieldSet' )
    call Clear ( FSC % Storage_FSC % Storage % Value )

    select type ( C  =>  FSC % Chart )
    class is ( Chart_GS_Form )

    associate ( nCB  =>  C % nCellsBrick )

    do iS  =  1, FSC % nFields
      iF  =  FSC % iaSelected ( iS )
      associate ( F  =>  FSC % Storage_FSC % Storage % Value ( :, iF ) )
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
      call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )

      end associate !-- F
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
      associate ( F  =>  FSC % Storage_FSC % Storage % Value ( :, iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchanges', CONSOLE % INFO_2 )
      call Show ( nGhostExchanges, 'nGhostExchanges', CONSOLE % INFO_2 )
      call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )
      end associate !-- F
    end do !-- iF

    if ( FSC % GhostExchange_FSC % DevicesCommunicate ) then
      call FSC % UpdateHost ( )
      do iS  =  1, FSC % nFields
        iF  =  FSC % iaSelected ( iS )
        associate ( F  =>  FSC % Storage_FSC % Storage % Value ( :, iF ) )
        call C % SetFieldPointer ( F, F_3D )
        call Show ( 'Field after update host', CONSOLE % INFO_2 )
        call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
        call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )
        end associate !-- F
      end do !-- iF
    end if !-- DevicesCommunicate

    end associate !-- nCB
    end select !-- G
    nullify ( F_3D )

  end subroutine SetFieldSet


  subroutine ShowField ( F_3D, nGL, nD )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_3D
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nGL  !-- nGhostLayers
    integer ( KDI ), intent ( in ) :: &
      nD  !-- nDimensions
    
    call Show ( 'Selected X strips', CONSOLE % INFO_2 )
    call Show ( F_3D ( :, nGL ( 2 ) + 1, nGL ( 3 ) + 1 ), &
                'F_3D ( :, 1, 1 )', CONSOLE % INFO_2 )
    if ( nD > 1 ) &
      call Show ( F_3D ( :, nGL ( 2 ), nGL ( 3 ) + 1 ), &
                  'F_3D ( :, 0, 1 )', CONSOLE % INFO_2 )
    if ( nD > 2 ) then
      call Show ( F_3D ( :, nGL ( 2 ) + 1, nGL ( 3 ) ), &
                  'F_3D ( :, 1, 0 )', CONSOLE % INFO_2 )
      call Show ( F_3D ( :, nGL ( 2 ),     nGL ( 3 ) ), &
                  'F_3D ( :, 0, 0 )', CONSOLE % INFO_2 )
    end if

    if ( nD > 1 ) then
      call Show ( 'Selected Y strips', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ) + 1, :, nGL ( 3 ) + 1 ), &
                  'F_3D ( 1, :, 1 )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ),     :, nGL ( 3 ) + 1 ), &
                  'F_3D ( 0, :, 1 )', CONSOLE % INFO_2 )
      if ( nD > 2 ) then
        call Show ( F_3D ( nGL ( 1 ) + 1, :, nGL ( 3 ) ), &
                    'F_3D ( 1, :, 0 )', CONSOLE % INFO_2 )
        call Show ( F_3D ( nGL ( 1 ) + 1, :, nGL ( 3 ) ), &
                    'F_3D ( 1, :, 0 )', CONSOLE % INFO_2 )
      end if
    end if

    if ( nD > 2 ) then
      call Show ( 'Selected Z strips', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ) + 1, nGL ( 2 ) + 1, : ), &
                  'F_3D ( 1, 1, : )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ),     nGL ( 2 ) + 1, : ), &
                  'F_3D ( 0, 1, : )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ) + 1, nGL ( 2 ),     : ), &
                  'F_3D ( 1, 0, : )', CONSOLE % INFO_2 )
      call Show ( F_3D ( nGL ( 1 ),     nGL ( 2 ),     : ), &
                  'F_3D ( 0, 0, : )', CONSOLE % INFO_2 )
    end if

  end subroutine ShowField


end program FieldSet_C__Form_Test
