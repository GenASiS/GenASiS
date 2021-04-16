program FieldSet_GS__Form_Test

  !-- FieldSet_GridStructured__Form_Test

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
  type ( Grid_S_Form ), allocatable :: &
    G
  type ( FieldSet_GS_Form ), allocatable :: &
    FSG, &
    FSG_234, &
    FSG_5

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_GS__Form_Test', DimensionalityOption = '2D' )

  allocate ( G )
  call G % Initialize &
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

  allocate ( FSG )
  allocate ( FSG_234 )
  allocate ( FSG_5 )
  call FSG % Initialize &
         ( G, &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FSG_234 % Clone &
         ( FSG, NameOption = 'Fields_234', iaSelectedOption = [ 2, 3, 4 ] )
  call FSG_5 % Clone &
         ( FSG, NameOption = 'Fields_5', iaSelectedOption = [ 5 ] )

  call   G     % Show ( )
  call FSG     % Show ( )
  call FSG_234 % Show ( )
  call FSG_5   % Show ( )

    nGhostExchanges  =  1000
    call PROGRAM_HEADER % GetParameter ( nGhostExchanges, 'nGhostExchanges' )

  call SetField ( FSG )
  call SetField ( FSG_234 )
  call SetField ( FSG_5 )

  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( FSG_5 )
  deallocate ( FSG_234 )
  deallocate ( FSG )
  deallocate ( G )  
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetField ( FSG )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG

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
    call Show ( FSG % Name, 'FieldSet' )
    call Clear ( FSG % FieldSet % Value )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate ( nCB  =>  G % nCellsBrick )

    do iS  =  1, FSG % nFields
      iF  =  FSG % iaSelected ( iS )
      associate ( F  =>  FSG % FieldSet % Value ( :, iF ) )
      call G % SetFieldPointer ( F, F_3D )

      oC  =  ( G % iaBrick  -  1 )  *  nCB
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
      call Show ( FSG % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, G % nGhostLayers, G % nDimensions )

      end associate !-- F
    end do !-- iF

    call FSG % UpdateDevice ( )

    do iGE  =  1, nGhostExchanges
      if ( .not. FSG % DevicesCommunicate ) &
        call FSG % UpdateHost ( )
      call FSG % ExchangeGhostData ( )
      if ( .not. FSG % DevicesCommunicate ) &
        call FSG % UpdateDevice ( )
    end do !-- iGE

    do iS  =  1, FSG % nFields
      iF  =  FSG % iaSelected ( iS )
      associate ( F  =>  FSG % FieldSet % Value ( :, iF ) )
      call G % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchanges', CONSOLE % INFO_2 )
      call Show ( nGhostExchanges, 'nGhostExchanges', CONSOLE % INFO_2 )
      call Show ( FSG % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, G % nGhostLayers, G % nDimensions )
      end associate !-- F
    end do !-- iF

    if ( FSG % DevicesCommunicate ) then
      call FSG % UpdateHost ( )
      do iS  =  1, FSG % nFields
        iF  =  FSG % iaSelected ( iS )
        associate ( F  =>  FSG % FieldSet % Value ( :, iF ) )
        call G % SetFieldPointer ( F, F_3D )
        call Show ( 'Field after update host', CONSOLE % INFO_2 )
        call Show ( FSG % Field ( iF ), 'Field', CONSOLE % INFO_2 )
        call ShowField ( F_3D, G % nGhostLayers, G % nDimensions )
        end associate !-- F
      end do !-- iF
    end if !-- DevicesCommunicate

    end associate !-- nCB
    end select !-- G
    nullify ( F_3D )

  end subroutine SetField


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


end program FieldSet_GS__Form_Test
