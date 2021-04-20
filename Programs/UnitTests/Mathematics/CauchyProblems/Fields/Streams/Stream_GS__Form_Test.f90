program Stream_GS__Form_Test

  !-- Stream_GridStream__Form_Test

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
  type ( Grid_S_Form ), allocatable :: &
    G
  type ( FieldSet_GS_Form ), allocatable :: &
    FSG,     FSG_R, &
    FSG_234, FSG_234_R, &
    FSG_5,   FSG_5_R
  type ( Stream_GS_Form ), allocatable :: &
    SG

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_GS__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( G )
  call G % Initialize &
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

  allocate ( FSG_R )
  allocate ( FSG_234_R )
  allocate ( FSG_5_R )
  call FSG_R % Initialize &
         ( G, &
           NameOption = 'Fields_R', &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate, &
           UnitOption = FieldUnit, &
           VectorIndicesOption = VectorIndices, &
           nFieldsOption = nFields )
  call FSG_234_R % Clone &
         ( FSG_R, NameOption = 'Fields_234_R', iaSelectedOption = [ 2, 3, 4 ] )
  call FSG_5_R % Clone &
         ( FSG_R, NameOption = 'Fields_5_R', iaSelectedOption = [ 5 ] )

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

  allocate ( SG )
  call SG % Initialize ( G, GIS )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SG % AddFieldSet ( FSG )
  call SG % AddFieldSet ( FSG_234 )
  call SG % AddFieldSet ( FSG_5 )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   G     % Show ( )
  call FSG     % Show ( )
  call FSG_234 % Show ( )
  call FSG_5   % Show ( )
  call  SG     % Show ( )

  nGhostExchanges  =  1
  call PROGRAM_HEADER % GetParameter ( nGhostExchanges, 'nGhostExchanges' )

  call TestReadWrite ( SG, FSG,     FSG_R )
  call TestReadWrite ( SG, FSG_234, FSG_234_R )
  call TestReadWrite ( SG, FSG_5,   FSG_5_R )

  deallocate ( SG )
  deallocate ( FSG_5 )
  deallocate ( FSG_234 )
  deallocate ( FSG )
  deallocate ( FSG_5_R )
  deallocate ( FSG_234_R )
  deallocate ( FSG_R )
  deallocate ( G )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestReadWrite ( SG, FSG, FSG_R )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG
    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG_R

  call SetFieldSet ( FSG_R )
  call SetFieldSet ( FSG )
  call WriteStream ( SG )
  call ClearFieldSet ( FSG )
  call PROGRAM_HEADER % Communicator % Synchronize ( )  !-- avoid dbopen errors
  call ReadStream ( SG, FSG )
  call CompareFieldSets ( FSG, FSG_R )

  end subroutine TestReadWrite


  subroutine SetFieldSet ( FSG )

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

    call Show ( 'Set FieldSet' )
    call Show ( FSG % Name, 'FieldSet' )
    call Clear ( FSG % Storage % Value )

    select type ( G  =>  FSG % Chart )
    class is ( Grid_S_Form )

    associate ( nCB  =>  G % nCellsBrick )

    do iS  =  1, FSG % nFields
      iF  =  FSG % iaSelected ( iS )
      associate &
        ( F    =>  FSG % Storage % Value ( :, iF ), &
          F_U  =>  FSG % Storage % Unit ( iF ) )
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
      call ShowField ( F_3D, F_U, G % nGhostLayers, G % nDimensions )

      end associate !-- F, etc.
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
      associate &
        ( F    =>  FSG % Storage % Value ( :, iF ), &
          F_U  =>  FSG % Storage % Unit ( iF ) )
      call G % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchanges', CONSOLE % INFO_2 )
      call Show ( nGhostExchanges, 'nGhostExchanges', CONSOLE % INFO_2 )
      call Show ( FSG % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, G % nGhostLayers, G % nDimensions )
      end associate !-- F
    end do !-- iF

    if ( FSG % DevicesCommunicate ) then
      call FSG % UpdateHost ( )
      do iS  =  1, FSG % nFields
        iF  =  FSG % iaSelected ( iS )
      associate &
        ( F    =>  FSG % Storage % Value ( :, iF ), &
          F_U  =>  FSG % Storage % Unit ( iF ) )
        call G % SetFieldPointer ( F, F_3D )
        call Show ( 'Field after update host', CONSOLE % INFO_2 )
        call Show ( FSG % Field ( iF ), 'Field', CONSOLE % INFO_2 )
        call ShowField ( F_3D, F_U, G % nGhostLayers, G % nDimensions )
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


  subroutine WriteStream ( SG )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG

    associate ( GIS  =>  SG % GridImageStream ) 
    call GIS % Open ( GIS % ACCESS_CREATE )
    call  SG % Write ( )
    call GIS % Close ( )
    end associate !-- GIS

  end subroutine WriteStream


  subroutine ClearFieldSet ( FSG )

    class ( FieldSet_GS_Form ), intent ( inout ) :: &
      FSG

    call Clear ( FSG % Storage % Value )

  end subroutine ClearFieldSet


  subroutine ReadStream ( SG, FSG )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
    class ( FieldSet_GS_Form ), intent ( in ) :: &
      FSG

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field

    associate ( GIS  =>  SG % GridImageStream ) 
    call GIS % Open ( GIS % ACCESS_READ, NumberOption = GIS % Number )
    call  SG % Read ( )
    call GIS % Close ( )
    end associate !-- GIS

    call Show ( 'FieldSet after reading', CONSOLE % INFO_2 )
    call Show ( FSG % Name, 'FieldSet', CONSOLE % INFO_2 )
    do iS  =  1, FSG % nFields
      iF  =  FSG % iaSelected ( iS )
      associate &
        ( F    =>  FSG % Storage % Value ( :, iF ), &
          F_U  =>  FSG % Storage % Unit ( iF ) )
      call G % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after reading', CONSOLE % INFO_2 )
      call Show ( FSG % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, F_U, G % nGhostLayers, G % nDimensions )
      end associate !-- F, etc.
    end do !-- iF

  end subroutine ReadStream


  subroutine CompareFieldSets ( FSG, FSG_R )

    class ( FieldSet_GS_Form ), intent ( in ) :: &
      FSG, &
      FSG_R

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField
    type ( CollectiveOperation_R_Form ) :: &
      CO 

    call Show ( 'Comparing FieldSet with reference' )
    call Show ( FSG % Name, 'FieldSet' )

    associate ( nF  =>  FSG % nFields )

    select type ( G => FSG % Chart )
    class is ( Grid_S_Form )

    associate ( Cmm  =>  G % Communicator )
    call CO % Initialize &
           ( Cmm, nOutgoing = [ 2 * nF ], nIncoming = [ 2 * nF ] )
    end associate !-- Cmm

    do iS  =  1, nF
      iF  =  FSG % iaSelected ( iS )
      associate &
        ( F_R  =>  FSG_R % Storage % Value ( :, iF ), &
          F    =>  FSG   % Storage % Value ( :, iF ) )

      !-- proper cells only
      CO % Outgoing % Value ( iS )  &
        =  sum ( pack ( abs ( F  -  F_R ), mask = G % ProperCell ) )
      CO % Outgoing % Value ( nF + iS )  &
        =  sum ( pack ( abs ( F_R ), mask = G % ProperCell ) )

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


end program Stream_GS__Form_Test
