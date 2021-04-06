program FieldSet_CB__Form_Test

  !-- FieldSet_ChartBase__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics
  use BaseCharts

  implicit none

  integer ( KDI ) :: &
    nFields = 1
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
!  type ( GridImageStreamForm ), allocatable :: &
!    GIS
  type ( Chart_BH_Form ), allocatable :: &
    C
  type ( FieldSet_CB_Form ), allocatable :: &
    FSC
!  type ( Stream_CH_Form ), allocatable :: &
!    SC
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    FSM
!  type ( Stream_MH_Form ), allocatable :: &
!    SM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_CB__Form_Test', DimensionalityOption = '2D' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize &
         ( M, 'Chart', Periodic )

  call CONSOLE % SetVerbosity ( 'INFO_3' )

  DeviceMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( PinnedMemory, 'PinnedMemory' )

  DevicesCommunicate  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize &
         ( M, 'Fields', &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate )
  call FSC % Initialize &
         ( C, FSM, nFields ) 

  ! allocate ( GIS )
  ! call GIS % Initialize &
  !        ( PROGRAM_HEADER % Name, CommunicatorOption = M % Communicator )

  ! allocate ( SM )
  ! allocate ( SC )
  ! call SM % Initialize ( M, GIS, 'Stream' )
  ! call SC % Initialize ( C, SM )

  ! call SM % AddFieldSet ( FSM )

  ! call FSC % AddStream ( SC )
  ! call FSC % AddStream ( SC )  !-- Test the prevention of duplication

  call M % Show ( )
  call Show ( M % nCharts,    'nCharts',    M % IGNORABILITY )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call FSM % Show ( )
  call Show ( FSM % nStreams,  'nStreams',  FSM % IGNORABILITY )

  call C % Show ( )
  call FSC % Show ( )

  call SetField ( FSC )

!  deallocate ( SC )
!  deallocate ( SM )
!  deallocate ( GIS )
  deallocate ( FSC )
  deallocate ( FSM )

  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

contains


  subroutine SetField ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    integer ( KDI ) :: &
      iF, &       !-- iField
      iC, jC, kC  !-- iCell, etc.
    integer ( KDI ), dimension ( 3 ) :: &
      oC
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      F_3D  !-- Field

    select type ( C  =>  FSC % Chart )
    class is ( Chart_BH_Form )

    associate ( nCB  =>  C % nCellsBrick )

    do iF  =  1, FSC % nFields
      associate ( F  =>  FSC % FieldSet % Value ( :, iF ) )
      call C % SetFieldPointer ( F, F_3D )

      oC  =  ( C % iaBrick  -  1 )  *  nCB
      do kC  =  1,  nCB ( 3 )
        do jC  =  1,  nCB ( 2 )
          do iC  =  1,  nCB ( 1 )
            F_3D ( iC, jC, kC )  &
              =  iF  *  (    1.e0  *  ( oC ( 1 )  +  iC  -  1 )  &
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

    call FSC % ExchangeGhostData ( )

    do iF  =  1, FSC % nFields
      associate ( F  =>  FSC % FieldSet % Value ( :, iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchange', CONSOLE % INFO_2 )
      call Show ( FSC % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )
      end associate !-- F
    end do !-- iF

    end associate !-- nCB
    end select !-- C
    nullify ( F_3D )

  end subroutine SetField


  subroutine ShowField ( F_3D, nGL, nD )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      F_3D
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nGL  !-- nGhostLayers
    integer ( KDI ), intent ( in ) :: &
      nD  !-- nDimensions
    
    call Show ( 'Selected X strips' )
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
      call Show ( 'Selected Y strips' )
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


end program FieldSet_CB__Form_Test

