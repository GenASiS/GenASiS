program FieldSet_Form_Test

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
  type ( MeasuredValueForm ), dimension ( :, : ), allocatable :: &
    FieldUnit
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( FieldSetForm ), allocatable :: &
    FS, &
    FS_234, &
    FS_5

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )
  call CONSOLE % SetVerbosity ( 'INFO_2' )

  nFields  =  5

  allocate ( FieldUnit ( nFields, A % nCharts ) )
  FieldUnit ( 1,     1 )  =  UNIT % MASS_DENSITY_MKS
  FieldUnit ( 2 : 4, 1 )  =  UNIT % SPEED_MKS
  FieldUnit ( 5,     1 )  =  UNIT % JOULE

  call VectorIndices ( 1 ) % Initialize ( [ 2, 3, 4 ] )

  DeviceMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( PinnedMemory, 'PinnedMemory' )

  DevicesCommunicate  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

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

  call  A     % Show ( )
  call FS     % Show ( )
  call FS_234 % Show ( )
  call FS_5   % Show ( )

  nGhostExchanges  =  1000
  call PROGRAM_HEADER % GetParameter ( nGhostExchanges, 'nGhostExchanges' )

  call SetFieldSet ( FS )
  call SetFieldSet ( FS_234 )
  call SetFieldSet ( FS_5 )

  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( FS_5 )
  deallocate ( FS_234 )
  deallocate ( FS )
  deallocate ( A )  
  deallocate ( PROGRAM_HEADER )


contains


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
      T_G, T_G_UH, T_G_EG, T_G_UD

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
      associate ( F  =>  FS % Storage_GS % Value ( :, iF ) )
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
      call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )

      end associate !-- F
    end do !-- iS

    call FS % UpdateDevice ( )

    do iGE  =  1, nGhostExchanges

      T_G  =>  FS % TimerGhost    ( LevelOption = 1 )

      if ( FS % DeviceMemory .and. .not. FS % DevicesCommunicate ) then
        T_G_UH  =>  FS % TimerGhost_UH ( LevelOption = T_G % Level + 1 )
        T_G_EG  =>  FS % TimerGhost_EG ( LevelOption = T_G % Level + 1 )
        T_G_UD  =>  FS % TimerGhost_UD ( LevelOption = T_G % Level + 1 )
      else
        T_G_UH  =>  null ( )
        T_G_EG  =>  null ( )
        T_G_UD  =>  null ( )
      end if

      call T_G % Start ( )

      if ( .not. FS % DevicesCommunicate ) then
        if ( associated ( T_G_UH ) ) call T_G_UH % Start ( ) 
        call FS % UpdateHost ( )
        if ( associated ( T_G_UH ) ) call T_G_UH % Stop ( ) 
      end if

      if ( associated ( T_G_EG ) ) call T_G_EG % Start ( ) 
      call FS % ExchangeGhostData ( )
      if ( associated ( T_G_EG ) ) call T_G_EG % Stop ( ) 

      if ( .not. FS % DevicesCommunicate ) then
        if ( associated ( T_G_UD ) ) call T_G_UD % Start ( ) 
        call FS % UpdateDevice ( )
        if ( associated ( T_G_UD ) ) call T_G_UD % Stop ( ) 
      end if

      call T_G % Stop ( )

    end do !-- iGE

    do iS  =  1, FS % nFields
      iF  =  FS % iaSelected ( iS )
      associate ( F  =>  FS % Storage_GS % Value ( :, iF ) )
      call C % SetFieldPointer ( F, F_3D )
      call Show ( 'Field after ghost exchanges', CONSOLE % INFO_2 )
      call Show ( FS % Field ( iF ), 'Field', CONSOLE % INFO_2 )
      call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )
      end associate !-- F
    end do !-- iF

    if ( FS % DevicesCommunicate ) then
      call FS % UpdateHost ( )
      do iS  =  1, FS % nFields
        iF  =  FS % iaSelected ( iS )
        associate ( F  =>  FS % Storage_GS % Value ( :, iF ) )
        call C % SetFieldPointer ( F, F_3D )
        call Show ( 'Field after update host', CONSOLE % INFO_2 )
        call Show ( FS % Field ( iF ), 'Field', CONSOLE % INFO_2 )
        call ShowField ( F_3D, C % nGhostLayers, C % nDimensions )
        end associate !-- F
      end do !-- iF
    end if !-- DevicesCommunicate

    end associate !-- nCB
    end associate !-- C
    end select !-- A

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


end program FieldSet_Form_Test
