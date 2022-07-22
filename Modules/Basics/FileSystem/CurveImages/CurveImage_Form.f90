!-- CurveImagesSiloForm provides a class to generate one-dimensional
!   grid data on XY-plane (curve) as columns in Silo format, suitable for
!   VisIt visualization tool.

module CurveImage_Form

  use Specifiers
  use DataManagement
  use Display
  use MessagePassing
  use GridImageBasics
        
  implicit none 
  private 
  
  include 'silo_f9x.inc'
  
  type, public, extends ( GridImageSiloTemplate ) :: CurveImageForm
  contains
    procedure, private, pass :: &
      SetGridWriteUnigrid
    procedure, private, pass :: &
      SetGridWriteRefinable
    generic, public :: &
      SetGridWrite => SetGridWriteUnigrid, SetGridWriteRefinable
    procedure, public, pass :: &
      SetGridRead
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    procedure, public, pass :: &
      ClearGrid
    final :: &
      Finalize
  end type CurveImageForm

contains

    
  subroutine SetGridWriteUnigrid &
               ( CI, Directory, Edge, nProperCells, oValue, &
                 CoordinateLabelOption, CoordinateUnitOption )
                 
    class ( CurveImageForm ), intent ( inout ) :: &
      CI
    character ( * ), intent ( in ) :: &
      Directory
    type ( Real_1D_Form ), intent ( in ) :: &
      Edge
    integer ( KDI ), intent ( in ) :: &
      nProperCells, &
      oValue
    character ( * ), intent ( in ), optional :: &
      CoordinateLabelOption
    type ( QuantityForm ), intent ( in ), optional :: &
      CoordinateUnitOption

    CI % oValue       = oValue
    CI % nTotalCells  = nProperCells
    CI % nGhostCells  = 0
    CI % lDirectory   = len_trim ( Directory )
    
    if ( size ( Edge % Value ) == nProperCells ) then
!-- FIXME: NAG 5.3.1 should support sourced allocation but doesn't compile
!      allocate ( CI % NodeCoordinate_1, source = NodeCoordinate )
      allocate ( CI % NodeCoordinate_1 ( nProperCells ) )
      CI % NodeCoordinate_1 = Edge % Value
    else
      associate &
        ( nP => nProperCells, &
          oV => oValue )
      allocate ( CI % NodeCoordinate_1 ( nProperCells ) )
      CI % NodeCoordinate_1 &
        = Edge % Value ( oV + 1 : oV + nP ) &
          + ( Edge % Value ( oV + 2 : oV + nP + 1 ) &
              - Edge % Value ( oV + 1 : oV + nP ) ) &
            * 0.5_KDR
      end associate !-- nP, etc.
    end if
    
    CI % Directory = Directory
    
    if ( present ( CoordinateUnitOption ) ) &
      CI % CoordinateUnit ( 1 ) = CoordinateUnitOption
    
    if ( present ( CoordinateLabelOption ) ) &
      CI % CoordinateLabel ( 1 ) &
        = CoordinateLabelOption

  end subroutine SetGridWriteUnigrid
  
  
  subroutine SetGridWriteRefinable &
               ( CI, Directory, NodeCoordinate, nProperCells, oValue, &
                 CoordinateUnitOption, CoordinateLabelOption )
                 
    class ( CurveImageForm ), intent ( inout ) :: &
      CI
    character ( * ), intent ( in ) :: &
      Directory
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      NodeCoordinate
    integer ( KDI ), intent ( in ) :: &
      nProperCells, &
      oValue
    type ( QuantityForm ), intent ( in ), optional :: &
      CoordinateUnitOption
    character ( * ), intent ( in ), optional :: &
      CoordinateLabelOption
    
    CI % oValue       = oValue
    CI % nTotalCells  = nProperCells
    CI % nGhostCells  = 0
    CI % lDirectory   = len_trim ( Directory )
    
    if ( size ( NodeCoordinate ) == nProperCells ) then
!-- FIXME: NAG 5.3.1 should support sourced allocation but doesn't compile
!      allocate ( CI % NodeCoordinate_1, source = NodeCoordinate )
      allocate ( CI % NodeCoordinate_1 ( size ( NodeCoordinate ) ) )
      CI % NodeCoordinate_1 = NodeCoordinate
    else
      allocate ( CI % NodeCoordinate_1 ( nProperCells ) )
      CI % NodeCoordinate_1 &
        = ( NodeCoordinate ( 2 : nProperCells + 1 ) &
              - NodeCoordinate ( 1 : nProperCells ) ) * 0.5_KDR &
          + NodeCoordinate ( 1 : nProperCells )
    end if
    
    CI % Directory = Directory
    
    if ( present ( CoordinateUnitOption ) ) &
      CI % CoordinateUnit ( 1 ) = CoordinateUnitOption
    
    if ( present ( CoordinateLabelOption ) ) &
      CI % CoordinateLabel ( 1 ) &
        = CoordinateLabelOption

  end subroutine SetGridWriteRefinable
  
  
  subroutine SetGridRead ( CI, Directory, nProperCells, oValue )
                 
    class ( CurveImageForm ), intent ( inout ) :: &
      CI
    character ( * ), intent ( in ) :: &
      Directory
    integer ( KDI ), intent ( in ) :: &
      nProperCells, &
      oValue

    CI % oValue       =  oValue
    CI % nTotalCells  =  nProperCells
    CI % nGhostCells  =  0
    CI % lDirectory   =  len_trim ( Directory )    
    CI % Directory    =  Directory

    allocate ( CI % NodeCoordinate_1 ( nProperCells ) )
    
  end subroutine SetGridRead
  
  
  subroutine Write ( GI, TimeOption, CycleNumberOption )
  
    class ( CurveImageForm ), intent ( inout ) :: &
      GI
    type ( QuantityForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iV, &      !-- iValue
      iStrg, &     !-- iStorage
      iS, &     !-- iSelected
      iVrbl, &  !-- iVariable
      nSiloOptions, &
      Error, &
      SiloOptionList
    real ( KDR ), dimension ( : ), allocatable :: &
      Coordinate, &
      Value, &
      Coordinate_MM, &  !-- CoordinateMultiMesh
      Value_MM          !-- ValueMultiMesh
    character ( LDF ) :: &
      WorkingDirectory
    type ( CollectiveOperation_R_Form ), target :: &
      CO_Coordinate, &
      CO_Variable

    if ( .not. GI % Stream % IsWritable ( ) ) return
    
    allocate ( Coordinate ( GI % nTotalCells ) )
    allocate ( Value ( GI % nTotalCells ) )

    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % MakeDirectory ( GI % Directory )
      
    SiloOptionList = DB_F77NULL

    call GI % WriteHeader ( TimeOption, CycleNumberOption )
    
    if ( GI % Stream % Parallel ) then

      call CO_Coordinate % Initialize &
             ( GI % Stream % Communicator, nOutgoing = [ GI % nTotalCells ], &
               nIncoming &
                 = [ GI % Stream % Communicator % Size * GI % nTotalCells ], &
               RootOption = 0 )

      CO_Coordinate % Outgoing % Value = GI % NodeCoordinate_1
      
      call CO_Coordinate % Gather ( )

      call CO_Variable % Initialize &
             ( GI % Stream % Communicator, nOutgoing = [ GI % nTotalCells ], & 
               nIncoming &
                 = [ GI % Stream % Communicator % Size * GI % nTotalCells ], &
               RootOption = 0 )
      
      if ( GI % Stream % IsWritable ( CheckMultiMeshOption = .true. ) ) then
        allocate ( Coordinate_MM ( size ( CO_Coordinate % Incoming % Value ) ) )
        allocate ( Value_MM ( size ( CO_Variable % Incoming % Value ) ) )
      end if !-- MultiMesh

      do iStrg = 1, GI % nStorages
        
        associate ( S => GI % Storage ( iStrg ) )

        call Show ( 'Writing a Storage (curve)', CONSOLE % INFO_5 )
        call Show ( iStrg, 'iStorage', CONSOLE % INFO_5 )
        call Show ( S % Name, 'Name', CONSOLE % INFO_5 )

        call GI % Stream % MakeDirectory ( S % Name )
        
        do iS = 1, S % nVariables

          iVrbl = S % iaSelected ( iS ) 
            
          call Show ( 'Writing a Variable (curve)', CONSOLE % INFO_6 )
          call Show ( iS, 'iSelected', CONSOLE % INFO_6 )

          nSiloOptions &
            = count &
                ( [ len_trim ( GI % CoordinateLabel ( 1 ) ) > 0, &
                    len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0, &
                    len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ] ) + 1
            
          Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
            
          if ( len_trim ( GI % CoordinateLabel ( 1 ) ) > 0 ) &
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_XLABEL, &
                        trim ( GI % CoordinateLabel ( 1 ) ), &
                        len_trim ( GI % CoordinateLabel ( 1 ) ) )
          Error = DBADDCOPT &
                    ( SiloOptionList, DBOPT_YLABEL, &
                      trim ( S % Variable ( iVrbl ) ), &
                      len_trim ( S % Variable ( iVrbl ) ) )
          if ( len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0 ) & 
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_XUNITS, &
                        trim ( GI % CoordinateUnit ( 1 ) % Label ), &
                        len_trim ( GI % CoordinateUnit ( 1 ) % Label ) )
          if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) &
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_YUNITS, &
                        trim ( S % Unit ( iVrbl ) % Label ), &
                        len_trim ( S % Unit ( iVrbl ) % Label ) )
              
          CO_Variable % Outgoing % Value &
            = S % Value &
                ( GI % oValue + 1 : GI % oValue + GI % nTotalCells, iVrbl )
          call CO_Variable % Gather ( )

          call Show ( trim ( S % Variable ( iVrbl ) ), 'Variable', &
                      CONSOLE % INFO_6 )
          call Show ( S % lVariable ( iVrbl ), 'lVariable', CONSOLE % INFO_6 )
          call Show ( nSiloOptions, 'nSiloOptions', CONSOLE % INFO_6 )
          if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) then
            call Show ( trim ( S % Unit ( iVrbl ) % Label ), 'Unit', &
                        CONSOLE % INFO_6 )
            call Show ( len_trim ( S % Unit ( iVrbl ) % Label ), 'lUnit', &
                        CONSOLE % INFO_6 )
          end if

          !-- FIXME: GCC 10.1 does not like some array operations
          do iV = 1, size ( Value )
            Coordinate ( iV )  =  GI % NodeCoordinate_1 ( iV ) &
                                  / GI % CoordinateUnit ( 1 ) % Number
            Value ( iV )  =  S % Value ( GI % oValue + iV, iVrbl ) &
                             / S % Unit ( iVrbl ) % Number
          end do !-- iV

          ! Error = DBPUTCURVE &
          !           ( GI % Stream % MeshBlockHandle, &
          !             trim ( S % Variable ( iVrbl ) ), &
          !             S % lVariable ( iVrbl ), &
          !             GI % NodeCoordinate_1 &
          !               / GI % CoordinateUnit ( 1 ) % Number, &
          !             S % Value ( GI % oValue + 1 &
          !                            : GI % oValue + GI % nTotalCells, &
          !                          iVrbl ) &
          !               / S % Unit ( iVrbl ) % Number, &
          !             DB_DOUBLE, GI % nTotalCells, SiloOptionList, Error )
          Error = DBPUTCURVE &
                    ( GI % Stream % MeshBlockHandle, &
                      trim ( S % Variable ( iVrbl ) ), &
                      S % lVariable ( iVrbl ), &
                      Coordinate, Value, &
                      DB_DOUBLE, GI % nTotalCells, SiloOptionList, Error )
          
          if ( GI % Stream % IsWritable ( CheckMultiMeshOption = .true. ) )&
          then

            !-- FIXME: GCC 10.1 does not like some array operations
            do iV = 1, size ( Value_MM )
              Coordinate_MM ( iV )  &
                =  CO_Coordinate % Incoming % Value ( iV )  &
                   / GI % CoordinateUnit ( 1 ) % Number
              Value_MM ( iV )  &
                =  CO_Variable % Incoming % Value ( iV )  & 
                   / S % Unit ( iVrbl ) % Number
            end do !-- iV

            ! associate &
            !   ( Coordinate    => CO_Coordinate % Incoming % Value, &
            !     VariableValue => CO_Variable % Incoming % Value )

            ! Error = DBPUTCURVE &
            !           ( GI % Stream % MultiMeshHandle, &
            !             trim ( S % Variable ( iVrbl ) ), &
            !             S % lVariable ( iVrbl ), &
            !             Coordinate / GI % CoordinateUnit ( 1 ) % Number, &
            !             VariableValue / S % Unit ( iVrbl ) % Number, &
            !             DB_DOUBLE, size ( VariableValue ), &
            !             SiloOptionList, Error )

            ! end associate
            
            Error = DBPUTCURVE &
                      ( GI % Stream % MultiMeshHandle, &
                        trim ( S % Variable ( iVrbl ) ), &
                        S % lVariable ( iVrbl ), &
                        Coordinate_MM, Value_MM, &
                        DB_DOUBLE, size ( Value_MM ), SiloOptionList, Error )

          end if
            
          if ( SiloOptionList /= DB_F77NULL ) then
            call Show ( SiloOptionList, 'SiloOptionList before free', &
                        CONSOLE % INFO_6 )        
            Error = DBFREEOPTLIST ( SiloOptionList )
            call Show ( SiloOptionList, 'SiloOptionList after free', &
                        CONSOLE % INFO_6 )        
            SiloOptionList = DB_F77NULL
            call Show ( SiloOptionList, 'SiloOptionList after nullification', &
                        CONSOLE % INFO_6 )        
            call Show ( DB_F77NULL, 'DB_F77NULL', CONSOLE % INFO_6 )
          end if
                
        end do
          
        call GI % Stream % ChangeDirectory ( '../' )
          
        end associate
        
      end do

    else !-- .not. Parallel
    
      do iStrg = 1, GI % nStorages
        
        associate ( S => GI % Storage ( iStrg ) )
        
        call Show ( 'Writing a Storage (curve)', CONSOLE % INFO_5 )
        call Show ( iStrg, 'iStorage', CONSOLE % INFO_5 )
        call Show ( S % Name, 'Name', CONSOLE % INFO_5 )

        call GI % Stream % MakeDirectory ( S % Name )
        
        do iS = 1 , S % nVariables

          iVrbl = S % iaSelected ( iS )
            
          call Show ( 'Writing a Variable (curve)', CONSOLE % INFO_6 )
          call Show ( iS, 'iSelected', CONSOLE % INFO_6 )

          nSiloOptions &
            = count &
                ( [ len_trim ( GI % CoordinateLabel ( 1 ) ) > 0, &
                    len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0, &
                    len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ] ) + 1
            
          Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
            
          if ( len_trim ( GI % CoordinateLabel ( 1 ) ) > 0 ) &
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_XLABEL, &
                        trim ( GI % CoordinateLabel ( 1 ) ), &
                        len_trim ( GI % CoordinateLabel ( 1 ) ) )
          Error = DBADDCOPT &
                    ( SiloOptionList, DBOPT_YLABEL, &
                      trim ( S % Variable ( iVrbl ) ), &
                      len_trim ( S % Variable ( iVrbl ) ) )
          if ( len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0 ) &
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_XUNITS, &
                        trim ( GI % CoordinateUnit ( 1 ) % Label ), &
                        len_trim ( GI % CoordinateUnit ( 1 ) % Label ) )
          if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) &
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_YUNITS, &
                        trim ( S % Unit ( iVrbl ) % Label ), &
                        len_trim ( S % Unit ( iVrbl ) % Label ) )

          call Show ( trim ( S % Variable ( iVrbl ) ), 'Variable', &
                      CONSOLE % INFO_6 )
          call Show ( S % lVariable ( iVrbl ), 'lVariable', CONSOLE % INFO_6 )
          call Show ( nSiloOptions, 'nSiloOptions', CONSOLE % INFO_6 )
          if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) then
            call Show ( trim ( S % Unit ( iVrbl ) % Label ), 'Unit', &
                        CONSOLE % INFO_6 )
            call Show ( len_trim ( S % Unit ( iVrbl ) % Label ), 'lUnit', &
                        CONSOLE % INFO_6 )
          end if

          !-- FIXME: GCC 10.1 does not like some array operations
          do iV = 1, size ( Value )
            Coordinate ( iV )  =  GI % NodeCoordinate_1 ( iV ) &
                                  / GI % CoordinateUnit ( 1 ) % Number
            Value ( iV )  =  S % Value ( GI % oValue + iV, iVrbl ) &
                             / S % Unit ( iVrbl ) % Number
          end do !-- iV

          ! Error = DBPUTCURVE &
          !           ( GI % Stream % MeshBlockHandle, &
          !             trim ( S % Variable ( iVrbl ) ), &
          !             S % lVariable ( iVrbl ), &
          !             GI % NodeCoordinate_1 &
          !               / GI % CoordinateUnit ( 1 ) % Number, &
          !             S % Value ( GI % oValue + 1 &
          !                            : GI % oValue + GI % nTotalCells, &
          !                          iVrbl ) &
          !               / S % Unit ( iVrbl ) % Number, &
          !             DB_DOUBLE, GI % nTotalCells, SiloOptionList, Error )
          Error = DBPUTCURVE &
                    ( GI % Stream % MeshBlockHandle, &
                      trim ( S % Variable ( iVrbl ) ), &
                      S % lVariable ( iVrbl ), &
                      Coordinate, Value, &
                      DB_DOUBLE, GI % nTotalCells, SiloOptionList, Error )
            
          if ( SiloOptionList /= DB_F77NULL ) then
            call Show ( SiloOptionList, 'SiloOptionList before free', &
                        CONSOLE % INFO_6 )        
            Error = DBFREEOPTLIST ( SiloOptionList )
            call Show ( SiloOptionList, 'SiloOptionList after free', &
                        CONSOLE % INFO_6 )        
            SiloOptionList = DB_F77NULL
            call Show ( SiloOptionList, 'SiloOptionList after nullification', &
                        CONSOLE % INFO_6 )        
            call Show ( DB_F77NULL, 'DB_F77NULL', CONSOLE % INFO_6 )
          end if
          
        end do
          
        call GI % Stream % ChangeDirectory ( '../' )
          
        end associate

      end do
      
    end if
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )

  end subroutine Write 
  
  
  subroutine Read ( GI, StorageOnlyOption, TimeOption, CycleNumberOption )
  
    class ( CurveImageForm ), intent ( inout ) :: &
      GI
    logical ( KDL ), intent ( in ), optional :: &
      StorageOnlyOption
    type ( QuantityForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iStrg, &     !-- iStorage
      iS, &     !-- iSelected
      iVrbl, &  !-- iVariable
      nSiloOptions, &
      nTotalCells, &
      nVariables, &
      Error, &
      DataType, &
      SiloOptionList
    real ( KDR ), dimension ( 1 ) :: &
      X_Scratch, &
      Y_Scratch
    logical ( KDL ) :: &
      StorageOnly
    character ( LDL ), dimension ( : ), allocatable :: &
      StorageName, &
      VariableName
    character ( LDF ) :: &
      WorkingDirectory
    
    if ( .not. GI % Stream % IsReadable ( ) ) return
    
    StorageOnly = .false.
    if ( present ( StorageOnlyOption ) ) &
      StorageOnly = StorageOnlyOption
    
    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % ChangeDirectory ( GI % Directory )
    
    SiloOptionList = DB_F77NULL
    
    call GI % ReadHeader ( TimeOption, CycleNumberOption )

    !-- prepare Storage to read into
    if ( GI % nStorages == 0 .and. .not. StorageOnly ) then
      call GI % Stream % ListContents ( ContentTypeOption = 'Directory' )
      GI % nStorages = size ( GI % Stream % ContentList )
!-- FIXME: NAG 5.3.1 should support sourced allocation
!      allocate ( StorageName, source = GI % Stream % ContentList )
      allocate ( StorageName ( size ( GI % Stream % ContentList ) ) )
      StorageName = GI % Stream % ContentList
      
      do iStrg = 1, GI % nStorages
      
        if ( len_trim ( StorageName ( iStrg ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( StorageName ( iStrg ) )
        call GI % Stream % ListContents &
               ( ContentTypeOption = 'Curve' )
        if ( allocated ( VariableName ) ) deallocate ( VariableName )
!        allocate ( VariableName, source = GI % Stream % ContentList )
        allocate ( VariableName ( size ( GI % Stream % ContentList ) ) )
        VariableName = GI % Stream % ContentList 
        nVariables = size ( GI % Stream % ContentList )
        if ( nVariables == 0 ) then
          GI % nStorages = 0
        else
          !-- read the first curve varriable to get the nTotalCells
          Error = DBGETCURVE &
                    ( GI % Stream % MeshBlockHandle, &
                      trim ( VariableName ( 1 ) ), &
                      len_trim ( VariableName ( 1 ) ), 1, &
                      X_Scratch, Y_Scratch, DataType, GI % nTotalCells )

          if ( allocated ( GI % NodeCoordinate_1 ) ) &
            deallocate ( GI % NodeCoordinate_1 )
          allocate ( GI % NodeCoordinate_1 ( GI % nTotalCells ) )

          call GI % Storage ( iStrg ) % Initialize &
                 ( [ GI % oValue + GI % nTotalCells, nVariables ], &
                     VariableOption = VariableName, &
                     NameOption = StorageName ( iStrg ) )
        end if 
        
        if ( len_trim ( StorageName ( iStrg ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( '..' )
          
      end do
    end if
    
    do iStrg = 1, GI % nStorages

      associate ( S => GI % Storage ( iStrg ) )

      call Show ( 'Reading a Storage (curve)', CONSOLE % INFO_5 )
      call Show ( iStrg, 'iStorage', CONSOLE % INFO_5 )
      call Show ( S % Name, 'Name', CONSOLE % INFO_5 )

      call GI % Stream % ChangeDirectory ( S % Name )
      do iS = 1, S % nVariables

        iVrbl = S % iaSelected ( iS )

        call Show ( 'Reading a Variable (curve)', CONSOLE % INFO_6 )
        call Show ( iS, 'iSelected', CONSOLE % INFO_6 )

        call Show ( trim ( S % Variable ( iVrbl ) ), 'Variable', &
                    CONSOLE % INFO_6 )
        call Show ( S % lVariable ( iVrbl ), 'lVariable', CONSOLE % INFO_6 )
        
        associate &
          ( S_V_Proper &
              => S % Value ( GI % oValue + 1 &
                               : GI % oValue + GI % nTotalCells, : ) )
         
        Error = DBGETCURVE &
                  ( GI % Stream % MeshBlockHandle, &
                    trim ( S % Variable ( iVrbl ) ), &
                    len_trim ( S % Variable ( iVrbl ) ), &
                    GI % nTotalCells, GI % NodeCoordinate_1, &
                    S_V_Proper ( :, iVrbl ), &
                    DataType, nTotalCells )

        !-- FIXME: An assumption is made that the unit used to write
        !          and read are the same. A better way would be to read
        !          the unit directly from Silo file.
        S_V_Proper ( :, iVrbl ) &
          = S_V_Proper ( :, iVrbl ) * S % Unit ( iVrbl ) % Number
          
        end associate !-- S_V_Proper

      end do

      call GI % Stream % ChangeDirectory ( '..' )

      end associate   !-- S

    end do
    
    !-- FIXME: Here we make the assumption that the CoordinateUnit for reading 
    !          is the same as the ones used to write. A better way would be to
    !          read the unit directly from Silo attribute. 
    GI % NodeCoordinate_1 &
      = GI % NodeCoordinate_1 * GI % CoordinateUnit ( 1 ) % Number
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )
     
  end subroutine Read
  
  
  subroutine ClearGrid ( CI )

    class ( CurveImageForm ), intent ( inout ) :: &
      CI

    if ( allocated ( CI % NodeCoordinate_3 ) ) &
      deallocate ( CI % NodeCoordinate_3 ) 
    if ( allocated ( CI % NodeCoordinate_2 ) ) &
      deallocate ( CI % NodeCoordinate_2 ) 
    if ( allocated ( CI % NodeCoordinate_1 ) ) &
      deallocate ( CI % NodeCoordinate_1 ) 
    
  end subroutine ClearGrid


  impure elemental subroutine Finalize ( CI )
    
    type ( CurveImageForm ), intent ( inout ) :: &
      CI
    
    nullify ( CI % Stream )
    
    if ( allocated ( CI % Storage ) ) &
      deallocate ( CI % Storage )

    call CI % ClearGrid ( )

  end subroutine Finalize
  
end module CurveImage_Form
