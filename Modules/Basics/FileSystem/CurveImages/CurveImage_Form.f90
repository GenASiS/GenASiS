!-- CurveImagesSiloForm provides a class to generate one-dimensional
!   grid data on XY-plane (curve) as columns in Silo format, suitable for
!   VisIt visualization tool.

module CurveImage_Form

  use VariableManagement
  use Display
  use MessagePassing
  use GridImageBasics
        
  implicit none 
  private 
  
  include 'silo_f9x.inc'
  
  type, public, extends ( GridImageSiloTemplate ) :: CurveImageForm
  contains
    procedure, public, pass :: &
      SetGridUnigrid
    procedure, public, pass :: &
      SetGridRefinable
    generic :: &
      SetGrid => SetGridUnigrid, SetGridRefinable
    procedure, public, pass :: &
      SetReadAttributes
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

    
  subroutine SetGridUnigrid &
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
    type ( MeasuredValueForm ), intent ( in ), optional :: &
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

  end subroutine SetGridUnigrid
  
  
  subroutine SetGridRefinable &
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
    type ( MeasuredValueForm ), intent ( in ), optional :: &
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

  end subroutine SetGridRefinable
  
  
  subroutine SetReadAttributes ( GI, Directory, oValue )

    class ( CurveImageForm ), intent ( inout ) :: &
      GI 
    character ( * ), intent ( in ) :: &
      Directory
    integer ( KDI ), intent ( in ) :: &
      oValue
      
    GI % oValue      = oValue
    GI % nGhostCells = 0
    
    GI % lDirectory  = len_trim ( Directory )

    GI % Directory   = Directory

  end subroutine SetReadAttributes
   
  
  subroutine Write ( GI, TimeOption, CycleNumberOption )
  
    class ( CurveImageForm ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iG, &     !-- iGroup
      iS, &     !-- iSelected
      iVrbl, &  !-- iVariable
      nSiloOptions, &
      Error, &
      SiloOptionList
    character ( LDF ) :: &
      WorkingDirectory
    type ( CollectiveOperation_R_Form ), target :: &
      CO_Coordinate, &
      CO_Variable

    if ( .not. GI % Stream % IsWritable ( ) ) return
    
    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % MakeDirectory ( GI % Directory )
      
    SiloOptionList = DB_F77NULL

    call GI % WriteHeader ( TimeOption, CycleNumberOption )
    
    if ( GI % Stream % Parallel ) then

      call CO_Coordinate % Initialize &
             ( GI % Stream % Communicator, nOutgoing = [ GI % nTotalCells ], &
               nIncoming &
                 = GI % Stream % Communicator % Size * [ GI % nTotalCells ], &
               RootOption = 0 )
      
      CO_Coordinate % Outgoing % Value = GI % NodeCoordinate_1
      
      call CO_Coordinate % Gather ( )

      call CO_Variable % Initialize &
             ( GI % Stream % Communicator, nOutgoing = [ GI % nTotalCells ], & 
               nIncoming &
                 = GI % Stream % Communicator % Size * [ GI % nTotalCells ], &
               RootOption = 0 )
      
      do iG = 1, GI % nVariableGroups
        
        associate ( VG => GI % VariableGroup ( iG ) )

          call GI % Stream % MakeDirectory ( VG % Name )
        
          do iS = 1, VG % nVariables

            iVrbl = VG % iaSelected ( iS ) 
            
            nSiloOptions &
              = count &
                  ( [ len_trim ( GI % CoordinateLabel ( 1 ) ) > 0, &
                      len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0, &
                      len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ] ) + 1
            
            Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
            
            if ( len_trim ( GI % CoordinateLabel ( 1 ) ) > 0 ) &
              Error = DBADDCOPT &
                        ( SiloOptionList, DBOPT_XLABEL, &
                          trim ( GI % CoordinateLabel ( 1 ) ), &
                          len_trim ( GI % CoordinateLabel ( 1 ) ) )
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_YLABEL, &
                        trim ( VG % Variable ( iVrbl ) ), &
                        len_trim ( VG % Variable ( iVrbl ) ) )
            if ( len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0 ) & 
              Error = DBADDCOPT &
                        ( SiloOptionList, DBOPT_XUNITS, &
                          trim ( GI % CoordinateUnit ( 1 ) % Label ), &
                          len_trim ( GI % CoordinateUnit ( 1 ) % Label ) )
            if ( len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ) &
              Error = DBADDCOPT &
                        ( SiloOptionList, DBOPT_YUNITS, &
                          trim ( VG % Unit ( iVrbl ) % Label ), &
                          len_trim ( VG % Unit ( iVrbl ) % Label ) )
              
            CO_Variable % Outgoing % Value &
              = VG % Value &
                  ( GI % oValue + 1 : GI % oValue + GI % nTotalCells, iVrbl )
            call CO_Variable % Gather ( )

            Error = DBPUTCURVE &
                      ( GI % Stream % MeshBlockHandle, &
                        trim ( VG % Variable ( iVrbl ) ), &
                        VG % lVariable ( iVrbl ), &
                        GI % NodeCoordinate_1 &
                          / GI % CoordinateUnit ( 1 ) % Number, &
                        VG % Value ( GI % oValue + 1 &
                                       : GI % oValue + GI % nTotalCells, &
                                     iVrbl ) &
                          / VG % Unit ( iVrbl ) % Number, &
                        DB_DOUBLE, GI % nTotalCells, SiloOptionList, Error )
          
            if ( GI % Stream % IsWritable ( CheckMultiMeshOption = .true. ) )&
            then

              associate &
                ( Coordinate    => CO_Coordinate % Incoming % Value, &
                  VariableValue => CO_Variable % Incoming % Value )

              Error = DBPUTCURVE &
                        ( GI % Stream % MultiMeshHandle, &
                          trim ( VG % Variable ( iVrbl ) ), &
                          VG % lVariable ( iVrbl ), &
                          Coordinate / GI % CoordinateUnit ( 1 ) % Number, &
                          VariableValue / VG % Unit ( iVrbl ) % Number, &
                          DB_DOUBLE, size ( VariableValue ), &
                          SiloOptionList, Error )

              end associate
            
            end if
            
            if ( SiloOptionList /= DB_F77NULL ) then
              Error = DBFREEOPTLIST ( SiloOptionList )
              SiloOptionList = DB_F77NULL
            end if
                
          end do
          
          call GI % Stream % ChangeDirectory ( '../' )
          
        end associate
        
      end do

    else !-- .not. Parallel
    
      do iG = 1, GI % nVariableGroups
        
        associate ( VG => GI % VariableGroup ( iG ) )
        
          call GI % Stream % MakeDirectory ( VG % Name )
        
          do iS = 1 , VG % nVariables

            iVrbl = VG % iaSelected ( iS )
            
            nSiloOptions &
              = count &
                  ( [ len_trim ( GI % CoordinateLabel ( 1 ) ) > 0, &
                      len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0, &
                      len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ] ) + 1
            
            Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
            
            if ( len_trim ( GI % CoordinateLabel ( 1 ) ) > 0 ) &
              Error = DBADDCOPT &
                        ( SiloOptionList, DBOPT_XLABEL, &
                          trim ( GI % CoordinateLabel ( 1 ) ), &
                          len_trim ( GI % CoordinateLabel ( 1 ) ) )
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_YLABEL, &
                        trim ( VG % Variable ( iVrbl ) ), &
                        len_trim ( VG % Variable ( iVrbl ) ) )
            if ( len_trim ( GI % CoordinateUnit ( 1 ) % Label ) > 0 ) &
              Error = DBADDCOPT &
                        ( SiloOptionList, DBOPT_XUNITS, &
                          trim ( GI % CoordinateUnit ( 1 ) % Label ), &
                          len_trim ( GI % CoordinateUnit ( 1 ) % Label ) )
            if ( len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ) &
              Error = DBADDCOPT &
                        ( SiloOptionList, DBOPT_YUNITS, &
                          trim ( VG % Unit ( iVrbl ) % Label ), &
                          len_trim ( VG % Unit ( iVrbl ) % Label ) )

            Error = DBPUTCURVE &
                      ( GI % Stream % MeshBlockHandle, &
                        trim ( VG % Variable ( iVrbl ) ), &
                        VG % lVariable ( iVrbl ), &
                        GI % NodeCoordinate_1 &
                          / GI % CoordinateUnit ( 1 ) % Number, &
                        VG % Value ( GI % oValue + 1 &
                                       : GI % oValue + GI % nTotalCells, &
                                     iVrbl ) &
                          / VG % Unit ( iVrbl ) % Number, &
                        DB_DOUBLE, GI % nTotalCells, SiloOptionList, Error )
            
            if ( SiloOptionList /= DB_F77NULL ) then
              Error = DBFREEOPTLIST ( SiloOptionList )
              SiloOptionList = DB_F77NULL
            end if
          
          end do
          
          call GI % Stream % ChangeDirectory ( '../' )
          
        end associate

      end do
      
    end if
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )

  end subroutine Write 
  
  
  subroutine Read ( GI, TimeOption, CycleNumberOption )
  
    class ( CurveImageForm ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iG, &     !-- iGroup
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
    character ( LDL ), dimension ( : ), allocatable :: &
      GroupName, &
      VariableName
    character ( LDF ) :: &
      WorkingDirectory
    
    if ( .not. GI % Stream % IsReadable ( ) ) return
    
    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % ChangeDirectory ( GI % Directory )
    
    SiloOptionList = DB_F77NULL
    
    call GI % ReadHeader ( TimeOption, CycleNumberOption )
    
    !-- prepare VariableGroup to read into
    if ( GI % nVariableGroups == 0 ) then
      call GI % Stream % ListContents ( ContentTypeOption = 'Directory' )
      GI % nVariableGroups = size ( GI % Stream % ContentList )
!-- FIXME: NAG 5.3.1 should support sourced allocation
!      allocate ( GroupName, source = GI % Stream % ContentList )
      allocate ( GroupName ( size ( GI % Stream % ContentList ) ) )
      GroupName = GI % Stream % ContentList
      do iG = 1, GI % nVariableGroups
        if ( len_trim ( GroupName ( iG ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( GroupName ( iG ) )
        call GI % Stream % ListContents &
               ( ContentTypeOption = 'Curve' )
        if ( allocated ( VariableName ) ) deallocate ( VariableName )
!        allocate ( VariableName, source = GI % Stream % ContentList )
        allocate ( VariableName ( size ( GI % Stream % ContentList ) ) )
        VariableName = GI % Stream % ContentList 
        nVariables = size ( GI % Stream % ContentList )
        if ( nVariables == 0 ) then
          GI % nVariableGroups = 0
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

          call GI % VariableGroup ( iG ) % Initialize &
                 ( [ GI % oValue + GI % nTotalCells, nVariables ], &
                     VariableOption = VariableName, &
                     NameOption = GroupName ( iG ) )
          
          associate ( VG => GI % VariableGroup ( iG ) )
          
          do iVrbl = 1, nVariables
            Error = DBGETCURVE &
                      ( GI % Stream % MeshBlockHandle, &
                        trim ( VariableName ( iVrbl ) ), &
                        len_trim ( VariableName ( iVrbl ) ), &
                        GI % nTotalCells, GI % NodeCoordinate_1, &
                        VG % Value ( GI % oValue + 1 &
                                       : GI % oValue + GI % nTotalCells, &
                                     iVrbl ), &
                        DataType, nTotalCells )
            
            !-- FIXME: An assumption is made that the unit used to write
            !          and read are the same. A better way would be to read
            !          the unit directly from Silo file.
            VG % Value ( :, iVrbl ) &
              = VG % Value ( :, iVrbl ) * VG % Unit ( iVrbl ) % Number

          end do
          
          end associate
        
        end if
        
        if ( len_trim ( GroupName ( iG ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( '..' )
      end do
    end if
    
    !-- FIXME: Here we make the assumption that the CoordinateUnit for reading 
    !          is the same as the ones used to write. A better way would be to
    !          read the unit directly from Silo attribute. 

    GI % NodeCoordinate_1 &
      = GI % NodeCoordinate_1 * GI % CoordinateUnit ( 1 ) % Number
    
    call GI % Stream % ChangeDirectory ( '../' )
     
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
    
    if ( allocated ( CI % VariableGroup ) ) &
      deallocate ( CI % VariableGroup )

    call CI % ClearGrid ( )

  end subroutine Finalize
  
end module CurveImage_Form
