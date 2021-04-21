module Geometry_F_GS__Form

  !-- Geometry_Flat_GridStructured_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometry_F_CH__Form

  implicit none
  private

  type, public, extends ( Geometry_F_CH_Form ) :: Geometry_F_GS_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Geometry_F_GS_Form

    private :: &
      SetCoordinates

contains


  subroutine Initialize &
               ( GG, G, NameOption, nFieldsOption, FieldOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption )

    class ( Geometry_F_GS_Form ), intent ( inout ) :: &
      GG
    class ( Grid_S_Form ), intent ( inout ), target :: &
      G
    character ( * ), intent ( inout ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    integer ( KDI ) :: &
      nFields
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    Name  =  ''
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    nFields  =  GG % N_FIELDS_F
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    if ( .not. allocated ( GG % FieldSet ) ) then
      allocate ( FieldSet_GS_Form :: GG % FieldSet )
      associate ( FS  =>  GG % FieldSet )
      FS % Type  =  'a Geometry_F_GS'    
      end associate !-- FS
    end if

    call GG % Geometry_F_CH_Form % Initialize_H &
           ( G, &
             NameOption = Name, &
             nFieldsOption = nFields, &
             FieldOption = Field, &
             UnitOption = Unit )

    select type ( FS  =>  GG % FieldSet )
    class is ( FieldSet_GS_Form )
    if ( FS % Type == 'a Geometry_F_GS' ) then
      call FS % Initialize &
             ( G, &
               FieldOption = Field, &
               NameOption = Name, &
               DeviceMemoryOption = DeviceMemoryOption, &
               PinnedMemoryOption = PinnedMemoryOption, &
               DevicesCommunicateOption = DevicesCommunicateOption, &
               UnitOption = Unit, &
               nFieldsOption = nFields )
    end if
    end select !-- FS

    call GG % Compute ( )

  end subroutine Initialize


  subroutine Compute ( GG )

    class ( Geometry_F_GS_Form ), intent ( inout ) :: &
      GG

    integer ( KDI ) :: &
      iD  !-- iDimension

    select type ( GFS  =>  GG % FieldSet )
    class is ( FieldSet_GS_Form )

    associate ( nD  =>  GFS % Chart % nDimensions )

    do iD = 1, nD
      call SetCoordinates ( GG, iD )
    end do !-- iD

    call GG % ComputeFromCoordinates ( GFS % Storage )

    end associate !-- nD
    end select !-- GFS

  end subroutine Compute


  impure elemental subroutine Finalize ( GG )

    type ( Geometry_F_GS_Form ), intent ( inout ) :: &
      GG

  end subroutine Finalize
  

  subroutine SetCoordinates ( GG, iD )

    class ( Geometry_F_GS_Form ), intent ( inout ) :: &
      GG
    integer ( KDI ), intent ( in ) :: &
      iD      !-- iDimension

    integer ( KDI ) :: &
      iaF, iaL, &  !-- iaFirst, iaLast
      iC, &        !-- iCell
      oC           !-- oCell
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Edge_I_3D, &
      Width_3D, &
      Center_3D

    select type ( GFS  =>  GG % FieldSet )
    class is ( FieldSet_GS_Form )

    select type ( G  =>  GFS % Chart )
    class is ( Grid_S_Form )

    associate ( GV  =>  GFS % Storage % Value )

    iaF  =  1  -  G % nGhostLayers ( iD ) 
    if ( G % Distributed ) then
      iaL  =  G % nCellsBrick ( iD )  +  G % nGhostLayers ( iD )
       oC  =  ( G % iaBrick ( iD )  -  1 )  *  G % nCellsBrick ( iD )
    else
      iaL  =  G % nCells ( iD )  +  G % nGhostLayers ( iD )
       oC  =  0
    end if

    call G % SetFieldPointer &
           ( GV ( :, GG % EDGE_I_U ( iD ) ), Edge_I_3D )
    call G % SetFieldPointer &
           ( GV ( :, GG % WIDTH_U ( iD ) ),  Width_3D )
    call G % SetFieldPointer &
           ( GV ( :, GG % CENTER_U ( iD ) ), Center_3D )

    associate &
      (   Edge_1D  =>  G %   Edge ( iD ) % Value, &
         Width_1D  =>  G %  Width ( iD ) % Value, &
        Center_1D  =>  G % Center ( iD ) % Value )
    do iC  =  iaF, iaL
      select case ( iD )
      case ( 1 )
        Edge_I_3D ( iC, :, : )  =    Edge_1D ( oC + iC )
         Width_3D ( iC, :, : )  =   Width_1D ( oC + iC )
        Center_3D ( iC, :, : )  =  Center_1D ( oC + iC )
      case ( 2 )
        Edge_I_3D ( :, iC, : )  =    Edge_1D ( oC + iC )
         Width_3D ( :, iC, : )  =   Width_1D ( oC + iC )
        Center_3D ( :, iC, : )  =  Center_1D ( oC + iC )
      case ( 3 )
        Edge_I_3D ( :, :, iC )  =    Edge_1D ( oC + iC )
         Width_3D ( :, :, iC )  =   Width_1D ( oC + iC )
        Center_3D ( :, :, iC )  =  Center_1D ( oC + iC )
      end select !-- iD
    end do !-- iC
    end associate !-- Edge_1D, etc.

    end associate !-- GV
    end select !-- G
    end select !-- GFS

  end subroutine SetCoordinates


end module Geometry_F_GS__Form
