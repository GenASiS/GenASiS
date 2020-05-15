module LaplacianMultipole_ASC__Form

  !-- LaplacianMultipole_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, extends ( LaplacianMultipoleTemplate ) :: &
    LaplacianMultipole_ASC_Form
      real ( KDR ), dimension ( :, :, : ), pointer :: &
        Rectangular_X => null ( ), &
        Rectangular_Y => null ( ), &
        Rectangular_Z => null ( ), &
        Radius        => null ( ), &
        RadiusSquared => null ( ), &
        Volume => null ( )
      real ( KDR ), dimension ( :, :, :, : ), pointer :: &
        SolidHarmonic_RC => null ( ), &
        SolidHarmonic_IC => null ( ), &
        SolidHarmonic_RS => null ( ), &
        SolidHarmonic_IS => null ( ), &
        Source => null ( )
      type ( StorageForm ), allocatable :: &
        RectangularCoordinates, &
        SolidHarmonics
      class ( ChartTemplate ), pointer :: &
        Chart => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetParametersAtlas
    procedure, private, pass :: &
      AllocateRectangularCoordinates
    procedure, private, pass :: &
      AllocateSolidHarmonics
    procedure, public, pass :: &
      ComputeSolidHarmonics_0_0
    procedure, public, pass :: &
      ComputeSolidHarmonics_iM_iM
    procedure, public, pass :: &
      ComputeSolidHarmonics_iL_iM_1
    procedure, public, pass :: &
      ComputeSolidHarmonics_iL_iM_2
    procedure, public, pass :: &
      ComputeMomentLocalAtlas
  end type LaplacianMultipole_ASC_Form

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      Compute_RC_CSL_S_Kernel, &
      Compute_SH_0_0_CSL_Kernel, &
      Compute_SH_iM_iM_CSL_Kernel, &
      Compute_SH_iL_iM_1_CSL_Kernel, &
      Compute_SH_iL_iM_2_CSL_Kernel, &
      ComputeMomentLocal_CSL_S_Kernel

    interface

      module subroutine Compute_RC_CSL_S_Kernel &
                          ( X, Y, Z, D, D_2, IsProperCell, X_1, X_2, X_3, &
                            nC, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          X, Y, Z, D, D_2
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          X_1, X_2, X_3
        integer ( KDI ), intent ( in ) :: &
          nC  !-- nCellsBrick
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_RC_CSL_S_Kernel

      module subroutine Compute_SH_0_0_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nC, oC, iSH_0, iSH_PD, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC
        integer ( KDI ), intent ( in ) :: &
          iSH_0, iSH_PD
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_0_0_CSL_Kernel

      module subroutine Compute_SH_iM_iM_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nC, oC, iM, iSH_0, iSH_PD, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC
        integer ( KDI ), intent ( in ) :: &
          iM, &
          iSH_0, iSH_PD
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_iM_iM_CSL_Kernel

      module subroutine Compute_SH_iL_iM_1_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nC, oC, iM, iSH_0, iSH_1, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, &
          R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC
        integer ( KDI ), intent ( in ) :: &
          iM, &
          iSH_0, iSH_1
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_iL_iM_1_CSL_Kernel

      module subroutine Compute_SH_iL_iM_2_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nC, oC, iL, iM, iSH_0, iSH_1, iSH_2, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, &
          R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC
        integer ( KDI ), intent ( in ) :: &
          iL, iM, &
          iSH_0, iSH_1, iSH_2
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_iL_iM_2_CSL_Kernel

      module subroutine ComputeMomentLocal_CSL_S_Kernel &
                          ( MyM_RC, MyM_IC, MyM_RS, MyM_IS, S, &
                             SH_RC,  SH_IC,  SH_RS,  SH_IS, dV, &
                            nC, oC, nE, oR, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          MyM_RC, MyM_IC, MyM_RS, MyM_IS  !-- MyMoments
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          S  !-- Source
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          SH_RC, SH_IC, SH_RS, SH_IS, &  !-- SolidHarmonics
          dV
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC
        integer ( KDI ), intent ( in ) :: &
          nE, &  !-- nEquations
          oR     !-- oRadius
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeMomentLocal_CSL_S_Kernel

    end interface


    private :: &
      AssignRectangularPointers, &
      AssignSolidHarmonicPointers, &
      AssignSourcePointer

contains


  subroutine Initialize ( L, A, MaxDegree, nEquations )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    if ( L % Type == '' ) &
      L % Type = 'a LaplacianMultipole_ASC' 

    call L % InitializeTemplate ( A, MaxDegree, nEquations )

  end subroutine Initialize


  impure elemental subroutine Finalize ( L )

    type ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    nullify ( L % Chart )

    if ( allocated ( L % SolidHarmonics ) ) &
      deallocate ( L % SolidHarmonics )
    if ( allocated ( L % RectangularCoordinates ) ) &
      deallocate ( L % RectangularCoordinates )

    nullify ( L % Source )
    nullify ( L % SolidHarmonic_IS )
    nullify ( L % SolidHarmonic_RS )
    nullify ( L % SolidHarmonic_IC )
    nullify ( L % SolidHarmonic_RC )

    nullify ( L % Volume )
    nullify ( L % RadiusSquared )
    nullify ( L % Radius )
    nullify ( L % Rectangular_Z )
    nullify ( L % Rectangular_Y )
    nullify ( L % Rectangular_X )

    call L % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetParametersAtlas ( L, A )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A

    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( A )
    class is ( Atlas_SC_Form )

      if ( A % nDimensions  <  3 ) &
        L % MaxOrder  =  0
      if ( A % nDimensions  <  2 ) &
        L % MaxDegree  =  0

      G  =>  A % Geometry ( )
      L % UseDevice  =  G % AllocatedDevice

      L % Chart  =>  A % Chart

    end select !-- A

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )

      L % nRadialCells  =  C % nCells ( 1 )

      allocate ( L % RadialEdges )
      associate &
        (  RE  =>  L % RadialEdges, &
          nRC  =>  L % nRadialCells )

      call RE % Initialize ( [ nRC + 1, 1 ] )
      call RE % AllocateDevice ( )

      !--- Note indexing of Edge value begins at -1 with ghost cells
      RE % Value ( :, 1 )  =  C % Edge ( 1 ) % Value ( 1 : nRC + 1 )
      call RE % UpdateDevice ( )

      end associate !-- RE, etc. 
    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C
    
    nullify ( G )

  end subroutine SetParametersAtlas


  subroutine AllocateRectangularCoordinates ( L )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    class ( GeometryFlatForm ), pointer :: &
      G

    allocate ( L % RectangularCoordinates )
    associate ( RC => L % RectangularCoordinates )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      associate ( nC => product ( C % nCellsBrick  +  2 * C % nGhostLayers ) )

      call RC % Initialize ( [ nC, 5 ], ClearOption = .true. )
             !-- 5: X, Y, Z, R, R ** 2
      if ( L % UseDevice ) then
        call RC % AllocateDevice ( )
        call RC % UpdateDevice ( )
      end if

      G  =>  C % Geometry ( )

      call AssignRectangularPointers &
             ( RC % Value ( :, 1 ), RC % Value ( :, 2 ), &
               RC % Value ( :, 3 ), RC % Value ( :, 4 ), &
               RC % Value ( :, 5 ), G % Value ( :, G % VOLUME ), &
               C % nCellsBrick, C % nGhostLayers, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z, &
               L % Radius, L % RadiusSquared, L % Volume )

      select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )
        call Compute_RC_CSL_S_Kernel &
               ( RC % Value ( :, 1 ), RC % Value ( :, 2 ), &
                 RC % Value ( :, 3 ), RC % Value ( :, 4 ), &
                 RC % Value ( :, 5 ), &
                 C % IsProperCell, &
                 G % Value ( :, G % CENTER_U ( 1 ) ), &
                 G % Value ( :, G % CENTER_U ( 2 ) ), &
                 G % Value ( :, G % CENTER_U ( 3 ) ), &
                 nC, &
                 UseDeviceOption = L % UseDevice )
      case default
        call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
        call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'AllocateRectangularCoordinates', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select

      nullify ( G )

      end associate !-- nC

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateRectangularCoordinates', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- RC

  end subroutine AllocateRectangularCoordinates


  subroutine AllocateSolidHarmonics ( L )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    allocate ( L % SolidHarmonics )
    associate ( SH => L % SolidHarmonics )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      associate ( nC_4 => product ( C % nCellsBrick  +  2 * C % nGhostLayers ) &
                          * 4 )
      call SH % Initialize ( [ nC_4, 4 ], ClearOption = .true. )
             !-- nC_4: nCellsBrick * ( Current, Previous_1, Previous_2, &
             !                         PreviousDiagonal )
             !-- 4: RegularCos, IrregularCos, RegularSin, IrregularSin
      if ( L % UseDevice ) then
        call SH % AllocateDevice ( )
        call SH % UpdateDevice ( )
      end if
      end associate !-- nC_4

      call AssignSolidHarmonicPointers &
             ( SH % Value ( :, L %   REGULAR_COS ), &
               SH % Value ( :, L % IRREGULAR_COS ), &
               SH % Value ( :, L %   REGULAR_SIN ), &
               SH % Value ( :, L % IRREGULAR_SIN ), &
               C % nCellsBrick, C % nGhostLayers, &
               L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
               L % SolidHarmonic_RS, L % SolidHarmonic_IS )

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'AllocateSolidHarmonics', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- SH

  end subroutine AllocateSolidHarmonics


  subroutine ComputeSolidHarmonics_0_0 ( L, iSH_0, iSH_PD )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      iSH_0, iSH_PD

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      call Compute_SH_0_0_CSL_Kernel &
             ( L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
               L % SolidHarmonic_RS, L % SolidHarmonic_IS, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z, &
               L % RadiusSquared, C % nCellsBrick, C % nGhostLayers, &
               iSH_0, iSH_PD, UseDeviceOption = L % UseDevice )

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_SH_0_0', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeSolidHarmonics_0_0


  subroutine ComputeSolidHarmonics_iM_iM ( L, iM, iSH_0, iSH_PD )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      iM, &
      iSH_0, iSH_PD

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      call Compute_SH_iM_iM_CSL_Kernel &
             ( L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
               L % SolidHarmonic_RS, L % SolidHarmonic_IS, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z, &
               L % RadiusSquared, C % nCellsBrick, C % nGhostLayers, &
               iM, iSH_0, iSH_PD, &
               UseDeviceOption = L % UseDevice )

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_SH_iM_iM', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeSolidHarmonics_iM_iM


  subroutine ComputeSolidHarmonics_iL_iM_1 ( L, iM, iSH_0, iSH_1 )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      iM, &
      iSH_0, iSH_1

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      call Compute_SH_iL_iM_1_CSL_Kernel &
             ( L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
               L % SolidHarmonic_RS, L % SolidHarmonic_IS, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z, &
               L % RadiusSquared, C % nCellsBrick, C % nGhostLayers, &
               iM, iSH_0, iSH_1, &
               UseDeviceOption = L % UseDevice )

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_SH_iL_iM_1', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeSolidHarmonics_iL_iM_1


  subroutine ComputeSolidHarmonics_iL_iM_2 ( L, iL, iM, iSH_0, iSH_1, iSH_2 )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    integer ( KDI ), intent ( in ) :: &
      iL, iM, &
      iSH_0, iSH_1, iSH_2

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      call Compute_SH_iL_iM_2_CSL_Kernel &
             ( L % SolidHarmonic_RC, L % SolidHarmonic_IC, &
               L % SolidHarmonic_RS, L % SolidHarmonic_IS, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z, &
               L % RadiusSquared, C % nCellsBrick, C % nGhostLayers, &
               iL, iM, iSH_0, iSH_1, iSH_2, &
               UseDeviceOption = L % UseDevice )

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_SH_iL_iM_2', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeSolidHarmonics_iL_iM_2


  subroutine ComputeMomentLocalAtlas ( L, Source, iA, iSH_0 )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    class ( FieldAtlasTemplate ), intent ( in ) :: &
      Source
    integer ( KDI ), intent ( in ) :: &
      iA, &  
      iSH_0  

    class ( StorageForm ), pointer :: &
      Source_S

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

    select type ( Source )
    class is ( Storage_ASC_Form )

    Source_S => Source % Storage ( )
    call AssignSourcePointer &
           ( Source_S % Value, C % nCellsBrick, C % nGhostLayers, &
             L % nEquations, L % Source )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call ComputeMomentLocal_CSL_S_Kernel &
             ( L % MyMoment_RC ( :, :, iA ), L % MyMoment_IC ( :, :, iA ), &
               L % MyMoment_RS ( :, :, iA ), L % MyMoment_IS ( :, :, iA ), &
               L % Source, &
               L % SolidHarmonic_RC ( :, :, :, iSH_0 ), &
               L % SolidHarmonic_IC ( :, :, :, iSH_0 ), &
               L % SolidHarmonic_RS ( :, :, :, iSH_0 ), &
               L % SolidHarmonic_IS ( :, :, :, iSH_0 ), &
               L % Volume, &
               C % nCellsBrick, C % nGhostLayers, L % nEquations, &
               oR = ( C % iaBrick ( 1 ) - 1 ) * C % nCellsBrick ( 1 ), &
               UseDeviceOption = L % UseDevice )
    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    class default
      call Show ( 'Source type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Source

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    nullify ( Source_S )

  end subroutine ComputeMomentLocalAtlas


  subroutine AssignRectangularPointers &
               ( X_1D, Y_1D, Z_1D, D_1D, D_2_1D, dV_1D, nC, nG, &
                 X_3D, Y_3D, Z_3D, D_3D, D_2_3D, dV_3D )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      X_1D, Y_1D, Z_1D, &
      D_1D, D_2_1D, &
      dV_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nC, &  !-- nCellsBrick
      nG     !-- nGhostLayers
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      X_3D, Y_3D, Z_3D, &
      D_3D, D_2_3D, &
      dV_3D

    associate &
      ( n1  =>  nC ( 1 )  +  2 * nG ( 1 ), &
        n2  =>  nC ( 2 )  +  2 * nG ( 2 ), &
        n3  =>  nC ( 3 )  +  2 * nG ( 3 ) )

      X_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    X_1D
      Y_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    Y_1D
      Z_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    Z_1D
      D_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    D_1D
    D_2_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>  D_2_1D
     dV_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>   dV_1D

    end associate !-- n1, etc.

  end subroutine AssignRectangularPointers


  subroutine AssignSolidHarmonicPointers &
               ( SH_RC_1D, SH_IC_1D, SH_RS_1D, SH_IS_1D, nC, nG, &
                 SH_RC_4D, SH_IC_4D, SH_RS_4D, SH_IS_4D )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      SH_RC_1D, SH_IC_1D, &
      SH_RS_1D, SH_IS_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nC, &  !-- nCellsBrick
      nG     !-- nGhostLayers
    real ( KDR ), dimension ( :, :, :, : ), intent ( out ), pointer :: &
      SH_RC_4D, SH_IC_4D, &
      SH_RS_4D, SH_IS_4D

    associate &
      ( n1  =>  nC ( 1 )  +  2 * nG ( 1 ), &
        n2  =>  nC ( 2 )  +  2 * nG ( 2 ), &
        n3  =>  nC ( 3 )  +  2 * nG ( 3 ) )

    SH_RC_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_RC_1D
    SH_IC_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_IC_1D
    SH_RS_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_RS_1D
    SH_IS_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_IS_1D
      !-- 1 : 4 Current, Previous_1, Previous_2, PreviousDiagonal 

    end associate !-- n1, etc.

  end subroutine AssignSolidHarmonicPointers


  subroutine AssignSourcePointer ( S_2D, nC, nG, nE, S_4D )

    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
      S_2D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nC, &  !-- nCellsBrick
      nG     !-- nGhostLayers
    integer ( KDI ), intent ( in ) :: &
      nE  !-- nEquations
    real ( KDR ), dimension ( :, :, :, : ), intent ( out ), pointer :: &
      S_4D

    associate &
      ( n1  =>  nC ( 1 )  +  2 * nG ( 1 ), &
        n2  =>  nC ( 2 )  +  2 * nG ( 2 ), &
        n3  =>  nC ( 3 )  +  2 * nG ( 3 ) )

    S_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : nE )  =>  S_2D

    end associate !-- n1, etc.

  end subroutine AssignSourcePointer


end module LaplacianMultipole_ASC__Form
