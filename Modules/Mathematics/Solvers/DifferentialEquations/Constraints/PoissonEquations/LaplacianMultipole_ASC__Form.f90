module LaplacianMultipole_ASC__Form

  !-- LaplacianMultipole_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, extends ( LaplacianMultipoleTemplate ) :: &
    LaplacianMultipole_ASC_Form
      integer ( KDI ) :: &
        COORDINATE_SYSTEM = 0
      real ( KDR ), dimension ( :, :, : ), pointer :: &
        Rectangular_X => null ( ), &
        Rectangular_Y => null ( ), &
        Rectangular_Z => null ( )
      real ( KDR ), dimension ( :, :, :, : ), pointer :: &
        SolidHarmonic_RC => null ( ), &
        SolidHarmonic_IC => null ( ), &
        SolidHarmonic_RS => null ( ), &
        SolidHarmonic_IS => null ( )
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
    procedure, private, pass :: &
      ComputeMomentContributions
  end type LaplacianMultipole_ASC_Form

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      ComputeRectangularCoordinates_CSL_Kernel, &
      ComputeMomentContributions_CSL_Kernel


    interface

      module subroutine ComputeRectangularCoordinates_CSL_Kernel &
                          ( X, Y, Z, IsProperCell, X_1, X_2, X_3, &
                            COORDINATE_SYSTEM, nC, CoordinateError, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          X, Y, Z
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          X_1, X_2, X_3
        integer ( KDI ), intent ( in ) :: &
          COORDINATE_SYSTEM, &
          nC  !-- nCells
        logical ( KDL ), intent ( out ) :: &
          CoordinateError
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeRectangularCoordinates_CSL_Kernel

      module subroutine ComputeMomentContributions_CSL_Kernel &
                          ( MyM_RC, MyM_RS, MyM_IC, MyM_IS, L, M, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MyM_RC, MyM_IC, &  !--MyMoments
          MyM_RS, MyM_IS
        integer ( KDI ), intent ( in ) :: &
          L, &
          M
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeMomentContributions_CSL_Kernel

    end interface


    private :: &
      AssignRectangularPointers, &
      AssignSolidHarmonicPointers

    integer ( KDI ), private, parameter :: &
      RECTANGULAR = 1, &
      CYLINDRICAL = 2, &
      SPHERICAL   = 3

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

    nullify ( L % SolidHarmonic_IS )
    nullify ( L % SolidHarmonic_RS )
    nullify ( L % SolidHarmonic_IC )
    nullify ( L % SolidHarmonic_RC )

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

    logical ( KDL ) :: &
      CoordinateError
    class ( GeometryFlatForm ), pointer :: &
      G

    allocate ( L % RectangularCoordinates )
    associate ( RC => L % RectangularCoordinates )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      associate ( nC => product ( C % nCells ) )

      call RC % Initialize ( [ nC, 3 ], ClearOption = .true. )
             !-- 3: X, Y, Z
      if ( L % UseDevice ) then
        call RC % AllocateDevice ( )
        call RC % UpdateDevice ( )
      end if

      call AssignRectangularPointers &
             ( RC % Value ( :, 1 ), RC % Value ( :, 2 ), RC % Value ( :, 3 ), &
               C % nCells, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z )

      select case ( trim ( C % CoordinateSystem ) )
      case ( 'RECTANGULAR' )
        L % COORDINATE_SYSTEM  =  RECTANGULAR
      case ( 'CYLINDRICAL' )
        L % COORDINATE_SYSTEM  =  CYLINDRICAL
      case ( 'SPHERICAL' )
        L % COORDINATE_SYSTEM  =  SPHERICAL
      end select

      G  =>  C % Geometry ( )
!      module subroutine ComputeRectangularCoordinates_CSL_Kernel &
!                          ( X, Y, Z, IsProperCell, X_1, X_2, X_3, &
!                            COORDINATE_SYSTEM, nC, CoordinateError, &
!                            UseDeviceOption )
      call ComputeRectangularCoordinates_CSL_Kernel &
             ( RC % Value ( :, 1 ), RC % Value ( :, 2 ), RC % Value ( :, 3 ), &
               C % IsProperCell, &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               G % Value ( :, G % CENTER_U ( 2 ) ), &
               G % Value ( :, G % CENTER_U ( 3 ) ), &
               L % COORDINATE_SYSTEM, nC, CoordinateError, &
               UseDeviceOption = L % UseDevice )
      if ( CoordinateError ) then
        call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
        call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'AllocateRectangularCoordinates', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

!call Show ( C % IsProperCell, '>>> IsProperCell' )
if ( L % UseDevice ) then
  RC % Value = huge ( 1.0_KDR )
  call RC % UpdateHost ( )
end if
!call Show ( C % IsProperCell, '>>> IsProperCell' )
call Show ( G % Value ( 9387 : 9396, G % CENTER_U ( 1 ) ), '>>> R ( 1 : 10 )' )
call Show ( G % Value ( 9387 : 9396, G % CENTER_U ( 2 ) ), '>>> Theta ( 1 : 10 )' )
call Show ( G % Value ( 9387 : 9396, G % CENTER_U ( 3 ) ), '>>> Phi ( 1 : 10 )' )
call Show ( RC % Value ( 9387 : 9396, 1 ), '>>> X ( 1 : 10 )' )
call Show ( RC % Value ( 9387 : 9396, 2 ), '>>> Y ( 1 : 10 )' )
call Show ( RC % Value ( 9387 : 9396, 3 ), '>>> Z ( 1 : 10 )' )
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

      associate ( nC_4 => product ( C % nCells ) )
      call SH % Initialize ( [ nC_4, 4 ] )
             !-- nC_4: nCells * ( Current, Previous_1, Previous_2, &
             !                    PreviousDiagonal )
             !-- 4: RegularCos, IrregularCos, RegularSin, IrregularSin
      if ( L % UseDevice ) &
        call SH % AllocateDevice ( )
      end associate !-- nC_4

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'AllocateSolidHarmonics', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- SH

  end subroutine AllocateSolidHarmonics


  subroutine ComputeMomentContributions ( L, Source )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L
    type ( StorageForm ), intent ( in ) :: &
      Source !-- array over levels    

    call ComputeMomentContributions_CSL_Kernel &
           ( L % MyM_RC, L % MyM_IC, L % MyM_RS, L % MyM_IS, &
             L % MaxDegree, L % MaxOrder, &
             UseDeviceOption = L % UseDevice )

  end subroutine ComputeMomentContributions


  subroutine AssignRectangularPointers &
               ( X_1D, Y_1D, Z_1D, nCells, X_3D, Y_3D, Z_3D )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      X_1D, Y_1D, Z_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCells
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      X_3D, Y_3D, Z_3D

    associate &
      ( n1  =>  nCells ( 1 ), &
        n2  =>  nCells ( 2 ), &
        n3  =>  nCells ( 3 ) )

    X_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>  X_1D
    Y_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>  Y_1D
    Z_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>  Z_1D

    end associate !-- n1, etc.

  end subroutine AssignRectangularPointers


  subroutine AssignSolidHarmonicPointers &
               ( SH_RC_1D, SH_IC_1D, SH_RS_1D, SH_IS_1D, nCells, &
                 SH_RC_4D, SH_IC_4D, SH_RS_4D, SH_IS_4D )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      SH_RC_1D, SH_IC_1D, &
      SH_RS_1D, SH_IS_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCells
    real ( KDR ), dimension ( :, :, :, : ), intent ( out ), pointer :: &
      SH_RC_4D, SH_IC_4D, &
      SH_RS_4D, SH_IS_4D

    associate &
      ( n1  =>  nCells ( 1 ), &
        n2  =>  nCells ( 2 ), &
        n3  =>  nCells ( 3 ) )

    SH_RC_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_RC_1D
    SH_IC_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_IC_1D
    SH_RS_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_RS_1D
    SH_IS_4D ( 1 : n1, 1 : n2, 1 : n3, 1 : 4 )  =>  SH_IS_1D
      !-- 1 : 4 Current, Previous_1, Previous_2, PreviousDiagonal 
    end associate !-- n1, etc.

  end subroutine AssignSolidHarmonicPointers


end module LaplacianMultipole_ASC__Form
