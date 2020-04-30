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
        RadiusSquared => null ( )
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
      ComputeSolidHarmonics_0_0
    procedure, private, pass :: &
      ComputeSolidHarmonics_iM_iM
    procedure, private, pass :: &
      ComputeSolidHarmonics_iL_iM_1
    procedure, private, pass :: &
      ComputeSolidHarmonics_iL_iM_2
  end type LaplacianMultipole_ASC_Form

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      Compute_RC_CSL_S_Kernel, &
      Compute_SH_0_0_CSL_Kernel, &
      Compute_SH_iM_iM_CSL_Kernel, &
      Compute_SH_iL_iM_1_CSL_Kernel, &
      Compute_SH_iL_iM_2_CSL_Kernel, &
      SumMomentContributions_CSL_SphericalKernel

    interface

      module subroutine Compute_RC_CSL_S_Kernel &
                          ( X, Y, Z, D_2, IsProperCell, X_1, X_2, X_3, &
                            nC, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          X, Y, Z, D_2
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          X_1, X_2, X_3
        integer ( KDI ), intent ( in ) :: &
          nC  !-- nCells
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_RC_CSL_S_Kernel

      module subroutine Compute_SH_0_0_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nCells, iSH_0, iSH_PD, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, &
          R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nCells
        integer ( KDI ), intent ( in ) :: &
          iSH_0, iSH_PD
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_0_0_CSL_Kernel

      module subroutine Compute_SH_iM_iM_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nCells, iM, iSH_0, iSH_PD, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, &
          R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nCells
        integer ( KDI ), intent ( in ) :: &
          iM, &
          iSH_0, iSH_PD
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_iM_iM_CSL_Kernel

      module subroutine Compute_SH_iL_iM_1_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nCells, iM, iSH_0, iSH_1, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, &
          R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nCells
        integer ( KDI ), intent ( in ) :: &
          iM, &
          iSH_0, iSH_1
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_iL_iM_1_CSL_Kernel

      module subroutine Compute_SH_iL_iM_2_CSL_Kernel &
                          ( R_C, I_C, R_S, I_S, X, Y, Z, D_2, &
                            nCells, iL, iM, iSH_0, iSH_1, iSH_2, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          R_C, I_C, &
          R_S, I_S
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          X, Y, Z, D_2
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nCells
        integer ( KDI ), intent ( in ) :: &
          iL, iM, &
          iSH_0, iSH_1, iSH_2
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SH_iL_iM_2_CSL_Kernel

      module subroutine SumMomentContributions_CSL_SphericalKernel &
                          ( )
        use Basics
        implicit none
      end subroutine SumMomentContributions_CSL_SphericalKernel

    end interface


    private :: &
      AssignRectangularPointers, &
      AssignSolidHarmonicPointers

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

    nullify ( L % RadiusSquared )
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

    class ( GeometryFlatForm ), pointer :: &
      G

    allocate ( L % RectangularCoordinates )
    associate ( RC => L % RectangularCoordinates )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      associate ( nC => product ( C % nCells ) )

      call RC % Initialize ( [ nC, 4 ], ClearOption = .true. )
             !-- 4: X, Y, Z, R ** 2
      if ( L % UseDevice ) then
        call RC % AllocateDevice ( )
        call RC % UpdateDevice ( )
      end if

      call AssignRectangularPointers &
             ( RC % Value ( :, 1 ), RC % Value ( :, 2 ), RC % Value ( :, 3 ), &
               RC % Value ( :, 4 ), C % nCells, &
               L % Rectangular_X, L % Rectangular_Y, L % Rectangular_Z, &
               L % RadiusSquared )

      G  =>  C % Geometry ( )
      select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )
        call Compute_RC_CSL_S_Kernel &
               ( RC % Value ( :, 1 ), RC % Value ( :, 2 ), &
                 RC % Value ( :, 3 ), RC % Value ( :, 4 ), &
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

      associate ( nC_4 => product ( C % nCells ) * 4 )
      call SH % Initialize ( [ nC_4, 4 ], ClearOption = .true. )
             !-- nC_4: nCells * ( Current, Previous_1, Previous_2, &
             !                    PreviousDiagonal )
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
               C % nCells, &
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
               L % RadiusSquared, C % nCells, iSH_0, iSH_PD, &
               UseDeviceOption = L % UseDevice )

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
               L % RadiusSquared, C % nCells, iM, iSH_0, iSH_PD, &
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
               L % RadiusSquared, C % nCells, iM, iSH_0, iSH_1, &
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
               L % RadiusSquared, C % nCells, iL, iM, iSH_0, iSH_1, iSH_2, &
               UseDeviceOption = L % UseDevice )

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_SH_iL_iM_2', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeSolidHarmonics_iL_iM_2


  subroutine AssignRectangularPointers &
               ( X_1D, Y_1D, Z_1D, D_2_1D, nCells, X_3D, Y_3D, Z_3D, D_2_3D )

    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      X_1D, Y_1D, Z_1D, D_2_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCells
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      X_3D, Y_3D, Z_3D, D_2_3D

    associate &
      ( n1  =>  nCells ( 1 ), &
        n2  =>  nCells ( 2 ), &
        n3  =>  nCells ( 3 ) )

      X_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    X_1D
      Y_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    Y_1D
      Z_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>    Z_1D
    D_2_3D ( 1 : n1, 1 : n2, 1 : n3 )  =>  D_2_1D

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
