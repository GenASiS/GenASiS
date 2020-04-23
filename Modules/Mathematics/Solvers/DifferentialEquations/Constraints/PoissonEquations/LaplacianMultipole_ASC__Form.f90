module LaplacianMultipole_ASC__Form

  !-- LaplacianMultipole_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, extends ( LaplacianMultipoleTemplate ) :: &
    LaplacianMultipole_ASC_Form
      class ( ChartTemplate ), pointer :: &
        Chart => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      SetParameters
    final :: &
      Finalize
    procedure, private, pass :: &
      ComputeMomentsLocal
  end type LaplacianMultipole_ASC_Form

    private :: &
      SetRadialEdgeSpherical

!-- FIXME: With GCC 6.1.0, must be public to trigger .smod generation
!    private :: &
    public :: &
      ComputeMomentsLocal_CSL_Kernel

    interface

      module subroutine ComputeMomentsLocal_CSL_Kernel &
                          ( MyM_RC, MyM_RS, MyM_IC, MyM_IS, &
                            CoordinateSystem, IsProperCell, Source, Origin, &
                            RadialEdge, Center_1, Center_2, Center_3, Volume, &
                            iaSource, MaxDegree, MaxOrder, nDimensions, &
                            nCells, nEquations, nAngularMomentCells, &
                            GridError, SH_RC, SH_IC, SH_RS, SH_IS, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MyM_RC, MyM_IC, &  !--MyMoments
          MyM_RS, MyM_IS
        character ( * ), intent ( in ) :: &
          CoordinateSystem
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          IsProperCell
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          Source
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Origin, &
          RadialEdge, &
          Center_1, Center_2, Center_3, &
          Volume
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaSource
        integer ( KDI ), intent ( in ) :: &
          MaxDegree, &
          MaxOrder, &
          nDimensions, &
          nCells, &
          nEquations, &
          nAngularMomentCells
        logical ( KDL ), intent ( out ) :: &
          GridError
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          SH_RC, SH_IC, &  !-- SolidHarmonics
          SH_RS, SH_IS  
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeMomentsLocal_CSL_Kernel

    end interface

contains


  subroutine Initialize ( LM, A, MaxDegree, nEquations )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    if ( LM % Type == '' ) &
      LM % Type = 'a LaplacianMultipole_ASC' 

    call LM % InitializeTemplate ( A, MaxDegree, nEquations )

    select case ( trim ( A % Chart % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call SetRadialEdgeSpherical ( LM )
      call LM % SetMomentStorage ( )
    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( A % Chart % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- CoordinateSystem

  end subroutine Initialize


  impure elemental subroutine Finalize ( LM )

    type ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM

    nullify ( LM % Chart )

    call LM % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetParameters ( LM, A, MaxDegree, nEquations )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    integer ( KDI ) :: &
      iL, &
      iM
    class ( GeometryFlatForm ), pointer :: &
      G

    LM % nEquations  =   nEquations
    LM % MaxDegree   =   MaxDegree
    LM % MaxOrder    =   MaxDegree

    select type ( A )
    class is ( Atlas_SC_Form )

      LM % Chart  =>  A % Chart

      G => A % Geometry ( )
      LM % UseDevice = G % AllocatedDevice

    end select !-- A

    select type ( C => LM % Chart )
    class is ( Chart_SL_Template )

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
       if ( C % nDimensions  <  3 ) &
         LM % MaxOrder  =  0
       if ( C % nDimensions  <  2 ) &
         LM % MaxDegree  =  0
    case default
      call Show ( 'CoordinateSystem not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- CoordinateSystem

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetParameters', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- C
    
    associate &
      ( L   => LM % MaxDegree, &
        M   => LM % MaxOrder, &
        nAM => LM % nAngularMomentCells )
    nAM = 0
    do iM = 0, M
      do iL = iM, L
        nAM = nAM + 1
      end do
    end do
    end associate !-- L, etc.

    call Show ( LM % MaxDegree, 'MaxDegree (l)', LM % IGNORABILITY )
    call Show ( LM % MaxOrder, 'MaxOrder (m)', LM % IGNORABILITY )
    call Show ( LM % nAngularMomentCells, 'nAngularMomentCells', &
                LM % IGNORABILITY )
    call Show ( LM % UseDevice, 'UseDevice', LM % IGNORABILITY )

    nullify ( G )

  end subroutine SetParameters


  subroutine ComputeMomentsLocal ( LM, Source )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM
    type ( StorageForm ), intent ( in ) :: &
      Source !-- array over levels    

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iR     !-- iRadius
    real ( KDR ) :: &
      R
    logical ( KDL ) :: &
      GridError
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( C => LM % Chart )
    class is ( Chart_SL_Template )

      G => C % Geometry ( )

      call ComputeMomentsLocal_CSL_Kernel &
             ( LM % MyM_RC, LM % MyM_RS, LM % MyM_IC, LM % MyM_IS, &
               C % CoordinateSystem, C % IsProperCell, &
               Source % Value, LM % Origin, LM % RadialEdge, &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               G % Value ( :, G % CENTER_U ( 2 ) ), &
               G % Value ( :, G % CENTER_U ( 3 ) ), &
               G % Value ( :, G % VOLUME ), Source % iaSelected, &
               LM % MaxDegree, LM % MaxOrder, C % nDimensions, G % nValues, &
               LM % nEquations, LM % nAngularMomentCells, GridError, &
               LM % SolidHarmonic_RC, LM % SolidHarmonic_IC, &
               LM % SolidHarmonic_RS, LM % SolidHarmonic_IS, &
               UseDeviceOption = LM % UseDevice )
      if ( GridError ) then
        call Show ( 'Radial grid not large enough', CONSOLE % ERROR )
      ! call Show ( R, 'R', CONSOLE % ERROR )
      ! call Show ( RadialEdge ( size ( RadialEdge ) ), 'R_Max', &
      !             CONSOLE % ERROR )
        call Show ( 'ComputeMomentsLocal', 'subroutine', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'ComputeMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
    end select !-- C
    
    nullify ( G )

  end subroutine ComputeMomentsLocal


  subroutine SetRadialEdgeSpherical ( LM )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      LM

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ), dimension ( : ), allocatable :: &
      Center, &
      Width

    select type ( C => LM % Chart )
    class is ( Chart_SL_Template )

      LM % Origin        =  0.0_KDR
      LM % nRadialCells  =  C % nCells ( 1 )

      allocate ( LM % RadialEdges )
      associate &
        ( nRC  =>  LM % nRadialCells, &
           RE  =>  LM % RadialEdges )

      call RE % Initialize ( [ nRC + 1, 1 ] )

      LM % RadialEdge  =>  RE % Value ( :, 1 )

      do iC = 1, nRC + 1
        LM % RadialEdge ( iC )  =  C % Edge ( 1 ) % Value ( iC )
      end do !-- iC

      end associate !-- nRC, etc.

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'SetRadialEdgeSpherical', 'subroutine', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC_Form', 'module', CONSOLE % ERROR )
    end select !-- C

  end subroutine SetRadialEdgeSpherical


end module LaplacianMultipole_ASC__Form
