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
        AngularFunction => null ( )
      type ( StorageForm ), allocatable :: &
        AngularFunctions
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
      SetAngularFunctions
  end type LaplacianMultipole_ASC_Form

    private :: &
      AssignAngularFunctionPointers, &
      ComputeAngularFunctions

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

    ! nullify ( L % Chart )

    ! if ( allocated ( L % SolidHarmonics ) ) &
    !   deallocate ( L % SolidHarmonics )
    ! if ( allocated ( L % RectangularCoordinates ) ) &
    !   deallocate ( L % RectangularCoordinates )

    ! nullify ( L % Source )
    ! nullify ( L % SolidHarmonic_IS )
    ! nullify ( L % SolidHarmonic_RS )
    ! nullify ( L % SolidHarmonic_IC )
    ! nullify ( L % SolidHarmonic_RC )

    ! nullify ( L % Volume )
    ! nullify ( L % RadiusSquared )
    ! nullify ( L % Radius )
    ! nullify ( L % Rectangular_Z )
    ! nullify ( L % Rectangular_Y )
    ! nullify ( L % Rectangular_X )
    nullify ( L % AngularFunction )

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
      
      select type ( C => A % Chart )
      class is ( Chart_SLD_Form )
        L % ReductionUseDevice = C % ExchangeGhostUseDevice
        call PROGRAM_HEADER % GetParameter &
               ( L % ReductionUseDevice, 'LaplacianReductionUseDevice', &
                 IgnorabilityOption = CONSOLE % INFO_2 )
      end select 

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


  subroutine SetAngularFunctions ( L )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      nAngularCells
real ( KDR ), dimension ( 100, 3 ), target :: &
  Test_2D
real ( KDR ), dimension ( :, :, : ), pointer :: &
  Test_3D
type ( StorageForm ) :: &
  TestStorage

    allocate ( L % AngularFunctions )
    associate ( AF => L % AngularFunctions )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )

        nAngularCells  =  product ( C % nCellsBrick ( 2 : 3 ) )
        call AF % Initialize &
               ( [ nAngularCells, L % nAngularMoments ], &
                 NameOption = 'AngularFunctions', &
                 VariableOption = L % AngularFunctionName )
        if ( L % UseDevice ) then
          call AF % AllocateDevice ( )
        end if

        call AssignAngularFunctionPointers &
               ( AF % Value, C % nCellsBrick, &
                 L % nAngularMoments, L % AngularFunction )

        call ComputeAngularFunctions &
               (     LM  =  L, &
                  Theta  =  C % Center ( 2 ) % Value, &
                    Phi  =  C % Center ( 3 ) % Value, &
                      L  =  L % MaxDegree, &
                      M  =  L % MaxOrder, &
                 nTheta  =  C % nCellsBrick ( 2 ), &
                   nPhi  =  C % nCellsBrick ( 3 ), &
                 oTheta  =  C % nGhostLayers ( 2 ), &
                   oPhi  =  C % nGhostLayers ( 3 ), &
                     AF  =  L % AngularFunction )

call Show ( size ( L % AngularFunction ), '>>> size AngularFunction' )
call Show ( size ( AF % Value ), '>>> size AF % Value' )
call Show ( L % AngularFunction ( :, 1, 1 ), '>>> AF (pointer)' )
call Show ( AF % Value ( 1 : 100, 1 ), '>>> AF (storage)' )

! Test_3D ( 1 : 10, 1 : 10, 1 : 3 )  =>  Test_2D
! Test_2D  =  0.0_KDR
! Test_3D  =  1.0_KDR
! call Show ( Test_3D ( :, :, 1 ), '>>> Test_3D' )
! call Show ( Test_2D ( :, 1 ), '>>> Test_2D' )

Test_3D  =>  null ( )
Test_2D  =   0.0_KDR 
call TestSetPointer ( Test_2D, Test_3D )
call TestSetValue ( 2.0_KDR, Test_3D )
call Show ( size ( Test_3D ), '>>> size Test_3D' )
call Show ( size ( Test_2D ), '>>> size Test_2D' )
call Show ( Test_3D ( :, :, 1 ), '>>> Test_3D 2' )
call Show ( Test_2D ( :, 1 ), '>>> Test_2D 2' )

Test_3D  =>  null ( )
call TestStorage % Initialize ( [ 100, 3 ] )
call TestSetPointer ( TestStorage % Value, Test_3D )
call TestSetValue ( 3.0_KDR, Test_3D )
call Show ( size ( Test_3D ), '>>> size Test_3D' )
call Show ( size ( Test_2D ), '>>> size Test_2D' )
call Show ( Test_3D ( :, :, 1 ), '>>> Test_3D 3' )
call Show ( TestStorage % Value ( :, 1 ), '>>> Test_2D 3' )

      case default
        call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
        call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetAngularFunctions', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetAngularFunctions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- AF

  end subroutine SetAngularFunctions


  subroutine AssignAngularFunctionPointers ( AF_2D, nCB, nAM, AF_3D )

    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
      AF_2D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCB  !-- nCellsBrick
    integer ( KDI ), intent ( in ) :: &
      nAM  !-- nAngularMoments
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      AF_3D

    AF_3D ( 1 : nCB ( 2 ), 1 : nCB ( 3 ), 1 : nAM )  =>  AF_2D

  end subroutine AssignAngularFunctionPointers


  subroutine ComputeAngularFunctions &
               ( LM, Theta, Phi, nTheta, L, M, nPhi, oTheta, oPhi, AF )

    class ( LaplacianMultipole_ASC_Form ), intent ( in ) :: &
      LM
    real ( KDR ), dimension ( -oTheta + 1 : ), intent ( in ) :: &
      Theta  !-- PolarAngle
    real ( KDR ), dimension ( -oPhi + 1 : ), intent ( in ) :: &
      Phi    !-- AzimuthalAngle
    integer ( KDI ), intent ( in ) :: &
      L, &   !-- MaxDegree
      M, &   !-- MaxOrder
      nTheta, nPhi, &
      oTheta, oPhi
    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      AF

    integer ( KDI ) :: &
      iL, &
      iM, &
      iA, &
      iTheta, &
      iPhi
    real ( KDR ) :: &
      P

    iA  =  1
    do iM  =  0, M
      do iL  =  iM, L
        do iTheta  =  1, nTheta
          P  =  LM % AssociatedLegendre ( cos ( Theta ( iTheta ) ), iL, iM )
          if ( nPhi > 1 ) then
            do iPhi  =  1, nPhi
              AF ( iTheta, iPhi, iA     )  =  P * cos ( iM * Phi ( iPhi ) )
              AF ( iTheta, iPhi, iA + 1 )  =  P * sin ( iM * Phi ( iPhi ) )
            end do !-- iPhi
          else
            AF ( iTheta, 1, iA     )  =  P
            AF ( iTheta, 1, iA + 1 )  =  0.0_KDR
          end if
        end do !-- iTheta
        iA  =  iA + 2 !-- Cos, Sin
      end do !-- iL
    end do !-- iM

  end subroutine ComputeAngularFunctions


subroutine TestSetPointer ( T_2D, T_3D )

  real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
    T_2D
  real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
    T_3D

  T_3D ( 1 : 10, 1 : 10, 1 : 3 )  =>  T_2D

end subroutine TestSetPointer


subroutine TestSetValue ( A, T_3D )

  real ( KDR ), intent ( in ) :: &
    A
  real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
    T_3D

  T_3D  =  A

end subroutine TestSetValue


end module LaplacianMultipole_ASC__Form
