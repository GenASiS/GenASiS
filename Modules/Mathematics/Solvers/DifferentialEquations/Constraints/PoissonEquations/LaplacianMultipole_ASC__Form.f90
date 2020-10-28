module LaplacianMultipole_ASC__Form

  !-- LaplacianMultipole_AtlasSingleChart__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Template

  implicit none
  private

  type, public, extends ( LaplacianMultipoleTemplate ) :: &
    LaplacianMultipole_ASC_Form
      real ( KDR ), dimension ( :, : ), pointer :: &
        SolidAngle => null ( )
      real ( KDR ), dimension ( :, :, : ), pointer :: &
        AngularFunction => null ( )
      type ( StorageForm ), allocatable :: &
        SolidAngles, &
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
    procedure, private, pass :: &
      ComputeMomentsLocal
  end type LaplacianMultipole_ASC_Form


    private :: &
      ComputeMomentsLocal_CSL_S_Kernel


    interface

      module subroutine ComputeMomentsLocal_CSL_S_Kernel &
                          ( UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeMomentsLocal_CSL_S_Kernel

    end interface


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

    nullify ( L % Chart )

    if ( allocated ( L % AngularFunctions ) ) &
      deallocate ( L % AngularFunctions )
    if ( allocated ( L % SolidAngles ) ) &
      deallocate ( L % SolidAngles )

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
    nullify ( L % SolidAngle )

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

    allocate ( L % SolidAngles )
    allocate ( L % AngularFunctions )
    associate &
      ( SA  =>  L % SolidAngles, &
        AF  =>  L % AngularFunctions )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )

        nAngularCells  =  product ( C % nCellsBrick ( 2 : 3 ) )
        call SA % Initialize &
               ( [ nAngularCells, 1 ] )
        call AF % Initialize &
               ( [ nAngularCells, L % nAngularMoments ], &
                 NameOption = 'AngularFunctions', &
                 VariableOption = L % AngularFunctionName )

        call AssignAngularFunctionPointers &
               ( AF % Value, SA % Value ( :, 1 ), C % nCellsBrick, &
                 L % nAngularMoments, L % AngularFunction, L % SolidAngle )

        associate &
          ( iaB  =>  C % iaBrick, &
            nCB  =>  C % nCellsBrick, &
            nGL  =>  C % nGhostLayers )
        call ComputeAngularFunctions &
               (        LM  =  L, &
                     Theta  =  C % Center ( 2 ) % Value, &
                  dTheta_L  =  C % WidthLeft ( 2 ) % Value, &
                  dTheta_R  =  C % WidthRight ( 2 ) % Value, &
                       Phi  =  C % Center ( 3 ) % Value, &
                    dPhi_L  =  C % WidthLeft ( 3 ) % Value, &
                    dPhi_R  =  C % WidthRight ( 3 ) % Value, &
                         L  =  L % MaxDegree, &
                         M  =  L % MaxOrder, &
                    nTheta  =  nCB ( 2 ), &
                      nPhi  =  nCB ( 3 ), &
                    oTheta  =  nGL ( 2 )  +  ( iaB ( 2 ) - 1 )  *  nCB ( 2 ), &
                      oPhi  =  nGL ( 3 )  +  ( iaB ( 3 ) - 1 )  *  nCB ( 3 ), &
                        AF  =  L % AngularFunction, &
                        SA  =  L % SolidAngle )
        end associate !-- iaB, etc.

      case default
        call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
        call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetAngularFunctions', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select

      if ( L % UseDevice ) then

        call SA % AllocateDevice ( )
        call SA % UpdateDevice ( )

        call AF % AllocateDevice ( )
        call AF % UpdateDevice ( )

      end if

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetAngularFunctions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- SA, etc.

  end subroutine SetAngularFunctions


  subroutine ComputeMomentsLocal ( L, Source )

      class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
        L
      class ( FieldAtlasTemplate ), intent ( in ) :: &
        Source

  end subroutine ComputeMomentsLocal


  subroutine AssignAngularFunctionPointers &
               ( AF_2D, SA_1D, nCB, nAM, AF_3D, SA_2D )

    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
      AF_2D
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      SA_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCB  !-- nCellsBrick
    integer ( KDI ), intent ( in ) :: &
      nAM  !-- nAngularMoments
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      AF_3D
    real ( KDR ), dimension ( :, : ), intent ( out ), pointer :: &
      SA_2D

    AF_3D ( 1 : nCB ( 2 ), 1 : nCB ( 3 ), 1 : nAM )  =>  AF_2D

    SA_2D ( 1 : nCB ( 2 ), 1 : nCB ( 3 ) )  =>  SA_1D

  end subroutine AssignAngularFunctionPointers


  subroutine ComputeAngularFunctions &
               ( LM, Theta, dTheta_L, dTheta_R, Phi, dPhi_L, dPhi_R, &
                 nTheta, L, M, nPhi, oTheta, oPhi, AF, SA )

    class ( LaplacianMultipole_ASC_Form ), intent ( in ) :: &
      LM
    real ( KDR ), dimension ( -oTheta + 1 : ), intent ( in ) :: &
      Theta, &  !-- PolarAngle
      dTheta_L, dTheta_R
    real ( KDR ), dimension ( -oPhi + 1 : ), intent ( in ) :: &
      Phi, &    !-- AzimuthalAngle
      dPhi_L, dPhi_R
    integer ( KDI ), intent ( in ) :: &
      L, &   !-- MaxDegree
      M, &   !-- MaxOrder
      nTheta, nPhi, &
      oTheta, oPhi
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      AF
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      SA

    integer ( KDI ) :: &
      iL, &
      iM, &
      iA, &
      iTheta, &
      iPhi
    real ( KDR ) :: &
      Pi, &
      Th_I, Th_O, &
      dPh, &
      P  !-- Normalized AssociatedLegendre polynomial

    Pi  =  CONSTANT % PI

    do iTheta  =  1, nTheta
      Th_I  =  Theta ( iTheta )  -  dTheta_L ( iTheta )
      Th_O  =  Theta ( iTheta )  +  dTheta_R ( iTheta )
      if ( nPhi > 1 ) then
        do iPhi  =  1, nPhi
          dPh  =  dPhi_L ( iPhi )  +  dPhi_R ( iPhi )
          SA ( iTheta, iPhi )  =  dPh * ( cos ( Th_I ) - cos ( Th_O ) )
        end do !-- iPhi
      else  !-- axisymmetry
        SA ( iTheta, 1 )  =  2.0_KDR * Pi *  ( cos ( Th_I ) - cos ( Th_O ) )
      end if  !-- nPhi > 1
    end do !-- iTheta

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
          else  !-- axisymmetry
            AF ( iTheta, 1, iA     )  =  P
            AF ( iTheta, 1, iA + 1 )  =  0.0_KDR
          end if  !-- nPhi > 1
        end do !-- iTheta
        iA  =  iA + 2 !-- Cos, Sin
      end do !-- iL
    end do !-- iM

  end subroutine ComputeAngularFunctions


end module LaplacianMultipole_ASC__Form
