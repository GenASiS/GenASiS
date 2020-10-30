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
        dSolidAngle => null ( )
      real ( KDR ), dimension ( :, :, : ), pointer :: &
        AngularFunction => null ( )
      real ( KDR ), dimension ( :, :, :, : ), pointer :: &
        Source => null ( )
      type ( StorageForm ), allocatable :: &
        dSolidAngles, &
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
                          ( MySM, S, AF, dSA, nC, oC, nE, nAM, oR, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MySM  !-- MyShellMoment_3D
        real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
          S  !-- Source
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          AF  !-- AngularFunction
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          dSA  !-- dSolidAngle
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, oC  !-- nCells, oCell
        integer ( KDI ), intent ( in ) :: &
          nE, &   !-- nEquations
          nAM, &  !-- nAngularMoments
          oR      !-- oRadius
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeMomentsLocal_CSL_S_Kernel

    end interface


    private :: &
      AssignAngularFunctionPointers, &
      ComputeAngularFunctions, &
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

    if ( allocated ( L % AngularFunctions ) ) &
      deallocate ( L % AngularFunctions )
    if ( allocated ( L % dSolidAngles ) ) &
      deallocate ( L % dSolidAngles )

    nullify ( L % Source )
    nullify ( L % AngularFunction )
    nullify ( L % dSolidAngle )

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

    allocate ( L % dSolidAngles )
    allocate ( L % AngularFunctions )
    associate &
      ( dSA  =>  L % dSolidAngles, &
         AF  =>  L % AngularFunctions )

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

      select case ( trim ( C % CoordinateSystem ) )
      case ( 'SPHERICAL' )

        nAngularCells  =  product ( C % nCellsBrick ( 2 : 3 ) )
        call dSA % Initialize &
               ( [ nAngularCells, 1 ] )
        call AF % Initialize &
               ( [ nAngularCells, L % nAngularMoments ], &
                 NameOption = 'AngularFunctions', &
                 VariableOption = L % AngularFunctionName )

        call AssignAngularFunctionPointers &
               ( AF % Value, dSA % Value ( :, 1 ), C % nCellsBrick, &
                 L % nAngularMoments, L % AngularFunction, L % dSolidAngle )

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
                       dSA  =  L % dSolidAngle )
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

        call dSA % AllocateDevice ( )
        call dSA % UpdateDevice ( )

        call AF % AllocateDevice ( )
        call AF % UpdateDevice ( )

      end if

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetAngularFunctions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- dSA, etc.

  end subroutine SetAngularFunctions


  subroutine ComputeMomentsLocal ( L, Source )

      class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
        L
      class ( FieldAtlasTemplate ), intent ( in ) :: &
        Source

    class ( StorageForm ), pointer :: &
      Source_S

    select type ( C => L % Chart )
    class is ( Chart_SL_Template )

    select type ( Source )
    class is ( Storage_ASC_Form )
    Source_S => Source % Storage ( )

    associate &
      (  nV => Source_S % nVariables, &
        iaS => Source_S % iaSelected )
 
    if ( nV /= L % nEquations ) then
      call Show ( 'Wrong number of variables in Solution', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iaS ( nV ) - iaS ( 1 ) + 1  /=  nV ) then
      call Show ( 'Solution variables must be contiguous', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call AssignSourcePointer &
           ( Source_S % Value ( :, iaS ( 1 ) : iaS ( nV ) ), &
             C % nCellsBrick, C % nGhostLayers, L % nEquations, L % Source )

    end associate !-- nV, etc.

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call ComputeMomentsLocal_CSL_S_Kernel &
             ( L % MyShellMoment_3D, L % Source, L % AngularFunction, &
               L % dSolidAngle, C % nCellsBrick, C % nGhostLayers, &
               L % nEquations, L % nAngularMoments, &
               oR = ( C % iaBrick ( 1 ) - 1 ) * C % nCellsBrick ( 1 ), &
               UseDeviceOption = L % UseDevice )
    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentAtlas', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    class default
      call Show ( 'Source type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Source

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    nullify ( Source_S )

  end subroutine ComputeMomentsLocal


  subroutine AssignAngularFunctionPointers &
               ( AF_2D, dSA_1D, nCB, nAM, AF_3D, dSA_2D )

    real ( KDR ), dimension ( :, : ), intent ( in ), target, contiguous :: &
      AF_2D
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      dSA_1D
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCB  !-- nCellsBrick
    integer ( KDI ), intent ( in ) :: &
      nAM  !-- nAngularMoments
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      AF_3D
    real ( KDR ), dimension ( :, : ), intent ( out ), pointer :: &
      dSA_2D

    AF_3D ( 1 : nCB ( 2 ), 1 : nCB ( 3 ), 1 : nAM )  =>  AF_2D

    dSA_2D ( 1 : nCB ( 2 ), 1 : nCB ( 3 ) )  =>  dSA_1D

  end subroutine AssignAngularFunctionPointers


  subroutine ComputeAngularFunctions &
               ( LM, Theta, dTheta_L, dTheta_R, Phi, dPhi_L, dPhi_R, &
                 nTheta, L, M, nPhi, oTheta, oPhi, AF, dSA )

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
      dSA

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
          dSA ( iTheta, iPhi )  =  dPh * ( cos ( Th_I ) - cos ( Th_O ) )
        end do !-- iPhi
      else  !-- axisymmetry
        dSA ( iTheta, 1 )  =  2.0_KDR * Pi *  ( cos ( Th_I ) - cos ( Th_O ) )
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
