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
      SetKernelFunctions
    procedure, private, pass :: &
      ComputeAngularMomentsLocal
  end type LaplacianMultipole_ASC_Form


    private :: &
      ComputeAngularMomentsLocal_CSL_S_Kernel


    interface

      module subroutine ComputeAngularMomentsLocal_CSL_S_Kernel &
                          ( MyAM, S, AF, dSA, nC, oC, nE, nAM, oR, &
                            UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
          MyAM  !-- MyAngularMoment_3D
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
      end subroutine ComputeAngularMomentsLocal_CSL_S_Kernel

    end interface


    private :: &
      AssignAngularFunctionPointers, &
      ComputeAngularFunctions, &
      ComputeRadialFunctions, &
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


  subroutine SetKernelFunctions ( L )

    class ( LaplacianMultipole_ASC_Form ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      nAngularCells

    allocate ( L % dSolidAngles )
    allocate ( L % AngularFunctions )
    allocate ( L % d_Radius_3_3 )
    allocate ( L % CellFraction )
    allocate ( L % RadialFunctions_R )
    allocate ( L % RadialFunctions_I )
    associate &
      (  dSA  =>  L % dSolidAngles, &
          AF  =>  L % AngularFunctions, &
        dR33  =>  L % d_Radius_3_3, &
          CF  =>  L % CellFraction, &
        RF_R  =>  L % RadialFunctions_R, &
        RF_I  =>  L % RadialFunctions_I )

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

        call dR33 % Initialize ( [ L % nRadialCells, 1 ] )
        call   CF % Initialize ( [ L % nRadialCells, 1 ] )
        call RF_R % Initialize ( [ L % nRadialCells, L % nAngularMoments ] )
        call RF_I % Initialize ( [ L % nRadialCells, L % nAngularMoments ] )

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

        call ComputeRadialFunctions &
                (    R  =  C % Center ( 1 ) % Value, &
                  dR_L  =  C % WidthLeft ( 1 ) % Value, &
                  dR_R  =  C % WidthRight ( 1 ) % Value, &
                     L  =  L % MaxDegree, &
                     M  =  L % MaxOrder, &
                    nR  =  L % nRadialCells, &
                    oR  =  nGL ( 1 ), &
                  dR33  =  L % d_Radius_3_3 % Value, &
                    CF  =  L % CellFraction % Value, &
                  RF_R  =  L % RadialFunctions_R % Value, &
                  RF_I  =  L % RadialFunctions_I % Value )

        end associate !-- iaB, etc.

      case default
        call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
        call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
        call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetKernelFunctions', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select

      if ( L % UseDevice ) then

        call  dSA % AllocateDevice ( AssociateVariablesOption = .false. )
        call   AF % AllocateDevice ( AssociateVariablesOption = .false. )
        call dR33 % AllocateDevice ( AssociateVariablesOption = .false. )
        call   CF % AllocateDevice ( AssociateVariablesOption = .false. )
        call RF_R % AllocateDevice ( AssociateVariablesOption = .false. )
        call RF_I % AllocateDevice ( AssociateVariablesOption = .false. )

        call  dSA % UpdateDevice ( )
        call   AF % UpdateDevice ( )
        call dR33 % UpdateDevice ( )
        call   CF % UpdateDevice ( )
        call RF_R % UpdateDevice ( )
        call RF_I % UpdateDevice ( )

      end if

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetKernelFunctions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    end associate !-- dSA, etc.

  end subroutine SetKernelFunctions


  subroutine ComputeAngularMomentsLocal ( L, Source )

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
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iaS ( nV ) - iaS ( 1 ) + 1  /=  nV ) then
      call Show ( 'Solution variables must be contiguous', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call AssignSourcePointer &
           ( Source_S % Value ( :, iaS ( 1 ) : iaS ( nV ) ), &
             C % nCellsBrick, C % nGhostLayers, L % nEquations, L % Source )

    end associate !-- nV, etc.

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call ComputeAngularMomentsLocal_CSL_S_Kernel &
             ( L % MyAngularMoment_3D, L % Source, L % AngularFunction, &
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
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Source

    class default
      call Show ( 'Chart type not supported', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    nullify ( Source_S )

  end subroutine ComputeAngularMomentsLocal


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
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Theta, &  !-- PolarAngle
      dTheta_L, dTheta_R
    real ( KDR ), dimension ( : ), intent ( in ) :: &
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

    if ( nTheta > 1 ) then
      do iTheta  =  1, nTheta
        Th_I  =  Theta ( oTheta + iTheta )  -  dTheta_L ( oTheta + iTheta )
        Th_O  =  Theta ( oTheta + iTheta )  +  dTheta_R ( oTheta + iTheta )
        if ( nPhi > 1 ) then
          do iPhi  =  1, nPhi
            dPh  =  dPhi_L ( oPhi + iPhi )  +  dPhi_R ( oPhi + iPhi )
            dSA ( iTheta, iPhi )  =  dPh * ( cos ( Th_I ) - cos ( Th_O ) )
          end do !-- iPhi
        else !-- axisymmetry
          dSA ( iTheta, 1 )  =  2.0_KDR * Pi *  ( cos ( Th_I ) - cos ( Th_O ) )
        end if  !-- nPhi > 1
      end do !-- iTheta
    else !-- spherical symmetry
      dSA ( 1, 1 )  =  4.0_KDR * Pi
    end if

    iA  =  1
    do iM  =  0, M
      do iL  =  iM, L
        if ( nTheta > 1 ) then
          do iTheta  =  1, nTheta
            P  =  LM % AssociatedLegendre &
                         ( cos ( Theta ( oTheta + iTheta ) ), iL, iM )
            if ( nPhi > 1 ) then
              do iPhi  =  1, nPhi
                AF ( iTheta, iPhi, iA     )  &
                  =  P * cos ( iM * Phi ( oPhi + iPhi ) )
                AF ( iTheta, iPhi, iA + 1 )  &
                  =  P * sin ( iM * Phi ( oPhi + iPhi ) )
              end do !-- iPhi
            else !-- axisymmetry
              AF ( iTheta, 1, iA     )  =  P
              AF ( iTheta, 1, iA + 1 )  =  0.0_KDR
            end if  !-- nPhi > 1
          end do !-- iTheta
        else !-- spherical symmetry
          P  =  LM % AssociatedLegendre ( cos ( 0.5_KDR * Pi ), iL, iM )
          AF ( 1, 1, iA     )  =  P
          AF ( 1, 1, iA + 1 )  =  0.0_KDR
        end if
        iA  =  iA + 2 !-- Cos, Sin    
      end do !-- iL
    end do !-- iM

  end subroutine ComputeAngularFunctions


  subroutine ComputeRadialFunctions &
               ( R, dR_L, dR_R, L, M, nR, oR, dR33, CF, RF_R, RF_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      dR_L, dR_R
    integer ( KDI ), intent ( in ) :: &
      L, &   !-- MaxDegree
      M, &   !-- MaxOrder
      nR, &  !-- nRadial
      oR     !-- oRadial
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      dR33, &
      CF, &
      RF_R, RF_I  !-- RadialFunction_Regular, _Irregular

    integer ( KDI ) :: &
      iL, &
      iM, &
      iA, &
      iR  !-- iRadial
    real ( KDR ) :: &
      R_C, &
      R_I, R_O

    do iR  =  1, nR
      R_C  =  R ( oR + iR )
      R_I  =  R ( oR + iR )  -  dR_L ( oR + iR )
      R_O  =  R ( oR + iR )  +  dR_R ( oR + iR )
      dR33 ( iR, 1 )  =  ( R_O ** 3  -  R_I ** 3 )  /  3.0_KDR
        CF ( iR, 1 )  =  ( R_C - R_I ) / ( R_O - R_I )
    end do !-- nRC

    iA  =  1
    do iM  =  0, M
      do iL  =  iM, L
        do iR  =  1, nR
          RF_R ( iR, iA     )  =  R ( oR + iR ) ** iL
          RF_R ( iR, iA + 1 )  =  R ( oR + iR ) ** iL
          RF_I ( iR, iA     )  =  R ( oR + iR ) ** ( - ( iL + 1 ) )
          RF_I ( iR, iA + 1 )  =  R ( oR + iR ) ** ( - ( iL + 1 ) ) 
        end do !-- iR
        iA  =  iA + 2 !-- Cos, Sin
      end do !-- iL
    end do !-- iM

  end subroutine ComputeRadialFunctions


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
