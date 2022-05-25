module Laplacian_M_ASCG__Form

  !-- Laplacian_Multipole_AtlasSingleChartGrid__Form

  use Basics
  use Manifolds
  use Fields
  use Laplacian_M_H__Form

  implicit none
  private

  type, public, extends ( Laplacian_M_H_Form ) :: Laplacian_M_ASCG_Form
    real ( KDR ), dimension ( :, : ), pointer :: &
      dSolidAngle_2D => null ( )
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      AngularFunction_3D => null ( )
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      Source_4D => null ( )
    type ( StorageForm ), allocatable :: &
      dSolidAngles, &
      AngularFunctions
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetParameters_A
    procedure, private, pass :: &
      SetKernelFunctions
    procedure, private, pass :: &
      ComputeAngularMomentsLocal
  end type Laplacian_M_ASCG_Form

    private :: &
      ComputeAngularMomentsLocal_CGS_S_Kernel

    interface

      module subroutine ComputeAngularMomentsLocal_CGS_S_Kernel &
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
      end subroutine ComputeAngularMomentsLocal_CGS_S_Kernel

    end interface


    private :: &
      AssignAngularFunctionPointers, &
      ComputeAngularFunctions, &
      ComputeRadialFunctions, &
      AssignSourcePointer


contains


  subroutine Initialize ( L, G, MaxDegree, nEquations )

    class ( Laplacian_M_ASCG_Form ), intent ( inout ) :: &
      L
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      MaxDegree, &
      nEquations

    if ( L % Type  ==  '' ) &
      L % Type  =  'a Laplacian_M_ASCG' 

    call L % Initialize_H ( G, MaxDegree, nEquations )

  end subroutine Initialize


  impure elemental subroutine Finalize ( L )

    type ( Laplacian_M_ASCG_Form ), intent ( inout ) :: &
      L

    nullify ( L % Geometry )

    if ( allocated ( L % Parameters_P ) ) &
      deallocate ( L % Parameters_P )
    if ( allocated ( L % Integral_P ) ) &
      deallocate ( L % Integral_P )
    if ( allocated ( L % AngularFunctions ) ) &
      deallocate ( L % AngularFunctions )
    if ( allocated ( L % dSolidAngles ) ) &
      deallocate ( L % dSolidAngles )

    nullify ( L % Source_4D )
    nullify ( L % AngularFunction_3D )
    nullify ( L % dSolidAngle_2D )

  end subroutine Finalize


  subroutine SetParameters_A ( L, G )

    class ( Laplacian_M_ASCG_Form ), intent ( inout ) :: &
      L
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    L % Geometry  =>  G

    select type ( A  =>  L % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    if ( C % nDimensions  <  3 ) &
      L % MaxOrder  =  0
    if ( C % nDimensions  <  2 ) &
      L % MaxDegree  =  0

    L % DeviceMemory        =  G % DeviceMemory
    L % PinnedMemory        =  G % PinnedMemory
    L % DevicesCommunicate  =  G % DevicesCommunicate

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )

      L % nRadialCells  =  C % nCells ( 1 )

    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetParameters_A', 'subroutine', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetParameters_A', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine SetParameters_A


  subroutine SetKernelFunctions ( L )

    class ( Laplacian_M_ASCG_Form ), intent ( inout ) :: &
      L

    integer ( KDI ) :: &
      nAngularCells

    allocate ( L % dSolidAngles )
    allocate ( L % AngularFunctions )
    allocate ( L % d_Radius_3_3 )
    allocate ( L % RadialFunctions_R )
    allocate ( L % RadialFunctions_I )
    associate &
      (  dSA  =>  L % dSolidAngles, &
          AF  =>  L % AngularFunctions, &
        dR33  =>  L % d_Radius_3_3, &
        RF_R  =>  L % RadialFunctions_R, &
        RF_I  =>  L % RadialFunctions_I )

    select type ( A  =>  L % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
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
    call RF_R % Initialize ( [ L % nRadialCells, L % nAngularMoments ] )
    call RF_I % Initialize ( [ L % nRadialCells, L % nAngularMoments ] )

    call AssignAngularFunctionPointers &
           ( AF % Value, dSA % Value ( :, 1 ), C % nCellsBrick, &
             L % nAngularMoments, L % AngularFunction_3D, L % dSolidAngle_2D )

    associate &
      ( iaB  =>  C % iaBrick, &
        nCB  =>  C % nCellsBrick, &
        nGL  =>  C % nGhostLayers )

    call ComputeAngularFunctions &
           (      LM  =  L, &
             Theta_E  =  C % Edge ( 2 ) % Value, &
               Phi_E  =  C % Edge ( 3 ) % Value, &
                   L  =  L % MaxDegree, &
                   M  =  L % MaxOrder, &
              nTheta  =  nCB ( 2 ), &
                nPhi  =  nCB ( 3 ), &
              oTheta  =  nGL ( 2 )  +  ( iaB ( 2 ) - 1 )  *  nCB ( 2 ), &
                oPhi  =  nGL ( 3 )  +  ( iaB ( 3 ) - 1 )  *  nCB ( 3 ), &
                  AF  =  L % AngularFunction_3D, &
                 dSA  =  L % dSolidAngle_2D )

    call ComputeRadialFunctions &
           (  R_E  =  C % Edge ( 1 ) % Value, &
                L  =  L % MaxDegree, &
                M  =  L % MaxOrder, &
               nR  =  L % nRadialCells, &
               oR  =  nGL ( 1 ), &
             dR33  =  L % d_Radius_3_3 % Value, &
             RF_R  =  L % RadialFunctions_R % Value, &
             RF_I  =  L % RadialFunctions_I % Value )

    end associate !-- iaB, etc.

    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetKernelFunctions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    if ( L % DeviceMemory ) then

      call  dSA % AllocateDevice ( AssociateVariablesOption = .false. )
      call   AF % AllocateDevice ( AssociateVariablesOption = .false. )
      call dR33 % AllocateDevice ( AssociateVariablesOption = .false. )
      call RF_R % AllocateDevice ( AssociateVariablesOption = .false. )
      call RF_I % AllocateDevice ( AssociateVariablesOption = .false. )

      call  dSA % UpdateDevice ( )
      call   AF % UpdateDevice ( )
      call dR33 % UpdateDevice ( )
      call RF_R % UpdateDevice ( )
      call RF_I % UpdateDevice ( )

    end if

    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetKernelFunctions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- dSA, etc.

  end subroutine SetKernelFunctions


  subroutine ComputeAngularMomentsLocal ( L, Source )

    class ( Laplacian_M_ASCG_Form ), intent ( inout ) :: &
      L
    class ( FieldSetForm ), intent ( inout ) :: &
      Source

    select type ( A  =>  L % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )
    associate &
      ( Source_S  =>  Source % Storage ( 1 ) )
    associate &
      (  nV => Source_S % nVariables, &
        iaS => Source_S % iaSelected )
 
    if ( nV /= L % nEquations ) then
      call Show ( 'Wrong number of variables in Solution', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iaS ( nV ) - iaS ( 1 ) + 1  /=  nV ) then
      call Show ( 'Solution variables must be contiguous', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    call AssignSourcePointer &
           ( Source_S % Value ( :, iaS ( 1 ) : iaS ( nV ) ), &
             C % nCellsBrick, C % nGhostLayers, L % nEquations, L % Source_4D )
    
    call Source_S % ReassociateHost ( AssociateVariablesOption = .false. )
    
    end associate !-- nV, etc.

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'SPHERICAL' )
      call ComputeAngularMomentsLocal_CGS_S_Kernel &
             ( L % MyAngularMoment_3D, L % Source_4D, L % AngularFunction_3D, &
               L % dSolidAngle_2D, C % nCellsBrick, C % nGhostLayers, &
               L % nEquations, L % nAngularMoments, &
               oR = ( C % iaBrick ( 1 ) - 1 ) * C % nCellsBrick ( 1 ), &
               UseDeviceOption = L % DeviceMemory )
    case default
      call Show ( 'Coordinate system not supported', CONSOLE % ERROR )
      call Show ( C % CoordinateSystem, 'CoordinateSystem', CONSOLE % ERROR )
      call Show ( 'LaplacianMultipole_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem
    
    call Source_S % ReassociateHost ( AssociateVariablesOption = .true. )
    
    end associate !-- Source_S
    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Laplacian_M_ASCG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeAngularMomentsLocal', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

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
               ( LM, Theta_E, Phi_E, L, M, nTheta, nPhi, oTheta, oPhi, AF, dSA )

    class ( Laplacian_M_ASCG_Form ), intent ( inout ) :: &
      LM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Theta_E  !-- PolarAngle
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Phi_E    !-- AzimuthalAngle
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
      Ph_I, Ph_O, &
      P  !-- Normalized AssociatedLegendre polynomial

    Pi  =  CONSTANT % PI

    associate &
      ( IP  =>  LM % Integral_P, &
        PP  =>  LM % Parameters_P )

    if ( nTheta > 1 ) then
      do iTheta  =  1, nTheta
        Th_I  =  Theta_E ( oTheta + iTheta )
        Th_O  =  Theta_E ( oTheta + iTheta + 1 )
        if ( nPhi > 1 ) then
          do iPhi  =  1, nPhi
            Ph_I  =  Phi_E ( oPhi + iPhi )
            Ph_O  =  Phi_E ( oPhi + iPhi + 1 )
            dSA ( iTheta, iPhi )  &
              =  ( cos ( Th_I ) - cos ( Th_O ) )  *  ( Ph_O - Ph_I )
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

        PP % iDegree  =  iL
        PP % iOrder   =  iM

        if ( nTheta > 1 ) then
          do iTheta  =  1, nTheta

!            P  =  LM % AssociatedLegendre &
!                         ( cos ( Theta_C ( oTheta + iTheta ) ), iL, iM )

            Th_I  =  Theta_E ( oTheta + iTheta )
            Th_O  =  Theta_E ( oTheta + iTheta + 1 )
            call IP % Compute ( cos ( Th_I ), cos ( Th_O ), P )
            P  =  P  /  ( cos ( Th_O )  -  cos ( Th_I ) )

            if ( nPhi > 1 ) then
              do iPhi  =  1, nPhi
                Ph_I  =  Phi_E ( oPhi + iPhi )
                Ph_O  =  Phi_E ( oPhi + iPhi + 1 )
                if ( iM  >  0 ) then
                  AF ( iTheta, iPhi, iA     )  &
                    =  P  *  ( sin ( iM * Ph_O )  -  sin ( iM * Ph_I ) ) &
                             / ( iM * ( Ph_O - Ph_I ) )
                  AF ( iTheta, iPhi, iA + 1 )  &
                    =  P  *  ( cos ( iM * Ph_I )  -  cos ( iM * Ph_O ) ) &
                             / ( iM * ( Ph_O - Ph_I ) )
                else
                  AF ( iTheta, iPhi, iA     )  =  P
                  AF ( iTheta, iPhi, iA + 1 )  =  0.0_KDR
                end if
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

    end associate !-- IP, etc.

  end subroutine ComputeAngularFunctions


  subroutine ComputeRadialFunctions &
               ( R_E, L, M, nR, oR, dR33, RF_R, RF_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_E
    integer ( KDI ), intent ( in ) :: &
      L, &   !-- MaxDegree
      M, &   !-- MaxOrder
      nR, &  !-- nRadial
      oR     !-- oRadial
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      dR33, &
      RF_R, RF_I  !-- RadialFunction_Regular, _Irregular

    integer ( KDI ) :: &
      iL, &
      iM, &
      iA, &
      iR  !-- iRadial
    real ( KDR ) :: &
      R_I, R_O, &
      R_C

    do iR  =  1, nR
      R_I  =  R_E ( oR + iR )
      R_O  =  R_E ( oR + iR + 1 )
      dR33 ( iR, 1 )  =  ( R_O ** 3  -  R_I ** 3 )  /  3.0_KDR
    end do !-- nRC

    ! iA  =  1
    ! do iM  =  0, M
    !   do iL  =  iM, L
    !     do iR  =  1, nR
    !       ! RF_R ( iR, iA     )  =  R_C ( oR + iR ) ** iL
    !       ! RF_R ( iR, iA + 1 )  =  R_C ( oR + iR ) ** iL
    !       ! RF_I ( iR, iA     )  =  R_C ( oR + iR ) ** ( - ( iL + 1 ) )
    !       ! RF_I ( iR, iA + 1 )  =  R_C ( oR + iR ) ** ( - ( iL + 1 ) ) 
    !       R_I  =  R_E ( oR + iR )
    !       R_O  =  R_E ( oR + iR + 1 )
    !       R_C  =  0.5_KDR * ( R_I + R_O )
    !       RF_R ( iR, iA     )  =  R_C ** iL
    !       RF_R ( iR, iA + 1 )  =  R_C ** iL
    !       RF_I ( iR, iA     )  =  R_C ** ( - ( iL + 1 ) )
    !       RF_I ( iR, iA + 1 )  =  R_C ** ( - ( iL + 1 ) ) 
    !     end do !-- iR
    !     iA  =  iA + 2 !-- Cos, Sin
    !   end do !-- iL
    ! end do !-- iM

    iA  =  1
    do iM  =  0, M
      do iL  =  iM, L
        do iR  =  1, nR

          if ( iR  ==  1  .and.  iL  >  1 ) then !-- avoid singularity below

            RF_R ( iR, iA )  =  0.0_KDR
            RF_I ( iR, iA )  =  0.0_KDR

          else

            R_I  =  R_E ( oR + iR )
            R_O  =  R_E ( oR + iR + 1 )

            RF_R ( iR, iA )  =  ( R_O ** ( iL + 3 )  -  R_I ** ( iL + 3 ) )  &
                              /  ( ( iL + 3 )  *  dR33 ( iR, 1 ) )

            if ( iL  ==  2 ) then
              RF_I ( iR, iA )  =  log ( R_O / R_I )  /  dR33 ( iR, 1 )
            else
              RF_I ( iR, iA )  =  ( R_O ** ( 2 - iL )  -  R_I ** ( 2 - iL ) )  &
                                  /  ( ( 2 - iL )  *  dR33 ( iR, 1 ) )
            end if !-- iL == 2

          end if !

          RF_R ( iR, iA + 1 )  =  RF_R ( iR, iA )
          RF_I ( iR, iA + 1 )  =  RF_I ( iR, iA )

        end do !-- iR
        iA  =  iA + 2  !-- Cos, Sin
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


end module Laplacian_M_ASCG__Form
