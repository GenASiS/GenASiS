module Coarsening_C__Form

  !-- Coarsening_Central_Form

  use Basics
  use Manifolds
  use FieldSets
  use Geometry_F__Form

  implicit none
  private

  type, public, extends ( FieldSetForm ) :: Coarsening_C_Form
    integer ( KDI ) :: &
      COARSENING_POLAR = 0, &
      COARSENING_AZIMUTHAL = 0
    integer ( KDI ) :: &
      nBlocksCoarsen
    integer ( KDI ), dimension ( : ), allocatable :: &
      iRadius
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      iTheta, &
      iPhi
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
  contains
    procedure, private, pass :: &
      Initialize_C
    generic, public :: &
      Initialize => Initialize_C
    procedure, public, pass ( C ) :: &
      Compute
    final :: &
      Finalize
  end type Coarsening_C_Form

    private :: &
      SetCoarseningPolar, &
      SetBlocks, &
      ComputeKernel
         
    interface
      
      module subroutine ComputeKernel &
               ( FS_4D, dV_3D, iTh, iPh, iR, iaS, oC, nBC )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          FS_4D
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          dV_3D
        integer ( KDI ), dimension ( :, : ), intent ( in ) :: &
          iTh, &
          iPh
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iR, &
          iaS, &
          oC
        integer ( KDI ), intent ( in ) :: &
          nBC
      end subroutine ComputeKernel

    end interface


contains


  subroutine Initialize_C ( C, G )

    class ( Coarsening_C_Form ), intent ( inout ) :: &
      C
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    if ( C % Type == '' ) &
      C % Type  =  'a Coarsening_C' 
    
    C % COARSENING_POLAR      =  1
    C % COARSENING_AZIMUTHAL  =  2

    call C % FieldSetForm % Initialize &
           ( G % Atlas, &
             FieldOption = [ 'CoarseningPolar    ', 'CoarseningAzimuthal' ], &
             NameOption = 'Coarsening', &
             DeviceMemoryOption = G % DeviceMemory, &
             PinnedMemoryOption = G % PinnedMemory, &
             DevicesCommunicateOption = G % DevicesCommunicate, &
             nFieldsOption = 2 )

    C % Geometry  =>  G

    select type ( A  =>  G % Atlas )
    class is ( Atlas_SCG_CC_Form )

      call SetCoarseningPolar &
             ( C  = A % Chart_GS_CC, &
               R  = G % Storage_GS % Value ( :, G % CENTER_U_1 ), &
               CP = C % Storage_GS % Value ( :, C % COARSENING_POLAR ) )

      call SetBlocks &
             (  C   = A % Chart_GS_CC, &
                CP  = C % Storage_GS % Value ( :, C % COARSENING_POLAR ), &
               iTh  = C % iTheta, &
               iPh  = C % iPhi, &
               iRad = C % iRadius, &
               nBC  = C % nBlocksCoarsen )

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_C', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Initialize_C


  subroutine SetCoarseningPolar ( C, R, CP )
    
    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      CP

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      dTheta

    dTheta  =  CONSTANT % PI  /  C % nCellsPolar

    do iV = 1, size ( CP )
      if ( .not. C % ProperCell ( iV ) ) &
        cycle
      CP ( iV )  =  1.0_KDR
      Coarsen_2: do
        if ( CP ( iV )  *  R ( iV )  *  dTheta  >  C % MinWidth ) &
          exit Coarsen_2
        CP ( iV )  =  2.0_KDR  *  CP ( iV )
      end do Coarsen_2
    end do !-- iV

  end subroutine SetCoarseningPolar


  subroutine SetBlocks ( C, CP, iTh, iPh, iRad, nBC )

    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      CP  !-- CoarsenPolar
    integer ( KDI ), dimension ( :, : ), intent ( out ), allocatable :: &
      iTh, &
      iPh
    integer ( KDI ), dimension ( : ), intent ( out ), allocatable :: &
      iRad
    integer ( KDI ), intent ( out ) :: &
      nBC  !-- nBlocksCoarsen

    integer ( KDI ) :: &
      iR, &   !-- iRadius
      iBC, &  !-- iBlockCoarsen
      iBP, &  !-- iBlockPolar
      oTh
    integer ( KDI ), dimension ( : ), allocatable :: &
      nCP, &  !-- nCoarsenPolar
      nBP     !-- nBlocksPolar
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      CP_3D

    call C % SetFieldPointer ( CP, CP_3D )

    associate ( nR  =>  C % nCellsBrick ( 1 ) )
    allocate ( nCP ( nR ), nBP ( nR ) )

    select case ( C % nDimensions )
    case ( 2 )

      do iR  =  1, nR
        nCP ( iR )  =  CP_3D ( iR, 1, 1 )  +  0.49_KDR
        nBP ( iR )  =  C % nCellsPolar  /  nCP ( iR )
      end do !-- iR

      where ( nCP == 1 )
        nBP  =  0
      end where

      nBC  =  sum ( nBP )

      allocate ( iRad ( nBC ) )
      allocate ( iTh ( 2, nBC ) )
      allocate ( iPh ( 2, nBC ) )

      iBC  =  0
      do iR  =  1, nR
        do iBP  =  1, nBP ( iR )

          iBC  =  iBC + 1

          iRad ( iBC )  =  iR

          oTh  =  ( iBP - 1 )  *  nCP ( iR )
          iTh ( 1 : 2, iBC )  =  [ oTh  +  1, oTh  +  nCP ( iR ) ]
          iPh ( 1 : 2, iBC )  =  [ 1, 1 ]

        end do !-- iBP
      end do !-- iR

    end select !-- nDimensions
    end associate !-- nR

  end subroutine SetBlocks


  subroutine Compute ( FS, C )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( Coarsening_C_Form ), intent ( in ) :: &
      C

    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dV_3D
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      FS_4D

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_CC_Form )
    associate &
      ( G        =>  C % Geometry, &
        C_GS_CC  =>  A % Chart_GS_CC )

    if ( C_GS_CC % nDimensions  ==  1 ) &
      return

    call C_GS_CC % SetFieldPointer &
           ( FS % Storage_GS % Value, FS_4D )
    call C_GS_CC % SetFieldPointer &
           ( G % Storage_GS % Value ( :, G % VOLUME ), dV_3D )
       
    call ComputeKernel &
           ( FS_4D, dV_3D, &
             iTh = C % iTheta, &
             iPh = C % iPhi, &
             iR  = C % iRadius, &
             iaS = FS % iaSelected, &
             oC  = C_GS_CC % nGhostLayers, &
             nBC = C % nBlocksCoarsen )

    end associate !-- G, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Compute


  impure elemental subroutine Finalize ( C )

    type ( Coarsening_C_Form ), intent ( inout ) :: &
      C

    nullify ( C % Geometry )

    if ( allocated ( C % iPhi ) ) &
      deallocate ( C % iPhi )
    if ( allocated ( C % iTheta ) ) &
      deallocate ( C % iTheta )
    if ( allocated ( C % iRadius ) ) &
      deallocate ( C % iRadius )

  end subroutine Finalize

  
end module Coarsening_C__Form
