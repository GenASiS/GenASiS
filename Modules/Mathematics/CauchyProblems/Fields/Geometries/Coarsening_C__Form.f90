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
      COARSENING_POLAR     = 0, &
      COARSENING_AZIMUTHAL = 0, &
      N_BLOCKS_POLAR       = 0, &
      N_BLOCKS_AZIMUTHAL   = 0, &
      BLOCK_LABEL          = 0  !-- Random label for visualization
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
      SetCoarseningAzimuthal, &
      SetBlocks

    private :: &
      ComputeKernel
         
    interface
      
      module subroutine ComputeKernel &
               ( FS_4D, dV_3D, iTh, iPh, iR, iaS, oC, nBC, UseDeviceOption )
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
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine Initialize_C ( C, G, NameOption )

    class ( Coarsening_C_Form ), intent ( inout ) :: &
      C
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    character ( * ), intent ( in ), optional :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( C % Type == '' ) &
      C % Type  =  'a Coarsening_C' 

    Name  =  'Coarsening'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    C % COARSENING_POLAR     =  1
    C % COARSENING_AZIMUTHAL =  2
    C % N_BLOCKS_POLAR       =  3
    C % N_BLOCKS_AZIMUTHAL   =  4
    C % BLOCK_LABEL          =  5

    call C % FieldSetForm % Initialize &
           ( G % Atlas, &
             FieldOption &
               = [ 'CoarseningPolar    ', &
                   'CoarseningAzimuthal', &
                   'nBlocksPolar       ', &
                   'nBlocksAzimuthal   ', &
                   'BlockLabel         ' ], &
             NameOption = Name, &
             DeviceMemoryOption = G % DeviceMemory, &
             PinnedMemoryOption = G % PinnedMemory, &
             DevicesCommunicateOption = G % DevicesCommunicate, &
             nFieldsOption = 5 )

    C % Geometry  =>  G

    select type ( A  =>  G % Atlas )
    class is ( Atlas_SCG_C_Form )

      call SetCoarseningPolar &
             ( C  = A % Chart_GS_C, &
               R  = G % Storage_GS % Value ( :, G % CENTER_U_1 ), &
               CP = C % Storage_GS % Value ( :, C % COARSENING_POLAR ) )
      call SetCoarseningAzimuthal &
             ( C  = A % Chart_GS_C, &
               R  = G % Storage_GS % Value ( :, G % CENTER_U_1 ), &
               Th = G % Storage_GS % Value ( :, G % CENTER_U_2 ), &
               CA = C % Storage_GS % Value ( :, C % COARSENING_AZIMUTHAL ) )

      call SetBlocks &
             (  BP  = C % Storage_GS % Value ( :, C % N_BLOCKS_POLAR ), &
                BA  = C % Storage_GS % Value ( :, C % N_BLOCKS_AZIMUTHAL ), &
                BL  = C % Storage_GS % Value ( :, C % BLOCK_LABEL ), &
                C   = A % Chart_GS_C, &
                CP  = C % Storage_GS % Value ( :, C % COARSENING_POLAR ), &
                CA  = C % Storage_GS % Value ( :, C % COARSENING_AZIMUTHAL ), &
               iTh  = C % iTheta, &
               iPh  = C % iPhi, &
               iRad = C % iRadius, &
               nBC  = C % nBlocksCoarsen )

      if ( C % DeviceMemory ) &
        call C % UpdateDevice ( )

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_C', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Initialize_C


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
      class is ( Atlas_SCG_C_Form )
    associate &
      ( G        =>  C % Geometry, &
        C_GS_C  =>  A % Chart_GS_C )

    if ( C_GS_C % nDimensions  ==  1 ) &
      return
      
    call FS % Storage_GS % ReassociateHost &
           ( AssociateVariablesOption = .false. )

    call C_GS_C % SetFieldPointer &
           ( FS % Storage_GS % Value, FS_4D )
    call C_GS_C % SetFieldPointer &
           ( G % Storage_GS % Value ( :, G % VOLUME ), dV_3D )
       
    call ComputeKernel &
           ( FS_4D, dV_3D, &
             iTh = C % iTheta, &
             iPh = C % iPhi, &
             iR  = C % iRadius, &
             iaS = FS % iaSelected, &
             oC  = C_GS_C % nGhostLayers, &
             nBC = C % nBlocksCoarsen, &
             UseDeviceOption = C % DeviceMemory )

    call FS % Storage_GS % ReassociateHost &
           ( AssociateVariablesOption = .true. )

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
      CoarsenPolar: do
        if ( CP ( iV )  *  R ( iV )  *  dTheta  >  C % MinWidth ) &
          exit CoarsenPolar
        CP ( iV )  =  2.0_KDR  *  CP ( iV )
      end do CoarsenPolar
    end do !-- iV

  end subroutine SetCoarseningPolar


  subroutine SetCoarseningAzimuthal ( C, R, Th, CA )
    
    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      Th
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      CA

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      dPhi

    dPhi  =  2.0_KDR * CONSTANT % PI  /  ( 2  *  C % nCellsPolar )

    do iV = 1, size ( CA )
      if ( .not. C % ProperCell ( iV ) ) &
        cycle
      CA ( iV )  =  1.0_KDR
      if ( C % nDimensions  ==  3 ) then
        CoarsenAzimuthal: do
          if ( CA ( iV )  *  R ( iV )  *  sin ( Th ( iV ) )  * dPhi  &
               >  C % MinWidth ) &
            exit CoarsenAzimuthal
          CA ( iV )  =  2.0_KDR  *  CA ( iV )
        end do CoarsenAzimuthal
      end if !-- nDimensions == 3
    end do !-- iV

  end subroutine SetCoarseningAzimuthal


  subroutine SetBlocks ( BP, BA, BL, C, CP, CA, iTh, iPh, iRad, nBC )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      BP, &  !-- nBlocksPolar
      BA, &  !-- nBlocksAzimuthal
      BL     !-- BlockLabel
    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      CP, &  !-- CoarsenPolar
      CA     !-- CoarsenAzimuthal
    integer ( KDI ), dimension ( :, : ), intent ( out ), allocatable :: &
      iTh, &
      iPh
    integer ( KDI ), dimension ( : ), intent ( out ), allocatable :: &
      iRad
    integer ( KDI ), intent ( out ) :: &
      nBC  !-- nBlocksCoarsen

    integer ( KDI ) :: &
      iR, &            !-- iRadius
      iTheta, &
      iTh_1, iTh_2, &  !-- iTheta_1, iTheta_2
      iPh_1, iPh_2, &  !-- iPhi_1, iPhi_2
      iBC, &           !-- iBlockCoarsen
      iBP, &           !-- iBlockPolar
      iBA, &           !-- iBlockAzimuthal
      oTh, &           !-- oTheta
      oPh, &           !-- oPhi
      CA_Max
    integer ( KDI ), dimension ( : ), allocatable :: &
      nCP, &  !-- nCoarsenPolar
      nBP     !-- nBlocksPolar
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      nCA, &  !-- nCoarsenAzimuthal
      nBA     !-- nBlocksAzimuthal
    real ( KDR ) :: &
      RandomLabel
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      CP_3D, &
      CA_3D, &
      BP_3D, &
      BA_3D, &
      BL_3D

    call C % SetFieldPointer ( CP, CP_3D )
    call C % SetFieldPointer ( CA, CA_3D )
    call C % SetFieldPointer ( BP, BP_3D )
    call C % SetFieldPointer ( BA, BA_3D )
    call C % SetFieldPointer ( BL, BL_3D )

    associate &
      ( nR   =>  C % nCellsBrick ( 1 ), &
        nTh  =>  C % nCellsBrick ( 2 ), &  !-- nCellsBrick ( 2 ) == nCells ( 2 )
        nPh  =>  C % nCellsBrick ( 3 ) )   !-- nCellsBrick ( 3 ) == nCells ( 3 )
    allocate &
      ( nCP ( nR ), nBP ( nR ), &
        nCA ( nR ), nBA ( nR ) )

    nBC  =  0

    select case ( C % nDimensions )
    case ( 2 )

      !-- Set nCoarsenPolar and nBlocksPolar
      do iR  =  1,  nR
        nCP ( iR )  =  CP_3D ( iR, 1, 1 )  +  0.5_KDR
        if ( nCP ( iR )  <=  1 ) then
          nBP ( iR )  =  0
        else if ( nCP ( iR )  >  1  .and.  nCP ( iR )  <=  32 ) then
          nBP ( iR )  =  nTh  /  nCP ( iR )
        else if ( nCP ( iR )  >  32 ) then !-- average over full polar range
          nBP ( iR )  =  1
        end if
      end do !-- iR

      !-- Set nBlocksCoarsen
      nBC  =  sum ( nBP )

      !-- Set block extents

      allocate ( iRad ( nBC ) )
      allocate ( iTh ( 2, nBC ) )
      allocate ( iPh ( 2, nBC ) )

      iBC  =  0
      do iR  =  1,  nR
        if ( nBP ( iR )  >  1 ) then
          do iBP  =  1,  nBP ( iR )

            iBC  =  iBC + 1

            iRad ( iBC )  =  iR

            oTh  =  ( iBP - 1 )  *  nCP ( iR )
            iTh ( 1 : 2, iBC )  =  [ oTh  +  1, oTh  +  nCP ( iR ) ]
            iPh ( 1 : 2, iBC )  =  [ 1, 1 ]

            !-- Visualization of nBlocksPolar
            iTh_1  =  iTh ( 1, iBC )
            iTh_2  =  iTh ( 2, iBC )
            iPh_1  =  iPh ( 1, iBC )
            iPh_2  =  iPh ( 2, iBC )
            BP_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  =  nTh  /  nCP ( iR )

          end do !-- iBP
        else if ( nBP ( iR )  ==  1 ) then  !-- average over full polar range

          iBC  =  iBC + 1

          iRad ( iBC )  =  iR

          iTh ( 1 : 2, iBC )  =  [ 1, nTh ]
          iPh ( 1 : 2, iBC )  =  [ 1, 1 ]          

          !-- Visualization of nBlocksPolar
          iTh_1  =  iTh ( 1, iBC )
          iTh_2  =  iTh ( 2, iBC )
          iPh_1  =  iPh ( 1, iBC )
          iPh_2  =  iPh ( 2, iBC )
          BP_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  =  1

        end if
      end do !-- iR

    case ( 3 )

      !-- Set nCoarsenPolar and nBlocksPolar
      do iR  =  1, nR
        nCP ( iR )  =  CP_3D ( iR, 1, 1 )  +  0.5_KDR
        nBP ( iR )  =  0
        if ( nCP ( iR )  <=  32 ) then
          do iTheta  =  1,  nTh,  nCP ( iR )
            CA_Max  &
              =  maxval ( CA_3D ( iR, iTheta : iTheta + nCP ( iR ) - 1, 1 ) ) &
                 +  0.5_KDR
            CA_Max  =  min ( CA_Max, nPh )
            if ( CA_Max  >  1 ) &
              nBP ( iR )  =  nBP ( iR )  +  1
          end do !-- iTheta
        else  !-- average over full polar range
          CA_Max  =  maxval ( CA_3D ( iR, 1 : nTh, 1 ) )  +  0.5_KDR
          CA_Max  =  min ( CA_Max, nPh )
          if ( CA_Max  >  1 ) &
            nBP ( iR )  =  nBP ( iR )  +  1
        end if
      end do !-- iR

      !-- Set nCoarsenAzimuthal and nBlocksAzimuthal
      do iR  =  1, nR
        call nCA ( iR ) % Initialize ( nBP ( iR ) )
        call nBA ( iR ) % Initialize ( nBP ( iR ) )
        iBP  =  0
        if ( nBP ( iR )  >  1 ) then
          do iTheta  =  1,  nTh,  nCP ( iR )
            CA_Max  &
              =  maxval ( CA_3D ( iR, iTheta : iTheta + nCP ( iR ) - 1, 1 ) ) &
                 +  0.5_KDR
            CA_Max  =  min ( CA_Max, nPh )
            if ( CA_Max  >  1 ) then
              iBP  =  iBP + 1
              if ( CA_Max  <=  2 * 32 ) then
                nCA ( iR ) % Value ( iBP )  =  CA_Max
                nBA ( iR ) % Value ( iBP )  =  nPh  /  CA_Max
              else  !-- average over full azimuthal range
                nCA ( iR ) % Value ( iBP )  =  nPh
                nBA ( iR ) % Value ( iBP )  =  1
              end if !-- CA_Max <= 2 * 32
            end if !-- CA_Max > 1
          end do !-- iTheta
        else if ( nBP ( iR )  ==  1 ) then  !-- average over full polar range
          CA_Max  =  maxval ( CA_3D ( iR, 1 : nTh, 1 ) )  +  0.5_KDR
          CA_Max  =  min ( CA_Max, nPh )
          if ( CA_Max  >  1 ) then
            iBP  =  iBP + 1
            if ( CA_Max  <=  2 * 32 ) then
              nCA ( iR ) % Value ( iBP )  =  CA_Max
              nBA ( iR ) % Value ( iBP )  =  nPh  /  CA_Max
            else  !-- average over full azimuthal range
              nCA ( iR ) % Value ( iBP )  =  nPh
              nBA ( iR ) % Value ( iBP )  =  1
            end if !-- CA_Max <= 2 * 32
          end if !-- CA_Max > 1
        end if !-- nBP ( iR ) > 1
      end do !-- iR

      !-- Set nBlocksCoarsen
      nBC  =  0
      do iR  =  1,  nR
        nBC  =  nBC  +  sum ( nBA ( iR ) % Value )
      end do !-- iR

      !-- Set block extents

      allocate ( iRad ( nBC ) )
      allocate ( iTh ( 2, nBC ) )
      allocate ( iPh ( 2, nBC ) )

      iBC  =  0
      do iR  =  1, nR
        iBP  =  0
        if ( nBP ( iR )  >  1 ) then
          do iTheta  =  1,  nTh,  nCP ( iR )
            CA_Max  &
              =  maxval ( CA_3D ( iR, iTheta : iTheta + nCP ( iR ) - 1, 1 ) ) &
                 +  0.5_KDR
            CA_Max  =  min ( CA_Max, nPh )
            if ( CA_Max  >  1 ) then
              iBP  =  iBP + 1
              if ( CA_Max  <=  2 * 32 ) then
                do iBA  =  1,  nBA ( iR ) % Value ( iBP )

                  iBC  =  iBC + 1

                  iRad ( iBC )  =  iR

                  oTh  =  iTheta - 1
                  iTh ( 1 : 2, iBC )  =  [ oTh  +  1, oTh  +  nCP ( iR ) ]

                  oPh  =  ( iBA - 1 )  *  nCA ( iR ) % Value ( iBP )
                  iPh ( 1 : 2, iBC )  =  [ oPh  +  1, &
                                           oPh  +  nCA ( iR ) % Value ( iBP ) ]

                  !-- Visualization of nBlocksPolar and nBlocksAzimuthal
                  iTh_1  =  iTh ( 1, iBC )
                  iTh_2  =  iTh ( 2, iBC )
                  iPh_1  =  iPh ( 1, iBC )
                  iPh_2  =  iPh ( 2, iBC )
                  BP_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                    =  nTh  /  nCP ( iR )
                  BA_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                    =  nPh  /  CA_Max

                end do !-- iBA
              else  !-- average over full azimuthal range

                  iBC  =  iBC + 1

                  iRad ( iBC )  =  iR

                  oTh  =  iTheta - 1
                  iTh ( 1 : 2, iBC )  =  [ oTh  +  1, oTh  +  nCP ( iR ) ]

                  iPh ( 1 : 2, iBC )  =  [ 1, nPh ]

                  !-- Visualization of nBlocksPolar and nBlocksAzimuthal
                  iTh_1  =  iTh ( 1, iBC )
                  iTh_2  =  iTh ( 2, iBC )
                  iPh_1  =  iPh ( 1, iBC )
                  iPh_2  =  iPh ( 2, iBC )
                  BP_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                    =  nTh  /  nCP ( iR )
                  BA_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                    =  1

              end if
            end if
          end do !-- iTheta
        else if ( nBP ( iR )  ==  1 ) then  !-- average over full polar range
          CA_Max  =  maxval ( CA_3D ( iR, 1 : nTh, 1 ) )  +  0.5_KDR
          CA_Max  =  min ( CA_Max, nPh )
          if ( CA_Max  >  1 ) then
            iBP  =  iBP + 1
            if ( CA_Max  <=  2 * 32 ) then
              do iBA  =  1,  nBA ( iR ) % Value ( iBP )

                iBC  =  iBC + 1

                iRad ( iBC )  =  iR

                iTh ( 1 : 2, iBC )  =  [ 1, nTh ]

                oPh  =  ( iBA - 1 )  *  nCA ( iR ) % Value ( iBP )
                iPh ( 1 : 2, iBC )  =  [ oPh  +  1, &
                                         oPh  +  nCA ( iR ) % Value ( iBP ) ]

                !-- Visualization of nBlocksPolar and nBlocksAzimuthal
                iTh_1  =  iTh ( 1, iBC )
                iTh_2  =  iTh ( 2, iBC )
                iPh_1  =  iPh ( 1, iBC )
                iPh_2  =  iPh ( 2, iBC )
                BP_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                  =  1
                BA_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                  =  nPh  /  CA_Max

              end do !-- iBA
            else  !-- average over full azimuthal range

              iBC  =  iBC + 1

              iRad ( iBC )  =  iR

              iTh ( 1 : 2, iBC )  =  [ 1, nTh ]

              iPh ( 1 : 2, iBC )  =  [ 1, nPh ]

              !-- Visualization of nBlocksPolar and nBlocksAzimuthal
              iTh_1  =  iTh ( 1, iBC )
              iTh_2  =  iTh ( 2, iBC )
              iPh_1  =  iPh ( 1, iBC )
              iPh_2  =  iPh ( 2, iBC )
              BP_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                =  1
              BA_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  &
                =  1

            end if
          end if
        end if !-- nBP ( iR )  >  1
      end do !-- iR

    end select !-- nDimensions

    do iBC  =  1, nBC

      iR     =  iRad ( iBC )
      iTh_1  =  iTh ( 1, iBC )
      iTh_2  =  iTh ( 2, iBC )
      iPh_1  =  iPh ( 1, iBC )
      iPh_2  =  iPh ( 2, iBC )

      call random_number ( RandomLabel )
      BL_3D ( iR, iTh_1 : iTh_2, iPh_1 : iPh_2 )  =  RandomLabel

    end do !-- iBC

    end associate !-- nR, etc.

  end subroutine SetBlocks


end module Coarsening_C__Form
