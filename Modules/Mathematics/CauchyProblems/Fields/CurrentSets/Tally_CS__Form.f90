module Tally_CS__Form

  !-- Tally_CurrentSet_Form
  
  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use CalculusFields

  implicit none
  private

  type, public :: Tally_CS_Form
    integer ( KDI ) :: &
      nIntegrals = 0, &
      nBalanced  = 0, &
      nSelected  = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaBalanced, &
      iaSelected
    real ( KDR ), dimension ( : ), allocatable :: &
      Value
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
    type ( VolumeIntegralForm ), allocatable :: &
      InteriorIntegral
  contains
    procedure, private, pass :: &
      InitializeBalanced
    generic, public :: &
      Initialize => InitializeBalanced
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      ComputeInterior
    procedure, private, pass :: &
      ComputeBoundary_SCG
    generic :: &
      ComputeBoundary => ComputeBoundary_SCG
    procedure, private, pass :: &
      Show_T
    generic :: &
      Show => Show_T
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand
!     procedure, public, pass :: &
!       ComputeBoundaryIntegrand_CSL
!     procedure, public, nopass :: &
!       ComputeFacePositions
  end type Tally_CS_Form

  type, public :: Tally_CS_Element
    class ( Tally_CS_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Tally_CS_Element

!   type, public :: Tally_C_PointerForm
!     class ( Tally_C_Form ), pointer :: &
!       Pointer => null ( )
!   end type Tally_C_PointerForm


contains


  subroutine InitializeBalanced &
               ( T, CS, G, iaBalanced, VariableOption, UnitOption )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaBalanced
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    integer ( KDI ) :: &
      iF, &  !-- iField
      iD     !-- iDimension

    T % nBalanced  =  size ( iaBalanced )
    allocate ( T % iaBalanced, source = iaBalanced )

    if ( T % nIntegrals  ==  0 ) &
      T % nIntegrals  =  T % nBalanced

    if ( .not. allocated ( T % Value ) ) then
      allocate ( T % Value ( T % nIntegrals ) )
      call Clear ( T % Value )
    end if

    if ( .not. allocated ( T % Variable ) ) &
      allocate ( T % Variable ( T % nIntegrals ) )

    if ( present ( VariableOption ) ) then
      T % Variable ( 1 : T % nBalanced ) &
        =  VariableOption ( 1 : T % nBalanced )
    else
      associate ( iaB  =>  T % iaBalanced )
      do iF  =  1,  T % nBalanced
        T % Variable ( iF ) &
          = trim ( CS % Field ( iaB ( iF ) ) ) // '_Integral'
      end do !-- iF
      end associate !-- iaB
    end if

    if ( .not. allocated ( T % Unit ) ) &
      allocate ( T % Unit ( T % nIntegrals ) )    

    if ( present ( UnitOption ) ) then

      T % Unit ( 1 : T % nBalanced )  =  UnitOption ( 1 : T % nBalanced )

    else 

      select type ( A  =>  G % Atlas )
        class is ( Atlas_SCG_Form )

      associate ( iaB => T % iaBalanced )
      do iF  =  1,  T % nBalanced
        T % Unit ( iF ) &
          =  CS % Unit ( iaB ( iF ), 1 )  *  G % Unit ( G % VOLUME, 1 )
      end do !-- iF
      end associate !-- iaB

      end select !-- A

    end if 

    call T % SelectVariables ( G % Atlas )

    T % Geometry  =>  G

  end subroutine InitializeBalanced


  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( Atlas_H_Form ), intent ( in ) :: &
      A

    integer ( KDI ) :: &
      iF  !-- iField

    T % nSelected = size ( T % Value )
    allocate ( T % iaSelected ( T % nSelected ) )
    T % iaSelected = [ ( iF, iF = 1, T % nSelected ) ]

  end subroutine SelectVariables


  subroutine ComputeInterior ( T, CS, ReduceOption )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    logical ( KDL ), intent ( in ), optional :: &
      ReduceOption

    integer ( KDI ) :: &
      iS  !-- iSelected

    select type ( A  =>  T % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( G  =>  T % Geometry, &
        C  =>  A % Chart_GS )

    if ( .not. allocated ( T % InteriorIntegral ) ) then
      allocate ( T % InteriorIntegral )
      associate ( II  =>  T % InteriorIntegral )
      call II % Initialize &
             ( G, nIntegrals = T % nSelected, NameOption = 'TallyInterior' )
      end associate !-- II
    end if

    call T % ComputeInteriorIntegrand &
           ( T % InteriorIntegral % Integrand, CS, G, C % nDimensions ) 

    associate ( II  =>  T % InteriorIntegral )
    call II % Compute ( ReduceOption )
    do iS  =  1,  T % nSelected
      T % Value ( T % iaSelected ( iS ) )  =  II % Output ( iS )
    end do !-- iS
    end associate !-- II, etc.

    end associate !-- G, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_CS__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeInterior', 'subroutine', CONSOLE % ERROR )
    end select !-- A

  end subroutine ComputeInterior


  subroutine ComputeBoundary_SCG ( T, CS, BoundaryFluence )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence  !-- boundary slab

!     integer ( KDI ) :: &
!       iD, &   !-- iDimension
!       iS      !-- iSelected
!     integer ( KDI ), dimension ( 3 ) :: &
!       nB  !-- nBoundary
!     real ( KDR ), dimension ( T % nSelected ) :: &
!       Integral
!     type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
!       Integrand
!     type ( SurfaceIntegralForm ) :: &
!       SI

    select type ( A  =>  T % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( G  =>  T % Geometry, &
        C  =>  A % Chart_GS )

!     associate ( Cnnct => CSL % Atlas % Connectivity )

!     allocate ( Integrand ( T % nSelected, Cnnct % nFaces ) )

!     do iD = 1, CSL % nDimensions
!       nB = shape ( BoundaryFluence ( 1, Cnnct % iaInner ( iD ) ) % Value ) 
!       do iS = 1, T % nSelected
!         call Integrand ( iS, Cnnct % iaInner ( iD ) ) % Initialize &
!                ( nB, ClearOption = .true. )
!         call Integrand ( iS, Cnnct % iaOuter ( iD ) ) % Initialize &
!                ( nB, ClearOption = .true. )
!       end do !-- iS
!     end do !-- iD

!     G => CSL % Geometry ( )
!     call T % ComputeBoundaryIntegrand_CSL &
!            ( Integrand, C, CSL, G, BoundaryFluence ) 

!     call SI % Compute ( CSL, Integrand, Integral, ReduceOption = .false. )

!     do iS = 1, T % nSelected
!       T % Value ( T % iaSelected ( iS ) ) &
!         = T % Value ( T % iaSelected ( iS ) ) + Integral ( iS )
!     end do !-- iS

!     end associate !-- Cnnct

!     end select !-- CSL
!     end select !-- C

    end associate !-- G, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_CS__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeBoundary_SCG', 'subroutine', CONSOLE % ERROR )
    end select !-- A

  end subroutine ComputeBoundary_SCG


  subroutine Show_T ( T, Description, IgnorabilityOption, &
                      nLeadingLinesOption, nTrailingLinesOption )
    
    class ( Tally_CS_Form ), intent ( in ) :: &
      T
    character ( * ), intent ( in ) :: &
      Description
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nLeadingLinesOption, &
      nTrailingLinesOption
      
    integer ( KDI ) :: &
      iV, &
      iS, &
      iLine, &
      Ignorability
    
    Ignorability  =  CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) &
      Ignorability = IgnorabilityOption
    
    if ( Ignorability  >  CONSOLE % WARNING ) then
      if ( CONSOLE % ProcessRank  /=  CONSOLE % DisplayRank &
           .or. Ignorability  >  CONSOLE % Verbosity ) &
        return
    end if

    if ( present ( nLeadingLinesOption )  ) then
      do iLine = 1, nLeadingLinesOption
        print *
      end do
    end if
    
    call Show ( trim ( Description ), Ignorability )
    
    do iV = 1, size ( T % iaSelected )
      iS = T % iaSelected ( iV )
      call Show ( T % Value ( iS ), T % Unit ( iS ), T % Variable ( iS ), &
                  Ignorability )
    end do
    
    if ( present ( nTrailingLinesOption )  ) then
      do iLine = 1, nTrailingLinesOption
        print *
      end do
    end if
  
  end subroutine Show_T
  
  
  impure elemental subroutine Finalize ( T )

    type ( Tally_CS_Form ), intent ( inout ) :: &
      T

    nullify ( T % Geometry )

    if ( allocated ( T % InteriorIntegral ) ) &
      deallocate ( T % InteriorIntegral )
    if ( allocated ( T % Unit ) ) &
      deallocate ( T % Unit ) 
    if ( allocated ( T % Variable ) ) &
      deallocate ( T % Variable )
    if ( allocated ( T % Value ) ) &
      deallocate ( T % Value )    
    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )    
    if ( allocated ( T % iaBalanced ) ) &
      deallocate ( T % iaBalanced )    

  end subroutine Finalize


  impure elemental subroutine Finalize_E ( TE )
    
    type ( Tally_CS_Element ), intent ( inout ) :: &
      TE

    if ( allocated ( TE % Element ) ) &
      deallocate ( TE % Element )

  end subroutine Finalize_E


  subroutine ComputeInteriorIntegrand ( T, I, CS, G, nD )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    type ( FieldSetForm ), intent ( inout ) :: &
      I
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nD

    integer ( KDI ) :: &
      iI  !-- iIntegral
    
    associate ( iaB  =>  T % iaBalanced )
    do iI  =  1,  T % nBalanced
      associate &
        ( CSV  =>  CS % Storage_GS % Value ( :, iaB ( iI ) ), &
           IV  =>   I % Storage_GS % Value ( :, iI ) )
      call Copy ( CSV, IV )
      end associate !-- CV, etc.
    end do !-- iI
    end associate !-- iaB

  end subroutine ComputeInteriorIntegrand

  
!   subroutine ComputeBoundaryIntegrand_CSL &
!                ( T, Integrand, C, CSL, G, BoundaryFluence )

!     class ( Tally_C_Form ), intent ( inout ) :: &
!       T
!     type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
!       Integrand
!     class ( CurrentTemplate ), intent ( in ) :: &
!       C
!     class ( Chart_SL_Template ), intent ( in ) :: &
!       CSL
!     class ( GeometryFlatForm ), intent ( in ) :: &
!       G
!     type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
!       BoundaryFluence

!     integer ( KDI ) :: &
!       iD, &   !-- iDimension
!       iF, &   !-- iFace
!       iC, &   !-- iConnectivity
!       iI      !-- iIntegral

!     associate ( Cnnct => CSL % Atlas % Connectivity )
!     do iD = 1, CSL % nDimensions
!       do iF = 1, 2

!         if ( iF == 1 ) then
!           iC = Cnnct % iaInner ( iD )
!         else if ( iF == 2 ) then
!           iC = Cnnct % iaOuter ( iD )
!         end if

!         associate ( iaC => C % iaConserved )
!         do iI = 1, C % N_CONSERVED
!           associate &
!             ( BFV => BoundaryFluence ( iI, iC ) % Value, &
!               IV => Integrand ( iI, iC ) % Value )
!           call Copy ( BFV, IV )
!           end associate !-- BFV, etc.          
!         end do !-- iI
!       end associate !-- iaC
      
!       end do !-- iF
!     end do !-- iD
!     end associate !-- Cnnct

!   end subroutine ComputeBoundaryIntegrand_CSL


!   subroutine ComputeFacePositions ( CSL, G, iD, iF, X_1, X_2, X_3 )

!     class ( Chart_SL_Template ), intent ( in ) :: &
!       CSL
!     class ( GeometryFlatForm ), intent ( in ) :: &
!       G
!     integer ( KDI ), intent ( in ) :: &
!       iD, &  !-- iDimension
!       iF     !-- iFace
!     real ( KDR ), dimension ( :, :, : ), intent ( out ), target :: &
!       X_1, X_2, X_3

!     integer ( KDI ), dimension ( 3 ) :: &
!       oB   !-- oBoundary
!     real ( KDR ), dimension ( :, :, : ), pointer :: &
!       XC_1, XC_2, XC_3, &
!       dXL_1, dXL_2, dXL_3, &
!       X_iD, dXL_iD

!     call CSL % SetVariablePointer &
!            ( G % Value ( :, G % CENTER_U ( 1 ) ), XC_1 )
!     call CSL % SetVariablePointer &
!            ( G % Value ( :, G % CENTER_U ( 2 ) ), XC_2 )
!     call CSL % SetVariablePointer &
!            ( G % Value ( :, G % CENTER_U ( 3 ) ), XC_3 )
!     call CSL % SetVariablePointer &
!            ( G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), dXL_1 )
!     call CSL % SetVariablePointer &
!            ( G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), dXL_2 )
!     call CSL % SetVariablePointer &
!            ( G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), dXL_3 )

!     select case ( iD )
!     case ( 1 )
!         X_iD =>   X_1
!       dXL_iD => dXL_1
!     case ( 2 ) 
!         X_iD =>   X_2
!       dXL_iD => dXL_2
!     case ( 3 ) 
!         X_iD =>   X_3
!       dXL_iD => dXL_3
!     end select !-- iD

!     !-- Geometry. Here proper cell indexing begins at 1
!     select case ( iF )
!       case ( 1 ) !-- inner
!         oB = 0
!       case ( 2 ) !-- outer
!         oB = 0
!         oB ( iD ) = oB ( iD ) + CSL % nCellsBrick ( iD )
!     end select !-- iF
!     ! X_1 = XC_1 ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
!     !              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
!     !              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
!     ! X_2 = XC_2 ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
!     !              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
!     !              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
!     ! X_3 = XC_3 ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
!     !              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
!     !              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
!     call CopyCollapse ( XC_1, X_1, oB + CSL % nGhostLayers )
!     call CopyCollapse ( XC_2, X_2, oB + CSL % nGhostLayers )
!     call CopyCollapse ( XC_3, X_3, oB + CSL % nGhostLayers )
        
!     !    X_iD  =  X_iD  -  dXL_iD ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
!     !                              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
!     !                              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
!     call MultiplyAddCollapse &
!            ( X_iD, dXL_iD, -1.0_KDR, oB + CSL % nGhostLayers )

!     nullify ( XC_1, XC_2, XC_3, dXL_1, dXL_2, dXL_3, X_iD, dXL_iD )

!   end subroutine ComputeFacePositions


end module Tally_CS__Form
