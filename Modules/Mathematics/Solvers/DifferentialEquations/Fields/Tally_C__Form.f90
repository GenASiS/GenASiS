!-- Tally_C handles volume and surface integrals related to the evolution of a
!   conserved current.

module Tally_C__Form

  !-- Tally_Current_Form
  
  use Basics
  use Manifolds
  use Operations
  use Current_Template

  implicit none
  private

  type, public :: Tally_C_Form
    integer ( KDI ) :: &
      N_INTEGRALS = 0
    integer ( KDI ) :: &
      nSelected = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected
    real ( KDR ), dimension ( : ), allocatable :: &
      Value
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    class ( AtlasHeaderForm ), pointer :: &
      Atlas => null ( )
  contains
    procedure, private, pass :: &
      InitializeConserved
    procedure, public, pass :: &
      SelectVariables
    generic, public :: &
      Initialize => InitializeConserved
    procedure, private, pass :: &
      ComputeInterior_CSL
    generic, public :: &
      ComputeInterior => ComputeInterior_CSL
    procedure, private, pass :: &
      ComputeBoundary_CSL
    generic :: &
      ComputeBoundary => ComputeBoundary_CSL
    procedure, private, pass :: &
      Show_T
    generic :: &
      Show => Show_T
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_CSL
    procedure, public, nopass :: &
      ComputeFacePositions
  end type Tally_C_Form

  type, public :: Tally_C_ElementForm
    class ( Tally_C_Form ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type Tally_C_ElementForm

  type, public :: Tally_C_PointerForm
    class ( Tally_C_Form ), pointer :: &
      Pointer => null ( )
  end type Tally_C_PointerForm

contains


  subroutine InitializeConserved ( T, C, A, VariableOption, UnitOption )

    class ( Tally_C_Form ), intent ( inout ) :: &
      T
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    integer ( KDI ) :: &
      iF, &  !-- iField
      iD     !-- iDimension
    class ( GeometryFlatForm ), pointer :: &
      G

    if ( T % N_INTEGRALS == 0 ) &
      T % N_INTEGRALS = C % N_CONSERVED

    if ( .not. allocated ( T % Value ) ) then
      allocate ( T % Value ( T % N_INTEGRALS ) )
      call Clear ( T % Value )
    end if

    if ( .not. allocated ( T % Variable ) ) &
      allocate ( T % Variable ( T % N_INTEGRALS ) )

    if ( present ( VariableOption ) ) then
      T % Variable ( 1 : C % N_CONSERVED ) &
        =  VariableOption ( 1 : C % N_CONSERVED )
    else
      associate ( iaC => C % iaConserved )
      do iF = 1, C % N_CONSERVED
        T % Variable ( iF ) &
          = trim ( C % Variable ( iaC ( iF ) ) ) // '_Integral'
      end do !-- iF
      end associate !-- iaC
    end if

    if ( .not. allocated ( T % Unit ) ) &
      allocate ( T % Unit ( T % N_INTEGRALS ) )    

    if ( present ( UnitOption ) ) then

      T % Unit ( 1 : C % N_CONSERVED ) = UnitOption ( 1 : C % N_CONSERVED )

    else 

      select type ( A )
      class is ( Atlas_SC_Form )
      associate ( CoordinateUnit => A % Chart % CoordinateUnit )   
      G => A % Geometry ( )

      associate ( iaC => C % iaConserved )
      do iF = 1, C % N_CONSERVED
        T % Unit ( iF ) &
          = C % Unit ( iaC ( iF ) )  *  G % Unit ( G % VOLUME )
      end do !-- iF
      end associate !-- iaC

      end associate !-- CoordinateUnit
      end select !-- A

    end if 

    call T % SelectVariables ( A )

    T % Atlas => A

    nullify ( G )

  end subroutine InitializeConserved


  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_C_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    integer ( KDI ) :: &
      iF  !-- iField

    T % nSelected = size ( T % Value )
    allocate ( T % iaSelected ( T % nSelected ) )
    T % iaSelected = [ ( iF, iF = 1, T % nSelected ) ]

  end subroutine SelectVariables


  subroutine ComputeInterior_CSL ( T, FC, ReduceOption )

    class ( Tally_C_Form ), intent ( inout ) :: &
      T
    class ( Field_CSL_Template ), intent ( in ) :: &
      FC
    logical ( KDL ), intent ( in ), optional :: &
      ReduceOption

    integer ( KDI ) :: &
      iS  !-- iSelected
    real ( KDR ), dimension ( T % nSelected ) :: &
      Integral
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      Integrand
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( VolumeIntegralForm ) :: &
      VI

    select type ( C => FC % Field )
    class is ( CurrentTemplate )

    select type ( CSL => FC % Chart )
    class is ( Chart_SL_Template )

    allocate ( Integrand ( T % nSelected ) )

    do iS = 1, T % nSelected
      call Integrand ( iS ) % Initialize ( C % nValues, ClearOption = .true. )
    end do !-- iF

    G => CSL % Geometry ( )
    call T % ComputeInteriorIntegrand &
           ( Integrand, C, G, CSL % nDimensions ) 
    
    call VI % Compute ( CSL, Integrand, Integral, ReduceOption )

    do iS = 1, T % nSelected
      T % Value ( T % iaSelected ( iS ) ) = Integral ( iS )
    end do !-- iS

    end select !-- CSL
    end select !-- C
    nullify ( G )

  end subroutine ComputeInterior_CSL


  subroutine ComputeBoundary_CSL ( T, FC, BoundaryFluence )

    class ( Tally_C_Form ), intent ( inout ) :: &
      T
    class ( Field_CSL_Template ), intent ( in ) :: &
      FC
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence  !-- boundary slab

    integer ( KDI ) :: &
      iD, &   !-- iDimension
      iS      !-- iSelected
    integer ( KDI ), dimension ( 3 ) :: &
      nB  !-- nBoundary
    real ( KDR ), dimension ( T % nSelected ) :: &
      Integral
    type ( Real_3D_Form ), dimension ( :, : ), allocatable :: &
      Integrand
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( SurfaceIntegralForm ) :: &
      SI

    select type ( C => FC % Field )
    class is ( CurrentTemplate )

    select type ( CSL => FC % Chart )
    class is ( Chart_SL_Template )

    associate ( Cnnct => CSL % Atlas % Connectivity )

    allocate ( Integrand ( T % nSelected, Cnnct % nFaces ) )

    do iD = 1, CSL % nDimensions
      nB = shape ( BoundaryFluence ( 1, Cnnct % iaInner ( iD ) ) % Value ) 
      do iS = 1, T % nSelected
        call Integrand ( iS, Cnnct % iaInner ( iD ) ) % Initialize &
               ( nB, ClearOption = .true. )
        call Integrand ( iS, Cnnct % iaOuter ( iD ) ) % Initialize &
               ( nB, ClearOption = .true. )
      end do !-- iS
    end do !-- iD

    G => CSL % Geometry ( )
    call T % ComputeBoundaryIntegrand_CSL &
           ( Integrand, C, CSL, G, BoundaryFluence ) 

    call SI % Compute ( CSL, Integrand, Integral, ReduceOption = .false. )

    do iS = 1, T % nSelected
      T % Value ( T % iaSelected ( iS ) ) &
        = T % Value ( T % iaSelected ( iS ) ) + Integral ( iS )
    end do !-- iS

    end associate !-- Cnnct

    end select !-- CSL
    end select !-- C
    nullify ( G )

  end subroutine ComputeBoundary_CSL


  subroutine Show_T &
               ( T, Description, IgnorabilityOption, nLeadingLinesOption, &
                 nTrailingLinesOption )
    
    class ( Tally_C_Form ), intent ( in ) :: &
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
    
    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption
    
    if ( Ignorability > CONSOLE % WARNING ) then
      if ( CONSOLE % ProcessRank /= CONSOLE % DisplayRank &
         .or. Ignorability > CONSOLE % Verbosity )  return
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

    type ( Tally_C_Form ), intent ( inout ) :: &
      T

    if ( allocated ( T % Unit ) ) &
      deallocate ( T % Unit ) 
    if ( allocated ( T % Variable ) ) &
      deallocate ( T % Variable )
    if ( allocated ( T % Value ) ) &
      deallocate ( T % Value )    
    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )    

  end subroutine Finalize


  impure elemental subroutine FinalizeElement ( TE )
    
    type ( Tally_C_ElementForm ), intent ( inout ) :: &
      TE

    if ( allocated ( TE % Element ) ) deallocate ( TE % Element )

  end subroutine FinalizeElement


  subroutine ComputeInteriorIntegrand &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_C_Form ), intent ( inout ) :: &
      T
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iI  !-- iIntegral
    
    associate ( iaC => C % iaConserved )
    do iI = 1, C % N_CONSERVED
      associate &
        ( CV => C % Value ( :, iaC ( iI ) ), &
          IV => Integrand ( iI ) % Value ( : ) )
      call Copy ( CV, IV )
      end associate !-- CV, etc.
    end do !-- iF
    end associate !-- iaC

  end subroutine ComputeInteriorIntegrand

  
  subroutine ComputeBoundaryIntegrand_CSL &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_C_Form ), intent ( inout ) :: &
      T
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence

    integer ( KDI ) :: &
      iD, &   !-- iDimension
      iF, &   !-- iFace
      iC, &   !-- iConnectivity
      iI      !-- iIntegral

    associate ( Cnnct => CSL % Atlas % Connectivity )
    do iD = 1, CSL % nDimensions
      do iF = 1, 2

        if ( iF == 1 ) then
          iC = Cnnct % iaInner ( iD )
        else if ( iF == 2 ) then
          iC = Cnnct % iaOuter ( iD )
        end if

        associate ( iaC => C % iaConserved )
        do iI = 1, C % N_CONSERVED
          associate &
            ( BFV => BoundaryFluence ( iI, iC ) % Value, &
              IV => Integrand ( iI, iC ) % Value )
          call Copy ( BFV, IV )
          end associate !-- BFV, etc.          
        end do !-- iI
      end associate !-- iaC
      
      end do !-- iF
    end do !-- iD
    end associate !-- Cnnct

  end subroutine ComputeBoundaryIntegrand_CSL


  subroutine ComputeFacePositions ( CSL, G, iD, iF, X_1, X_2, X_3 )

    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iF     !-- iFace
    real ( KDR ), dimension ( :, :, : ), intent ( out ), target :: &
      X_1, X_2, X_3

    integer ( KDI ), dimension ( 3 ) :: &
      oB   !-- oBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      XC_1, XC_2, XC_3, &
      dXL_1, dXL_2, dXL_3, &
      X_iD, dXL_iD

    call CSL % SetVariablePointer &
           ( G % Value ( :, G % CENTER_U ( 1 ) ), XC_1 )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % CENTER_U ( 2 ) ), XC_2 )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % CENTER_U ( 3 ) ), XC_3 )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), dXL_1 )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), dXL_2 )
    call CSL % SetVariablePointer &
           ( G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), dXL_3 )

    select case ( iD )
    case ( 1 )
        X_iD =>   X_1
      dXL_iD => dXL_1
    case ( 2 ) 
        X_iD =>   X_2
      dXL_iD => dXL_2
    case ( 3 ) 
        X_iD =>   X_3
      dXL_iD => dXL_3
    end select !-- iD

    !-- Geometry. Here proper cell indexing begins at 1
    select case ( iF )
      case ( 1 ) !-- inner
        oB = 0
      case ( 2 ) !-- outer
        oB = 0
        oB ( iD ) = oB ( iD ) + CSL % nCellsBrick ( iD )
    end select !-- iF
    ! X_1 = XC_1 ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
    !              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
    !              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
    ! X_2 = XC_2 ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
    !              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
    !              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
    ! X_3 = XC_3 ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
    !              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
    !              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
    call CopyCollapse ( XC_1, X_1, oB + CSL % nGhostLayers )
    call CopyCollapse ( XC_2, X_2, oB + CSL % nGhostLayers )
    call CopyCollapse ( XC_3, X_3, oB + CSL % nGhostLayers )
        
    !    X_iD  =  X_iD  -  dXL_iD ( oB ( 1 ) + 1 : oB ( 1 ) + nB ( 1 ), &
    !                              oB ( 2 ) + 1 : oB ( 2 ) + nB ( 2 ), &
    !                              oB ( 3 ) + 1 : oB ( 3 ) + nB ( 3 ) )
    call MultiplyAddCollapse &
           ( X_iD, dXL_iD, -1.0_KDR, oB + CSL % nGhostLayers )

    nullify ( XC_1, XC_2, XC_3, dXL_1, dXL_2, dXL_3, X_iD, dXL_iD )

  end subroutine ComputeFacePositions


end module Tally_C__Form
