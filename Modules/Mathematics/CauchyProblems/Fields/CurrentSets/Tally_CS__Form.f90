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
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Unit
    class ( Geometry_F_Form ), pointer :: &
      Geometry => null ( )
    type ( VolumeIntegralForm ), allocatable :: &
      InteriorIntegral
    type ( SurfaceIntegralForm ), allocatable :: &
      BoundaryIntegral
  contains
    procedure, private, pass :: &
      InitializeBalanced
    generic, public :: &
      Initialize => InitializeBalanced
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      ComputeInterior
    procedure, public, pass :: &
      ComputeBoundary
    procedure, private, pass :: &
      Show_T
    generic :: &
      Show => Show_T
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand
    procedure, public, pass :: &
      ComputeBoundaryIntegrand
    procedure, public, nopass :: &
      ComputeFacePositions
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
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
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

    T % Geometry  =>  G

    call T % SelectVariables ( )

  end subroutine InitializeBalanced


  subroutine SelectVariables ( T ) 
    
    class ( Tally_CS_Form ), intent ( inout ) :: &
      T

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

    call T % ComputeInteriorIntegrand ( CS ) 

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


  subroutine ComputeBoundary ( T, CS, BoundaryFluence )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence  !-- boundary slab

    integer ( KDI ) :: &
      iS  !-- iSelected

    select type ( A  =>  T % Geometry % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( G  =>  T % Geometry, &
        C  =>  A % Chart_GS )

    if ( .not. allocated ( T % BoundaryIntegral ) ) then
      allocate ( T % BoundaryIntegral )
      associate ( BI  =>  T % BoundaryIntegral )
      call BI % Initialize &
             ( G, nIntegrals = T % nSelected, NameOption = 'TallyBoundary' )
      end associate !-- BI
    end if

    call T % ComputeBoundaryIntegrand ( CS, C, BoundaryFluence ) 

    associate ( BI  =>  T % BoundaryIntegral )
    call BI % Compute ( )
    do iS  =  1,  T % nSelected
      T % Value ( T % iaSelected ( iS ) ) &
        =  T % Value ( T % iaSelected ( iS ) )  +  BI % Output ( iS )
    end do !-- iS
    end associate !-- BI

    end associate !-- G, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_CS__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeBoundary_SCG', 'subroutine', CONSOLE % ERROR )
    end select !-- A

  end subroutine ComputeBoundary


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

    if ( allocated ( T % BoundaryIntegral ) ) &
      deallocate ( T % BoundaryIntegral )
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


  subroutine ComputeInteriorIntegrand ( T, CS )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iI  !-- iIntegral
    
    associate &
      (   G  =>  T % Geometry, &
          I  =>  T % InteriorIntegral % Integrand, &
        iaB  =>  T % iaBalanced )
    do iI  =  1,  T % nBalanced
      associate &
        ( CSV  =>  CS % Storage_GS % Value ( :, iaB ( iI ) ), &
           IV  =>   I % Storage_GS % Value ( :, iI ) )
      call Copy ( CSV, IV )
      end associate !-- CV, etc.
    end do !-- iI
    end associate !-- I, etc.

  end subroutine ComputeInteriorIntegrand

  
  subroutine ComputeBoundaryIntegrand ( T, CS, C, BF )

    class ( Tally_CS_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BF

    integer ( KDI ) :: &
      iD, &   !-- iDimension
      iF, &   !-- iFace
      iC, &   !-- iConnectivity
      iI      !-- iIntegral

    associate &
      ( I   =>  T % BoundaryIntegral % Integrand, &
        Cy  =>  C % Connectivity )

    do iD  =  1,  C % nDimensions
      do iF = 1, 2

        if ( iF == 1 ) then
          iC  =  Cy % iaInner ( iD )
        else if ( iF == 2 ) then
          iC  =  Cy % iaOuter ( iD )
        end if

        do iI  =  1,  T % nSelected
          associate &
            ( BFV  =>  BF ( iI, iC ) % Value, &
               IV  =>   I ( iI, iC ) % Value )
          call Copy ( BFV, IV )
          end associate !-- BFV, etc.          
        end do !-- iI
      
      end do !-- iF
    end do !-- iD

    end associate !-- Cy, etc.

  end subroutine ComputeBoundaryIntegrand


  subroutine ComputeFacePositions ( G, C, iD, iF, X_1, X_2, X_3 )

    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iD, &  !-- iDimension
      iF     !-- iFace
    real ( KDR ), dimension ( :, :, : ), intent ( out ), target :: &
      X_1, X_2, X_3

    integer ( KDI ), dimension ( 3 ) :: &
      oB, &   !-- oBoundary
      nB      !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      XC_1, XC_2, XC_3, &
      XE_I_1, XE_I_2, XE_I_3, &
      X_iD, &
      XE_I_iD

    associate ( GV  =>  G % Storage_GS % Value )
    call C % SetFieldPointer ( GV ( :, G % CENTER_U_1 ), XC_1 )
    call C % SetFieldPointer ( GV ( :, G % CENTER_U_2 ), XC_2 )
    call C % SetFieldPointer ( GV ( :, G % CENTER_U_3 ), XC_3 )
    call C % SetFieldPointer ( GV ( :, G % EDGE_I_U_1 ), XE_I_1 )
    call C % SetFieldPointer ( GV ( :, G % EDGE_I_U_2 ), XE_I_2 )
    call C % SetFieldPointer ( GV ( :, G % EDGE_I_U_3 ), XE_I_3 )
    end associate !-- GSV

    select case ( iD )
    case ( 1 )
         X_iD  =>  X_1
      XE_I_iD  =>  XE_I_1
    case ( 2 ) 
         X_iD  =>  X_2
      XE_I_iD  =>  XE_I_2
    case ( 3 ) 
         X_iD  =>  X_3
      XE_I_iD  =>  XE_I_3
    end select !-- iD

    !-- Geometry. Here proper cell indexing begins at 1
    select case ( iF )
    case ( 1 ) !-- inner
      oB  =  0
    case ( 2 ) !-- outer
      oB  =  0
      oB ( iD )  =  C % nCellsBrick ( iD )
    end select !-- iF

    nB  =  shape ( X_1 )

    X_1 = XC_1 ( oB ( 1 )  +  1  :  oB ( 1 )  +  nB ( 1 ), &
                 oB ( 2 )  +  1  :  oB ( 2 )  +  nB ( 2 ), &
                 oB ( 3 )  +  1  :  oB ( 3 )  +  nB ( 3 ) )
    X_2 = XC_2 ( oB ( 1 )  +  1  :  oB ( 1 )  +  nB ( 1 ), &
                 oB ( 2 )  +  1  :  oB ( 2 )  +  nB ( 2 ), &
                 oB ( 3 )  +  1  :  oB ( 3 )  +  nB ( 3 ) )
    X_3 = XC_3 ( oB ( 1 )  +  1  :  oB ( 1 )  +  nB ( 1 ), &
                 oB ( 2 )  +  1  :  oB ( 2 )  +  nB ( 2 ), &
                 oB ( 3 )  +  1  :  oB ( 3 )  +  nB ( 3 ) )
!     call CopyCollapse ( XC_1, X_1, oB + C % nGhostLayers )
!     call CopyCollapse ( XC_2, X_2, oB + C % nGhostLayers )
!     call CopyCollapse ( XC_3, X_3, oB + C % nGhostLayers )

    X_iD  =  XE_I_iD ( oB ( 1 )  +  1  :  oB ( 1 )  +  nB ( 1 ), &
                       oB ( 2 )  +  1  :  oB ( 2 )  +  nB ( 2 ), &
                       oB ( 3 )  +  1  :  oB ( 3 )  +  nB ( 3 ) )

  end subroutine ComputeFacePositions


end module Tally_CS__Form
