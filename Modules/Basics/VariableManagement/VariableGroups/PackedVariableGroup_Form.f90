!-- PackedVariableGroupForm loads and stores selected rows and columns of a
!   VariableGroup data array into a contiguous data array. 

module PackedVariableGroup_Form

  use Specifiers
  use ArrayOperations
  use VariableGroup_Form
  
  implicit none
  private
  
  type, public :: PackedVariableGroupForm
    integer ( KDI ) :: &
      nValues = 0, &
      nValuesUnpacked = 0, &
      nVariables = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaUnpacked
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Value
  contains
    procedure, private, pass :: &
      InitializeEmpty
    procedure, private, pass :: &
      InitializeCopy
    generic :: &
      Initialize => InitializeEmpty, InitializeCopy
    procedure, public, pass :: &
      Load
    procedure, public, pass :: &
      Add
    procedure, public, pass :: &
      Store
    final :: &
      Finalize
  end type PackedVariableGroupForm

  type, public :: PackedVariableGroupElementForm
    class ( PackedVariableGroupForm ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type PackedVariableGroupElementForm

    private :: &
      LoadVariable, &
      AddVariable, &
      StoreVariable

contains


  subroutine InitializeEmpty &
               ( PVG, iaUnpacked, nValuesUnpacked, nVariables, ClearOption )

    class ( PackedVariableGroupForm ), intent ( inout ) :: &
      PVG
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaUnpacked
    integer ( KDI ), intent ( in ) :: &
      nValuesUnpacked, &
      nVariables
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption

    logical ( KDL ) :: &
      ClearRequested

    ClearRequested = .false.
    if ( present ( ClearOption ) ) ClearRequested = ClearOption

    PVG % nVariables = nVariables    
    PVG % nValues = size ( iaUnpacked )
    PVG % nValuesUnpacked = nValuesUnpacked

    allocate ( PVG % iaUnpacked ( PVG % nValues ) )
    allocate ( PVG % Value ( PVG % nValues, PVG % nVariables ) )
    if ( ClearRequested ) call Clear ( PVG % Value )

    if ( size ( iaUnpacked ) > 0 ) &
      call Copy ( iaUnpacked, PVG % iaUnpacked )

  end subroutine InitializeEmpty
  
  
  subroutine InitializeCopy ( PVG, PVG_Source )

    class ( PackedVariableGroupForm ), intent ( inout ) :: &
      PVG
    class ( PackedVariableGroupForm ), intent ( in ) :: &
      PVG_Source

    PVG % nVariables = PVG_Source % nVariables    
    PVG % nValues = PVG_Source % nValues
    PVG % nValuesUnpacked = PVG_Source % nValuesUnpacked

    allocate ( PVG % iaUnpacked ( PVG % nValues ) )
    allocate ( PVG % Value ( PVG % nValues, PVG % nVariables ) )

    if ( size ( PVG % iaUnpacked ) > 0 ) then
      call Copy ( PVG_Source % iaUnpacked, PVG % iaUnpacked )
      call Copy ( PVG_Source % Value, PVG % Value )
    end if

  end subroutine InitializeCopy


  subroutine Load ( PVG, VG )
    
    class ( PackedVariableGroupForm ), intent ( inout ) :: &
      PVG
    class ( VariableGroupForm ), intent ( in ) :: &
      VG
      
    integer ( KDI ) :: &
      iV  !-- iVariable
    
    if ( PVG % nValues == 0 ) return
    
    do iV = 1, PVG % nVariables
      call LoadVariable &
             ( PVG % Value ( :, iV ), &
               VG % Value ( :, VG % iaSelected ( iV ) ), &
               PVG % iaUnpacked, PVG % nValues )
    end do

  end subroutine Load
  
    
  subroutine Add ( PVG, VG )
    
    class ( PackedVariableGroupForm ), intent ( inout ) :: &
      PVG
    class ( VariableGroupForm ), intent ( in ) :: &
      VG
      
    integer ( KDI ) :: &
      iV  !-- iVariable
    
    if ( PVG % nValues == 0 ) return
    
    do iV = 1, PVG % nVariables
      call AddVariable &
             ( PVG % Value ( :, iV ), &
               VG % Value ( :, VG % iaSelected ( iV ) ), &
               PVG % iaUnpacked, PVG % nValues )
    end do

  end subroutine Add
  
    
  subroutine Store ( PVG, VG )
    
    class ( PackedVariableGroupForm ), intent ( inout ) :: &
      PVG
    class ( VariableGroupForm ), intent ( inout ) :: &
      VG
    
    integer ( KDI ) :: &
      iV  !-- iVariable
    
    if ( PVG % nValues == 0 ) return
    
    do iV = 1, PVG % nVariables
      call StoreVariable &
             ( VG % Value ( :, VG % iaSelected ( iV ) ), &
               PVG % Value ( :, iV ), &
               PVG % iaUnpacked, PVG % nValues )
    end do
    
  end subroutine Store


  elemental subroutine Finalize ( PVG )

    type ( PackedVariableGroupForm ), intent ( inout ) :: &
      PVG

    if ( allocated ( PVG % Value ) ) deallocate ( PVG % Value )
    if ( allocated ( PVG % iaUnpacked ) ) deallocate ( PVG % iaUnpacked )

  end subroutine Finalize
  

  subroutine LoadVariable &
               ( PackedValue, UnpackedValue, iaUnpacked, nValues )
      
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      PackedValue
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      UnpackedValue
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaUnpacked
    integer ( KDI ), intent ( in ) :: &
      nValues
    
    integer ( KDI ) :: &
      iV  !-- iValue
    
    !-- $OMP parallel do default ( shared ) private ( iV ) schedule ( static )
    
    do iV = 1, nValues
      PackedValue ( iV ) = UnpackedValue ( iaUnpacked ( iV ) )
    end do
    
    !-- $OMP end parallel do

  end subroutine LoadVariable


  subroutine AddVariable &
               ( PackedValue, UnpackedValue, iaUnpacked, nValues )
      
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      PackedValue
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      UnpackedValue
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaUnpacked
    integer ( KDI ), intent ( in ) :: &
      nValues
    
    integer ( KDI ) :: &
      iV  !-- iValue
    
    !-- $OMP parallel do default ( shared ) private ( iV ) schedule ( static )
    
    do iV = 1, nValues
      PackedValue ( iV ) &
        = PackedValue ( iV ) + UnpackedValue ( iaUnpacked ( iV ) )
    end do
    
    !-- $OMP end parallel do

  end subroutine AddVariable


  subroutine StoreVariable &
               ( UnpackedValue, PackedValue, iaUnpacked, nValues )
      
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      UnpackedValue
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      PackedValue
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaUnpacked
    integer ( KDI ), intent ( in ) :: &
      nValues
    
    integer ( KDI ) :: &
      iV  !-- iValue
    
    !-- $OMP parallel do default ( shared ) private ( iV ) schedule ( static )

    do iV = 1, nValues
      UnpackedValue ( iaUnpacked ( iV ) ) = PackedValue ( iV )
    end do
    
    !-- $OMP end parallel do

  end subroutine StoreVariable
  

  impure elemental subroutine FinalizeElement ( PVGE )
    
    type ( PackedVariableGroupElementForm ), intent ( inout ) :: &
      PVGE

    if ( allocated ( PVGE % Element ) ) deallocate ( PVGE % Element )

  end subroutine FinalizeElement


end module PackedVariableGroup_Form
