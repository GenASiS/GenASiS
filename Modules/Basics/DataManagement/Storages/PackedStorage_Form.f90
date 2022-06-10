!-- PackedStorageForm loads and stores selected rows and columns of a
!   Storage data array into a contiguous data array. 

module PackedStorage_Form

  use Specifiers
  use ArrayOperations
  use Storage_Form
  
  implicit none
  private
  
  type, public :: PackedStorageForm
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
      Add => Add_PS
    procedure, public, pass :: &
      Store
    final :: &
      Finalize
  end type PackedStorageForm

  type, public :: PackedStorageElementForm
    class ( PackedStorageForm ), allocatable :: &
      Element
  contains
    final :: &
      FinalizeElement
  end type PackedStorageElementForm

    private :: &
      LoadVariable, &
      AddVariable, &
      StoreVariable

contains


  subroutine InitializeEmpty &
               ( PS, iaUnpacked, nValuesUnpacked, nVariables, ClearOption )

    class ( PackedStorageForm ), intent ( inout ) :: &
      PS
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

    PS % nVariables = nVariables    
    PS % nValues = size ( iaUnpacked )
    PS % nValuesUnpacked = nValuesUnpacked

    allocate ( PS % iaUnpacked ( PS % nValues ) )
    allocate ( PS % Value ( PS % nValues, PS % nVariables ) )
    if ( ClearRequested ) call Clear ( PS % Value )

    if ( size ( iaUnpacked ) > 0 ) &
      call Copy ( iaUnpacked, PS % iaUnpacked )

  end subroutine InitializeEmpty
  
  
  subroutine InitializeCopy ( PS, PS_Source )

    class ( PackedStorageForm ), intent ( inout ) :: &
      PS
    class ( PackedStorageForm ), intent ( in ) :: &
      PS_Source

    PS % nVariables = PS_Source % nVariables    
    PS % nValues = PS_Source % nValues
    PS % nValuesUnpacked = PS_Source % nValuesUnpacked

    allocate ( PS % iaUnpacked ( PS % nValues ) )
    allocate ( PS % Value ( PS % nValues, PS % nVariables ) )

    if ( size ( PS % iaUnpacked ) > 0 ) then
      call Copy ( PS_Source % iaUnpacked, PS % iaUnpacked )
      call Copy ( PS_Source % Value, PS % Value )
    end if

  end subroutine InitializeCopy


  subroutine Load ( PS, S )
    
    class ( PackedStorageForm ), intent ( inout ) :: &
      PS
    class ( StorageForm ), intent ( in ) :: &
      S
      
    integer ( KDI ) :: &
      iV  !-- iVariable
    
    if ( PS % nValues == 0 ) return
    
    do iV = 1, PS % nVariables
      call LoadVariable &
             ( PS % Value ( :, iV ), &
               S % Value ( :, S % iaSelected ( iV ) ), &
               PS % iaUnpacked, PS % nValues )
    end do

  end subroutine Load
  
    
  subroutine Add_PS ( PS, S )
    
    class ( PackedStorageForm ), intent ( inout ) :: &
      PS
    class ( StorageForm ), intent ( in ) :: &
      S
      
    integer ( KDI ) :: &
      iV  !-- iVariable
    
    if ( PS % nValues == 0 ) return
    
    do iV = 1, PS % nVariables
      call AddVariable &
             ( PS % Value ( :, iV ), &
               S % Value ( :, S % iaSelected ( iV ) ), &
               PS % iaUnpacked, PS % nValues )
    end do

  end subroutine Add_PS
  
    
  subroutine Store ( PS, S )
    
    class ( PackedStorageForm ), intent ( inout ) :: &
      PS
    class ( StorageForm ), intent ( inout ) :: &
      S
    
    integer ( KDI ) :: &
      iV  !-- iVariable
    
    if ( PS % nValues == 0 ) return
    
    do iV = 1, PS % nVariables
      call StoreVariable &
             ( S % Value ( :, S % iaSelected ( iV ) ), &
               PS % Value ( :, iV ), &
               PS % iaUnpacked, PS % nValues )
    end do
    
  end subroutine Store


  elemental subroutine Finalize ( PS )

    type ( PackedStorageForm ), intent ( inout ) :: &
      PS

    if ( allocated ( PS % Value ) ) deallocate ( PS % Value )
    if ( allocated ( PS % iaUnpacked ) ) deallocate ( PS % iaUnpacked )

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
  

  impure elemental subroutine FinalizeElement ( PSE )
    
    type ( PackedStorageElementForm ), intent ( inout ) :: &
      PSE

    if ( allocated ( PSE % Element ) ) deallocate ( PSE % Element )

  end subroutine FinalizeElement


end module PackedStorage_Form
