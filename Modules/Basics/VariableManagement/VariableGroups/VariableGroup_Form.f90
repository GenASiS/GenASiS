!-- VariableGroupForm provides infrastructure in handling collection of
!   variables, typically sets of related physical fields. It includes both the
!   metadata about the variables (names, units, etc) and storage for the
!   variable data itself.

module VariableGroup_Form
  
  use Specifiers
  use ArrayOperations
  use ArrayArrays
    
  implicit none
  private
  
  type, public :: VariableGroupForm
    integer ( KDI ) :: &
      nValues    = 0, &
      nVariables = 0, &
      nVectors   = 0, &
      lName      = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected, &
      lVariable, &
      lVector
    real ( KDR ), dimension ( :, : ), pointer :: &
      Value => null (  )
    logical ( KDL ) :: &
      AllocatedValue = .false.
    character ( LDF ) :: &
      Name = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
  contains
    !-- FIXME: Changed "private" to "public" since marking this private
    !          cause problem in the overriding subroutine in extension
    !          with CCE-8.x (possibly related to bug CCS # 121803)
    procedure, public, pass :: &
      InitializeAllocate
    procedure, public, pass, non_overridable :: &
      InitializeAssociate
    procedure, public, pass :: &
      InitializeClone
    generic :: &
      Initialize => InitializeAllocate, InitializeAssociate, InitializeClone
    final :: &
      Finalize
  end type VariableGroupForm

    private :: &
      InitializeOptionalMembers

contains


  subroutine InitializeAllocate &
               ( VG, ValueShape, VectorIndicesOption, UnitOption, &
                 VectorOption, VariableOption, NameOption, ClearOption )
    
    class ( VariableGroupForm ), intent ( inout ) :: &
      VG
    integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
      ValueShape
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption, &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    
    integer ( KDI ) :: &
      iVrbl
    logical ( KDL ) :: &
      ClearRequested
    
    VG % nValues = ValueShape ( 1 )
    
    VG % nVariables = ValueShape ( 2 )

    allocate ( VG % iaSelected ( VG % nVariables ) )
    VG % iaSelected = [ ( iVrbl, iVrbl = 1, VG % nVariables ) ]

    allocate ( VG % Value ( VG % nValues, ValueShape ( 2 ) ) )
    VG % AllocatedValue = .true.
    
    ClearRequested = .false.
    if ( present ( ClearOption ) ) ClearRequested = ClearOption
    if ( ClearRequested ) call Clear ( VG % Value )  
    
    call InitializeOptionalMembers &
           ( VG, VectorIndicesOption, UnitOption, VectorOption, &
             VariableOption, NameOption )
  
  end subroutine InitializeAllocate
  
  
  subroutine InitializeAssociate &
               ( VG, Value, VectorIndicesOption, UnitOption, VectorOption, &
                 VariableOption, NameOption, iaSelectedOption )
    
    class ( VariableGroupForm ), intent ( inout ) :: &
      VG
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      Value
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption, &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption
    
    integer ( KDI ) :: &
      iVrbl

    VG % nValues = size ( Value, dim = 1 )
    
    if ( present ( iaSelectedOption ) ) then
      VG % nVariables = size ( iaSelectedOption )
    else
      VG % nVariables = size ( Value, dim = 2 )
    end if
    
    allocate ( VG % iaSelected ( VG % nVariables ) )
    if ( present ( iaSelectedOption ) ) then
      VG % iaSelected = iaSelectedOption
    else
      VG % iaSelected = [ ( iVrbl, iVrbl = 1, VG % nVariables ) ]
    end if

    VG % Value => Value
    VG % AllocatedValue = .false.
    
    call InitializeOptionalMembers &
           ( VG, VectorIndicesOption, UnitOption, VectorOption, &
             VariableOption, NameOption )
  
  end subroutine InitializeAssociate
  
  
  subroutine InitializeClone (  &
               VG_Target, VG_Source, VectorIndicesOption, VectorOption, &
               NameOption, iaSelectedOption )

    class ( VariableGroupForm ), intent ( inout ) :: &
      VG_Target
    class ( VariableGroupForm ), intent ( in ) :: &
      VG_Source
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ) :: &
      iV     !-- iVector

    VG_Target % nValues = VG_Source % nValues
    
    if ( present ( iaSelectedOption ) ) then
      VG_Target % nVariables = size ( iaSelectedOption )
    else
      VG_Target % nVariables = VG_Source % nVariables
    end if

    if ( .not. present ( VectorIndicesOption ) ) &
      VG_Target % nVectors = VG_Source % nVectors

    VG_Target % lName = VG_Source % lName

    if ( allocated ( VG_Source % lVariable ) ) then
      allocate ( VG_Target % lVariable ( size ( VG_Source % lVariable ) ) )
      VG_Target % lVariable = VG_Source % lVariable
    end if
      
    if ( .not. present ( VectorOption ) ) then
      allocate ( VG_Target % lVector ( size ( VG_Source % lVector ) ) )
      VG_Target % lVector = VG_Source % lVector 
    end if

    allocate ( VG_Target % iaSelected ( VG_Target % nVariables ) )
    if ( present ( iaSelectedOption ) ) then
      VG_Target % iaSelected = iaSelectedOption
    else
      VG_Target % iaSelected = VG_Source % iaSelected
    end if
  
    VG_Target % Value => VG_Source % Value
    VG_Target % AllocatedValue = .false.
    
    if ( .not. present ( NameOption ) ) &
      VG_Target % Name = trim ( VG_Source % Name )

    if ( allocated ( VG_Source % Variable ) ) then
      allocate ( VG_Target % Variable ( size (  VG_Source % Variable ) ) )
      VG_Target % Variable = VG_Source % Variable
    end if
      
    if ( .not. present ( VectorOption ) ) then
      allocate ( VG_Target % Vector ( size ( VG_Source % Vector ) ) )
      VG_Target % Vector = VG_Source % Vector 
    end if
      
    if ( allocated ( VG_Source % Unit ) ) then
      allocate ( VG_Target % Unit ( size ( VG_Source % Unit ) ) )
      VG_Target % Unit = VG_Source % Unit
    end if

    if ( .not. present ( VectorIndicesOption ) ) then
      allocate ( VG_Target % VectorIndices ( VG_Target % nVectors ) )
      do iV = 1, VG_Target % nVectors
        call VG_Target % VectorIndices ( iV ) % Initialize &
               ( VG_Source % VectorIndices ( iV ) )
      end do
    end if
      
    call InitializeOptionalMembers &
           ( VG_Target, VectorIndicesOption = VectorIndicesOption, &
             VectorOption = VectorOption, NameOption = NameOption )

  end subroutine InitializeClone


  impure elemental subroutine Finalize ( VG )

    type ( VariableGroupForm ), intent ( inout ) :: &
      VG

!-- FIXME: this deallocation in a cloned variable group with no vectors
!          caused trouble with Intel 12.1.2 
!    if ( allocated ( VG % VectorIndices ) ) deallocate ( VG % VectorIndices )
    if ( allocated ( VG % VectorIndices ) .and. VG % nVectors > 0 ) &
      deallocate ( VG % VectorIndices )

    if ( allocated ( VG % Unit ) )     deallocate ( VG % Unit )

    if ( allocated ( VG % Vector ) )   deallocate ( VG % Vector )
    if ( allocated ( VG % Variable ) ) deallocate ( VG % Variable )

    if ( VG % AllocatedValue ) then
      if ( associated ( VG % Value ) ) deallocate ( VG % Value )
    end if
    nullify ( VG % Value )

    if ( allocated ( VG % iaSelected ) )  deallocate ( VG % iaSelected )
    if ( allocated ( VG % lVector ) )   deallocate ( VG % lVector )
    if ( allocated ( VG % lVariable ) ) deallocate ( VG % lVariable )

  end subroutine Finalize
  

  subroutine InitializeOptionalMembers &
               ( VG, VectorIndicesOption, UnitOption, VectorOption, &
                 VariableOption, NameOption )
                 
    class ( VariableGroupForm ), intent ( inout ) :: &
      VG
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption, &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
      
    integer ( KDI ) :: &
      iS, &      !-- iSelected
      iVrbl, &   !-- iVariable
      iVctr      !-- iVector
      
    if ( present ( VectorIndicesOption ) ) then
      VG % nVectors = size ( VectorIndicesOption )
      allocate ( VG % VectorIndices ( VG % nVectors ) )
      do iVctr = 1, VG % nVectors
        call VG % VectorIndices ( iVctr ) % Initialize &
               ( VectorIndicesOption ( iVctr ) )
      end do
    end if
    
    if ( .not. allocated ( VG % VectorIndices ) ) &
      allocate ( VG % VectorIndices ( 0 ) )

    if ( present ( NameOption ) ) then
      VG % lName = len_trim ( NameOption )
      VG % Name = NameOption
    end if 
    
    if ( .not. allocated ( VG % lVariable ) ) then
      allocate ( VG % lVariable ( size ( VG % Value, dim = 2 ) ) )
      VG % lVariable = 0
    end if
    
    if ( .not. allocated ( VG % Variable ) ) then
      allocate ( VG % Variable ( size ( VG % Value, dim = 2 ) ) )
      VG % Variable = ''
    end if
    
    if ( present ( VariableOption ) ) then
      do iS = 1, VG % nVariables
        iVrbl = VG % iaSelected ( iS )
        VG % Variable ( iVrbl ) = VariableOption ( iS )
        VG % lVariable ( iVrbl ) = len_trim ( VariableOption ( iS ) )
      end do
    end if
      
    if ( present ( VectorOption ) ) then
      allocate ( VG % lVector ( size ( VectorOption ) ) )
      VG % lVector = len_trim ( VectorOption )
      allocate ( VG % Vector ( size ( VectorOption ) ) )
      VG % Vector = VectorOption
    end if
    
    if ( .not. allocated ( VG % lVector ) ) &
      allocate ( VG % lVector ( 0 ) )
    if ( .not. allocated ( VG % Vector ) ) &
      allocate ( VG % Vector ( 0 ) )

    if ( .not. allocated ( VG % Unit ) ) &
      allocate ( VG % Unit ( size ( VG % Value, dim = 2 ) ) )
    if ( present ( UnitOption ) ) VG % Unit = UnitOption
    
  end subroutine InitializeOptionalMembers


end module VariableGroup_Form
