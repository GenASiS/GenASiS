!-- StorageForm provides infrastructure in handling collection of
!   variables, typically sets of related physical fields. It includes both the
!   metadata about the variables (names, units, etc) and storage for the
!   variable data itself.

module Storage_Form
  
  use iso_c_binding
  use Specifiers
  use Devices
  use ArrayOperations
  use ArrayArrays
    
  implicit none
  private
  
  type, public :: StorageForm
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
      AllocatedValue = .false., &
      Pinned = .false.
    character ( LDF ) :: &
      Name = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector
    type ( c_ptr ), dimension ( : ), allocatable :: &
      D_Selected     !-- Device pointer for Selected Value
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
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_S
    procedure, private, pass :: &
      UpdateDeviceAll
    procedure, private, pass :: &
      UpdateDeviceSingle
    generic :: &
      UpdateDevice => UpdateDeviceAll, UpdateDeviceSingle
    procedure, private, pass :: &
      UpdateHostAll
    procedure, private, pass :: &
      UpdateHostSingle
    generic :: &
      UpdateHost => UpdateHostAll, UpdateHostSingle
    final :: &
      Finalize
  end type StorageForm

    private :: &
      InitializeOptionalMembers

contains


  subroutine InitializeAllocate &
               ( S, ValueShape, VectorIndicesOption, UnitOption, &
                 VectorOption, VariableOption, NameOption, ClearOption, &
                 PinnedOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      S
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
      ClearOption, &
      PinnedOption
    
    integer ( KDI ) :: &
      iVrbl
    logical ( KDL ) :: &
      ClearRequested

    S % nValues = ValueShape ( 1 )
    
    S % nVariables = ValueShape ( 2 )

    allocate ( S % iaSelected ( S % nVariables ) )
    S % iaSelected = [ ( iVrbl, iVrbl = 1, S % nVariables ) ]

    S % Pinned = .false.
    if ( present ( PinnedOption ) ) S % Pinned = PinnedOption
    if ( S % Pinned ) then
      call AllocateHost ( S % Value, [ S % nValues, ValueShape ( 2 ) ] )
    else
      allocate ( S % Value ( S % nValues, ValueShape ( 2 ) ) )
    end if
    S % AllocatedValue = .true.
    
    ClearRequested = .false.
    if ( present ( ClearOption ) ) ClearRequested = ClearOption
    if ( ClearRequested ) call Clear ( S % Value )  
    
    allocate ( S % D_Selected ( S % nVariables ) )
    S % D_Selected = c_null_ptr
    
    call InitializeOptionalMembers &
           ( S, VectorIndicesOption, UnitOption, VectorOption, &
             VariableOption, NameOption )
  
  end subroutine InitializeAllocate
  
  
  subroutine InitializeAssociate &
               ( S, Value, VectorIndicesOption, UnitOption, VectorOption, &
                 VariableOption, NameOption, iaSelectedOption )
    
    class ( StorageForm ), intent ( inout ) :: &
      S
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

    S % nValues = size ( Value, dim = 1 )
    
    if ( present ( iaSelectedOption ) ) then
      S % nVariables = size ( iaSelectedOption )
    else
      S % nVariables = size ( Value, dim = 2 )
    end if
    
    allocate ( S % iaSelected ( S % nVariables ) )
    if ( present ( iaSelectedOption ) ) then
      S % iaSelected = iaSelectedOption
    else
      S % iaSelected = [ ( iVrbl, iVrbl = 1, S % nVariables ) ]
    end if

    S % Value => Value
    S % AllocatedValue = .false.
    
    allocate ( S % D_Selected ( S % nVariables ) )
    S % D_Selected = c_null_ptr
    
    call InitializeOptionalMembers &
           ( S, VectorIndicesOption, UnitOption, VectorOption, &
             VariableOption, NameOption )
  
  end subroutine InitializeAssociate
  
  
  subroutine InitializeClone (  &
               S_Target, S_Source, VectorIndicesOption, VectorOption, &
               NameOption, iaSelectedOption )

    class ( StorageForm ), intent ( inout ) :: &
      S_Target
    class ( StorageForm ), intent ( in ) :: &
      S_Source
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ) :: &
      iV, &     !-- iVector
      iS_S, iS_T
    type ( c_ptr ), dimension ( : ), allocatable :: &
      Scratch

    S_Target % nValues = S_Source % nValues
    
    if ( present ( iaSelectedOption ) ) then
      S_Target % nVariables = size ( iaSelectedOption )
    else
      S_Target % nVariables = S_Source % nVariables
    end if

    if ( .not. present ( VectorIndicesOption ) ) &
      S_Target % nVectors = S_Source % nVectors

    S_Target % lName = S_Source % lName

    if ( allocated ( S_Source % lVariable ) ) then
      allocate ( S_Target % lVariable ( size ( S_Source % lVariable ) ) )
      S_Target % lVariable = S_Source % lVariable
    end if
      
    if ( .not. present ( VectorOption ) ) then
      allocate ( S_Target % lVector ( size ( S_Source % lVector ) ) )
      S_Target % lVector = S_Source % lVector 
    end if

    allocate ( S_Target % iaSelected ( S_Target % nVariables ) )
    if ( present ( iaSelectedOption ) ) then
      S_Target % iaSelected = iaSelectedOption
    else
      S_Target % iaSelected = S_Source % iaSelected
    end if
  
    S_Target % Value => S_Source % Value
    S_Target % AllocatedValue = .false.
    
    if ( .not. present ( NameOption ) ) &
      S_Target % Name = trim ( S_Source % Name )

    if ( allocated ( S_Source % Variable ) ) then
      allocate ( S_Target % Variable ( size (  S_Source % Variable ) ) )
      S_Target % Variable = S_Source % Variable
    end if
      
    if ( .not. present ( VectorOption ) ) then
      allocate ( S_Target % Vector ( size ( S_Source % Vector ) ) )
      S_Target % Vector = S_Source % Vector 
    end if
      
    if ( allocated ( S_Source % Unit ) ) then
      allocate ( S_Target % Unit ( size ( S_Source % Unit ) ) )
      S_Target % Unit = S_Source % Unit
    end if

    if ( .not. present ( VectorIndicesOption ) ) then
      allocate ( S_Target % VectorIndices ( S_Target % nVectors ) )
      do iV = 1, S_Target % nVectors
        call S_Target % VectorIndices ( iV ) % Initialize &
               ( S_Source % VectorIndices ( iV ) )
      end do
    end if
    
    allocate ( S_Target % D_Selected ( S_Target % nVariables ) )
    allocate ( Scratch, source = S_Source % D_Selected )
    if ( .not. present ( iaSelectedOption ) ) then
      S_Target % D_Selected = S_Source % D_Selected
    else
      S_Target % D_Selected = c_null_ptr
      do iS_T = 1, S_Target % nVariables
        iV = S_Target % iaSelected ( iS_T )
        do iS_S = 1, S_Source % nVariables
          if ( iV == S_Source % iaSelected ( iS_S ) ) then
             S_Target % D_Selected ( iS_T ) = Scratch ( iS_S )
!            S_Target % D_Selected ( iS_T ) = S_Source % D_Selected ( iS_S )
            exit
          end if
        end do
      end do
    end if
    
    call InitializeOptionalMembers &
           ( S_Target, VectorIndicesOption = VectorIndicesOption, &
             VectorOption = VectorOption, NameOption = NameOption )

  end subroutine InitializeClone
  
  
  subroutine AllocateDevice_S ( S )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
      
    integer ( KDI ) :: &
      iV
      
    do iV = 1, S % nVariables
      if ( c_associated ( S % D_Selected ( iV ) ) ) cycle
      call AllocateDevice &
             ( S % Value ( :, S % iaSelected ( iV ) ), S % D_Selected ( iV ) )
    end do
  
  end subroutine AllocateDevice_S
  
  
  subroutine UpdateDeviceAll ( S )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
      
    integer ( KDI ) :: &
      iS
      
    do iS = 1, S % nVariables
      call UpdateDevice &
             ( S % Value ( :, S % iaSelected ( iS ) ), S % D_Selected ( iS ) )
    end do
  
  end subroutine UpdateDeviceAll


  subroutine UpdateDeviceSingle ( S, iV )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iV
    
    integer ( KDI ) :: &
      iS
    
    call Search ( S % iaSelected, iV, iS )
    call UpdateDevice ( S % Value ( :, iV ), S % D_Selected ( iS ) )
  
  end subroutine UpdateDeviceSingle


  subroutine UpdateHostAll ( S )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
      
    integer ( KDI ) :: &
      iS
      
    do iS = 1, S % nVariables
      call UpdateHost &
             ( S % D_Selected ( iS ), S % Value ( :, S % iaSelected ( iS ) ) )
    end do
  
  end subroutine UpdateHostAll
  
  
  subroutine UpdateHostSingle ( S, iV )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iV
    
    integer ( KDI ) :: &
      iS
    
    call Search ( S % iaSelected, iV, iS )
    call UpdateHost ( S % D_Selected ( iS ), S % Value ( :, iV ) )
  
  end subroutine UpdateHostSingle
  
  
  impure elemental subroutine Finalize ( S )

    type ( StorageForm ), intent ( inout ) :: &
      S
      
    integer ( KDI ) :: &
      iV

!-- FIXME: this deallocation in a cloned variable group with no vectors
!          caused trouble with Intel 12.1.2 
!    if ( allocated ( S % VectorIndices ) ) deallocate ( S % VectorIndices )
    if ( allocated ( S % VectorIndices ) .and. S % nVectors > 0 ) &
      deallocate ( S % VectorIndices )

    if ( allocated ( S % Unit ) )       deallocate ( S % Unit )
    
    if ( allocated ( S % D_Selected ) ) then
      do iV = 1, S % nVariables
        call DeallocateDevice ( S % D_Selected ( iV ) )
      end do
      deallocate ( S % D_Selected )
    end if

    if ( allocated ( S % Vector ) )     deallocate ( S % Vector )
    if ( allocated ( S % Variable ) )   deallocate ( S % Variable )

    if ( S % AllocatedValue ) then
      if ( associated ( S % Value ) ) then
        if ( S % Pinned ) then
          call DeallocateHost ( S % Value )
        else
         deallocate ( S % Value )
        end if
      end if
    end if
    nullify ( S % Value )

    if ( allocated ( S % iaSelected ) ) deallocate ( S % iaSelected )
    if ( allocated ( S % lVector ) )    deallocate ( S % lVector )
    if ( allocated ( S % lVariable ) )  deallocate ( S % lVariable )

  end subroutine Finalize
  

  subroutine InitializeOptionalMembers &
               ( S, VectorIndicesOption, UnitOption, VectorOption, &
                 VariableOption, NameOption )
                 
    class ( StorageForm ), intent ( inout ) :: &
      S
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
      S % nVectors = size ( VectorIndicesOption )
      allocate ( S % VectorIndices ( S % nVectors ) )
      do iVctr = 1, S % nVectors
        call S % VectorIndices ( iVctr ) % Initialize &
               ( VectorIndicesOption ( iVctr ) )
      end do
    end if
    
    if ( .not. allocated ( S % VectorIndices ) ) &
      allocate ( S % VectorIndices ( 0 ) )

    if ( present ( NameOption ) ) then
      S % lName = len_trim ( NameOption )
      S % Name = NameOption
    end if 
    
    if ( .not. allocated ( S % lVariable ) ) then
      allocate ( S % lVariable ( size ( S % Value, dim = 2 ) ) )
      S % lVariable = 0
    end if
    
    if ( .not. allocated ( S % Variable ) ) then
      allocate ( S % Variable ( size ( S % Value, dim = 2 ) ) )
      S % Variable = ''
    end if
    
    if ( present ( VariableOption ) ) then
      do iS = 1, S % nVariables
        iVrbl = S % iaSelected ( iS )
        S % Variable ( iVrbl ) = VariableOption ( iS )
        S % lVariable ( iVrbl ) = len_trim ( VariableOption ( iS ) )
      end do
    end if
      
    if ( present ( VectorOption ) ) then
      allocate ( S % lVector ( size ( VectorOption ) ) )
      S % lVector = len_trim ( VectorOption )
      allocate ( S % Vector ( size ( VectorOption ) ) )
      S % Vector = VectorOption
    end if
    
    if ( .not. allocated ( S % lVector ) ) &
      allocate ( S % lVector ( 0 ) )
    if ( .not. allocated ( S % Vector ) ) &
      allocate ( S % Vector ( 0 ) )

    if ( .not. allocated ( S % Unit ) ) &
      allocate ( S % Unit ( size ( S % Value, dim = 2 ) ) )
    if ( present ( UnitOption ) ) S % Unit = UnitOption
    
  end subroutine InitializeOptionalMembers


end module Storage_Form
