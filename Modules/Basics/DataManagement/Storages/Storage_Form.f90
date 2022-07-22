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
    type ( c_ptr ), private :: &
      D_Value = c_null_ptr  !-- Device pointer to Value
    type ( c_ptr ), dimension ( : ), allocatable, private :: &
      D_Selected     !-- Device pointer for Selected Value
    integer ( KDI ) :: &
      nValues     = 0, &
      nVariables  = 0, &
      nVectors    = 0, &
      lName       = 0, &
      ErrorDevice = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected, &
      lVariable, &
      lVector
    real ( KDR ), dimension ( :, : ), pointer, contiguous :: &
      Value => null (  )
    logical ( KDL ) :: &
      AllocatedValue  = .false., &
      AllocatedDevice = .false., &
      ClearRequested  = .false., &
      Pinned          = .false.
    character ( LDF ) :: &
      Name = ''
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Unit
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( StorageForm ), pointer, private :: &
      Primary => null ( )
  contains
    !-- FIXME: Changed "private" to "public" since marking this private
    !          cause problem in the overriding subroutine in extension
    !          with CCE-8.x (possibly related to bug CCS # 121803)
    procedure, public, pass :: &
      InitializeAllocate
    !procedure, public, pass, non_overridable :: &
    !  InitializeAssociate
    procedure, public, pass :: &
      InitializeClone
    generic :: &
      Initialize => InitializeAllocate, InitializeClone
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_S
    procedure, public, pass :: &
      Clear => Clear_S
    procedure, public, pass :: &
      ReassociateHost
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
    procedure, public, pass :: &
      ShowAssociation
    procedure, private, pass :: &
      AssociateHost_S
    procedure, private, pass :: &
      DisassociateHost_S
    final :: &
      Finalize
  end type StorageForm

    private :: &
      InitializeOptionalMembers, &
      AdjustNonPrimary_D_Selected

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
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
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
    
    if ( present ( ClearOption ) ) S % ClearRequested = ClearOption
    if ( S % ClearRequested ) call Clear ( S % Value )  
    
    allocate ( S % D_Selected ( S % nVariables ) )
    S % D_Selected = c_null_ptr
    
    call InitializeOptionalMembers &
           ( S, VectorIndicesOption, UnitOption, VectorOption, &
             VariableOption, NameOption )
  
  end subroutine InitializeAllocate

  
!  subroutine InitializeAssociate &
!               ( S, Value, VectorIndicesOption, UnitOption, VectorOption, &
!                 VariableOption, NameOption, iaSelectedOption )
!    
!    class ( StorageForm ), intent ( inout ) :: &
!      S
!    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
!      Value
!    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
!      VectorIndicesOption
!    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
!      UnitOption
!    character ( * ), dimension ( : ), intent ( in ), optional :: &
!      VectorOption, &
!      VariableOption
!    character ( * ), intent ( in ), optional :: &
!      NameOption
!    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
!      iaSelectedOption
!    
!    integer ( KDI ) :: &
!      iVrbl
!
!    S % nValues = size ( Value, dim = 1 )
!    
!    if ( present ( iaSelectedOption ) ) then
!      S % nVariables = size ( iaSelectedOption )
!    else
!      S % nVariables = size ( Value, dim = 2 )
!    end if
!    
!    allocate ( S % iaSelected ( S % nVariables ) )
!    if ( present ( iaSelectedOption ) ) then
!      S % iaSelected = iaSelectedOption
!    else
!      S % iaSelected = [ ( iVrbl, iVrbl = 1, S % nVariables ) ]
!    end if
!
!    S % Value => Value
!    S % AllocatedValue = .false.
!    
!    allocate ( S % D_Selected ( S % nVariables ) )
!    S % D_Selected = c_null_ptr
!    
!    call InitializeOptionalMembers &
!           ( S, VectorIndicesOption, UnitOption, VectorOption, &
!             VariableOption, NameOption )
!  
!  end subroutine InitializeAssociate
  
  
  subroutine InitializeClone (  &
               S_Target, S_Source, VectorIndicesOption, VectorOption, &
               NameOption, iaSelectedOption )

    class ( StorageForm ), intent ( inout ) :: &
      S_Target
    class ( StorageForm ), intent ( in ), target :: &
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
    S_Target % AllocatedValue   = .false.
    S_Target % AllocatedDevice  = S_Source % AllocatedDevice
    S_Target % Pinned           = S_Source % Pinned
    
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
    
    if ( S_Source % AllocatedValue ) then
      S_Target % Primary => S_Source
    else
      S_Target % Primary => S_Source % Primary
    end if
    
    allocate ( S_Target % D_Selected ( S_Target % nVariables ) )
    if ( .not. present ( iaSelectedOption ) ) then
      S_Target % D_Selected = S_Source % D_Selected
    else
      call AdjustNonPrimary_D_Selected ( S_Target )
    end if
    
    call InitializeOptionalMembers &
           ( S_Target, VectorIndicesOption = VectorIndicesOption, &
             VectorOption = VectorOption, NameOption = NameOption )

  end subroutine InitializeClone
  
  
  subroutine AllocateDevice_S ( S, AssociateVariablesOption )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    logical ( KDL ), intent ( in ), optional :: &
      AssociateVariablesOption
      
    logical ( KDL ) :: &
      AssociateVariables
    
    if ( S % AllocatedValue ) then
      
      AssociateVariables = .true.
      if ( present ( AssociateVariablesOption ) ) &
        AssociateVariables = AssociateVariablesOption

      call AllocateDevice ( S % nValues * S % nVariables, S % D_Value )
      
      if ( .not. c_associated ( S % D_Value ) ) &
        return
      
      call S % AssociateHost_S ( AssociateVariables )
      
      S % AllocatedDevice = .true.
      
      if ( S % ClearRequested ) &
        call Clear ( S % Value, UseDeviceOption = .true. )
      
    else
      
      call AdjustNonPrimary_D_Selected ( S )
      
    end if !-- S % AllocatedValue
  
  end subroutine AllocateDevice_S
  
  
  subroutine ReassociateHost ( S, AssociateVariablesOption )
  
    class ( StorageForm ), intent ( inout ), target :: &
      S
    logical ( KDL ), intent ( in ), optional :: &
      AssociateVariablesOption
    
    logical ( KDL ) :: &
      AssociateVariables
    class ( StorageForm ), pointer :: &
      S_Primary
      
    if ( .not. S % AllocatedDevice ) &
      return
    
    AssociateVariables = .true.
    if ( present ( AssociateVariablesOption ) ) &
      AssociateVariables = AssociateVariablesOption
      
    if ( S % AllocatedValue ) then
      call S % DisassociateHost_S ( )
      call S % AssociateHost_S ( AssociateVariables )
    else
      S_Primary => S % Primary
      call S_Primary % DisassociateHost_S ( )
      call S_Primary % AssociateHost_S ( AssociateVariables )
      call AdjustNonPrimary_D_Selected ( S )
    end if
    
  end subroutine ReassociateHost


  subroutine Clear_S ( S )

    class ( StorageForm ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSelected

    do iS = 1, S % nVariables
      call Clear ( S % Value ( :, S % iaSelected ( iS ) ), &
                   UseDeviceOption = S % AllocatedDevice )
    end do !-- iS

  end subroutine Clear_S


  subroutine UpdateDeviceAll ( S )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    integer ( KDI ) :: &
      iS
      
    if ( .not. S % AllocatedDevice ) &
      return
      
    if ( S % AllocatedValue ) then
      call UpdateDevice &
             ( S % Value, S % D_Selected ( 1 ), &
               ErrorOption = S % ErrorDevice )
    else
      do iS = 1, S % nVariables
        call UpdateDevice &
               ( S % Value ( :, S % iaSelected ( iS ) ), &
                 S % D_Selected ( iS ), ErrorOption = S % ErrorDevice )
      end do
    end if
  
  end subroutine UpdateDeviceAll


  subroutine UpdateDeviceSingle ( S, iV )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iV
    
    integer ( KDI ) :: &
      iS
      
    if ( .not. S % AllocatedDevice ) &
      return
    
    call Search ( S % iaSelected, iV, iS )
    call UpdateDevice &
           ( S % Value ( :, iV ), S % D_Selected ( iS ), &
             ErrorOption = S % ErrorDevice )
  
  end subroutine UpdateDeviceSingle


  subroutine UpdateHostAll ( S )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
      
    integer ( KDI ) :: &
      iS
    
    if ( .not. S % AllocatedDevice ) &
      return
      
    if ( S % AllocatedValue ) then
      call UpdateHost ( S % D_Selected ( 1 ), S % Value )
    else
      do iS = 1, S % nVariables
        call UpdateHost &
               ( S % D_Selected ( iS ), S % Value ( :, S % iaSelected ( iS ) ) )
      end do
    end if
    
  end subroutine UpdateHostAll
  
  
  subroutine UpdateHostSingle ( S, iV )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iV
    
    integer ( KDI ) :: &
      iS
      
    if ( .not. S % AllocatedDevice ) &
      return
    
    call Search ( S % iaSelected, iV, iS )
    call UpdateHost ( S % D_Selected ( iS ), S % Value ( :, iV ) )
  
  end subroutine UpdateHostSingle
  
  
  subroutine ShowAssociation ( S, Description )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ) :: &
      Description
    
    integer ( KDI ) :: &
      iV
    integer ( KBI ) :: &
      D_Address, &
      H_Address
    character ( LDN ) :: &
      IndexLabel
    character ( LDB ) :: &
      D_AddressStr, &
      H_AddressStr
      
    print '(a35)', trim ( Description )

    print '(a20, a20, a20)', &
      'Index   ', &
      'Host Address        ', &
      'Device Address      '

    do iV = 1, S % nVariables 
      write ( IndexLabel, fmt = '( i7 )' ) iV
      
      D_Address = transfer ( S % D_Selected ( iV ), 1_KBI )
      write ( D_AddressStr, fmt = ' ( z64 ) ' ) D_Address
      D_AddressStr = '0x' //  adjustl ( D_AddressStr )
      
      H_Address &
        = transfer ( c_loc ( S % Value ( :, S % iaSelected ( iV ) ) ), &
                     1_KBI )
      write ( H_AddressStr, fmt = ' ( z64 ) ' ) H_Address
      H_AddressStr = '0x' //  adjustl ( H_AddressStr )
      
      print &
        '(a20, a20, a20)', &
        '( ' // trim ( adjustl ( IndexLabel ) ) // ' ) = ', &
        H_AddressStr, D_AddressStr
    end do
  
  end subroutine ShowAssociation
  
  
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
      if ( S % AllocatedValue .and. S % AllocatedDevice ) then
        call DeallocateDevice ( S % D_Selected ( 1 ) )
      end if
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
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
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
  
  
  subroutine AdjustNonPrimary_D_Selected ( S )
  
    class ( StorageForm ), intent ( inout ) :: &
      S
    
    integer ( KDI ) :: &
      iS_T, iS_S, &
      iV
    
    if ( .not. associated ( S % Primary ) ) &
      return
    
    S % AllocatedDevice = S % Primary % AllocatedDevice
    
    if ( .not. S % AllocatedDevice ) &
      return 
    
    S % D_Selected = c_null_ptr
    do iS_T = 1, S % nVariables
      iV = S % iaSelected ( iS_T )
      do iS_S = 1, S % Primary % nVariables
        if ( iV == S % Primary % iaSelected ( iS_S ) ) then
          S % D_Selected ( iS_T ) &
            = S % Primary % D_Selected ( iS_S )
          exit
        end if
      end do
    end do
  
  end subroutine AdjustNonPrimary_D_Selected
  
  
  subroutine AssociateHost_S ( S, AssociateVariables )
    
    class ( StorageForm ), intent ( inout ) :: &
      S
    logical ( KDL ), intent ( in ) :: &
      AssociateVariables
      
    integer ( KDI ) :: &
      iV
    real ( KDR ), dimension ( : ), pointer :: &
      Variable
    real ( KDR ), dimension ( :, : ), pointer :: &
      D_Scratch
      
    !-- AssociateHost_S assumes S is primary
    
    if ( AssociateVariables ) then
      
      !-- Associate individual variables (columns of S % Value ) on host
      !   to locations on device.

      call c_f_pointer &
             ( S % D_Value, D_Scratch, &
               [ S % nValues, S % nVariables ] )

      do iV = 1, S % nVariables
        S % D_Selected ( iV ) = c_loc ( D_Scratch ( :, iV ) )
        Variable => S % Value ( :, iV )
        call AssociateHost ( S % D_Selected ( iV ), Variable )
      end do

    else

      !-- Associate S % Value (as an entire block) on host to the head
      !   location on the device.
      
      S % D_Selected ( 1 ) = S % D_Value
      call AssociateHost ( S % D_Selected ( 1 ), S % Value )

    end if !-- AssociateVariables

  end subroutine AssociateHost_S
  
  
  subroutine DisassociateHost_S ( S )
    
    class ( StorageForm ), intent ( inout ) :: &
      S
    
    integer ( KDI ) :: &
      iV
    real ( KDR ), dimension ( : ), pointer :: &
      Variable
    
    do iV = 1, S % nVariables
      if ( c_associated ( S % D_Selected ( iV ) ) ) then
        Variable => S % Value ( :, iV )
        call DisassociateHost ( Variable )
        S % D_Selected ( iV ) = c_null_ptr
      end if 
    end do
        
  end subroutine DisassociateHost_S


end module Storage_Form
