!-- ParameterStreamForm is a class to read a parameter file with multiple
!   label=value lines.

module ParameterStream_Form
  
  use Specifiers
  use Display
  use DelayFileAccess_Command
  use FindParameter_Command

  implicit none
  private
  
    integer ( KDI ), parameter, private :: &
      HANDLE_UNINITIALIZED = - huge ( 1_KDI )

  type, public :: ParameterStreamForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      Handle       = HANDLE_UNINITIALIZED, &
      ProcessRank  = 0
    character ( LDF ) :: &
      Filename = '', &
      Path     = ''
    character ( LDB ), dimension ( : ), allocatable :: &
      Buffer
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Read_0D_Integer
    procedure, private, pass :: &
      Read_0D_Real
    procedure, private, pass :: &
      Read_0D_Quantity
    procedure, private, pass :: &
      Read_0D_Logical
    procedure, private, pass :: &
      Read_0D_Character
    procedure, private, pass :: &
      Read_1D_Integer
    procedure, private, pass :: &
      Read_1D_Real
    procedure, private, pass :: &
      Read_1D_Quantity
    procedure, private, pass :: &
      Read_1D_Logical
    procedure, private, pass :: &
      Read_1D_Character
    generic :: &
      Read &
        => Read_0D_Integer, Read_0D_Real, Read_0D_Quantity, &
           Read_0D_Logical, Read_0D_Character, &
           Read_1D_Integer, Read_1D_Real, Read_1D_Quantity, &
           Read_1D_Logical, Read_1D_Character
    final :: &
      Finalize
  end type ParameterStreamForm 

    character ( 7 ), private, parameter :: &
      FORMAT_BUFFER = '(a1023)'  !-- must correspond to 
                                 !   LDB = LEN_DEFAULT % BUFFER

contains
  
  
  subroutine Initialize &
               ( PS, Filename, ProcessRank, PathOption, IgnorabilityOption, &
                 nLinesOption )
    
    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    character ( * ), intent ( in ) :: &
      Filename
    integer ( KDI ), intent ( in ) :: &
      ProcessRank
    character ( * ), intent ( in ), optional :: &
      PathOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nLinesOption
    
    integer ( KDI ) :: &
      iL, &  !-- iLine
      nLines, &
      Status
    character ( LDB ) :: &
      ErrorMessage, &
      Dummy

    PS % IGNORABILITY = CONSOLE % INFO_2
    if ( present ( IgnorabilityOption ) ) &
      PS % IGNORABILITY = IgnorabilityOption

    PS % ProcessRank = ProcessRank

    PS % Filename = trim ( Filename )
    
    PS % Path = '../Parameters/'
    if ( present ( PathOption ) ) PS % Path = trim ( PathOption )
    
    call Show ( 'Opening a parameter file', PS % IGNORABILITY )
    call Show ( PS % Filename, 'Name', PS % IGNORABILITY )
    
    call DelayFileAccess ( ProcessRank )
    open &
      ( newunit = PS % Handle, &
        file = trim ( PS % Path ) // trim ( PS % Filename ), &
        action = 'read', status = 'old', iostat = Status, &
        iomsg = ErrorMessage )

    if ( Status /= 0 ) then
      PS % Handle = HANDLE_UNINITIALIZED
      allocate ( PS % Buffer ( 0 ) )
      call Show ( ErrorMessage, PS % IGNORABILITY )
      return
    end if

    if ( present ( nLinesOption ) ) then
      nLines = nLinesOption
    else
      nLines = 0
      do 
        read ( unit = PS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) &
          Dummy
        if ( Status /= 0 ) exit
        nLines = nLines + 1
      end do
    end if

    allocate ( PS % Buffer ( nLines ) )

    rewind PS % Handle
    do iL = 1, nLines
      read ( unit = PS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) &
        PS % Buffer ( iL )      
    end do

  end subroutine Initialize
  
  
  subroutine Read_0D_Integer &
               ( PS, Value, Name, IgnorabilityOption, SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    integer ( KDI ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, IgnorabilityOption, &
             SuccessOption )

  end subroutine Read_0D_Integer


  subroutine Read_0D_Real &
               ( PS, Value, Name, InputUnitOption, IgnorabilityOption, &
                 SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    real ( KDR ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( QuantityForm ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, InputUnitOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_0D_Real


  subroutine Read_0D_Quantity &
               ( PS, Value, Name, InputUnitOption, IgnorabilityOption, &
                 ConvertOption, SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    type ( QuantityForm ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( QuantityForm ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( in ), optional :: &
      ConvertOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, InputUnitOption, &
             IgnorabilityOption, ConvertOption, SuccessOption )

  end subroutine Read_0D_Quantity


  subroutine Read_0D_Logical &
               ( PS, Value, Name, IgnorabilityOption, SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    logical ( KDL ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, IgnorabilityOption, &
             SuccessOption )

  end subroutine Read_0D_Logical


  subroutine Read_0D_Character &
               ( PS, Value, Name, IgnorabilityOption, SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    character ( * ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, IgnorabilityOption, &
             SuccessOption )

  end subroutine Read_0D_Character


  subroutine Read_1D_Integer &
               ( PS, Value, Name, nValuesOption, IgnorabilityOption, &
                 SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, nValuesOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Integer


  subroutine Read_1D_Real &
               ( PS, Value, Name, InputUnitOption, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( QuantityForm ), dimension ( : ), intent ( inout ), &
      optional :: &
        InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, InputUnitOption, &
             nValuesOption, IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Real


  subroutine Read_1D_Quantity &
               ( PS, Value, Name, InputUnitOption, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    type ( QuantityForm ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( QuantityForm ), dimension ( : ), intent ( inout ), &
      optional :: &
        InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, InputUnitOption, &
             nValuesOption, IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Quantity


  subroutine Read_1D_Logical &
               ( PS, Value, Name, nValuesOption, IgnorabilityOption, &
                 SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    logical ( KDL ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, nValuesOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Logical


  subroutine Read_1D_Character &
               ( PS, Value, Name, nValuesOption, IgnorabilityOption, &
                 SuccessOption )

    class ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    character ( * ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, PS % Buffer, PS % Filename, Name, nValuesOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Character


  impure elemental subroutine Finalize ( PS )
  
    type ( ParameterStreamForm ), intent ( inout ) :: &
      PS
    
    if ( PS % Handle == HANDLE_UNINITIALIZED ) return
    
    call Show ( 'Closing a parameter file', PS % IGNORABILITY )
    call Show ( PS % Filename, 'Name', PS % IGNORABILITY )

    if ( allocated ( PS % Buffer ) ) deallocate ( PS % Buffer )
    
    call DelayFileAccess ( PS % ProcessRank )
    close ( unit = PS % Handle )
    
  end subroutine Finalize
  
  
end module ParameterStream_Form
