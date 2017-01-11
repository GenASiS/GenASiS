!-- CommandLineOptions is a class to handle reading and parsing command-line
!   options.

module CommandLineOptions_Form
  
  use VariableManagement
  use FileSystem

  implicit none
  private
  
  type, public :: CommandLineOptionsForm
    integer ( KDI ), private :: &
      nOptions
    character ( LDF ), dimension ( : ), allocatable, private :: &
      Option
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Read_0D_Integer
    procedure, private, pass :: &
      Read_0D_Real
    procedure, private, pass :: &
      Read_0D_MeasuredValue
    procedure, private, pass :: &
      Read_0D_Logical
    procedure, private, pass :: &
      Read_0D_Character
    procedure, private, pass :: &
      Read_1D_Integer
    procedure, private, pass :: &
      Read_1D_Real
    procedure, private, pass :: &
      Read_1D_MeasuredValue
    procedure, private, pass :: &
      Read_1D_Logical
    procedure, private, pass :: &
      Read_1D_Character
    generic :: &
      Read &
        => Read_0D_Integer, Read_0D_Real, Read_0D_MeasuredValue, &
           Read_0D_Logical, Read_0D_Character, &
           Read_1D_Integer, Read_1D_Real, Read_1D_MeasuredValue, &
           Read_1D_Logical, Read_1D_Character
    final :: &
      Finalize
  end type CommandLineOptionsForm
    
contains

  
  subroutine Initialize ( CLO )
    
    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO

    integer ( KDI ) :: &
      iO, &   !-- iOption
      Status
    
    CLO % nOptions = command_argument_count ( ) 
    
    allocate ( CLO % Option ( CLO % nOptions ) )
    
    do iO = 1, CLO % nOptions
      call get_command_argument ( iO, CLO % Option ( iO ), Status )
    end do
    
  end subroutine Initialize
  
  
  subroutine Read_0D_Integer &
               ( CLO, Value, Name, IgnorabilityOption, SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    integer ( KDI ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, IgnorabilityOption, &
             SuccessOption )

  end subroutine Read_0D_Integer


  subroutine Read_0D_Real &
               ( CLO, Value, Name, InputUnitOption, IgnorabilityOption, &
                 SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    real ( KDR ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, InputUnitOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_0D_Real


  subroutine Read_0D_MeasuredValue &
               ( CLO, Value, Name, InputUnitOption, IgnorabilityOption, &
                 ConvertOption, SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    type ( MeasuredValueForm ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( in ), optional :: &
      ConvertOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, InputUnitOption, &
             IgnorabilityOption, ConvertOption, SuccessOption )

  end subroutine Read_0D_MeasuredValue


  subroutine Read_0D_Logical &
               ( CLO, Value, Name, IgnorabilityOption, SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    logical ( KDL ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, IgnorabilityOption, &
             SuccessOption )

  end subroutine Read_0D_Logical


  subroutine Read_0D_Character &
               ( CLO, Value, Name, IgnorabilityOption, SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    character ( * ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, IgnorabilityOption, &
             SuccessOption )

  end subroutine Read_0D_Character


  subroutine Read_1D_Integer &
               ( CLO, Value, Name, nValuesOption, IgnorabilityOption, &
                 SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
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
           ( Value, CLO % Option, 'command line', Name, nValuesOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Integer


  subroutine Read_1D_Real &
               ( CLO, Value, Name, InputUnitOption, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, InputUnitOption, &
             nValuesOption, IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Real


  subroutine Read_1D_MeasuredValue &
               ( CLO, Value, Name, InputUnitOption, nValuesOption, &
                 IgnorabilityOption, SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), &
      optional :: &
        InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    call FindParameter &
           ( Value, CLO % Option, 'command line', Name, InputUnitOption, &
             nValuesOption, IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_MeasuredValue


  subroutine Read_1D_Logical &
               ( CLO, Value, Name, nValuesOption, IgnorabilityOption, &
                 SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
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
           ( Value, CLO % Option, 'command line', Name, nValuesOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Logical


  subroutine Read_1D_Character &
               ( CLO, Value, Name, nValuesOption, IgnorabilityOption, &
                 SuccessOption )

    class ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
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
           ( Value, CLO % Option, 'command line', Name, nValuesOption, &
             IgnorabilityOption, SuccessOption )

  end subroutine Read_1D_Character


  elemental subroutine Finalize ( CLO )
  
    type ( CommandLineOptionsForm ), intent ( inout ) :: &
      CLO
      
    if ( allocated ( CLO % Option ) ) deallocate ( CLO % Option )
      
  end subroutine Finalize  

  
end module CommandLineOptions_Form
