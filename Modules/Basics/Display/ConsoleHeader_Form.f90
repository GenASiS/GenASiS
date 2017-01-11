!-- ConsoleHeaderForm defines message verbosity flags and labels

module ConsoleHeader_Form

  use VariableManagement

  implicit none
  private

  type, public :: ConsoleHeaderForm
    integer ( KDI ) :: &
      ERROR   = 1, &
      WARNING = 2, &
      INFO_1  = 3, &
      INFO_2  = 4, &
      INFO_3  = 5, &
      INFO_4  = 6, &
      INFO_5  = 7, &
      INFO_6  = 8, &
      INFO_7  = 9
    character ( LDL ), dimension ( 9 ) :: &
      LABEL &
        = [ 'ERROR                          ', &
            'WARNING                        ', &
            'INFO_1                         ', &
            'INFO_2                         ', &
            'INFO_3                         ', &
            'INFO_4                         ', &
            'INFO_5                         ', &
            'INFO_6                         ', &
            'INFO_7                         ' ]
  end type ConsoleHeaderForm

end module ConsoleHeader_Form
