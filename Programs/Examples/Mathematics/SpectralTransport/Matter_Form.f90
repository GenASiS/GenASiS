module Matter_Form

  use Basics

  implicit none
  private

  type, public, extends ( VariableGroupForm ) :: MatterForm
    integer ( KDI ) :: &
      TEMPERATURE = 1, &
      CHEMICAL_POTENTIAL = 2
  contains
    procedure, public, pass :: &
      InitializeAllocate_M
    generic, public :: &
      Initialize => InitializeAllocate_M
  end type MatterForm

contains


  subroutine InitializeAllocate_M &
               ( M, TemperatureUnit, ChemicalPotentialUnit, nValues, &
                 NameOption )

    class ( MatterForm ), intent ( inout ) :: &
      M
    type ( MeasuredValueForm ), intent ( in ) :: &
      TemperatureUnit, &
      ChemicalPotentialUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOption

    character ( LDF ) :: &
      Name 

    Name = 'Matter'
    if ( present ( NameOption ) ) &
      Name = NameOption

    call M % Initialize &
           ( [ nValues, 2 ], &
             UnitOption = [ TemperatureUnit, ChemicalPotentialUnit ], &
             VariableOption = [ 'Temperature      ', 'ChemicalPotential' ], &
             NameOption = Name, ClearOption = .true. )

  end subroutine InitializeAllocate_M


end module Matter_Form
