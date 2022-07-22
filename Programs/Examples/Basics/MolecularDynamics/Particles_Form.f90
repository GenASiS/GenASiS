module Particles_Form

  use Basics

  implicit none
  private

  type, public, extends ( StorageForm ) :: ParticlesForm
    integer, dimension ( 3 ) :: &
      POSITION     = [ 1, 2, 3 ], &
      POSITION_BOX = [ 4, 5, 6 ], &
      VELOCITY     = [ 7, 8, 9 ]
  contains
    procedure, public, pass :: &
      InitializeParticles
    generic, public :: &
      Initialize => InitializeParticles
  end type ParticlesForm

contains


  subroutine InitializeParticles &
               ( P, nParticles, LengthUnitOption, TimeUnitOption )
    
    class ( ParticlesForm ), intent ( inout ) :: &
      P
    integer ( KDI ), intent ( in ) :: &
      nParticles
    type ( QuantityForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      TimeUnitOption

    type ( Integer_1D_Form ), dimension ( 1 ) :: &
      VectorIndices
    type ( QuantityForm ), dimension ( 9 ) :: &
      VariableUnit

    call VectorIndices ( 1 ) % Initialize ( P % VELOCITY )

    VariableUnit = UNIT % IDENTITY
    if ( present ( LengthUnitOption ) ) then
      VariableUnit ( P % POSITION ( 1 ) : P % POSITION ( 3 ) ) &
        = LengthUnitOption
      VariableUnit ( P % POSITION_BOX ( 1 ) : P % POSITION_BOX ( 3 ) ) &
        = LengthUnitOption
    end if
    if ( present ( LengthUnitOption ) .and. present ( TimeUnitOption ) ) &
      VariableUnit ( P % VELOCITY ( 1 ) : P % VELOCITY ( 3 ) ) &
        = LengthUnitOption / TimeUnitOption
    
    call P % Initialize &
           ( [ nParticles, 9 ], UnitOption = VariableUnit, &
             VariableOption = [ 'Position_1                     ', &
                                'Position_2                     ', &
                                'Position_3                     ', &
                                'PositionBox_1                  ', &
                                'PositionBox_2                  ', &
                                'PositionBox_3                  ', &
                                'Velocity_1                     ', &
                                'Velocity_2                     ', &
                                'Velocity_3                     ' ], &
             NameOption = 'Particles', VectorIndicesOption = VectorIndices, &
             VectorOption = [ 'Velocity                       ' ] )

  end subroutine InitializeParticles


end module Particles_Form
