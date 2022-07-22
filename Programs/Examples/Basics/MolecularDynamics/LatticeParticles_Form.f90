module LatticeParticles_Form

  use Basics
  use NormalRandomNumber_Function
  use DistributedParticles_Form

  implicit none
  private

  type, public, extends ( DistributedParticlesForm ) :: LatticeParticlesForm
    integer ( KDI ) :: &  !-- Extensive
      SQUARE_DISPLACEMENT
    integer ( KDI ) :: &
      N_EXTENSIVE_LP = 1
    integer ( KDI ) :: &  !-- Intensive
      TEMPERATURE, &
      COMPRESSIBILITY_FACTOR, &
      MEAN_SQUARE_DISPLACEMENT
    integer ( KDI ) :: &
      N_INTENSIVE_LP = 3
    real ( KDR ) :: &
      NumberDensity, &
      TargetTemperature
    real ( KDR ), dimension ( :, : ), allocatable :: &
      InitialPosition
    type ( QuantityForm ) :: &
      TemperatureUnit
  contains
    procedure, private, pass :: &
      Initialize_LP
    generic, public :: &
      Initialize => Initialize_LP
    procedure, public, pass :: &
      EnforceTargetTemperature
    procedure, public, pass :: &
      TemperatureExpression
    procedure, public, pass :: &
      CompressibilityFactorExpression
    procedure, public, pass :: &
      MySquareDisplacement
  end type LatticeParticlesForm

    private :: &
      SetParticles_FCC

contains


  subroutine Initialize_LP ( LP, C, TimeUnitOption )

    class ( LatticeParticlesForm ), intent ( inout ) :: &
      LP
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    type ( QuantityForm ), intent ( in ), optional :: &
      TimeUnitOption

    integer ( KDI ) :: &
      nUnitCellsRoot, &
      nParticlesPerUnitCell, &
      nParticles
    real ( KDR ) :: &
      BoxLength, &
      ParticleMass
    type ( QuantityForm ) :: &
      NumberDensityUnit, &
      LengthUnit
    character ( LDL ) :: &
      LatticeType

    LatticeType = 'FCC'
    call PROGRAM_HEADER % GetParameter ( LatticeType, 'LatticeType' )

    select case ( trim ( LatticeType ) )
    case ( 'FCC' )       
      LP % Type = 'FCC_LatticeParticles'
      nParticlesPerUnitCell = 4
    case default
      call Show ( 'Lattice type not implemented', CONSOLE % ERROR )
      call Show ( LatticeType, 'LatticeType', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    call Show ( nParticlesPerUnitCell, 'nParticlesPerUnitCell' )

    nUnitCellsRoot = 4
    call PROGRAM_HEADER % GetParameter &
           ( nUnitCellsRoot, 'nUnitCellsRoot' )

    associate ( nUnitCells => nUnitCellsRoot ** 3 ) 
    if ( mod ( nUnitCells, C % Size ) /= 0 ) then
      call Show ( 'The number of MPI processes must divide evenly ' &
                  // 'into nUnitCells', CONSOLE % ERROR )
      call Show ( C % Size, 'nProcesses', CONSOLE % ERROR )
      call Show ( mod ( nUnitCells, C % Size ), &
                  'mod ( nUnitCells, nProcesses )', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if
    end associate !-- nUnitCells

    nParticles = nParticlesPerUnitCell * ( nUnitCellsRoot ** 3 )
    call Show ( 'nParticles = nParticlesPerUnitCell ' &
                // '* ( nUnitCellsRoot ** 3 )', CONSOLE % INFO_2 )
    call Show ( nParticles, 'nParticles', CONSOLE % INFO_2 )

    LP % NumberDensity = nParticles
    call PROGRAM_HEADER % GetParameter &
           ( LP % NumberDensity, 'NumberDensity', &
             InputUnitOption = NumberDensityUnit )

    if ( NumberDensityUnit == UNIT % NUMBER_DENSITY_ANGSTROM ) then
      LengthUnit = UNIT % ANGSTROM
    end if

    associate ( OneThird => 1.0_KDR / 3.0_KDR )
    BoxLength = ( nParticles / LP % NumberDensity ) ** OneThird
    call Show ( BoxLength, LengthUnit, 'BoxLength', CONSOLE % INFO_2 )
    end associate !-- OneThird

    ParticleMass = 1.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( ParticleMass, 'ParticleMass', InputUnitOption = LP % MassUnit )

    call LP % DistributedParticlesForm % Initialize &
           ( C, BoxLength, ParticleMass, nParticles, &
             LengthUnitOption = LengthUnit, TimeUnitOption = TimeUnitOption )
    LP % IsPeriodic = .true.

    LP % TargetTemperature = 1.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( LP % TargetTemperature, 'TargetTemperature', &
             InputUnitOption = LP % TemperatureUnit )

    call InitializeRandomSeed ( C )
    select case ( trim ( LatticeType ) )
    case ( 'FCC' )       
      call SetParticles_FCC ( LP, nUnitCellsRoot )
    end select

    associate ( MP => LP % MyParticles )
    allocate ( LP % InitialPosition ( nParticles, 3 ) )
    LP % InitialPosition = MP % Value ( :, MP % POSITION )
    end associate !-- MP

  end subroutine Initialize_LP


  subroutine EnforceTargetTemperature ( LP, T )

    class ( LatticeParticlesForm ), intent ( inout ) :: &
      LP
    real ( KDR ), intent ( in ) :: &
      T

    call Show ( 'Enforcing TargetTemperature', CONSOLE % INFO_2 )
    call Show ( sqrt ( LP % TargetTemperature / T ), 'Velocity rescaling', &
                CONSOLE % INFO_2 )

    associate ( MP => LP % MyParticles )
    associate ( V => MP % Value &
                      ( :, MP % VELOCITY ( 1 ) : MP % VELOCITY ( 3 ) ) )
    V = V * sqrt ( LP % TargetTemperature / T )
    end associate !-- V
    end associate !-- MP

  end subroutine EnforceTargetTemperature


  function TemperatureExpression ( LP, KE ) result ( T )

    class ( LatticeParticlesForm ), intent ( in ) :: &
      LP
    real ( KDR ), intent ( in ) :: &
      KE
    real ( KDR ) :: &
      T

    T = KE / ( 1.5_KDR * ( LP % nParticles - 1 ) * CONSTANT % BOLTZMANN )

  end function TemperatureExpression


  function CompressibilityFactorExpression ( LP, Virial, T ) result ( CF )

    class ( LatticeParticlesForm ), intent ( in ) :: &
      LP
    real ( KDR ), intent ( in ) :: &
      Virial, &
      T
    real ( KDR ) :: &
      CF 

    CF = 1.0_KDR &
         + Virial / ( 3.0_KDR * LP % nParticles * CONSTANT % BOLTZMANN * T )

  end function CompressibilityFactorExpression


  function MySquareDisplacement ( LP ) result ( MSD )

    class ( LatticeParticlesForm ), intent ( in ) :: &
      LP
    real ( KDR ) :: &
      MSD

    integer ( KDI ) :: &
      iP  !-- iParticle
    real ( KDR ), dimension ( LP % nMyParticles, 3 ) :: &
      D

    associate ( MP => LP % MyParticles )
    associate &
      ( R   => MP % Value ( :, MP % POSITION ( 1 ) : MP % POSITION ( 3 ) ), &
        R_0 => LP % InitialPosition )

    D = R - R_0
    MSD = 0.0_KDR
    do iP = 1, LP % nMyParticles
      MSD = MSD + dot_product ( D ( iP, : ), D ( iP, : ) )
    end do

    end associate !-- R, R_0
    end associate !-- MP

  end function MySquareDisplacement


  subroutine SetParticles_FCC ( LP, nUnitCellsRoot )

    class ( LatticeParticlesForm ), intent ( inout ) :: &
      LP
    integer ( KDI ), intent ( in ) :: &
      nUnitCellsRoot

    integer ( KDI ) :: &
      iC, jC, kC, &  !-- iCell, etc.
      iCP, &  !-- iCellParticle
      iD, & !-- iDimension
      iParticle, &
      iMyParticle
    real ( KDR ), dimension ( 3 ) :: &
      CellOrigin, &
      AverageSpeed
    type ( CollectiveOperation_R_Form ) :: &
      CO

    call Show ( 'Setting particle initial conditions', CONSOLE % INFO_3 )

    associate &
      ( Rank => LP % Communicator % Rank, &
        LatticeParameter => LP % BoxLength / nUnitCellsRoot, &
        MP => LP % MyParticles, &
        Variance => CONSTANT % BOLTZMANN * LP % TargetTemperature &
                    / LP % ParticleMass )

    call Show ( LatticeParameter, LP % LengthUnit, 'LatticeParameter', &
                CONSOLE % INFO_3 )

    iParticle = 0
    iMyParticle = 0
    do kC = 1, nUnitCellsRoot
      do jC = 1, nUnitCellsRoot
        do iC = 1, nUnitCellsRoot

          iParticle = iParticle + 4
          if ( iParticle <= Rank * LP % nMyParticles &
               .or. iParticle > ( Rank + 1 ) * LP % nMyParticles ) cycle

          CellOrigin &
            = - 0.5_KDR * LP % BoxLength &
              + [ iC - 1, jC - 1, kC - 1 ] * LatticeParameter &
              + 0.25_KDR * LatticeParameter

          MP % Value ( iMyParticle + 1, MP % POSITION ) &
            = CellOrigin
          MP % Value ( iMyParticle + 2, MP % POSITION ) &
            = CellOrigin  +  0.5_KDR * LatticeParameter * [ 1, 1, 0 ]
          MP % Value ( iMyParticle + 3, MP % POSITION ) &
            = CellOrigin  +  0.5_KDR * LatticeParameter * [ 1, 0, 1 ]
          MP % Value ( iMyParticle + 4, MP % POSITION ) &
            = CellOrigin  +  0.5_KDR * LatticeParameter * [ 0, 1, 1 ]

          do iD = 1, 3
            do iCP = 1, 4
              MP % Value ( iMyParticle + iCP, MP % VELOCITY ( iD ) ) &
                = NormalRandomNumber ( Variance )
            end do
          end do

          iMyParticle = iMyParticle + 4

        end do
      end do
    end do

    !-- Positions in box

    MP % Value ( :, MP % POSITION_BOX ) = MP % Value ( :, MP % POSITION )

    !-- Ensure zero total momentum

    call CO % Initialize ( LP % Communicator, [ 3 ], [ 3 ] )

    do iD = 1, 3
      associate ( V_Dim => MP % Value ( :, MP % VELOCITY ( iD ) ) )
      CO % Outgoing % Value ( iD ) = sum ( V_Dim )
      end associate !-- V_Dim
    end do

    call CO % Reduce ( REDUCTION % SUM )
    AverageSpeed = CO % Incoming % Value / LP % nParticles

    do iD = 1, 3
      associate ( V_Dim => MP % Value ( :, MP % VELOCITY ( iD ) ) )
      V_Dim = V_Dim - AverageSpeed ( iD )
      end associate !-- V_Dim
    end do

    end associate !-- Rank, etc.

  end subroutine SetParticles_FCC


end module LatticeParticles_Form
