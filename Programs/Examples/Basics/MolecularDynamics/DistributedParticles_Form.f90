module DistributedParticles_Form

  use Basics
  use Particles_Form

  implicit none
  private

  type, public :: DistributedParticlesForm
    integer ( KDI ) :: &
      N_EXTENSIVE, &
      N_INTENSIVE
    integer ( KDI ) :: &  !-- Extensive
      KINETIC_ENERGY, &
      POTENTIAL_ENERGY, &
      TOTAL_ENERGY, &
      VIRIAL
    integer ( KDI ), dimension ( 3 ) :: &
      MASS_MOMENT, &
      MOMENTUM, &
      ANGULAR_MOMENTUM
    integer ( KDI ) :: &
      N_EXTENSIVE_DP = 13
    integer ( KDI ), dimension ( 3 ) :: &  !-- Intensive  
      CENTER_OF_MASS
    integer ( KDI ) :: &
      N_INTENSIVE_DP = 3
    integer ( KDI ) :: &
      nParticles, &
      nMyParticles, &
      nCorrelationBins
    integer ( KDI ), private :: &
      iTimer_IO
    real ( KDR ) :: &
      BoxLength, &
      ParticleMass
    real ( KDR ), dimension ( : ), allocatable :: &
      CorrelationBinEdge, &
      MyPairCount, &
      PairCount
    type ( MeasuredValueForm ) :: &
      LengthUnit, &
      MassUnit
    real ( KDR ), dimension ( :, : ), allocatable :: &
      GuestPosition
    logical ( KDL ) :: &
      IsPeriodic = .true.
    character ( LDL ) :: &
      Type = ''
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( MessageIncoming_R_Form ) :: &
      Incoming
    type ( MessageOutgoing_R_Form ) :: &
      Outgoing
    type ( GridImageStreamForm ) :: &
      GridImageStream
    type ( PointGridImageForm ) :: &
      GridImagePosition, &
      GridImagePositionBox, &
      GridImageVelocity
    type ( ParticlesForm ) :: &
      MyParticles
  contains
    procedure, private, pass :: &
      Initialize_DP
    generic, public :: &
      Initialize => Initialize_DP
    procedure, public, pass :: &
      MyKineticEnergy
    procedure, public, pass :: &
      MyPotentialEnergy
    procedure, public, pass :: &
      MyVirial
    procedure, public, pass :: &
      MyMassMoment
    procedure, public, pass :: &
      MyMomentum
    procedure, public, pass :: &
      MyAngularMomentum
    procedure, public, pass :: &
      SetImage
    procedure, public, pass :: &
      Write
    final :: &
      Finalize
  end type DistributedParticlesForm

contains


  subroutine Initialize_DP &
               ( DP, C, BoxLength, ParticleMass, nParticles, &
                 LengthUnitOption, TimeUnitOption )

    class ( DistributedParticlesForm ), intent ( inout ) :: &
      DP
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    real ( KDR ), intent ( in ) :: &
      BoxLength, &
      ParticleMass
    integer ( KDI ), intent ( in ) :: &
      nParticles
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      LengthUnitOption, &
      TimeUnitOption

    integer ( KDI ) :: &
      iE  !-- iEdge

    if ( DP % Type == '' ) DP % Type = 'DistributedParticles'
    call Show ( 'Initializing ' // trim ( DP % Type ), CONSOLE % INFO_2 )

    DP % nParticles = nParticles
    call Show ( DP % nParticles, 'nParticles', CONSOLE % INFO_2 )

    if ( mod ( DP % nParticles, C % Size ) /= 0 ) then
      call Show ( 'The number of MPI processes must divide evenly ' &
                  // 'into nParticles', CONSOLE % ERROR )
      call Show ( C % Size, 'nProcesses', CONSOLE % ERROR )
      call Show ( mod ( DP % nParticles, C % Size ), &
                  'mod ( nParticles, nProcesses )', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    DP % nMyParticles = DP % nParticles / C % Size
    call Show ( DP % nMyParticles, 'nMyParticles' )

    DP % LengthUnit = UNIT % IDENTITY
    if ( present ( LengthUnitOption ) ) DP % LengthUnit = LengthUnitOption

    DP % BoxLength = BoxLength
    call Show ( DP % BoxLength, DP % LengthUnit, 'BoxLength', CONSOLE % INFO_2 )

    allocate ( DP % GuestPosition ( DP % nMyParticles, 3 ) )

    call DP % MyParticles % Initialize &
           ( DP % nMyParticles, LengthUnitOption = DP % LengthUnit, &
             TimeUnitOption = TimeUnitOption )

    DP % ParticleMass = ParticleMass

    DP % Communicator => C
    associate &
      ( SourceRank => modulo ( C % Rank - 1, C % Size ), &
        TargetRank => modulo ( C % Rank + 1, C % Size ), &
        Tag => 1, &
        nSendReceive => 3 * DP % nMyParticles )
    call DP % Incoming % Initialize ( C, Tag, SourceRank, nSendReceive )
    call DP % Outgoing % Initialize ( C, Tag, TargetRank, nSendReceive )
    end associate !-- SourceRank, etc.

    DP % nCorrelationBins = 256
    call PROGRAM_HEADER % GetParameter &
           ( DP % nCorrelationBins, 'nCorrelationBins' )
    allocate ( DP % CorrelationBinEdge ( DP % nCorrelationBins + 1 ) )
    allocate ( DP % MyPairCount ( DP % nCorrelationBins ) )
    allocate ( DP % PairCount ( DP % nCorrelationBins ) )

    associate ( dr => 0.5_KDR * DP % BoxLength / DP % nCorrelationBins )
    DP % CorrelationBinEdge = [ ( iE * dr, iE = 0, DP % nCorrelationBins ) ]
    end associate !-- dr

  end subroutine Initialize_DP

  
  function MyKineticEnergy ( DP ) result ( MKE )

    class ( DistributedParticlesForm ), intent ( in ) :: &
      DP
    real ( KDR ) :: &
      MKE

    integer ( KDI ) :: &
      iP  !-- iParticle

    associate ( MP => DP % MyParticles )
    associate &
      ( M => DP % ParticleMass, &
        V => MP % Value ( :, MP % VELOCITY ( 1 ) : MP % VELOCITY ( 3 ) ) )

    MKE = 0.0_KDR
    do iP = 1, DP % nMyParticles
      MKE = MKE + 0.5_KDR * M * dot_product ( V ( iP, : ), V ( iP, : ) )
    end do

    end associate !-- M, V
    end associate !-- MP

  end function MyKineticEnergy


  function MyPotentialEnergy ( DP, PotentialValue ) result ( MPE )

    class ( DistributedParticlesForm ), intent ( in ) :: &
      DP
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      PotentialValue
    real ( KDR ) :: &
      MPE

    MPE = 0.5_KDR * sum ( PotentialValue )

  end function MyPotentialEnergy


  function MyVirial ( DP, VirialValue ) result ( MV )

    class ( DistributedParticlesForm ), intent ( in ) :: &
      DP
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      VirialValue
    real ( KDR ) :: &
      MV

    MV = 0.5_KDR * sum ( VirialValue )

  end function MyVirial


  function MyMassMoment ( DP ) result ( MMM )

    class ( DistributedParticlesForm ), intent ( in ) :: &
      DP
    real ( KDR ), dimension ( 3 ) :: &
      MMM

    integer ( KDI ) :: &
      iP  !-- iParticle

    associate ( MP => DP % MyParticles )
    associate &
      ( M => DP % ParticleMass, &
        R => MP % Value ( :, MP % POSITION ( 1 ) : MP % POSITION ( 3 ) ) )

    MMM = 0.0_KDR
    do iP = 1, DP % nMyParticles
      MMM = MMM + M * R ( iP, : )
    end do

    end associate !-- M, V
    end associate !-- MP

  end function MyMassMoment


  function MyMomentum ( DP ) result ( MM )

    class ( DistributedParticlesForm ), intent ( in ) :: &
      DP
    real ( KDR ), dimension ( 3 ) :: &
      MM

    integer ( KDI ) :: &
      iP  !-- iParticle

    associate ( MP => DP % MyParticles )
    associate &
      ( M => DP % ParticleMass, &
        V => MP % Value ( :, MP % VELOCITY ( 1 ) : MP % VELOCITY ( 3 ) ) )

    MM = 0.0_KDR
    do iP = 1, DP % nMyParticles
      MM = MM + M * V ( iP, : )
    end do

    end associate !-- M, V
    end associate !-- MP

  end function MyMomentum


  function MyAngularMomentum ( DP, CenterOfMass ) result ( MAM )

    class ( DistributedParticlesForm ), intent ( in ) :: &
      DP
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      CenterOfMass
    real ( KDR ), dimension ( 3 ) :: &
      MAM

    integer ( KDI ) :: &
      iP  !-- iParticle

    associate ( MP => DP % MyParticles )
    associate &
      ( M    => DP % ParticleMass, &
        CM   => CenterOfMass, &
        R    => MP % Value ( :, MP % POSITION ( 1 ) : MP % POSITION ( 3 ) ), &
        V    => MP % Value ( :, MP % VELOCITY ( 1 ) : MP % VELOCITY ( 3 ) ) )

    MAM = 0.0_KDR
    do iP = 1, DP % nMyParticles
      MAM ( 1 ) = MAM ( 1 ) &
                  + M * (    ( R ( iP, 2 ) - CM ( 2 ) ) * V ( iP, 3 )  &
                          -  ( R ( iP, 3 ) - CM ( 3 ) ) * V ( iP, 2 ) ) 
      MAM ( 2 ) = MAM ( 2 ) &
                  + M * (    ( R ( iP, 3 ) - CM ( 3 ) ) * V ( iP, 1 )  &
                          -  ( R ( iP, 1 ) - CM ( 1 ) ) * V ( iP, 3 ) ) 
      MAM ( 3 ) = MAM ( 3 ) &
                  + M * (    ( R ( iP, 1 ) - CM ( 1 ) ) * V ( iP, 2 )  &
                          -  ( R ( iP, 2 ) - CM ( 2 ) ) * V ( iP, 1 ) ) 
    end do

    end associate !-- M, V
    end associate !-- MP

  end function MyAngularMomentum


  subroutine SetImage ( DP, VG_P, VG_V, Name )

    class ( DistributedParticlesForm ), intent ( inout ) :: &
      DP
    class ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      VG_P, &
      VG_V
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iVG  !-- iVariableGroup
    character ( LDF ) :: &
      OutputDirectory

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )    
    call PROGRAM_HEADER % AddTimer ( 'InputOutput', DP % iTimer_IO )
    
    call PROGRAM_HEADER % Timer ( DP % iTimer_IO ) % Start ( )

    associate ( GIS => DP % GridImageStream )
    call GIS % Initialize &
           ( Name, CommunicatorOption = DP % Communicator, &
             WorkingDirectoryOption = OutputDirectory )
    call GIS % Open ( GIS % ACCESS_SET_GRID )

    associate ( GIP => DP % GridImagePosition )
    call GIP % Initialize ( GIS ) 
    do iVG = 1, size ( VG_P )
      call GIP % AddVariableGroup ( VG_P ( iVG ) )
    end do
    end associate !-- GIP

    associate ( GIPB => DP % GridImagePositionBox )
    call GIPB % Initialize ( GIS ) 
    do iVG = 1, size ( VG_P )
      call GIPB % AddVariableGroup ( VG_P ( iVG ) )
    end do
    end associate !-- GIPB

    associate ( GIV => DP % GridImageVelocity )
    call GIV % Initialize ( GIS ) 
    do iVG = 1, size ( VG_V )
      call GIV % AddVariableGroup ( VG_V ( iVG ) )
    end do
    end associate !-- GIV

    call GIS % Close ( )

    end associate !-- GIS

    call PROGRAM_HEADER % Timer ( DP % iTimer_IO ) % Stop ( )

  end subroutine SetImage


  subroutine Write ( DP, TimeOption, CycleNumberOption )

    class ( DistributedParticlesForm ), intent ( inout ) :: &
      DP
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
      
    call PROGRAM_HEADER % Timer ( DP % iTimer_IO ) % Start ( )
    
    associate ( GIS => DP % GridImageStream )
    call GIS % Open ( GIS % ACCESS_CREATE )

    associate &
      ( GIP => DP % GridImagePosition, &
        MP => DP % MyParticles )
    call GIP % SetGrid &
           ( 'PositionSpace', MP % Value ( :, MP % POSITION ), &
             DP % nMyParticles, 3, 0, &
             CoordinateUnitOption = MP % Unit ( MP % POSITION ) )
    call GIP % Write &
           ( TimeOption = TimeOption, CycleNumberOption = CycleNumberOption )
    end associate !-- GIP

    associate &
      ( GIPB => DP % GridImagePositionBox, &
        MP => DP % MyParticles )
    call GIPB % SetGrid &
           ( 'PositionSpaceBox', MP % Value ( :, MP % POSITION_BOX ), &
             DP % nMyParticles, 3, 0, &
             CoordinateUnitOption = MP % Unit ( MP % POSITION_BOX ) )
    call GIPB % Write &
           ( TimeOption = TimeOption, CycleNumberOption = CycleNumberOption )
    end associate !-- GIPB

    associate &
      ( GIV => DP % GridImageVelocity, &
        MP => DP % MyParticles )
    call GIV % SetGrid &
           ( 'VelocitySpace', MP % Value ( :, MP % VELOCITY ), &
             DP % nMyParticles, 3, 0, &
             CoordinateUnitOption = MP % Unit ( MP % VELOCITY ) )
    call GIV % Write &
           ( TimeOption = TimeOption, CycleNumberOption = CycleNumberOption )
    end associate !-- GV

    call GIS % Close ( )
    end associate !-- GIS
    
    call PROGRAM_HEADER % Timer ( DP % iTimer_IO ) % Stop ( )
    
  end subroutine Write


  subroutine Finalize ( DP )

    type ( DistributedParticlesForm ), intent ( inout ) :: &
      DP

    nullify ( DP % Communicator )

    if ( allocated ( DP % GuestPosition ) ) deallocate ( DP % GuestPosition )

  end subroutine Finalize


end module DistributedParticles_Form
