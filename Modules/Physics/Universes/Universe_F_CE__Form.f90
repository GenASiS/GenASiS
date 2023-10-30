module Universe_F_CE__Form

  !-- Universe_Fluid_CentralExcision__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_F_C__Form

  implicit none
  private

  type, public, extends ( Universe_F_C_Form ) :: Universe_F_CE_Form
  contains
    procedure, private, pass :: &
      Initialize_F_CE
    generic, public :: &
      Initialize => Initialize_F_CE
    final :: &
      Finalize
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, private, pass :: &
      InitializeAtlas
  end type Universe_F_CE_Form

    private :: &
      Compute_dT_Local

contains


  subroutine Initialize_F_CE &
               ( U, FluidType, GravitationType, Name, &
                 DimensionlessOption, FinishTimeOption, RadiusMaxOption, &
                 RadiusExcisionOption, RadialRatioOption, CentralMassOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_F_CE_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType, &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusExcisionOption, &
      RadialRatioOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( U % Type == '' ) &
      U % Type = 'a Universe_F_CE'

    call U % Initialize_F_C &
           ( FluidType, GravitationType, Name, &
             DimensionlessOption = DimensionlessOption, &
             FinishTimeOption = FinishTimeOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusExcisionOption = RadiusExcisionOption, &
             RadialRatioOption = RadialRatioOption, &
             CentralMassOption = CentralMassOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nWriteOption = nWriteOption ) 

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % Compute_dT_Local  =>  Compute_dT_Local
    end associate !-- I

  end subroutine Initialize_F_CE


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_CE_Form ), intent ( inout ) :: &
      U

  end subroutine Finalize


  subroutine SetBoundaryConditions ( U )

    class ( Universe_F_CE_Form ), intent ( inout ) :: &
      U
    
    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    associate &
      ( F  =>  I % CurrentSet_X )
    call F % SetBoundaryConditionsFace &
           ( [ 'OUTFLOW', 'INFLOW ' ], iC = 1, iD = 1 )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    call F % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end associate !-- F

    end select !-- I

  end subroutine SetBoundaryConditions


  subroutine InitializeAtlas &
               ( U, CommunicatorOption, RadiusMaxOption, RadiusCoreOption, &
                 RadiusExcisionOption, RadialRatioOption, nCellsPolarOption )

    class ( Universe_F_CE_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusExcision, &
      RadialRatio
    type ( CommunicatorForm ), pointer :: &
      Communicator

    if ( present ( CommunicatorOption ) ) then
      Communicator  =>  CommunicatorOption
    else
      Communicator  =>  U % Communicator
    end if

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_CE_Form :: I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_CE_Form )

    if ( U % Dimensionless ) then

      RadiusMax       =  10.0_KDR
      RadiusExcision  =  0.45_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax  =  RadiusMaxOption
      if ( present ( RadiusExcisionOption ) ) &
        RadiusExcision  =  RadiusExcisionOption

      call PS % Initialize &
             ( RadiusMax = RadiusMax, &
               RadiusExcision = RadiusExcision, &
               CommunicatorOption = Communicator, &
               NameOption = 'PositionSpace', &
               DeviceMemoryOption = U % DeviceMemory, &
               nCellsPolarOption = nCellsPolarOption )

    else

      RadiusMax       =  1.0e3_KDR  *  UNIT % KILOMETER
      RadiusExcision  =   40.0_KDR  *  UNIT % KILOMETER
      RadialRatio     =    1.0_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax = RadiusMaxOption
      if ( present ( RadiusExcisionOption ) ) &
        RadiusExcision = RadiusExcisionOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio  =  RadialRatioOption

      call PS % Initialize &
             ( RadiusMax = RadiusMax, &
               RadiusExcision = RadiusExcision, &
               CommunicatorOption = PROGRAM_HEADER % Communicator, &
               NameOption = 'PositionSpace', &
               DeviceMemoryOption = U % DeviceMemory, &
               CoordinateUnitOption = U % Units_F ( 1 ) % Coordinate_PS, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    end if !-- Dimensionless

    !-- Azimuthal average
    if ( PS % Chart_GS_CE % nDimensions  >  2 ) then
      allocate ( Atlas_SCG_CE_Form :: U % PositionSpace_AA )
      select type ( PS_SA  =>  U % PositionSpace_AA )
        class is ( Atlas_SCG_CE_Form )
      call PS_SA % Initialize ( PS, nDimensions = 2 )
      end select !-- PS_SA
    end if !-- nDimensions

    !-- Spherical average
    if ( PS % Chart_GS_CE % nDimensions  >  1 ) then
      allocate ( Atlas_SCG_CE_Form :: U % PositionSpace_SA )
      select type ( PS_SA  =>  U % PositionSpace_SA )
        class is ( Atlas_SCG_CE_Form )
      call PS_SA % Initialize ( PS, nDimensions = 1 )
      end select !-- PS_SA
    end if !-- nDimensions

    end select !-- PS
    end associate !-- I

  end subroutine InitializeAtlas


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )
    select type ( I )
      class is ( Integrator_CS_Form )
    associate &
      ( dT_1  =>  dT_Candidate ( 1 ) )

    !-- Advection step

    if ( U % Coarsen ) then
      call U % Compute_dT_CS_CGS_C ( dT_1, iC, T_Option )
    else !-- .not. Coarsen
      call I % Compute_dT_CS_CGS &
             ( I % EigenspeedSet_X, dT_1, iC, T_Option )
    end if !-- Coarsen
    dT_1  =  I % CourantFactor  *  dT_1
    
    end associate !-- dT_1
    end select !-- I
    end select !-- U

  end subroutine Compute_dT_Local


end module Universe_F_CE__Form
