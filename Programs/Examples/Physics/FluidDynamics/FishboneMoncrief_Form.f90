module FishboneMoncrief_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_CE_Form ) :: FishboneMoncriefForm
    real ( KDR ), private :: &
      CentralMass, &
      AngularMomentumParameter, &
      RadiusInner, &
      RadiusOuter, &
      SpecificAngularMomentum
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      Show => Show_U
    final :: &
      Finalize
  end type FishboneMoncriefForm

    private :: &
      InitializeUniverse, &
      SetInitial


contains


  subroutine Initialize_H ( U, NameOption )

    class ( FishboneMoncriefForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a FishboneMoncrief'

    Name  =  'FishboneMoncrief'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )

  end subroutine Initialize_H


  subroutine Show_U ( U )

    class ( FishboneMoncriefForm ), intent ( in ) :: &
      U

    call U % Universe_H_Form % Show ( )

    call Show ( U % CentralMass, &
                UNIT % SOLAR_MASS, &
                'M' )
    call Show ( U % AngularMomentumParameter, &
                'Kappa' )
    call Show ( U % RadiusInner, &
                UNIT % KILOMETER,  &
                'R_In' )
    call Show ( U % RadiusOuter, &
                UNIT % KILOMETER,  &
                'R_Out' )
    call Show ( U % SpecificAngularMomentum, &
                UNIT % KILOMETER * UNIT % SPEED_OF_LIGHT, &
                'L' )

  end subroutine Show_U


  impure elemental subroutine Finalize ( FM )
    
    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

  end subroutine Finalize


  subroutine InitializeUniverse ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ), target :: &
      FM
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      G, &
      c, &
      R_Min, &
      T_Finish

     G  =  CONSTANT % GRAVITATIONAL
     c  =  CONSTANT % SPEED_OF_LIGHT

    associate &
      ( M      =>  FM % CentralMass, &
        Kappa  =>  FM % AngularMomentumParameter, &  !-- Between 1 and 2
        R_In   =>  FM % RadiusInner, &
        R_Out  =>  FM % RadiusOuter, &
        L      =>  FM % SpecificAngularMomentum )

    M      =  3.0_KDR * UNIT % SOLAR_MASS
    Kappa  =  1.85_KDR
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )
    call PROGRAM_HEADER % GetParameter ( Kappa, 'AngularMomentumParameter' )

    R_In   =  6.0_KDR * G * M  /  c ** 2
    call PROGRAM_HEADER % GetParameter ( R_In,  'RadiusInner' )

    R_Out  =  R_In * Kappa / ( 2.0_KDR - Kappa )
    L      =  sqrt ( Kappa * G * M * R_In )

    R_Min  =  ( 2.0_KDR / 3.0_KDR ) * R_In
    call PROGRAM_HEADER % GetParameter ( R_Min, 'RadiusMin' )

    T_Finish  =  1.0e-3_KDR  *  UNIT % SECOND

    call FM % Initialize &
           ( FluidType = 'IDEAL', &
             GravitationType = 'NEWTON_CM', &
             NameOption = Name, &
             FinishTimeOption = T_Finish, &
             RadiusMaxOption = R_Out, &
             RadiusExcisionOption = R_Min, &
             CentralMassOption = M )

    end associate !-- M, etc.

    !-- Modify from default OUTFLOW to INFLOW outer radial boundary condition
    select type ( I  =>  FM % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'INFLOW    ' ], iC = 1, iD = 1 )
    end associate !-- F
    end select !-- I

    FM % Integrator % SetInitial    =>  SetInitial
    FM % Integrator % System        =>  FM

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

  end subroutine SetInitial


end module FishboneMoncrief_Form
