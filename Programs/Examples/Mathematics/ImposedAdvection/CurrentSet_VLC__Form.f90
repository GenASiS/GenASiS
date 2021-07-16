module CurrentSet_VLC__Form
  
  !-- CurrentSet_VelocityLinearConstant_Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public, extends ( CurrentSetForm ) :: CurrentSet_VLC_Form
    real ( KDR ) :: &
      SpeedInitial, &
      LengthInitial
  contains
    procedure, private, pass :: &
      Initialize_CS_VL
    generic, public :: &
      Initialize => Initialize_CS_VL
    procedure, public, pass :: &
      SetVelocityLinear
    final :: &
      Finalize
  end type CurrentSet_VLC_Form


contains


  subroutine Initialize_CS_VL ( CS, G, SpeedInitial, LengthInitial )

    class ( CurrentSet_VLC_Form ), intent ( inout ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    real ( KDR ) :: &
      SpeedInitial, &
      LengthInitial

    CS % SpeedInitial   =  SpeedInitial
    CS % LengthInitial  =  LengthInitial

    call CS % CurrentSetForm % Initialize ( G )

  end subroutine Initialize_CS_VL


  subroutine SetVelocityLinear ( CS, Speed, Length )

    class ( CurrentSet_VLC_Form ), intent ( inout ) :: &
      CS
    real ( KDR ), intent ( in ) :: &
      Speed, &
      Length

    associate &
      ( G  =>  CS % Geometry )
    associate &
      ( CSV  =>  CS % Storage_GS % Value, &
         GV  =>   G % Storage_GS % Value )
    associate &
      ( V_1  =>  CSV ( :, CS % VELOCITY_CS_U_1 ), &
        V_2  =>  CSV ( :, CS % VELOCITY_CS_U_2 ), &
        V_3  =>  CSV ( :, CS % VELOCITY_CS_U_3 ), &
        X_1  =>   GV ( :,  G % CENTER_U_1 ) )

    V_1  =  Speed  *  ( X_1 / Length )
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR
    
    end associate !-- V_1, etc.
    end associate !-- CV, etc.
    end associate !-- G

  end subroutine SetVelocityLinear


  impure elemental subroutine Finalize ( CS )

    type ( CurrentSet_VLC_Form ), intent ( inout ) :: &
      CS

  end subroutine Finalize


end module CurrentSet_VLC__Form
