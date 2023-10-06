module SedovTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_SC_Form ) :: SedovTaylorForm
    real ( KDR ), private :: &
      AdiabaticIndex, &
      Density, &
      BlastEnergy, &
      BlastRadiusRatio, &
      ContinuousBlastVolume, &
      DiscreteBlastVolume, &
      EnergyDensity
  contains
    procedure, private, pass :: &
      Initialize_H
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type SedovTaylorForm

    private :: &
      InitializeUniverse, &
      SetInitial

      private :: &
        SetFluid

        private :: &
          SetBlastVolumeFraction, &
          SetFluidKernel
    

contains


  subroutine Initialize_H ( U, Name, CommunicatorOption )

    class ( SedovTaylorForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a SedovTaylor'

    call InitializeUniverse ( U, Name )

  end subroutine Initialize_H


  subroutine Finalize ( ST )

    type ( SedovTaylorForm ), intent ( inout ) :: &
      ST

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( SedovTaylorForm ), intent ( in ) :: &
      U

    call U % Universe_F_SC_Form % ShowParameters ( )

    call Show ( U % Density,               'Density' )
    call Show ( U % AdiabaticIndex,        'AdiabaticIndex' )
    call Show ( U % BlastEnergy,           'BlastEnergy' )
    call Show ( U % BlastRadiusRatio,      'BlastRadiusRatio' )
    call Show ( U % ContinuousBlastVolume, 'ContinuousBlastVolume' )
    call Show ( U % DiscreteBlastVolume,   'DiscreteBlastVolume' )
    call Show ( U % EnergyDensity,         'EnergyDensity' )

  end subroutine ShowParameters


  subroutine InitializeUniverse ( ST, Name )

    class ( SedovTaylorForm ), intent ( inout ), target :: &
      ST
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      RadiusMax, &
      FinishTime

    RadiusMax   =  0.35_KDR
    FinishTime  =  0.05_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'RadiusMax' )
    call PROGRAM_HEADER % GetParameter ( FinishTime, 'FinishTime' )

    call ST % Initialize &
           ( FluidType = 'IDEAL', &
             RadiusMax = RadiusMax, &
             Name = Name, &
             FinishTimeOption = FinishTime, &
             nCellsRadiusOption = 64 )

    ST % Integrator % SetInitial  =>  SetInitial
    ST % Integrator % System      =>  ST

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( ST  =>  I % System )
      class is ( SedovTaylorForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    ST % AdiabaticIndex   = 1.4_KDR
    ST % Density          = 1.0_KDR
    ST % BlastEnergy      = 1.0_KDR
    ST % BlastRadiusRatio = 0.03_KDR

    call PROGRAM_HEADER % GetParameter &
           ( ST % AdiabaticIndex, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter &
           ( ST % Density, 'Density' )
    call PROGRAM_HEADER % GetParameter &
           ( ST % BlastEnergy, 'BlastEnergy' )
    call PROGRAM_HEADER % GetParameter &
           ( ST % BlastRadiusRatio, 'BlastRadiusRatio' )
  
    call SetFluid ( ST, F )

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- ST

  end subroutine SetInitial


  subroutine SetFluid ( ST, F )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    real ( KDR ), dimension ( : ), allocatable :: &
      BVF  !-- BlastVolumeFraction
    type ( CollectiveOperation_R_Form ) :: &
      CO

    call F % SetAdiabaticIndex ( ST % AdiabaticIndex )

    associate &
      ( G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    !-- Find energy density

    allocate ( BVF ( F % Storage_GS % nValues ) )
    call SetBlastVolumeFraction ( ST, G, BVF )

    associate ( dV => GV ( :, G % VOLUME ) )
    call CO % Initialize &
           ( PROGRAM_HEADER % Communicator, &
             nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    CO % Outgoing % Value ( 1 ) = sum ( dV * BVF )
    call CO % Reduce ( REDUCTION % SUM )
    ST % DiscreteBlastVolume  =  CO % Incoming % Value ( 1 )
    end associate !-- dV

    ST % EnergyDensity  =  ST % BlastEnergy  /  ST % DiscreteBlastVolume 

    call SetFluidKernel &
           ( BVF, ST % Density, ST % EnergyDensity, &
               N = FV ( :, F % BARYON_DENSITY_C ), &
               E = FV ( :, F % ENERGY_DENSITY_C ), &
             V_1 = FV ( :, F % VELOCITY_U_1 ), &
             V_2 = FV ( :, F % VELOCITY_U_2 ), &
             V_3 = FV ( :, F % VELOCITY_U_3 ) )

    end associate !-- FV, etc.
    end associate !-- G

  end subroutine SetFluid


  subroutine SetBlastVolumeFraction ( ST, G, BVF )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      BVF

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iS, jS, kS     !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nS  !-- nSubcells
    real ( KDR ) :: &
      Pi, &
      BlastRadius, &
      dVS, &  !-- dVolumeSubcell
      VS      !-- VolumeSubcell
    real ( KDR ), dimension ( 3 ) :: &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS      !--  X_Subcell

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        GV  =>  G % Storage_GS % Value )

    call Clear ( BVF )

    nS  =  1
    nS ( : C % nDimensions )  =  20

    Pi = CONSTANT % PI

    BlastRadius  =  ST % BlastRadiusRatio  *  C % MaxCoordinate ( 1 )

    do iC  =  1,  G % Storage_GS % nValues
      associate &
        (  X    =>  GV ( iC, G % CENTER_U_1 : G % CENTER_U_3 ), &
           X_I  =>  GV ( iC, G % EDGE_I_U_1 : G % EDGE_I_U_3 ), &
          dX    =>  GV ( iC, G %  WIDTH_U_1 : G %  WIDTH_U_3 ) )

      X_O  =  X_I  +  dX

      if ( .not. C % ProperCell ( iC ) ) &
        cycle

      if ( sqrt ( dot_product ( X, X ) ) &
           - 0.5_KDR * sqrt ( dot_product ( dX, dX ) )  >  BlastRadius ) &
      then
        BVF ( iC ) = 0.0_KDR
        cycle
      else if ( sqrt ( dot_product ( X, X ) ) &
                + 0.5_KDR * sqrt ( dot_product ( dX, dX ) )  <=  BlastRadius ) &
      then
        BVF ( iC ) = 1.0_KDR
        cycle
      end if

      dXS  =  ( X_O  -  X_I )  /  nS

      VS  =  0.0_KDR
      do kS  =  1,  nS ( 3 )
        do jS  =  1,  nS ( 2 )
          do iS  =  1,  nS ( 1 )
            XS  =  X_I  +  ( [ iS, jS, kS ] - 0.5_KDR ) * dXS
            select case ( C % nDimensions )
            case ( 1 ) !-- Spherical coordinates
              dVS  =  4.0_KDR * Pi  *  XS ( 1 ) ** 2  *  dXS ( 1 )
            case ( 2 ) !-- Cylindrical coordinates
              dVS  =  2.0_KDR * Pi * XS ( 1 ) * dXS ( 1 ) * dXS ( 2 )
            case ( 3 ) !-- Rectangular coordinates
              dVS  =  product ( dXS )
            end select !-- nD
            VS = VS + dVS
            !-- All cases here have Rectangular distances!
            if ( sqrt ( dot_product ( XS, XS ) ) <= BlastRadius ) &
              BVF ( iC )  =  BVF ( iC ) + dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iC ) = BVF ( iC ) / VS

      end associate !-- X_I, etc.
    end do !-- iC

    ST % ContinuousBlastVolume &
      =  4.0_KDR / 3.0_KDR * CONSTANT % Pi  *  BlastRadius ** 3

    end associate !-- C
    end select !-- A

  end subroutine SetBlastVolumeFraction


  subroutine SetFluidKernel ( BVF, N0, E0, N, E, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      BVF
    real ( KDR ), intent ( in ) :: &
      N0, &
      E0
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      E, &
      V_1, V_2, V_3

      N  =  N0
      E  =  max ( E0  *  BVF, E0  *  1.0e-10 )
    V_1  =  0.0_KDR
    V_2  =  0.0_KDR
    V_3  =  0.0_KDR

  end subroutine SetFluidKernel


end module SedovTaylor_Form
