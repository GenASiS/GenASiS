module SedovTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SedovTaylorForm
    real ( KDR ), private :: &
      AdiabaticIndex, &
      Density, &
      BlastEnergy, &
      BlastRadiusRatio
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type SedovTaylorForm

    private :: &
      SetFluid

contains


  subroutine Initialize ( ST, Name )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      RadiusMax, &
      FinishTime

    if ( ST % Type == '' ) &
      ST % Type = 'a SedovTaylor'

    call ST % InitializeTemplate ( Name )

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


    !-- Integrator

    RadiusMax  = 0.35_KDR
    FinishTime = 0.05_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'MaxRadius' )
    call PROGRAM_HEADER % GetParameter ( FinishTime, 'FinishTime' )

    allocate ( FluidSymmetricCurvilinearForm :: ST % Integrator )
    select type ( FSC => ST % Integrator )
    type is ( FluidSymmetricCurvilinearForm )
    call FSC % Initialize &
           ( Name, FluidType = 'IDEAL', &
             FinishTimeOption = FinishTime, &
             RadiusMaxOption = RadiusMax, &
             nCellsRadiusOption = 64 )
    end select !-- FSC


    !-- Initial conditions

    call SetFluid ( ST )


  end subroutine Initialize


  impure elemental subroutine Finalize ( ST )
    
    type ( SedovTaylorForm ), intent ( inout ) :: &
      ST

    call ST % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( ST )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iS, jS, kS     !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ) :: &
      BlastRadius, &
      BlastVolume, &
      EnergyDensity, &
      dVS, &  !-- dVolumeSubcell
      VS      !-- VolumeSubcell
    real ( KDR ), dimension ( 3 ) :: &
      X_I, &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS      !--  X_Subcell
    real ( KDR ), dimension ( : ), allocatable :: &
      BVF  !-- VolumeFraction
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( FSC => ST % Integrator )
    class is ( FluidSymmetricCurvilinearForm )

    select type ( FA => FSC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    call F % SetAdiabaticIndex ( ST % AdiabaticIndex )

    select type ( PS => FSC % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )

    !-- Find energy density

    BlastRadius  =  ST % BlastRadiusRatio  *  C % MaxCoordinate ( 1 )
    call Show ( BlastRadius, 'BlastRadius' )

    associate &
      (   N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
        V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
          E => F % Value ( :, F % INTERNAL_ENERGY ), &
         dV => G % Value ( :, G % VOLUME ) )

    allocate ( BVF ( size ( dV ) ) )
    call Clear ( BVF )

    nSubCells = 1
    nSubcells ( : C % nDimensions ) = 20

    do iC = 1, size ( dV )
      associate &
        (  X   => G % Value ( iC, G % CENTER_U ( 1 ) &
                                  : G % CENTER_U ( 3 ) ), &
          dX_L => G % Value ( iC, G % WIDTH_LEFT_U ( 1 ) &
                                  : G % WIDTH_LEFT_U  ( 3 ) ), &
          dX_R => G % Value ( iC, G % WIDTH_RIGHT_U ( 1 ) &
                                  : G % WIDTH_RIGHT_U  ( 3 ) ), &
          Pi   => CONSTANT % PI )
      associate ( dX => dX_L + dX_R )

      if ( .not. C % IsProperCell ( iC ) ) cycle

      if ( sqrt ( dot_product ( X, X ) ) &
           - 0.5_KDR * sqrt ( dot_product ( dX, dX ) ) > BlastRadius ) &
      then
        BVF ( iC ) = 0.0_KDR
        cycle
      else if ( sqrt ( dot_product ( X, X ) ) &
                + 0.5_KDR * sqrt ( dot_product ( dX, dX ) ) <= BlastRadius ) &
      then
        BVF ( iC ) = 1.0_KDR
        cycle
      end if

      X_I  =  X  -  dX_L
      X_O  =  X  +  dX_R
      dXS  =  ( X_O  -  X_I ) / nSubcells

      VS = 0.0_KDR
      do kS = 1, nSubcells ( 3 )
        do jS = 1, nSubcells ( 2 )
          do iS = 1, nSubcells ( 1 )
            XS  =  X_I  +  ( [ iS, jS, kS ] - 0.5_KDR ) * dXS
            select case ( C % nDimensions )
            case ( 1 ) !-- Spherical coordinates
              dVS = 4.0_KDR * Pi  *  XS ( 1 ) ** 2  *  dXS ( 1 )
            case ( 2 ) !-- Cylindrical coordinates
              dVS = 2.0_KDR * Pi * XS ( 1 ) * dXS ( 1 ) * dXS ( 2 )
            case ( 3 ) !-- Cartesian coordinates
              dVS = product ( dXS )
            end select !-- nD
            VS = VS + dVS
            !-- All cases here have Cartesian distances!
            if ( sqrt ( dot_product ( XS, XS ) ) <= BlastRadius ) &
              BVF ( iC ) = BVF ( iC ) + dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iC ) = BVF ( iC ) / VS

      end associate !-- dX
      end associate !-- X, etc.
    end do !-- iC

    call CO % Initialize &
           ( PS % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    CO % Outgoing % Value ( 1 ) = sum ( dV * BVF )
    call CO % Reduce ( REDUCTION % SUM )
    BlastVolume = CO % Incoming % Value ( 1 )
    call Show ( BlastVolume, 'Discrete BlastVolume' )
    call Show ( 4. / 3. * CONSTANT % Pi  *  BlastRadius ** 3, &
                'Continuous BlastVolume' )

    EnergyDensity = ST % BlastEnergy / BlastVolume 
    call Show ( EnergyDensity, 'EnergyDensity' )

    !-- Set fluid values

      N = ST % Density
    V_1 = 0.0_KDR
    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    E = EnergyDensity * BVF

    call F % ComputeFromPrimitive ( G )
    call C % ExchangeGhostData ( F )

    end associate !-- N, etc.
    end select !-- C
    end select !-- PS
    end select !-- FA
    end select !-- FSC
    nullify ( F, G )

  end subroutine SetFluid

    
end module SedovTaylor_Form
