module SedovTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidSymmetricCurvilinearForm ) :: SedovTaylorForm
    real ( KDR ), private :: &
      AdiabaticIndex, &
      Density, &
      BlastEnergy, &
      BlastRadiusRatio
  contains
    procedure, private, pass :: &
      Initialize_ST
    generic, public :: &
      Initialize => Initialize_ST
    final :: &
      Finalize
  end type SedovTaylorForm

    private :: &
      InitializeFluidSymmetricCurvilinear, &
      SetProblem

      private :: &
        SetFluid

        private :: &
          SetBlastVolumeFraction, &
          SetFluidKernel

contains


  subroutine Initialize_ST ( ST, Name )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    character ( * ), intent ( in ) :: &
      Name

    if ( ST % Type == '' ) &
      ST % Type = 'a SedovTaylor'

    call InitializeFluidSymmetricCurvilinear ( ST, Name )
    call SetProblem ( ST )

  end subroutine Initialize_ST


  impure elemental subroutine Finalize ( ST )
    
    type ( SedovTaylorForm ), intent ( inout ) :: &
      ST

  end subroutine Finalize


  subroutine InitializeFluidSymmetricCurvilinear ( ST, Name )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      RadiusMax, &
      FinishTime

    RadiusMax  = 0.35_KDR
    FinishTime = 0.05_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'MaxRadius' )
    call PROGRAM_HEADER % GetParameter ( FinishTime, 'FinishTime' )

    call ST % Initialize &
           ( FluidType = 'IDEAL', Name = Name, &
             FinishTimeOption = FinishTime, &
             RadiusMaxOption = RadiusMax, &
             nCellsRadiusOption = 64 )

  end subroutine InitializeFluidSymmetricCurvilinear


  subroutine SetProblem ( ST )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( I => ST % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

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

    G => PS % Geometry ( )
    F => FA % Fluid_P_I ( )
    call SetFluid ( ST, F, G )

    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine SetFluid ( ST, F, G )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G

    real ( KDR ) :: &
      BlastVolume, &
      EnergyDensity
    real ( KDR ), dimension ( : ), allocatable :: &
      BVF  !-- BlastVolumeFraction
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( I => ST % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( PSC => PS % Chart )
    class is ( Chart_SLD_Form )

    call F % SetAdiabaticIndex ( ST % AdiabaticIndex )

    !-- Find energy density

    allocate ( BVF ( F % nValues ) )
    call SetBlastVolumeFraction ( ST, PSC, G, BVF )

    associate ( dV => G % Value ( :, G % VOLUME ) )
    call CO % Initialize &
           ( PS % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    CO % Outgoing % Value ( 1 ) = sum ( dV * BVF )
    call CO % Reduce ( REDUCTION % SUM )
    BlastVolume = CO % Incoming % Value ( 1 )
    call Show ( BlastVolume, 'Discrete BlastVolume' )
    end associate !-- dV

    EnergyDensity = ST % BlastEnergy / BlastVolume 
    call Show ( EnergyDensity, 'EnergyDensity' )

    call SetFluidKernel &
           ( BVF, ST % Density, EnergyDensity, &
               N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
               E = F % Value ( :, F % INTERNAL_ENERGY ), &
             V_1 = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             V_2 = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             V_3 = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

!     end associate !-- N, etc.
    end select !-- PSC
    end select !-- PS
    end select !-- I

  end subroutine SetFluid


  subroutine SetBlastVolumeFraction ( ST, PSC, G, BVF )

    class ( SedovTaylorForm ), intent ( in ) :: &
      ST
    class ( Chart_SL_Template ), intent ( in ) :: &
      PSC
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      BVF

    integer ( KDI ) :: &
      iC, &  !-- iCell
      iS, jS, kS     !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ) :: &
      Pi, &
      BlastRadius, &
      dVS, &  !-- dVolumeSubcell
      VS      !-- VolumeSubcell
    real ( KDR ), dimension ( 3 ) :: &
      X_I, &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS      !--  X_Subcell

    call Clear ( BVF )

    nSubCells = 1
    nSubcells ( : PSC % nDimensions ) = 20

    Pi = CONSTANT % PI

    BlastRadius  =  ST % BlastRadiusRatio  *  PSC % MaxCoordinate ( 1 )
    call Show ( BlastRadius, 'BlastRadius' )

    do iC = 1, G % nValues
      associate &
        (  X   => G % Value ( iC, G % CENTER_U ( 1 ) &
                                  : G % CENTER_U ( 3 ) ), &
          dX_L => G % Value ( iC, G % WIDTH_LEFT_U ( 1 ) &
                                  : G % WIDTH_LEFT_U  ( 3 ) ), &
          dX_R => G % Value ( iC, G % WIDTH_RIGHT_U ( 1 ) &
                                  : G % WIDTH_RIGHT_U  ( 3 ) ) )
      associate ( dX => dX_L + dX_R )

      if ( .not. PSC % IsProperCell ( iC ) ) cycle

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
            select case ( PSC % nDimensions )
            case ( 1 ) !-- Spherical coordinates
              dVS = 4.0_KDR * Pi  *  XS ( 1 ) ** 2  *  dXS ( 1 )
            case ( 2 ) !-- Cylindrical coordinates
              dVS = 2.0_KDR * Pi * XS ( 1 ) * dXS ( 1 ) * dXS ( 2 )
            case ( 3 ) !-- Rectangular coordinates
              dVS = product ( dXS )
            end select !-- nD
            VS = VS + dVS
            !-- All cases here have Rectangular distances!
            if ( sqrt ( dot_product ( XS, XS ) ) <= BlastRadius ) &
              BVF ( iC ) = BVF ( iC ) + dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iC ) = BVF ( iC ) / VS

      end associate !-- dX
      end associate !-- X, etc.
    end do !-- iC

    call Show ( 4. / 3. * CONSTANT % Pi  *  BlastRadius ** 3, &
                'Continuous BlastVolume' )

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

      N = N0
      E = E0 * BVF
    V_1 = 0.0_KDR
    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

  end subroutine SetFluidKernel


end module SedovTaylor_Form
