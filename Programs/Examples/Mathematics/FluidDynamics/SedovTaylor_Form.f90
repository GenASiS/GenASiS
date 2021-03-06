module SedovTaylor_Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  use Fluid_P_P__Form
  use Fluid_ASC__Form
  use ApplyCurvilinear_F__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Form ) :: SedovTaylorForm
  contains
    procedure, private, pass :: &
      Initialize_ST
    generic, public :: &
      Initialize => Initialize_ST
    final :: &
      Finalize  
  end type SedovTaylorForm

    private :: &
      SetFluid

    real ( KDR ), private :: &
      AdiabaticIndex, &
      Density, &
      BlastEnergy, &
      BlastRadiusRatio

    type ( UniverseHeaderForm ), private :: &
      Universe  !-- Non-functional dummy argument

contains


  subroutine Initialize_ST ( ST, Name )

    class ( SedovTaylorForm ), intent ( inout ) :: &
      ST
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nCellsRadius
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      MaxRadius, &
      FinishTime
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDL ) :: &
      CoordinateSystem

    MaxRadius        = 0.35_KDR
    AdiabaticIndex   = 1.4_KDR
    Density          = 1.0_KDR
    BlastEnergy      = 1.0_KDR
    BlastRadiusRatio = 0.03_KDR
    FinishTime       = 0.05_KDR

    call PROGRAM_HEADER % GetParameter ( MaxRadius, 'MaxRadius' )
    call PROGRAM_HEADER % GetParameter ( AdiabaticIndex, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )
    call PROGRAM_HEADER % GetParameter ( BlastEnergy, 'BlastEnergy' )
    call PROGRAM_HEADER % GetParameter ( BlastRadiusRatio, 'BlastRadiusRatio' )
    call PROGRAM_HEADER % GetParameter ( FinishTime, 'FinishTime' )

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: ST % PositionSpace )
    select type ( PS => ST % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    nCellsRadius = 64
    call PROGRAM_HEADER % GetParameter ( nCellsRadius, 'nCellsRadius' )

    select case ( PS % nDimensions )
    case ( 1 )  !-- Spherical coordinates

      CoordinateSystem = 'SPHERICAL' 

      call PS % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )

      MinCoordinate = [ 0.0_KDR,   0.0_KDR, 0.0_KDR ]
      MaxCoordinate = [ MaxRadius, 0.0_KDR, 0.0_KDR ]

      nCells = [ nCellsRadius, 1, 1 ]

    case ( 2 )  !-- Cylindrical coordinates

      CoordinateSystem = 'CYLINDRICAL' 

      call PS % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
      call PS % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW   ', 'OUTFLOW   ' ], iDimension = 2 )

      MinCoordinate = [ 0.0_KDR,   - MaxRadius, 0.0_KDR ]
      MaxCoordinate = [ MaxRadius, + MaxRadius, 0.0_KDR ]

      nCells = [ nCellsRadius, 2 * nCellsRadius, 1 ]

    case ( 3 )  !-- Rectangular coordinates
 
      CoordinateSystem = 'RECTANGULAR' 

      call PS % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW   ', 'OUTFLOW   ' ], iDimension = 1 )
      call PS % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW   ', 'OUTFLOW   ' ], iDimension = 2 )
      call PS % SetBoundaryConditionsFace &
             ( [ 'OUTFLOW   ', 'OUTFLOW   ' ], iDimension = 3 )

      MinCoordinate  =  - MaxRadius
      MaxCoordinate  =  + MaxRadius

      nCells  =  2 * nCellsRadius
      
    end select !-- nDimensions

    call PS % CreateChart &
           ( CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nCells )

    call PS % SetGeometry ( )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: ST % Current_ASC )
    select type ( FA => ST % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'POLYTROPIC' )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: ST % Step )
    select type ( S => ST % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( ST, FA, Name )
    S % ApplySources % Pointer => ApplyCurvilinear_F
    end select !-- S

    !-- Set fluid and initialize Integrator template

    call SetFluid ( ST )
    call ST % Initialize ( Universe, Name, FinishTimeOption = FinishTime )

    !-- Cleanup

    end select !-- FA
    end select !-- PS

  end subroutine Initialize_ST


  subroutine Finalize ( ST )

    type ( SedovTaylorForm ), intent ( inout ) :: &
      ST

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
    class ( Fluid_P_P_Form ), pointer :: &
      F
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( FA => ST % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_P ( )

    call F % SetAdiabaticIndex ( AdiabaticIndex )

    select type ( PS => ST % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )

    !-- Find energy density

    BlastRadius = BlastRadiusRatio * C % MaxCoordinate ( 1 )
    call Show ( BlastRadius, 'BlastRadius' )

    associate &
      (   N => F % Value ( :, F % COMOVING_DENSITY ), &
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

    call CO % Initialize &
           ( PS % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    CO % Outgoing % Value ( 1 ) = sum ( dV * BVF )
    call CO % Reduce ( REDUCTION % SUM )
    BlastVolume = CO % Incoming % Value ( 1 )
    call Show ( BlastVolume, 'Discrete BlastVolume' )
    call Show ( 4. / 3. * CONSTANT % Pi  *  BlastRadius ** 3, &
                'Continuous BlastVolume' )

    EnergyDensity = BlastEnergy / BlastVolume 
    call Show ( EnergyDensity, 'EnergyDensity' )

    !-- Set fluid values

      N = Density
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
    nullify ( F, G )

  end subroutine SetFluid

    
end module SedovTaylor_Form
