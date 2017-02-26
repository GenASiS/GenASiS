module SedovTaylor_Form

  use Basics
  use Mathematics
  use Fluid_P__Template
  use Fluid_P_P__Form
  use Fluid_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: SedovTaylorForm
  contains
    procedure, public, pass :: &
      Initialize
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

contains


  subroutine Initialize ( ST, Name )

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

    if ( ST % Type == '' ) &
      ST % Type = 'a SedovTaylor' 

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
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

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

    case ( 3 )  !-- Cartesian coordinates
 
      CoordinateSystem = 'CARTESIAN' 

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

    !-- Geometry of PositionSpace

    allocate ( ST % Geometry_ASC )
    associate ( GA => ST % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Fluid

    allocate ( Fluid_ASC_Form :: ST % Current_ASC )
    select type ( FA => ST % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'POLYTROPIC' )
    end select !-- FA

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: ST % Step )
    select type ( S => ST % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( Name )
    S % ApplySources % Pointer => ApplySourcesCurvilinear_Fluid_P
    end select !-- S

    !-- Set fluid and initialize Integrator template

    call SetFluid ( ST )
    call ST % InitializeTemplate_C_PS ( Name, FinishTimeOption = FinishTime )

    !-- Cleanup

    end select !-- PS

  end subroutine Initialize


  subroutine Finalize ( ST )

    type ( SedovTaylorForm ), intent ( inout ) :: &
      ST

    call ST % FinalizeTemplate_C_PS ( )

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
      dV, &
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
         VJ => G % Value ( :, G % VOLUME_JACOBIAN ) )

    allocate ( dV ( size ( VJ ) ) )
    allocate ( BVF ( size ( VJ ) ) )
    call Clear ( dV )
    call Clear ( BVF )

    nSubCells = 1
    nSubcells ( : C % nDimensions ) = 20

    do iC = 1, size ( VJ )
      associate &
        (  X => G % Value ( iC, G % CENTER ( 1 ) : G % CENTER ( 3 ) ), &
          dX => G % Value ( iC, G % WIDTH  ( 1 ) : G % WIDTH  ( 3 ) ), &
          Pi => CONSTANT % PI )

      if ( .not. C % IsProperCell ( iC ) ) cycle

      dV ( iC )  =  VJ ( iC ) * product ( dX ( : C % nDimensions ) )

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

      X_I  =  X  -  0.5_KDR * dX
      X_O  =  X  +  0.5_KDR * dX
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
