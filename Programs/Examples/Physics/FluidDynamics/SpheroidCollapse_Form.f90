module SpheroidCollapse_Form

! From Lin, Mestel, and Shu 1965
! "The Gravitational Collapse of a Uniform Spheroid"

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SpheroidCollapseForm
    Real ( KDR ) :: &
      SemiMajor, &
      SemiMinor, &
      Eccentricity, &
      InitialDensity
  contains
    procedure, public, pass :: &
      Initialize 
    final :: &
      Finalize
  end type SpheroidCollapseForm

    private :: &
      SetFluid!, &
      !SetReference

contains


  subroutine Initialize ( SC, Name )

    class ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC
    character ( * ), intent ( in ) :: &
      Name

    class ( Fluid_D_Form ), pointer :: &
      F

    real ( KDR ) :: &
     TimeScale

    if ( SC % Type == '' ) &
      SC % Type = ' a SpheroidCollapse'

    call SC % InitializeTemplate ( Name )

    !-- Integrator

    allocate ( FluidCentralCoreForm :: SC % Integrator )
    select type ( FCC => SC % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'DUST' ,&
             GeometryType = 'NEWTONIAN', &
             DimensionlessOption = .true., &
             GravityFactorOption = 0.01_KDR, &
             LimiterParameterOption = 1.0_KDR )

   ! FCC % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )

    !-- Initial conditions

    SC % SemiMajor = PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR
    SC % SemiMinor = SC % SemiMajor * 0.8
    SC % InitialDensity = 1.0_KDR

    TimeScale &
      = sqrt ( 3.0_KDR / ( 8.0_KDR * CONSTANT % PI * SC % InitialDensity ) )

    FCC % FinishTime = 1.2 * TimeScale
    call Show ( FCC % FinishTime, 'Reset FinishTime' )

    F => FA % Fluid_D ( )
    call SetFluid ( SC, F, Time = 0.0_KDR )

    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SC )

    type ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC 
    
    call SC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( SC, F, Time )

    class ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => SC % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetFluidKernel &
           (  PS, &
              a_1 = SC % SemiMajor, &
              a_3 = SC % SemiMinor, &
              D0  = SC % InitialDensity, &
              N   = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             V_1 = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             V_2 = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             V_3 = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetFluid

  
   subroutine SetFluidKernel ( A, a_1, a_3, D0, N, V_1, V_2, V_3 )

    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      a_1, &
      a_3, &
      D0
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      V_1, &
      V_2, &
      V_3

    class ( GeometryFlatForm ), pointer :: &
      G
    integer ( KDI ) :: &
      i, &
      iC, &  !-- iCell
      iS, jS, kS     !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ), dimension ( 3 ) :: &
      X_I, &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS      !--  X_Subcell
    real ( KDR ), dimension ( : ), allocatable :: &
      BVF, &  !-- VolumeFraction
      rho_sq, & !-- X^2 + Y^2 
      Z_sq
    real ( KDR ) :: &
      dVS, &  !-- dVolumeSubcell
      VS, &   !-- VolumeSubcell
      rho_sq_in_in, &
      rho_sq_in_out, &
      rho_sq_out_in, &
      rho_sq_out_out, &
      Z_sq_in_in, &
      Z_sq_in_out, &
      Z_sq_out_in, &
      Z_sq_out_out, &
      rho_S_sq, &
      Z_s_sq

    N   = 0.0_KDR
    V_1 = 0.0_KDR
    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )

    associate ( dV => G % Value ( :, G % VOLUME ) )

    allocate ( BVF ( size ( dV ) ) )
    call Clear ( BVF )

    nSubCells = 1
    nSubcells ( : C % nDimensions ) = 20

    do iC = 1, size ( dV )
      associate &
        (  X => &
             G % Value ( iC, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ), &
          dX_L => &
            G % Value ( iC, G % WIDTH_LEFT_U ( 1 ) &
                              : G % WIDTH_LEFT_U ( 3 ) ), &
          dX_R => &
            G % Value ( iC, G % WIDTH_RIGHT_U ( 1 ) &
                              : G % WIDTH_RIGHT_U ( 3 ) ), &
          Pi => CONSTANT % PI )

      if ( .not. C % IsProperCell ( iC ) ) cycle

       X_I  =  X  - dX_L
       X_O  =  X  + dX_R

       rho_sq_in_in   = ( X_I ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
       rho_sq_in_out  = ( X_O ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
       rho_sq_out_in  = ( X_I ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2
       rho_sq_out_out = ( X_O ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2

       Z_sq_in_in   = ( X_I ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
       Z_sq_in_out  = ( X_O ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
       Z_sq_out_in  = ( X_I ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2
       Z_sq_out_out = ( X_O ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2

       if ( rho_sq_in_in / a_1 ** 2 + Z_sq_in_in / a_3 ** 2 <= 1.0_KDR &
            .and. rho_sq_in_out / a_1 ** 2 &
                    + Z_sq_in_out / a_3 ** 2 <= 1.0_KDR &
            .and. rho_sq_out_in / a_1 ** 2 &
                    + Z_sq_out_in / a_3 ** 2 <= 1.0_KDR &
            .and. rho_sq_out_out / a_1 ** 2 &
                    + Z_sq_out_out / a_3 ** 2 <= 1.0_KDR ) then 
         BVF ( iC ) = 1.0_KDR
         cycle
       end if

      if ( rho_sq_in_in / a_1 ** 2 + Z_sq_in_in / a_3 ** 2 > 1.0_KDR &
           .and. rho_sq_in_out / a_1 ** 2 &
                   + Z_sq_in_out / a_3 ** 2 > 1.0_KDR &
           .and. rho_sq_out_in / a_1 ** 2 &
                   + Z_sq_out_in / a_3 ** 2 > 1.0_KDR &
           .and. rho_sq_out_out / a_1 ** 2 &
                   + Z_sq_out_out / a_3 ** 2 > 1.0_KDR ) then 
        BVF ( iC ) = 0.0_KDR
        cycle
      end if
      
      dXS  =  ( X_O  -  X_I ) / nSubcells

      VS = 0.0_KDR
      do kS = 1, nSubcells ( 3 )
        do jS = 1, nSubcells ( 2 )
          do iS = 1, nSubcells ( 1 )
            XS  =  X_I  +  ( [ iS, jS, kS ] - 0.5_KDR ) * dXS
            rho_S_sq = ( XS ( 1 ) * sin ( XS ( 2 ) ) ) ** 2
            Z_S_sq = ( XS ( 1 ) * cos ( XS ( 2 ) ) ) ** 2
            select case ( C % nDimensions )
            case ( 2 )
              dVS = 2 * Pi * XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                      * dXS ( 1 ) * dXS ( 2 )
            case ( 3 )
              dVS = XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                     * dXS ( 1 ) * dXS ( 2 ) * dXS ( 3 )
            end select 
            VS = VS + dVS
            if ( rho_S_sq / a_1 ** 2 + Z_S_sq / a_3 ** 2 <= 1.0_KDR ) &
              BVF ( iC ) = BVF ( iC ) + dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iC ) = BVF ( iC ) / VS

      end associate !-- X, etc.
    end do !-- iC

    end associate !-- dV
    end select !-- C

    nullify ( G )

    N = D0 * BVF

  end subroutine SetFluidKernel

end module SpheroidCollapse_Form


