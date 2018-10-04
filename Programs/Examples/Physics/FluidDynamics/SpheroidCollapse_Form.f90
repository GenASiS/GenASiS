module SpheroidCollapse_Form

! From Lin, Mestel, and Shu 1965
! "The Gravitational Collapse of a Uniform Spheroid"

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SpheroidCollapseForm
    Real ( KDR ) :: &
      Eccentricity
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
    
    real ( KDR ) :: &
      Radius, &
      Density
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => SC % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetFluidKernel &
           (    R = G % Value ( :, G % CENTER_U ( 1 ) ), &
             dR_L = G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
             dR_R = G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
             Density = 1.0_KDR, & !OS % Density, &
             RadiusDensity = 1.0_KDR, & !OS % Radius, &
              N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetFluid

  
    subroutine SetFluidKernel &
               ( R, dR_L, dR_R, Density, RadiusDensity, N, VX, VY, VZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      dR_L, &
      dR_R
    real ( KDR ), intent ( in ) :: &
      Density, &
      RadiusDensity
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      VX, VY, VZ

    N  = 0.0_KDR
    VX = 0.0_KDR
    VY = 0.0_KDR
    VZ = 0.0_KDR

    associate &
      ( R_In  => R - dR_L, &
        R_Out => R + dR_R )

    where ( R_Out <= RadiusDensity )
      N = Density
    end where
    where ( R_In < RadiusDensity .and. R_Out > RadiusDensity )
      N = Density * ( RadiusDensity ** 3  -  R_In ** 3 ) &
                    / ( R_Out ** 3  -  R_In ** 3 )
    end where

!    N = Density / ( 1 + exp ( ( R - RadiusDensity ) &
!                              / ( 3 * ( dR_L + dR_R ) ) ) ) 
 
    end associate !-- R_In, etc.

  end subroutine SetFluidKernel

end module SpheroidCollapse_Form


