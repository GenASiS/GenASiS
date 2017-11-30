module DustCollapse_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: DustCollapseForm
    real ( KDR ) :: &
      RadiusInitial, &
      DensityInitial, &
      Radius, &
      Density
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type DustCollapseForm

    private :: &
      SetFluid!, &
!      SetReference

      private :: &
        SetFluidKernel

contains


  subroutine Initialize ( DC, Name )

    class ( DustCollapseForm ), intent ( inout ) :: &
      DC
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      TimeScale
    class ( Fluid_D_Form ), pointer :: &
      F

    if ( DC % Type == '' ) &
      DC % Type = 'a DustCollapse'

    call DC % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidCentralCoreForm :: DC % Integrator )
    select type ( FCC => DC % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'DUST', GeometryType = 'NEWTONIAN', &
             GravitySolverTypeOption = 'MULTIPOLE' )
!    FB % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )


    !-- Initial conditions

    DC % RadiusInitial &
      =  PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR
    DC % DensityInitial &
      =  CONSTANT % SOLAR_BARYON_NUMBER &
         / ( 4.0_KDR / 3.0_KDR *  CONSTANT % PI  *  DC % RadiusInitial ** 3 )
    call PROGRAM_HEADER % GetParameter &
           ( DC % RadiusInitial, 'RadiusInitial' )
    call PROGRAM_HEADER % GetParameter &
           ( DC % DensityInitial, 'DensityInitial' )

    F => FA % Fluid_D ( )
    call SetFluid ( DC, F, Time = 0.0_KDR )

    TimeScale  &
      =  ( CONSTANT % GRAVITATIONAL  *  DC % DensityInitial ) ** ( -0.5_KDR )
    FCC % WriteTimeInterval &
      =  TimeScale / FCC % nWrite
    call Show ( 'Time Scales', DC % IGNORABILITY )
    call Show ( TimeScale, FCC % TimeUnit, 'TimeScaleDensityAve', &
                DC % IGNORABILITY )
    

    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F )


  end subroutine Initialize


  impure elemental subroutine Finalize ( DC )
    
    type ( DustCollapseForm ), intent ( inout ) :: &
      DC

    call DC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( DC, F, Time )

    class ( DustCollapseForm ), intent ( inout ) :: &
      DC
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time
    
    real ( KDR ) :: &
      Radius, &
      Density
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => DC % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    if ( Time == 0.0_KDR ) then
      DC % Radius   =  DC % RadiusInitial
      DC % Density  =  DC % DensityInitial
    end if

    call SetFluidKernel &
           (  R = G % Value ( :, G % CENTER ( 1 ) ), &
             dR = G % Value ( :, G % WIDTH ( 1 ) ), &
             Density = DC % Density, &
             RadiusDensity = DC % Radius, &
              N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetFluid


  subroutine SetFluidKernel ( R, dR, Density, RadiusDensity, N, VX, VY, VZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      dR
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
      ( R_In  => R - 0.5_KDR * dR, &
        R_Out => R + 0.5_KDR * dR )

    where ( R_Out <= RadiusDensity )
      N = Density
    end where
    where ( R_In < RadiusDensity .and. R_Out > RadiusDensity )
      N = Density * ( RadiusDensity - R_In ) / dR
    end where

    end associate !-- R_In, etc.

  end subroutine SetFluidKernel


end module DustCollapse_Form
