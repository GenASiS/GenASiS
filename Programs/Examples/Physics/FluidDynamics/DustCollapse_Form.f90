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
      TimeScale, &
      RadiusFraction, &
      Beta
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
           ( Name, FluidType = 'DUST', &
             GeometryType = 'NEWTONIAN', &
             GravitySolverTypeOption = 'MULTIPOLE', &
             DimensionlessOption = .true. )
!    FB % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )


    !-- Initial conditions

    associate &
      ( R0 => DC % RadiusInitial, &
        D0 => DC % DensityInitial, &
        RF => RadiusFraction, &
        PI => CONSTANT % PI )

    R0  =  PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR
    D0  =  1.0_KDR
    RF  =  0.5_KDR
    call PROGRAM_HEADER % GetParameter ( R0, 'RadiusInitial' )
    call PROGRAM_HEADER % GetParameter ( D0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter ( RF, 'RadiusFraction' )

    TimeScale         =  sqrt ( 3.0 / ( 8.0 * PI * D0 ) )
    Beta              =  acos ( sqrt ( RF ) )
    FCC % FinishTime  =  ( Beta  +  0.5 * sin ( 2 * Beta ) )  *  TimeScale

    call Show ( 'DustCollapse parameters' )
    call Show ( DC % RadiusInitial, 'RadiusInitial' )
    call Show ( DC % DensityInitial, 'DensityInitial' )
    call Show ( PI / 2  *  TimeScale, 'CollapseTime' )
    call Show ( RadiusFraction, 'RadiusFraction' )
    call Show ( FCC % FinishTime, 'Reset FinishTime' )

    end associate !-- R0, etc.

!    FCC % WriteTimeInterval &
!      =  TimeScale / FCC % nWrite
!    call Show ( 'Time Scales', DC % IGNORABILITY )
!    call Show ( TimeScale, FCC % TimeUnit, 'TimeScaleDensityAve', &
!                DC % IGNORABILITY )
    
    F => FA % Fluid_D ( )
    call SetFluid ( DC, F, Time = 0.0_KDR )


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
