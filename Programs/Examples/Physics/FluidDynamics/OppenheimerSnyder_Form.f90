module OppenheimerSnyder_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: OppenheimerSnyderForm
    real ( KDR ) :: &
      DensityInitial, &
      RadiusInitial, &
      Density, &
      Radius
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type OppenheimerSnyderForm

    private :: &
      SetFluid!, &
!      SetReference

      private :: &
        SetFluidKernel

contains


  subroutine Initialize ( OS, Name )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      TimeScale, &
      DensityFactor, &
      RadiusFactor, &
      Beta
    class ( Fluid_D_Form ), pointer :: &
      F

    if ( OS % Type == '' ) &
      OS % Type = 'a OppenheimerSnyder'

    call OS % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidCentralCoreForm :: OS % Integrator )
    select type ( FCC => OS % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'DUST', &
             GeometryType = 'NEWTONIAN', &
             DimensionlessOption = .true., &
             GravityFactorOption = 0.01_KDR, &
             LimiterParameterOption = 1.0_KDR )
!    FB % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )


    !-- Initial conditions

    associate &
      ( R0 => OS % RadiusInitial, &
        D0 => OS % DensityInitial, &
        DF => DensityFactor, &
        RF => RadiusFactor, &
        PI => CONSTANT % PI )

    R0  =  PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR
    D0  =  1.0_KDR
    DF  =  10.0_KDR
    call PROGRAM_HEADER % GetParameter ( R0, 'RadiusInitial' )
    call PROGRAM_HEADER % GetParameter ( D0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter ( DF, 'DensityFactor' )

    RF = DF ** ( - 1.0_KDR / 3.0_KDR ) 

    TimeScale         =  sqrt ( 3.0 / ( 8.0 * PI * D0 ) )
    Beta              =  acos ( sqrt ( RF ) )
    FCC % FinishTime  =  ( Beta  +  0.5 * sin ( 2 * Beta ) )  *  TimeScale

    call Show ( 'OppenheimerSnyder parameters' )
    call Show ( OS % DensityInitial, 'DensityInitial' )
    call Show ( OS % RadiusInitial, 'RadiusInitial' )
    call Show ( DensityFactor, 'DensityFactor' )
    call Show ( RadiusFactor, 'RadiusFactor' )
    call Show ( PI / 2  *  TimeScale, 'CollapseTime' )
    call Show ( FCC % FinishTime, 'Reset FinishTime' )

    end associate !-- R0, etc.

    F => FA % Fluid_D ( )
    call SetFluid ( OS, F, Time = 0.0_KDR )


    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F )


  end subroutine Initialize


  impure elemental subroutine Finalize ( OS )
    
    type ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    call OS % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( OS, F, Time )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time
    
    real ( KDR ) :: &
      Radius, &
      Density
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => OS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    if ( Time == 0.0_KDR ) then
      OS % Radius   =  OS % RadiusInitial
      OS % Density  =  OS % DensityInitial
    end if

    call SetFluidKernel &
           (    R = G % Value ( :, G % CENTER_U ( 1 ) ), &
             dR_L = G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
             dR_R = G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
             Density = OS % Density, &
             RadiusDensity = OS % Radius, &
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


end module OppenheimerSnyder_Form