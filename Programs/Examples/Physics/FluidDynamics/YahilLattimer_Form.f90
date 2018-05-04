module YahilLattimer_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: YahilLattimerForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type YahilLattimerForm

    private :: &
      SetFluid, &
      PrepareInterpolation
!      SetReference

!      private :: &
!        SetFluidKernel

    integer ( KDI ), private, parameter :: &
      iX_TS  = 1, &  !-- must match the profile file columns
      iD_TS  = 2, &
      iV_TS  = 3, &
      iM_TS  = 4
    integer ( KDI ), private, parameter :: &
      iX_SI  = 1, & !-- spline interpolation
      iD_SI  = 2, & 
      iV_SI  = 3, &
      iM_SI  = 4

contains


  subroutine Initialize ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    character ( * ), intent ( in ) :: &
      Name

    class ( Fluid_P_I_Form ), pointer :: &
      F

    if ( YL % Type == '' ) &
      YL % Type = 'a YahilLattimer'

    call YL % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidCentralCoreForm :: YL % Integrator )
    select type ( FCC => YL % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'IDEAL', &
             GeometryType = 'NEWTONIAN' )
!    FB % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )


    !-- Initial Conditions

!-- FIXME: Set default parameters, read parameters

    F => FA % Fluid_P_I ( )
    call SetFluid ( YL, F, Time = 0.0_KDR )


    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F )

  end subroutine Initialize


  impure elemental subroutine Finalize ( YL )
    
    type ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    call YL % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( YL, F, Time )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

  end subroutine SetFluid


    subroutine PrepareInterpolation ( SI )

    type ( SplineInterpolationForm ), dimension ( 4 ), intent ( inout ) :: &
      SI

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Slope_X, &
      Slope_D, &
      Slope_V,&
      Slope_M
    real ( KDR ), dimension ( : ), allocatable :: &
      XC, &   !-- XCenter
      dXC, &  !-- WidthCenter
      X_edge, &
      D_edge, &
      V_edge
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Profile
    character ( LDF ) :: &
      Path, &
      Filename
    type ( TableStreamForm ) :: &
      TS

    Path = '../Parameters/'
    Filename = 'YahilHomologousCollapse_Gm_130.dat'
    
    call PROGRAM_HEADER % GetParameter &
           ( Filename, 'Filename' )

    call TS % Initialize &
           ( Filename, PROGRAM_HEADER % Communicator % Rank, &
             PathOption = Path )
    call TS % Read ( Profile, oRowOption = 1 )

    !-- Set "edge" values

    associate &
      ( X => Profile ( :, iX_TS ), &  !-- cell outer edge
        D => Profile ( :, iD_TS ), &  !-- cell center 
        V => Profile ( :, iV_TS ), &  !-- cell center
        M => Profile ( :, iM_TS ), &  !-- cell center
        nProfile => size ( Profile, dim = 1 ) )

   ! R   = ( sqrt ( Kappa ) * G ** ( ( 1.0_KDR - Gamma ) / 2 ) &
   !            * t ** ( 2.0_KDR - Gamma ) ) * X
   ! Rho = D_P / ( G * t * t )
   ! v   = ( sqrt ( Kappa ) * G ** ( ( 1.0 - Gamma ) / 2 ) &
   !           * t ** ( 1.0_KDR - Gamma ) ) * V
    

    allocate ( X_edge ( nProfile + 1 ) )
    allocate ( dXC ( nProfile ) )
    allocate ( XC ( nProfile ) )
    X_edge ( 1 )          =  0.0_KDR
    do iV = 2, nProfile + 1
      X_edge  ( iV )  =  X ( iV - 1 )
      dXC ( iV - 1 )  =  X_edge ( iV )  -  X_edge ( iV - 1 )
      XC  ( iV - 1 )  =  X_edge ( iV - 1 )  +  0.5_KDR * dXC ( iV - 1 )
    end do

    allocate ( Density ( nProfile + 1 ) )
    allocate ( Velocity ( nProfile + 1 ) )
   ! allocate ( M ( nProfile + 1 ) )

    !-- First edge extrapolated
    Slope_D      =  ( D ( 2 )  -  D ( 1 ) )  &
                    /  ( 0.5_KDR * ( dXC ( 1 )  +  dXC ( 2 ) ) )
    Slope_V     =  ( V ( 2 )  -  V ( 1 ) )  &
                    /  ( 0.5_KDR * ( dVC ( 1 )  +  dVC ( 2 ) ) )
    !Slope_M    =  ( M_P ( 2 )  -  M_P ( 1 ) )  &
    !                /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )

    D_edge ( 1 )  =  Rho ( 1 ) +  Slope_D  * ( Radius ( 1 )  -  RC ( 1 ) )
    V_edge ( 1 ) =  v   ( 1 ) +  Slope_V  * ( Radius ( 1 )  -  RC ( 1 ) )
  !  M ( 1 )  =  M_P ( 1 )  +  Slope_SF   * ( Radius ( 1 )  -  RC ( 1 ) )

    do iV = 2, nProfile + 1

      if ( iV <= nProfile ) then
        Slope_D      =  ( Rho ( iV )  - Rho ( iV - 1 ) )  &
                        /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_V     =  ( v ( iV )  -  v ( iV - 1 ) )  &
                        /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
       ! Slope_M    =  ( M_P ( iV )  -  M_P ( iV - 1 ) )  &
       !                 /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
      else
        !-- Last edge extrapolated with same slope
      end if

      Density ( iV )  =  Rho ( iV - 1 )  &
                         +  Slope_D * ( Radius ( iV )  -  RC ( iV - 1 ) )
      Velocity ( iV ) =  v ( iV - 1 ) &
                         +  Slope_V   * ( Radius ( iV )  -  RC ( iV - 1 ) )
    !  M ( iV ) =  M_P ( iV -1 ) &
    !              +  Slope_SF   * ( Radius ( iV )  -  RC ( iV - 1 ) )
    end do

    end associate !-- R, etc.

    call Show ( 'First few values' )
    call Show ( Profile ( 1 : 5, iX_TS ), 'XTable' )
    call Show ( Radius ( 1 : 5 ), 'RadiusEdge' )
    call Show ( Profile ( 1 : 5, iD_TS ), 'DimensionlessDensityTable' )
    call Show ( Density ( 1 : 5 ), 'DimensionlessDensityEdge' )

    !-- SplineInterpolation initialization

    call SI ( iD_SI ) % Initialize &
           ( Radius, Density, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iV_SI ) % Initialize &
           ( Radius, Velocity, VerbosityOption = CONSOLE % INFO_3 )
   ! call SI ( iM_SI ) % Initialize &
   !        ( Radius, M, VerbosityOption = CONSOLE % INFO_3 )

  end subroutine PrepareInterpolation


end module YahilLattimer_Form
