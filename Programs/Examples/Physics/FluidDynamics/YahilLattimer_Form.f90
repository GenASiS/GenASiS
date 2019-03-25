module YahilLattimer_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: YahilLattimerForm
  !   real ( KDR ) :: &
  !     AdiabaticIndex, &
  !     CentralDensity, &
  !     CentralPressure, &
  !     Kappa, &
  !     t_collapse, &
  !     t_end
  !   real ( KDR ), dimension ( :, : ), allocatable :: &
  !     AnalyticProfile
  !   type ( Fluid_ASC_Form ), allocatable :: &
  !     Reference, &
  !     Difference
  contains
    procedure, private, pass :: &
      Initialize_YL
    generic, public :: &
      Initialize => Initialize_YL
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type YahilLattimerForm

    private :: &
      SetReference, &
      InitializeFluidCentralCore

!     private :: &
!       SetFluid, &
!       ReadTable, &
!       PrepareInterpolation, &

! !      private :: &
! !        SetFluidKernel

!     integer ( KDI ), private, parameter :: &
!       iX_TS  = 1, &  !-- must match the profile file columns
!       iD_TS  = 2, &
!       iV_TS  = 3, &
!       iM_TS  = 4
!     integer ( KDI ), private, parameter :: &
!       iRho_SI  = 1, & !-- spline interpolation
!       iV_SI    = 2

contains


  subroutine Initialize_YL ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ), target :: &
      YL
    character ( * ), intent ( in ) :: &
      Name

!     class ( Fluid_P_I_Form ), pointer :: &
!       F, &
!       F_R, &  !-- F_Reference
!       F_D     !-- F_Difference
!     real ( KDR ) :: &
!       Rho_final, &
!       FinishTime
!     type ( MeasuredValueForm ), dimension ( 3 ) :: &
!       CoordinateUnit
!     class ( GeometryFlatForm ), pointer :: &
!       G

    if ( YL % Type == '' ) &
      YL % Type = 'a YahilLattimer'

!     !--Set default parameters, read parameters
    
!     YL % AdiabaticIndex  = 1.3_KDR
!     YL % CentralDensity  = 7e9_KDR * UNIT % MASS_DENSITY_CGS
!     YL % CentralPressure = 6e27_KDR * UNIT % BARYE
    
!     call PROGRAM_HEADER % GetParameter ( YL % AdiabaticIndex, &
!                                          'AdiabaticIndex' )
!     call PROGRAM_HEADER % GetParameter ( YL % CentralDensity, &
!                                          'CentralDensity' )
!     call PROGRAM_HEADER % GetParameter ( YL % CentralPressure, &
!                                          'CentralPressure' )

!     call Show &
!          ( YL % CentralDensity, UNIT % MASS_DENSITY_CGS, 'CentralDensity' )
!     call Show &
!          ( YL % CentralPressure, UNIT % BARYE, 'CentralPressure' )

!     YL % Kappa= YL % CentralPressure &
!                   / ( YL % CentralDensity ** YL % AdiabaticIndex )

!     YL % t_collapse = 0.2_KDR * UNIT % SECOND
!     call PROGRAM_HEADER % GetParameter ( YL % t_collapse, 't_collapse' )
!     call Show ( YL % t_collapse, UNIT % SECOND, 't_collapse' )

!     Rho_final = 1e14_KDR * UNIT % MASS_DENSITY_CGS
!     call PROGRAM_HEADER % GetParameter ( Rho_final, 'Rho_final' )
!     call Show ( Rho_final, UNIT % MASS_DENSITY_CGS, 'Rho_final' )

!     call ReadTable ( YL, YL % AnalyticProfile, Rho_final )
!     call PROGRAM_HEADER % GetParameter ( YL % t_end, 't_end' )

!     FinishTime = YL % t_collapse - YL % t_end
!     Call Show ( FinishTime, UNIT % SECOND, 'FinishTime')

!     !-- Initial Conditions

!     F => FA % Fluid_P_I ( )
!     call F % SetAdiabaticIndex ( YL % AdiabaticIndex )

!     call SetFluid ( YL, F, Time = 0.0_KDR )


!     !-- Diagnostics

!     CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

!     allocate ( YL % Reference )
!     allocate ( YL % Difference )
!     call YL % Reference % Initialize &
!            ( PS, 'IDEAL', NameShortOption = 'Reference', &
!                Velocity_U_UnitOption &
!                  =  CoordinateUnit / UNIT % SECOND, &
!                BaryonMassUnitOption &
!                  =  UNIT % ATOMIC_MASS_UNIT, &
!                NumberDensityUnitOption &
!                  =  UNIT % FEMTOMETER ** ( -3 ), &
!                EnergyDensityUnitOption &
!                  =  UNIT % MEGA_ELECTRON_VOLT  &
!                     *  UNIT % FEMTOMETER ** ( -3 ), &
!                BaryonMassReferenceOption = CONSTANT % ATOMIC_MASS_UNIT, &
!              AllocateSourcesOption = .false., &
!              IgnorabilityOption = CONSOLE % INFO_2 )
!     call YL % Difference % Initialize &
!            ( PS, 'IDEAL', NameShortOption = 'Difference', &
!               ! Velocity_U_UnitOption &
!               !   =  CoordinateUnit / UNIT % SECOND, &
!                BaryonMassUnitOption &
!                  =  UNIT % ATOMIC_MASS_UNIT, &
!               ! NumberDensityUnitOption &
!               !   =  UNIT % FEMTOMETER ** ( -3 ), &
!               ! EnergyDensityUnitOption &
!               !   =  UNIT % MEGA_ELECTRON_VOLT  &
!               !      *  UNIT % FEMTOMETER ** ( -3 ), &
!                BaryonMassReferenceOption = CONSTANT % ATOMIC_MASS_UNIT, &
!              AllocateSourcesOption = .false., &
!              IgnorabilityOption = CONSOLE % INFO_2 )

!     F_R => YahilLattimer % Reference % Fluid_P_I ( )
!     call F_R % SetAdiabaticIndex ( YL % AdiabaticIndex )

!     !-- Cleanup

!     end select !-- FA
!     end select !-- PS
!     end select !-- FCC
!     nullify ( F, F_R )

  end subroutine Initialize_YL


  impure elemental subroutine Finalize ( YL )
    
    type ( YahilLattimerForm ), intent ( inout ) :: &
      YL

  end subroutine Finalize


  subroutine SetReference ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

!     class ( Fluid_P_I_Form ), pointer :: &
!       F, &
!       F_R, &  !-- F_Reference
!       F_D     !-- F_Difference
!     integer ( KDI ) :: &
!       iV
!     real ( KDR ) :: &
!       R_Max

!     select type ( FCC )
!     class is ( FluidCentralCoreForm )
!     select type ( FA => FCC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_P_I ( )

!     F_R => YahilLattimer % Reference % Fluid_P_I ( )
!     call SetFluid ( YahilLattimer, F_R, FCC % Time )

!     F_D => YahilLattimer % Difference % Fluid_P_I ( )
!     call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

!     F_D % Value = abs ( F_D % Value &
!                           / max ( abs ( F_R % Value ), &
!                                     sqrt ( tiny ( 0.0_KDR ) ) ) )

!     end select !-- FA
!     end select !-- FCC
!     nullify ( F, F_R, F_D )

  end subroutine SetReference


  subroutine InitializeFluidCentralCore ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    character ( * ), intent ( in )  :: &
      Name

    call YL % Initialize &
           ( FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', Name = Name )
    YL % Integrator % SetReference => SetReference

  end subroutine InitializeFluidCentralCore


!   subroutine SetFluid ( YL, F, Time )

!     class ( YahilLattimerForm ), intent ( inout ) :: &
!       YL
!     class ( Fluid_P_I_Form ), intent ( inout ) :: &
!       F
!     real ( KDR ), intent ( in ) :: &
!       Time

!     class ( GeometryFlatForm ), pointer :: &
!       G
!     integer ( KDI ) :: &
!       iV
!     real ( KDR ) :: &
!       T, &
!       R_TM
!     real ( KDR ), dimension ( : ), allocatable :: &
!       R, &
!       Rho, &
!       V
!     type ( SplineInterpolationForm ), dimension ( 2 ) :: &
!       SI

!     !-- Interpolate self similar solution

!     T = YL % t_collapse - Time

!     call PrepareInterpolation &
!            ( SI, YL % AnalyticProfile, YL % Kappa, &
!              YL % AdiabaticIndex, T )

!     select type ( PS => YL % Integrator % PositionSpace )
!     class is ( Atlas_SC_Form )
!     G => PS % Geometry ( )

!     select type ( Grid => PS % Chart )
!     class is ( Chart_SL_Template )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_Form )

!     associate &
!       ( R    => G % Value ( :, G % CENTER_U ( 1 ) ), &
!         Rho  => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
!         V_1  => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
!         V_2  => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
!         V_3  => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
!         E    => F % Value ( :, F % INTERNAL_ENERGY ), &
!         P    => F % Value ( :, F % PRESSURE ) )

!     V_2 = 0.0_KDR
!     V_3 = 0.0_KDR

!       do iV = 1, size ( R )
!         if ( R ( iV ) < 0.0_KDR ) &
!              cycle
!         call SI ( iRHO_SI ) &
!                 % Evaluate ( R ( iV ), Rho ( iV ) )
!         call SI ( iV_SI ) &
!                 % Evaluate ( R ( iV ), V_1 ( iV ) )
!         P ( iV ) = YL % Kappa * Rho ( iV ) ** ( YL % AdiabaticIndex )
!         E ( iV ) = P ( iV ) / ( YL % AdiabaticIndex - 1.0_KDR )
!       end do
      
!       Rho = Rho / CONSTANT % ATOMIC_MASS_UNIT

!       call F % ComputeFromPrimitive ( G )

!     end associate !-- R, etc.
!     end select    !-- C
!     end select    !-- Grid
!     end select    !-- PS

!     nullify ( G )

!   end subroutine SetFluid


!   subroutine ReadTable ( YL, AP, Rho_final )
    
!     class ( YahilLattimerForm ), intent ( inout ) :: &
!       YL
!     real ( KDR ), dimension ( :, : ), allocatable, intent ( inout ) :: &
!       AP
!     real ( KDR ), intent ( in ) :: &
!       Rho_final
   
!     character ( LDF ) :: &
!       Path, &
!       FileName
!     type ( TableStreamForm ) :: &
!       TS

!     Path = '../Parameters/'
!     FileName = 'YahilLattimerCollapse_Gm_130.dat'
    
!     call PROGRAM_HEADER % GetParameter &
!            ( FileName, 'FileName' )
!     call Show ( FileName, 'FileName' )

!     call TS % Initialize &
!            ( Filename, PROGRAM_HEADER % Communicator % Rank, &
!              PathOption = Path )
!     call TS % Read ( AP, oRowOption = 1 )

!     YL % t_end &
!          = sqrt ( AP ( 1, iD_TS ) / ( Rho_final * CONSTANT % GRAVITATIONAL ) )

!   end subroutine ReadTable


!   subroutine PrepareInterpolation ( SI, AP, Kappa, Gamma, T )

!     type ( SplineInterpolationForm ), dimension ( 2 ), intent ( inout ) :: &
!       SI
!     real ( KDR ), dimension ( :, : ), intent ( in ) :: &
!       AP
!     real ( KDR ), intent ( in ) :: &
!       Kappa, &
!       Gamma, &
!       T

!     integer ( KDI ) :: &
!       iV
!     real ( KDR ) :: &
!       Slope_Rho, &
!       Slope_V
!     real ( KDR ), dimension ( : ), allocatable :: &
!       R, &   !-- RCenter
!       R_Edge, &
!       dR, &  !-- CellWidth
!       Rho, &
!       V

!     associate &
!       ( X   => AP ( :, iX_TS ), &  !-- cell center
!         D   => AP ( :, iD_TS ), &  !-- cell center 
!         V_P => AP ( :, iV_TS ), &  !-- cell center
!         M   => AP ( :, iM_TS ), &  !-- cell center
!         nProfile => size ( AP, dim = 1 ) )
    
!     allocate ( R      ( nProfile + 1 ), &
!                Rho    ( nProfile + 1 ), &
!                V      ( nProfile + 1 ) )
    
!     R ( 2 : nProfile + 1 ) &
!       = ( sqrt ( Kappa ) * CONSTANT % GRAVITATIONAL &
!                              ** ( ( 1.0_KDR - Gamma ) / 2 ) &
!          * T ** ( 2.0_KDR - Gamma ) ) * X 

!     Rho ( 2 : nProfile + 1 ) &
!       = D / ( CONSTANT % GRAVITATIONAL * T * T ) 

!     V ( 2 : nProfile + 1 ) &
!       = ( sqrt ( Kappa ) * CONSTANT % GRAVITATIONAL &
!                              ** ( ( 1.0 - Gamma ) / 2 ) &
!           * T ** ( 1.0_KDR - Gamma ) ) * V_P 

!     R ( 1 ) = 0.0_KDR
    
!     Slope_Rho = ( Rho ( 3 ) - Rho ( 2 ) ) / ( R ( 3 ) - R ( 2 ) )
!     Slope_V   = ( V   ( 3 ) - V   ( 2 ) ) / ( R ( 3 ) - R ( 2 ) )

!      Rho ( 1 ) = Rho ( 2 ) - Slope_Rho * R ( 2 )
!      V   ( 1 ) = V   ( 2 ) - Slope_V   * R ( 2 ) 

!    !-- SplineInterpolation initialization

!     call SI ( iRho_SI ) % Initialize &
!            ( R, Rho, VerbosityOption = CONSOLE % INFO_3 )
!     call SI ( iV_SI ) % Initialize &
!            ( R, V, VerbosityOption = CONSOLE % INFO_3 )

!     end associate !-- R, etc

!   end subroutine PrepareInterpolation


  subroutine ComputeError ( YL )

    class ( YahilLattimerForm ) , intent ( in ) :: &
      YL
    
!     real ( KDR ) :: &
!       L1_Rho, &
!       L1_V, &
!       L1_P
!     class ( Fluid_P_I_Form ), pointer :: &
!       F, &
!       F_D, &
!       F_R
!     type ( CollectiveOperation_R_Form ) :: &
!       CO

!     select type ( FCC => YL % Integrator )
!     type is ( FluidCentralCoreForm )
!     select type ( PS => FCC % PositionSpace )
!     class is ( Atlas_SC_Form ) 
!     select type ( C => PS % Chart )
!     class is ( Chart_SL_Template )
!     select type ( FA => FCC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_P_I ( )

!     F_R => YL % Reference % Fluid_P_I ( )
!     call SetFluid ( YahilLattimer, F_R, FCC % Time )

!     F_D => YL % Difference % Fluid_P_I ( )
!     call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )
       
    
!     associate &
!       ( Difference_Rho &
!           => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ), &
!         Reference_Rho &
!           => F_R % Value ( :, F_R % COMOVING_BARYON_DENSITY ), &
!         Difference_V   &
!           => F_D % Value ( :, F_D % VELOCITY_U ( 1 ) ), &
!         Reference_V   &
!           => F_R % Value ( :, F_R % VELOCITY_U ( 1 ) ), &
!         Difference_P   &
!           => F_D % Value ( :, F_D % PRESSURE ), &
!         Reference_P   &
!           => F_R % Value ( :, F_R % PRESSURE ) )

!     call CO % Initialize ( PS % Communicator, [ 6 ], [ 6 ] )

!     CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference_Rho ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 2 ) = sum ( abs ( Reference_Rho ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 3 ) = sum ( abs ( Difference_V ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 4 ) = sum ( abs ( Reference_V ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 5 ) = sum ( abs ( Difference_P ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 6 ) = sum ( abs ( Reference_P ), &
!                                         mask = C % IsProperCell )

!     call CO % Reduce ( REDUCTION % SUM )

!     end associate !-- Difference_Rho, etc.
    
!     associate &
!       ( DifferenceSum_Rho  => CO % Incoming % Value ( 1 ), &
!         ReferenceSum_Rho   => CO % Incoming % Value ( 2 ), &
!         DifferenceSum_V    => CO % Incoming % Value ( 3 ), &
!         ReferenceSum_V     => CO % Incoming % Value ( 4 ), &
!         DifferenceSum_P    => CO % Incoming % Value ( 5 ), &
!         ReferenceSum_P     => CO % Incoming % Value ( 6 ) )

!     L1_Rho = DifferenceSum_Rho / ReferenceSum_Rho
!     L1_V   = DifferenceSum_V   / ReferenceSum_V
!     L1_P   = DifferenceSum_P   / ReferenceSum_P
!     end associate

!     call Show ( L1_Rho, '*** L1_Rho error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )
!     call Show ( L1_V, '*** L1_V error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )
!     call Show ( L1_P, '*** L1_P error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )

!     end select !-- FA
!     end select !-- C
!     end select !-- PS
!     end select !--FCC

!     nullify ( F, F_R, F_D )

  end subroutine ComputeError


end module YahilLattimer_Form
