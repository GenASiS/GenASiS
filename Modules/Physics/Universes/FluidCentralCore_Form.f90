module FluidCentralCore_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use FluidCentral_Template

  implicit none
  private

  type, public, extends ( FluidCentralTemplate ) :: FluidCentralCoreForm
  !   real ( KDR ) :: &
  !     GravityFactor
  contains
  !   procedure, public, pass :: &
  !     Initialize
  !   final :: &
  !     Finalize
    procedure, private, pass :: &
      InitializeAtlas
  !   procedure, private, pass :: &
  !     InitializeGeometry
  !   procedure, private, pass :: &
  !     SetCoarsening
  !   procedure, public, nopass :: &
  !     CoarsenSingularities
  !   procedure, private, pass :: &
  !     SetWriteTimeInterval
  !   procedure, public, pass :: &
  !     PrepareCycle
  !   procedure, public, pass :: &
  !     ComputeTimeStep_G_ASC
  end type FluidCentralCoreForm

    ! class ( FluidCentralCoreForm ), private, pointer :: &
    !   FluidCentralCore => null ( )

    !   private :: &
    !     ComputeTimeStepLocal, &
    !     LocalMax, &
    !     ComputeTimeStep_G_CSL

contains


!   subroutine Initialize &
!                ( FCC, Name, FluidType, GeometryType, &
!                  DimensionlessOption, TimeUnitOption, FinishTimeOption, &
!                  CourantFactorOption, GravityFactorOption, &
!                  LimiterParameterOption, ShockThresholdOption, &
!                  RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
!                  nWriteOption, nCellsPolarOption )

!     class ( FluidCentralCoreForm ), intent ( inout ), target :: &
!       FCC
!     character ( * ), intent ( in )  :: &
!       Name, &
!       FluidType, &
!       GeometryType
!     logical ( KDL ), intent ( in ), optional :: &
!       DimensionlessOption
!     type ( MeasuredValueForm ), intent ( in ), optional :: &
!       TimeUnitOption
!     real ( KDR ), intent ( in ), optional :: &
!       FinishTimeOption, &
!       CourantFactorOption, &
!       GravityFactorOption, &
!       LimiterParameterOption, &
!       ShockThresholdOption, &
!       RadiusMaxOption, &
!       RadiusCoreOption, &
!       RadialRatioOption
!     integer ( KDI ), intent ( in ), optional :: &
!       nCellsPolarOption, &
!       nWriteOption

!     if ( FCC % Type == '' ) &
!       FCC % Type = 'a FluidCentralCore'

!     FluidCentralCore => FCC

!     if ( any ( trim ( GeometryType ) &
!                  == [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
!     then

!       if ( .not. allocated ( FCC % TimeStepLabel ) ) &
!         allocate ( FCC % TimeStepLabel ( 2 ) )
!       FCC % TimeStepLabel ( 1 )  =  'Fluid'
!       FCC % TimeStepLabel ( 2 )  =  'Gravity'

!       FCC % GravityFactor = 0.7_KDR
!       if ( present ( GravityFactorOption ) ) &
!         FCC % GravityFactor = GravityFactorOption
!       call PROGRAM_HEADER % GetParameter &
!              ( FCC % GravityFactor, 'GravityFactor' )

!     end if

!     call FCC % InitializeTemplate_FC &
!                ( Name, FluidType, GeometryType, &
!                  DimensionlessOption = DimensionlessOption, &
!                  TimeUnitOption = TimeUnitOption, &
!                  FinishTimeOption = FinishTimeOption, &
!                  CourantFactorOption = CourantFactorOption, &
!                  LimiterParameterOption = LimiterParameterOption, &
!                  ShockThresholdOption = ShockThresholdOption, &
!                  RadiusMaxOption = RadiusMaxOption, &
!                  RadiusCoreOption = RadiusCoreOption, &
!                  RadialRatioOption = RadialRatioOption, &
!                  nWriteOption = nWriteOption, &
!                  nCellsPolarOption = nCellsPolarOption )

!     FCC % ComputeTimeStepLocal => ComputeTimeStepLocal

!     select type ( S => FCC % Step )
!     class is ( Step_RK2_C_ASC_Form )
!       S % CoarsenSingularities => CoarsenSingularities
!     end select !-- S

!     if ( any ( trim ( GeometryType ) &
!                  == [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
!       call Show ( FCC % GravityFactor, 'GravityFactor' )

!   end subroutine Initialize


  subroutine InitializeAtlas &
               ( FC, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

    class ( FluidCentralCoreForm ), intent ( inout ) :: &
      FC
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusCore, &
      RadialRatio

    associate ( I => FC % Integrator )

    allocate ( Atlas_SC_CC_Form :: I % PositionSpace )
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_CC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( FC % Dimensionless ) then

      call PS % CreateChart_CC ( )

    else

      RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
      RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
      RadialRatio  =  2.4_KDR

      call PS % CreateChart_CC &
             ( CoordinateUnitOption = FC % Units % Coordinate_PS, &
               RadiusCoreOption = RadiusCore, &
               RadiusMaxOption = RadiusMax, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

      FC % RadiusPolarMomentum  =  8.0_KDR  *  UNIT % KILOMETER
      call PROGRAM_HEADER % GetParameter &
             ( FC % RadiusPolarMomentum, 'RadiusPolarMomentum' )

    end if !-- Dimensionless

    end select !-- PS
    end associate !-- I

  end subroutine InitializeAtlas


!   subroutine InitializeGeometry &
!                ( FC, GA, PS, GeometryType, CentralMassOption )

!     class ( FluidCentralCoreForm ), intent ( inout ) :: &
!       FC
!     type ( Geometry_ASC_Form ), intent ( inout ) :: &
!       GA
!     class ( Atlas_SC_Form ), intent ( in ) :: &
!       PS
!     character ( * ), intent ( in )  :: &
!       GeometryType
!     real ( KDR ), intent ( in ), optional :: &
!       CentralMassOption

!     if ( FC % Dimensionless ) then
!       call GA % Initialize &
!              ( PS, GeometryType, GravitySolverTypeOption = 'MULTIPOLE', &
!                GravitationalConstantOption = 1.0_KDR )
!     else
!       call GA % Initialize &
!              ( PS, GeometryType, GravitySolverTypeOption = 'MULTIPOLE' )
!     end if !-- Dimensionless

!   end subroutine InitializeGeometry


!   impure elemental subroutine Finalize ( FCC )

!     type ( FluidCentralCoreForm ), intent ( inout ) :: &
!       FCC

!     call FCC % FinalizeTemplate_FC ( )

!   end subroutine Finalize


!   subroutine SetCoarsening ( FC )

!     class ( FluidCentralCoreForm ), intent ( inout ) :: &
!       FC

!     select type ( PS => FC % PositionSpace )
!     class is ( Atlas_SC_CC_Form )

!     !-- Polar coarsening

!     if ( PS % nDimensions < 2 ) &
!       return
!     call Show ( 'SetCoarsening_2', FC % IGNORABILITY )
!     call Show ( FC % Name, 'Universe', FC % IGNORABILITY )
!     call FC % SetCoarseningTemplate ( iAngular = 2 )

!     !-- Azimuthal coarsening

!     if ( PS % nDimensions < 3 ) &
!       return
!     call Show ( 'SetCoarsening_3', FC % IGNORABILITY )
!     call Show ( FC % Name, 'Universe', FC % IGNORABILITY )
!     call FC % SetCoarseningTemplate ( iAngular = 3 )

!     end select !-- PS

!   end subroutine SetCoarsening


!   subroutine CoarsenSingularities ( S )

!     class ( StorageForm ), intent ( inout ) :: &
!       S

!     associate ( FCC => FluidCentralCore )
!     select type ( PS => FCC % PositionSpace )
!     class is ( Atlas_SC_CC_Form )

!     select case ( mod ( FCC % iCycle, 2 ) )
!     case ( 0 )
!       if ( PS % nDimensions > 2 ) &
!         call FCC % CoarsenSingularityTemplate ( S, iAngular = 3 )
!       if ( PS % nDimensions > 1 ) &
!         call FCC % CoarsenSingularityTemplate ( S, iAngular = 2 )
!     case ( 1 )
!       if ( PS % nDimensions > 1 ) &
!         call FCC % CoarsenSingularityTemplate ( S, iAngular = 2 )
!       if ( PS % nDimensions > 2 ) &
!         call FCC % CoarsenSingularityTemplate ( S, iAngular = 3 )
!     end select
  
!     end select !-- PS
!     end associate !-- FCC

!   end subroutine CoarsenSingularities


!   subroutine SetWriteTimeInterval ( I )

!     class ( FluidCentralCoreForm ), intent ( inout ) :: &
!       I

!     integer ( KDI ) :: &
!       iProcess, &
!       iRadius
!     real ( KDR ) :: &
!       VelocityMax, &
!       VelocityMaxRadius, &
!       DensityAve, &
!       TimeScaleDensityAve, &
!       TimeScaleVelocityMax
!     type ( CollectiveOperation_R_Form ), allocatable :: &
!       CO
!     class ( GeometryFlatForm ), pointer :: &
!       G
!     class ( Fluid_D_Form ), pointer :: &
!       F

!     select type ( FA => I % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_D ( )

!     select type ( PS => I % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( GA => PS % Geometry_ASC )
!     class is ( Geometry_ASC_Form )

!     G => PS % Geometry ( )

!     select type ( Chart => PS % Chart )
!     class is ( Chart_SL_Template )

!     associate ( C => PS % Communicator ) 

!     !-- Find max velocity
!     allocate ( CO )
!     call CO % Initialize ( C, nOutgoing = [ 1 ], nIncoming = [ C % Size ] )
!     CO % Outgoing % Value ( 1 ) &
!       =  LocalMax ( Chart % IsProperCell, &
!                     abs ( F % Value ( :, F % VELOCITY_U ( 1 ) ) ) ) 
!     call CO % Gather ( )
!     VelocityMax  =  maxval ( CO % Incoming % Value )
!     iProcess     =  maxloc ( CO % Incoming % Value, dim = 1 )  -  1
!     deallocate ( CO )

!     if ( VelocityMax == 0.0_KDR ) &
!       return

!     !-- Find radius of max velocity
!     allocate ( CO )
!     call CO % Initialize &
!            ( C, nOutgoing = [ 1 ], nIncoming = [ 1 ], RootOption = iProcess )
!     if ( C % Rank == iProcess ) then
!       iRadius  =  maxloc ( abs ( F % Value ( :, F % VELOCITY_U ( 1 ) ) ), &
!                            dim = 1, mask = Chart % IsProperCell )
!       CO % Outgoing % Value ( 1 )  =  G % Value ( iRadius, G % CENTER_U ( 1 ) )
!     end if
!     call CO % Broadcast ( )
!     VelocityMaxRadius = CO % Incoming % Value ( 1 )
!     deallocate ( CO )

!     !-- Compute average density
!     select type ( TI => FA % TallyInterior )
!     class is ( Tally_F_D_Form )
!     DensityAve  =  F % BaryonMassReference &
!                    * TI % Value ( TI % BARYON_NUMBER ) &
!                    / ( 4.0 / 3.0  *  CONSTANT % PI  *  VelocityMaxRadius ** 3 )

!     !-- Time scales
!     TimeScaleVelocityMax &
!       =  VelocityMaxRadius  /  VelocityMax
!     TimeScaleDensityAve &
!       =  ( GA % GravitationalConstant  *  DensityAve ) ** ( -0.5_KDR )

!     I % WriteTimeInterval  &
!       =  min ( TimeScaleDensityAve, TimeScaleVelocityMax )  /  I % nWrite

!     call Show ( 'Time Scales', I % IGNORABILITY )
!     call Show ( VelocityMax, Chart % CoordinateUnit ( 1 ) / I % TimeUnit, &
!                 'VelocityMax', I % IGNORABILITY )
!     call Show ( VelocityMaxRadius, Chart % CoordinateUnit ( 1 ), &
!                 'VelocityMaxRadius', I % IGNORABILITY )
!     call Show ( DensityAve, UNIT % IDENTITY, 'DensityAve', &
!                 I % IGNORABILITY )
!     call Show ( TimeScaleDensityAve, I % TimeUnit, 'TimeScaleDensityAve', &
!                 I % IGNORABILITY )
!     call Show ( TimeScaleVelocityMax, I % TimeUnit, 'TimeScaleVelocityMax', &
!                 I % IGNORABILITY )

!     !-- Cleanup
!     end select !-- TI
!     end associate !-- C
!     end select !-- Chart
!     end select !-- GA
!     end select !-- PS
!     end select !-- FA
!     nullify ( G, F )

!   end subroutine SetWriteTimeInterval


!   subroutine PrepareCycle ( I )

!     class ( FluidCentralCoreForm ), intent ( inout ) :: &
!       I

!     ! class ( Fluid_D_Form ), pointer :: &
!     !   F

!     ! select type ( FA => I % Current_ASC )
!     ! class is ( Fluid_ASC_Form )

!     ! select type ( PS => I % PositionSpace )
!     ! class is ( Atlas_SC_Form )

!     ! select type ( GA => PS % Geometry_ASC )
!     ! class is ( Geometry_ASC_Form )

!     ! F => FA % Fluid_D ( )

!     ! call GA % ComputeGravity &
!     !       ( I % Current_ASC, &
!     !         iBaryonMass = F % BARYON_MASS, &
!     !         iBaryonDensity = F % COMOVING_BARYON_DENSITY )

!     ! end select !-- GA
!     ! end select !-- PS
!     ! end select !-- FA
!     ! nullify ( F )

!   end subroutine PrepareCycle


!   subroutine ComputeTimeStep_G_ASC ( FCC, TimeStepCandidate )

!     class ( FluidCentralCoreForm ), intent ( inout ), target :: &
!       FCC
!     real ( KDR ), intent ( inout ) :: &
!       TimeStepCandidate

!     class ( Geometry_N_Form ), pointer :: &
!       G

!     select type ( PS => FCC % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( GA => PS % Geometry_ASC )
!     class is ( Geometry_ASC_Form )

!     if ( all ( trim ( GA % GeometryType ) &
!                  /= [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
!       return

!     select type ( CSL => PS % Chart )
!     class is ( Chart_SL_Template )

!     G  =>  GA % Geometry_N ( )

!     call ComputeTimeStep_G_CSL &
!            ( CSL % IsProperCell, &
!              G % Value ( :, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
!              G % Value ( :, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
!              G % Value ( :, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
!              G % Value ( :, G % METRIC_UU_22 ), &
!              G % Value ( :, G % METRIC_UU_33 ), &
!              G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
!              G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), & 
!              G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
!              G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
!              G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), & 
!              G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
!              CSL % nDimensions, TimeStepCandidate )

!     TimeStepCandidate  =  FCC % GravityFactor * TimeStepCandidate

!     end select !-- CSL
!     end select !-- GA
!     end select !-- PS
!     nullify ( G )

!   end subroutine ComputeTimeStep_G_ASC


!   subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

!     class ( IntegratorTemplate ), intent ( inout ), target :: &
!       I
!     real ( KDR ), dimension ( : ), intent ( inout ) :: &
!       TimeStepCandidate

!     select type ( I )
!     class is ( FluidCentralCoreForm )

!     call I % ComputeTimeStep_C_ASC &
!            ( TimeStepCandidate ( 1 ), I % Current_ASC )

!     call I % ComputeTimeStep_G_ASC &
!            ( TimeStepCandidate ( 2 ) )

!     end select !-- I

!   end subroutine ComputeTimeStepLocal


!   function LocalMax ( IsProperCell, V ) result ( ML ) 

!     logical ( KDL ), dimension ( : ), intent ( in ) :: &
!       IsProperCell
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       V
!     real ( KDR ) :: &
!       ML

!     integer ( KDI ) :: &
!       iV

!     ML = - huge ( 0.0_KDR )
!     !$OMP parallel do private ( iV ) &
!     !$OMP reduction ( max : ML )
!     do iV = 1, size ( V )
!       if ( IsProperCell ( iV ) ) &
!         ML  =  max ( ML, V ( iV ) )
!     end do !-- iV
!     !$OMP end parallel do
 
!   end function LocalMax


!   subroutine ComputeTimeStep_G_CSL &
!                ( IsProperCell, GradPhi_1, GradPhi_2, GradPhi_3, &
!                  M_UU_22, M_UU_33, dXL_1, dXL_2, dXL_3, dXR_1, dXR_2, dXR_3, &
!                  nDimensions, TimeStep )

!     logical ( KDL ), dimension ( : ), intent ( in ) :: &
!       IsProperCell
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       GradPhi_1, GradPhi_2, GradPhi_3, &
!       M_UU_22, M_UU_33, &
!       dXL_1, dXL_2, dXL_3, &
!       dXR_1, dXR_2, dXR_3
!     integer ( KDI ), intent ( in ) :: &
!       nDimensions
!     real ( KDR ), intent ( inout ) :: &
!       TimeStep

!     integer ( KDI ) :: &
!       iV, &
!       nV
!     real ( KDR ) :: &
!       TimeStepInverse

!     nV = size ( dXL_1 )

!     select case ( nDimensions )
!     case ( 1 )
!       TimeStepInverse &
!         = maxval ( sqrt ( abs ( GradPhi_1 ) / ( dXL_1 + dXR_1 ) ), &
!                    mask = IsProperCell )
!     case ( 2 )
!       TimeStepInverse &
!         = maxval ( sqrt (    abs ( GradPhi_1 ) &
!                              / ( dXL_1 + dXR_1 ) &
!                           +  abs ( M_UU_22 * GradPhi_2 ) &
!                              / ( dXL_2 + dXR_2 ) ), &
!                    mask = IsProperCell )
!     case ( 3 )
!       ! TimeStepInverse &
!       !   = maxval ( sqrt (   abs ( GradPhi_1 ) / dX_1 &
!       !                     + abs ( M_UU_22 * GradPhi_2 ) / dX_2 &
!       !                     + abs ( M_UU_33 * GradPhi_3 ) / dX_3 ) ), &
!       !              mask = IsProperCell )
!       TimeStepInverse = - huge ( 0.0_KDR )
!       !$OMP parallel do private ( iV ) &
!       !$OMP reduction ( max : TimeStepInverse )
!       do iV = 1, nV
!         if ( IsProperCell ( iV ) ) &
!           TimeStepInverse &
!             = max ( TimeStepInverse, &
!                     sqrt (   abs ( GradPhi_1 ( iV ) ) &
!                                    / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
!                            + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
!                                    / ( dXL_2 ( iV ) + dXR_2 ( iV ) ) &
!                            + abs ( M_UU_33 ( iV ) * GradPhi_3 ( iV ) ) &
!                                    / ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
!       end do
!       !$OMP end parallel do
!     end select !-- nDimensions

!     TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
!     TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

!   end subroutine ComputeTimeStep_G_CSL


end module FluidCentralCore_Form
