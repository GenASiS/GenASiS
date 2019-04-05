module RadiationCentralCore_Form

  use Basics
  use Mathematics
  use FluidCentralCore_Form

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: RadiationCentralCoreForm
    character ( LDL ) :: &
      MomentsType = ''
  contains
    procedure, private, pass :: &
      Initialize_RCC
    generic, public :: &
      Initialize => Initialize_RCC
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_RCC
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_RCC
  end type RadiationCentralCoreForm

contains


  subroutine Initialize_RCC &
               ( RCC, RadiationName, RadiationType, MomentsType, GeometryType, &
                 Name, FinishTimeOption, CourantFactorOption, &
                 GravityFactorOption, LimiterParameterOption, &
                 ShockThresholdOption, RadiusMaxOption, RadiusCoreOption, &
                 RadialRatioOption, nCellsPolarOption, nWriteOption )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      MomentsType, &
      GeometryType, &
      Name
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      GravityFactorOption, &
      LimiterParameterOption, &
      ShockThresholdOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( RCC % Type == '' ) &
      RCC % Type = 'a RadiationCentralCore'
    
    call RCC % InitializeTemplate ( Name )

    RCC % MomentsType  =  MomentsType

    call RCC % AllocateIntegrator &
           ( RadiationName )
    call RCC % InitializePositionSpace &
           ( GeometryType, RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )

  end subroutine Initialize_RCC


  impure elemental subroutine Finalize ( RCC )

    type ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC

  end subroutine Finalize


  subroutine AllocateIntegrator_RCC ( RCC, RadiationName )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName

    integer ( KDI ) :: &
      iC  !-- iCurrent

    select case ( trim ( RCC % MomentsType ) )
    case ( 'NONE' )
      allocate ( Integrator_C_PS_Form :: RCC % Integrator )
    case ( 'GREY' )
      allocate ( Integrator_C_1D_PS_C_PS_Form :: RCC % Integrator )
    case ( 'SPECTRAL' )
      allocate ( Integrator_C_1D_MS_C_PS_Form :: RCC % Integrator )
    case default
      call Show ( 'MomentsType not recognized', CONSOLE % ERROR )
      call Show ( RCC % MomentsType, 'MomentsType', CONSOLE % ERROR )
      call Show ( 'RadiationCentralCore_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator_RCC', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- MomentsType

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

      I % N_CURRENTS_1D  =  size ( RadiationName )
      allocate ( I % TimeStepLabel &
                   ( I % N_CURRENTS_1D  +  1  +  1  +  I % N_CURRENTS_1D ) )

      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( iC )  &
          =  trim ( RadiationName ( iC ) ) // ' Streaming'
      end do !-- iC

      I % TimeStepLabel ( I % N_CURRENTS_1D  +  1 )  &
          =  'Fluid Advection'

      I % TimeStepLabel ( I % N_CURRENTS_1D  +  1  +  1 )  &
          =  'Gravity'

      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( I % N_CURRENTS_1D  +  1  +  iC )  &
          =  trim ( RadiationName ( iC ) ) // ' Interactions'
      end do !-- iC

    end select !-- I

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )
      allocate ( I % Current_ASC_1D ( I % N_CURRENTS_1D ) )
    class is ( Integrator_C_1D_MS_C_PS_Form )
      allocate ( I % Current_BSLL_ASC_CSLD_1D ( I % N_CURRENTS_1D ) )
    end select !-- I

  end subroutine AllocateIntegrator_RCC


end module RadiationCentralCore_Form
