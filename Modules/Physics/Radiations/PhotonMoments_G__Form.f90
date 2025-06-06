module PhotonMoments_G__Form

  !-- PhotonMoments_Grey__Form

  use Basics
  use Mathematics
  use Fluids
  use Units_R__Form
  use RadiationMoments_BM__Form
  use Interactions_BM__Form

  implicit none
  private

  integer ( KDI ), private, parameter :: &
      N_FIELDS_PM = 1

  type, public, extends ( RadiationMoments_BM_Form ) :: PhotonMoments_G_Form
    integer ( KDI ) :: &
      N_FIELDS_PM = N_FIELDS_PM
    integer ( KDI ) :: &
      TEMPERATURE_GREY = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    final :: &
      Finalize
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      ComputeSpectralParametersAll
    procedure, public, pass :: &
      ComputeSpectralParametersSingle
    generic, public :: &
      ComputeSpectralParameters => ComputeSpectralParametersAll, &
                                   ComputeSpectralParametersSingle
    procedure, public, pass :: &
      ComputeEquilibriumAll
    procedure, public, pass :: &
      ComputeEquilibriumSingle
  end type PhotonMoments_G_Form

    private :: &
      Compute_SP_A_Kernel, &
      Compute_SP_S_Kernel, &
      Compute_Eq_A_Kernel, &
      Compute_Eq_S_Kernel
    
    interface

      module subroutine Compute_SP_A_Kernel ( T_R, J, UseDeviceOption )
        !-- Compute_SpectralParameters_All_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          T_R
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SP_A_Kernel
 
      module subroutine Compute_SP_S_Kernel ( T_R, J, iV )
        !-- Compute_SpectralParameters_Single_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          T_R
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_SP_S_Kernel
 
      module subroutine Compute_Eq_A_Kernel &
               ( J_Eq, J_RD, J, T, UseDeviceOption )
        !-- Compute_Equilibrium_All_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J_Eq, J_RD
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J, &
          T
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_Eq_A_Kernel

      module subroutine Compute_Eq_S_Kernel &
               ( J_Eq, J_RD, J, T, iV )
        !-- Compute_Equilibrium_Single_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J_Eq, J_RD
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J, &
          T
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_Eq_S_Kernel

    end interface


contains


  subroutine InitializeAllocate_RM &
               ( RM, F, Units_R, RadiationType, FieldOption, VectorOption, &
                 TallyVariableOption, NameOption, UnitOption, TallyUnitOption, &
                 VectorIndicesOption, iaPrimitiveOption, iaBalancedOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      RM
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    character ( * ), intent ( in ) :: &
      RadiationType
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption, &
      TallyVariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      TallyUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      oF, &  !-- oField
      nFields
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RM % Type  ==  '' ) &
      RM % Type  =  'a PhotonMoments_G' 
    
    Name  =  'PhotonMoments'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  RM % N_FIELDS_CS  +  RM % N_FIELDS_RM

    RM % TEMPERATURE_GREY  =  oF + 1

    nFields  =  oF  +  RM % N_FIELDS_PM
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RM % N_FIELDS_PM ) &
      = [ 'TemperatureGrey' ]

    !-- Units

    associate ( nC  =>  F % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( RM % TEMPERATURE_GREY, iC ) &
        =  Units_R ( iC ) % Temperature
    end do !-- iC

    end associate !-- nC

    !-- RadiationMoments_BM

    call RM % RadiationMoments_BM_Form % Initialize &
           ( F, Units_R, RadiationType, &             
             FieldOption = Field, &
             VectorOption = VectorOption, &
             TallyVariableOption = TallyVariableOption, &
             NameOption = Name, &
             UnitOption = FieldUnit, &
             TallyUnitOption = TallyUnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitiveOption, &
             iaBalancedOption = iaBalancedOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_RM


  impure elemental subroutine Finalize ( PM )

    type ( PhotonMoments_G_Form ), intent ( inout ) :: &
      PM

  end subroutine Finalize


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( PhotonMoments_G_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % ENERGY_DENSITY_C, &
                    CS % ENERGY_DENSITY_C_EQ, &
                    CS % ENERGY_DENSITY_C_RD, &
                    CS % MOMENTUM_DENSITY_C_U, &
                    CS % FLUX_FACTOR, &
                    CS % STRESS_FACTOR, &
                    CS % STRESS_FACTOR_RD, &
                    CS % DIFFUSION_INDICATOR, &
                    CS % TEMPERATURE_GREY ] )

  end subroutine SetStream


  subroutine ComputeEquilibriumAll ( RM )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeEquilibriumAll', CONSOLE % INFO_6 )
    call Show ( RM % Name, 'PhotonMoments', CONSOLE % INFO_6 )

    select type ( I  =>  RM % Interactions )
      class is ( Interactions_BM_Form )
    associate &
      ( F  =>  I % Fluid )

    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RMV  =>  RM % Storage ( iC ) % Value, &
           FV  =>   F % Storage ( iC ) % Value )
      associate &
        ( J_Eq  =>  RMV ( :, RM % ENERGY_DENSITY_C_EQ ), &
          J_RD  =>  RMV ( :, RM % ENERGY_DENSITY_C_RD ), &
          J     =>  RMV ( :, RM % ENERGY_DENSITY_C ), &
          T     =>   FV ( :,  F % TEMPERATURE ) )
               
      call Compute_Eq_A_Kernel &
             ( J_Eq, J_RD, J, T, &
               UseDeviceOption = RM % DeviceMemory )

      end associate !-- T_R, etc.
      end associate !-- RV, etc.
    end do !-- iC

    end associate !-- F
    end select !-- I

  end subroutine ComputeEquilibriumAll


  subroutine ComputeEquilibriumSingle ( RM, iC, iV )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      RM
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

!    call Show ( 'ComputeEquilibriumAll', CONSOLE % INFO_6 )
!    call Show ( RM % Name, 'PhotonMoments', CONSOLE % INFO_6 )

    select type ( I  =>  RM % Interactions )
      class is ( Interactions_BM_Form )
    associate &
      ( F  =>  I % Fluid )

    associate &
      ( RMV  =>  RM % Storage ( iC ) % Value, &
         FV  =>   F % Storage ( iC ) % Value )
    associate &
      ( J_Eq  =>  RMV ( :, RM % ENERGY_DENSITY_C_EQ ), &
        J_RD  =>  RMV ( :, RM % ENERGY_DENSITY_C_RD ), &
        J     =>  RMV ( :, RM % ENERGY_DENSITY_C ), &
        T     =>   FV ( :,  F % TEMPERATURE ) )
             
    call Compute_Eq_S_Kernel &
           ( J_Eq, J_RD, J, T, iV )

    end associate !-- T_R, etc.
    end associate !-- RV, etc.

    end associate !-- F
    end select !-- I

  end subroutine ComputeEquilibriumSingle


  subroutine ComputeSpectralParametersAll ( RM )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeSpectralParametersAll', CONSOLE % INFO_6 )
    call Show ( RM % Name, 'PhotonMoments', CONSOLE % INFO_6 )

    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RMV  =>  RM % Storage ( iC ) % Value )
      associate &
        ( T_R  =>  RMV ( :, RM % TEMPERATURE_GREY ), &
          J    =>  RMV ( :, RM % ENERGY_DENSITY_C ) )
               
      call Compute_SP_A_Kernel &
             ( T_R, J, &
               UseDeviceOption  =  RM % DeviceMemory )

      end associate !-- T_R, etc.
      end associate !-- RV, etc.
    end do !-- iC

  end subroutine ComputeSpectralParametersAll


  subroutine ComputeSpectralParametersSingle ( RM, iC, iV )

    class ( PhotonMoments_G_Form ), intent ( inout ) :: &
      RM
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

!    call Show ( 'ComputeSpectralParametersAll', CONSOLE % INFO_6 )
!    call Show ( RM % Name, 'PhotonMoments', CONSOLE % INFO_6 )

    associate &
      ( RMV  =>  RM % Storage ( iC ) % Value )
    associate &
      ( T_R  =>  RMV ( :, RM % TEMPERATURE_GREY ), &
        J    =>  RMV ( :, RM % ENERGY_DENSITY_C ) )
             
    call Compute_SP_S_Kernel ( T_R, J, iV )

    end associate !-- T_R, etc.
    end associate !-- RV, etc.

  end subroutine ComputeSpectralParametersSingle


end module PhotonMoments_G__Form
