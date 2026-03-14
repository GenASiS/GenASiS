module Fluid_P_I__Form

  !-- Fluid_Perfect_Ideal__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_P__Form
  use Features_F_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_I = 0, &
      N_CONSERVED_I = 0, &
      N_FIELDS_I    = 0, &
      N_VECTORS_I   = 0

  type, public, extends ( Fluid_P_Form ) :: Fluid_P_I_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_I = N_PRIMITIVE_I, &
      N_CONSERVED_I = N_CONSERVED_I, &
      N_FIELDS_I    = N_FIELDS_I, &
      N_VECTORS_I   = N_VECTORS_I
    real ( KDR ) :: &
      BoltzmannConstant, &
      AdiabaticIndex, &
      MeanMolecularWeight, &
      SpecificHeatVolume, &  !-- per baryon
      FiducialBaryonDensity, &
      FiducialPressure
  contains
    procedure, public, pass :: &
      InitializeAllocate_F
    procedure, public, pass :: &
      SetAdiabaticIndex
    procedure, public, pass :: &
      SetMeanMolecularWeight
    procedure, public, pass :: &
      SetSpecificHeatVolume
    procedure, public, pass :: &
      SetFiducialParameters
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      ComputeFromTemperature
    procedure, public, pass ( CS ) :: &
      ComputeFromPrimitive
    procedure, public, pass :: &
      ComputeFromBalancedAll
    procedure, public, pass :: &
      ComputeFromBalancedSingle
    final :: &
      Finalize
  end type Fluid_P_I_Form

    private :: &
      Apply_EOS_I_T_Kernel, &
      Apply_EOS_I_E_A_Kernel, &
      Apply_EOS_I_E_S_Kernel, &
      Apply_EOS_I_E_SB_A_Kernel
      
  interface
  
    module subroutine Apply_EOS_I_T_Kernel &
             ( M, N, E, P, T, SB, SS, Gmm, &
               M_Ref, N_Min, T_Min, Gamma, C_V, N0, P0, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        M, &
        N, &
        E, &
        P, &
        T, &
        SB, &
        SS, &
        Gmm
      real ( KDR ), intent ( in ) :: &
        M_Ref, &
        N_Min, &
        T_Min, &
        Gamma, &
        C_V, &
        N0, &
        P0
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_I_T_Kernel

    module subroutine Apply_EOS_I_E_A_Kernel &
             ( M, N, E, P, T, SB, SS, Gmm, &
               M_Ref, N_Min, E_Min, Gamma, C_V, N0, P0, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        M, &
        N, &
        E, &
        P, &
        T, &
        SB, &
        SS, &
        Gmm
      real ( KDR ), intent ( in ) :: &
        M_Ref, &
        N_Min, &
        E_Min, &
        Gamma, &
        C_V, &
        N0, &
        P0
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_I_E_A_Kernel

    module subroutine Apply_EOS_I_E_S_Kernel &
             ( M, N, E, P, T, SB, SS, &
               M_Ref, N_Min, E_Min, Gamma, C_V, N0, P0, &
               iV )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        M, &
        N, &
        E, &
        P, &
        T, &
        SB, &
        SS
      real ( KDR ), intent ( in ) :: &
        M_Ref, &
        N_Min, &
        E_Min, &
        Gamma, &
        C_V, &
        N0, &
        P0
      integer ( KDI ), intent ( in ) :: &
        iV
    end subroutine Apply_EOS_I_E_S_Kernel

    module subroutine Apply_EOS_I_E_SB_A_Kernel &
             ( M, N, E, P, T, SB, SS, Gmm, Shock, Threshold, &
               M_Ref, N_Min, E_Min, Gamma, C_V, N0, P0, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        M, &
        N, &
        E, &
        P, &
        T, &
        SB, &
        SS, &
        Gmm
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        Shock
      real ( KDR ), intent ( in ) :: &
        Threshold, &
        M_Ref, &
        N_Min, &
        E_Min, &
        Gamma, &
        C_V, &
        N0, &
        P0
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_I_E_SB_A_Kernel

  end interface


contains


  subroutine InitializeAllocate_F &
               ( F, G, Units_F, FieldOption, VectorOption, NameOption, &
                 UnitOption, VectorIndicesOption, iaPrimitiveOption, &
                 iaBalancedOption, nFieldsOption, IgnorabilityOption )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Units_F_Form ), dimension ( : ), intent ( in ) :: &
      Units_F
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    real ( KDR ) :: &
      ShockThreshold

    if ( F % Type  ==  '' ) &
      F % Type  =  'a Fluid_P_I' 
    
    !-- Field indices: no additional fields
    
    !-- Field names: no additional fields

    !-- Units: no additional fields

    !-- Vector indices: no additional vectors

    !-- Vector names: no additional vectors

    !-- Primitive fields: no additional fields

    !-- Balanced fields: no additional fields

    !-- Fluid_P

    call F % Fluid_P_Form % Initialize &
           ( G, Units_F, &
             FieldOption = FieldOption, &
             VectorOption = VectorOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitiveOption, &
             iaBalancedOption = iaBalancedOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

    !-- Features

    if ( .not. allocated ( F % Features ) ) then

      ShockThreshold  =  0.1_KDR
      call PROGRAM_HEADER % GetParameter ( ShockThreshold, 'ShockThreshold' )

      allocate ( Features_F_P_Form :: F % Features )
      select type ( FFP  =>  F % Features )
      type is ( Features_F_P_Form )
        call FFP % Initialize ( F, ShockThreshold = ShockThreshold )
      end select !-- FFP

    end if

    !-- Parameters

    associate &
      ( k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume, &
        n_0   => F % FiducialBaryonDensity, &
        p_0   => F % FiducialPressure )

    if ( Units_F ( 1 ) % Temperature  ==  UNIT % IDENTITY ) then
      k  =  1.0_KDR
    else
      k  =  CONSTANT % BOLTZMANN
    end if

    gamma  =  1.4_KDR
    mu     =  1.0_KDR
    c_v    =  k / ( mu * ( gamma - 1.0_KDR ) )
    n_0    =  1.0_KDR
    p_0    =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( mu,    'MeanMolecularWeight' )
    call PROGRAM_HEADER % GetParameter ( c_v,   'SpecificHeatVolume' )
    call PROGRAM_HEADER % GetParameter ( n_0,   'FiducialBaryonDensity' )
    call PROGRAM_HEADER % GetParameter ( p_0,   'FiducialPressure' )

    end associate !-- k, etc.

  end subroutine InitializeAllocate_F

  
  subroutine SetAdiabaticIndex ( F, AdiabaticIndex )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      AdiabaticIndex

    F % AdiabaticIndex  =  AdiabaticIndex

    call Show ( 'Setting AdiabaticIndex of a Fluid_P_I', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % AdiabaticIndex, 'AdiabaticIndex', F % IGNORABILITY + 1 )

  end subroutine SetAdiabaticIndex


  subroutine SetMeanMolecularWeight ( F, MeanMolecularWeight )

    !-- Assumes AdiabaticIndex already set.

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      MeanMolecularWeight

    associate &
      ( k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume )

    mu  =  MeanMolecularWeight

    c_v  =  k / ( mu * ( gamma - 1.0_KDR ) )

    call Show ( 'Setting MeanMolecularWeight of a Fluid', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( mu, 'MeanMolecularWeight', &
                F % IGNORABILITY + 1 )
    call Show ( c_v, F % Unit ( F % ENTROPY_PER_BARYON, 1 ), &
                'SpecificHeatVolume', F % IGNORABILITY + 1 )

    end associate !-- k, etc.

  end subroutine SetMeanMolecularWeight


  subroutine SetSpecificHeatVolume ( F, SpecificHeatVolume )

    !-- Assumes AdiabaticIndex already set.

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      SpecificHeatVolume

    associate &
      ( k     => F % BoltzmannConstant, &
        gamma => F % AdiabaticIndex, &
        mu    => F % MeanMolecularWeight, &
        c_v   => F % SpecificHeatVolume )

    c_v  =  SpecificHeatVolume

    mu  =  k / ( c_v * ( gamma - 1.0_KDR ) )

    call Show ( 'Setting SpecificHeatVolume of a Fluid', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( c_v, F % Unit ( F % ENTROPY_PER_BARYON, 1 ), &
                'SpecificHeatVolume', F % IGNORABILITY + 1 )
    call Show ( F % MeanMolecularWeight, 'MeanMolecularWeight', &
                F % IGNORABILITY + 1 )

    end associate !-- amu, etc.

  end subroutine SetSpecificHeatVolume


  subroutine SetFiducialParameters &
               ( F, FiducialBaryonDensity, FiducialPressure )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      FiducialBaryonDensity, &
      FiducialPressure

    F % FiducialBaryonDensity = FiducialBaryonDensity
    F % FiducialPressure = FiducialPressure

    call Show ( 'Setting fiducial parameters of a Fluid', F % IGNORABILITY + 1 )
    call Show ( F % Name, 'Name', F % IGNORABILITY + 1 )
    call Show ( F % FiducialBaryonDensity, &
                F % Unit ( F % BARYON_DENSITY_C, 1 ), &
                'FiducialBaryonDensity', F % IGNORABILITY + 1 )
    call Show ( F % FiducialPressure, F % Unit ( F % PRESSURE, 1 ), &
                'FiducialPressure', F % IGNORABILITY + 1 )

  end subroutine SetFiducialParameters


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % BARYON_DENSITY_C, CS % VELOCITY_U, &
                    CS % ENERGY_DENSITY_C, CS % PRESSURE, CS % TEMPERATURE, &
                    CS % ENTROPY_PER_BARYON ] )

    if ( allocated ( CS % Features ) ) &
      call CS % Features % SetStream ( S )

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( Fluid_P_I_Form ), intent ( in ) :: &
      FS

    call FS % Fluid_P_Form % Show ( )

    call Show ( FS % BoltzmannConstant, &
                FS % Unit ( FS % ENTROPY_PER_BARYON, 1 ), &
                'BoltzmannConstant', FS % IGNORABILITY )
    call Show ( FS % AdiabaticIndex, 'AdiabaticIndex', FS % IGNORABILITY )
    call Show ( FS % MeanMolecularWeight, 'MeanMolecularWeight', &
                FS % IGNORABILITY )
    call Show ( FS % SpecificHeatVolume, &
                FS % Unit ( FS % ENTROPY_PER_BARYON, 1 ), &
                'SpecificHeatVolume', FS % IGNORABILITY )
    call Show ( FS % FiducialBaryonDensity, &
                FS % Unit ( FS % BARYON_DENSITY_C, 1 ), &
                'FiducialBaryonDensity', FS % IGNORABILITY )
    call Show ( FS % FiducialPressure, &
                FS % Unit ( FS % PRESSURE, 1 ), &
                'FiducialPressure', FS % IGNORABILITY )

  end subroutine Show_FS


  subroutine ComputeFromTemperature ( F )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromTemperature', CONSOLE % INFO_6 )
    call Show ( F % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, F % Atlas % nCharts

      associate &
        (    FV  =>  F % Storage ( iC ) % Value, &
          M_Ref  =>  F % BaryonMass, &
          N_Min  =>  F % BaryonDensityMin, &
          E_Min  =>  F % EnergyDensityMin, &
          T_Min  =>  F % TemperatureMin, &
          Gamma  =>  F % AdiabaticIndex, &
          C_V    =>  F % SpecificHeatVolume, &
          N_0    =>  F % FiducialBaryonDensity, &
          P_0    =>  F % FiducialPressure )
      associate &
        ( M    =>  FV ( :, F % BARYON_MASS ), &
          N    =>  FV ( :, F % BARYON_DENSITY_C ), &
          V_1  =>  FV ( :, F % VELOCITY_U_1 ), &
          V_2  =>  FV ( :, F % VELOCITY_U_2 ), &
          V_3  =>  FV ( :, F % VELOCITY_U_3 ), &
          D    =>  FV ( :, F % BARYON_DENSITY_B ), &
          S_1  =>  FV ( :, F % MOMENTUM_DENSITY_D_1 ), &
          S_2  =>  FV ( :, F % MOMENTUM_DENSITY_D_2 ), &
          S_3  =>  FV ( :, F % MOMENTUM_DENSITY_D_3 ), &
          E    =>  FV ( :, F % ENERGY_DENSITY_C ), &
          G    =>  FV ( :, F % ENERGY_DENSITY_B ), &
          P    =>  FV ( :, F % PRESSURE ), &
          T    =>  FV ( :, F % TEMPERATURE ), &
          SB   =>  FV ( :, F % ENTROPY_PER_BARYON ), &
          DS   =>  FV ( :, F % ENTROPY_DENSITY_B ), &
          SS   =>  FV ( :, F % SOUND_SPEED ), &
          Gmm  =>  FV ( :, F % ADIABATIC_INDEX ) )
   
      call Apply_EOS_I_T_Kernel &
             ( M, N, E, P, T, SB, Gmm, SS, M_Ref, N_Min, T_Min, Gamma, C_V, &
               N_0, P_0, UseDeviceOption = F % DeviceMemory )

      select type ( Gn  =>  F % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Gn % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, Gn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Gn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Gn % METRIC_F_DD_33 ) )

        call F % Compute_D_S_G_DS_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, E, SB, M, SS, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, D, S_1, S_2, S_3, G, DS, &
                 UseDeviceOption = F % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_P_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromTemperature', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Gn

      end associate !-- M, etc.
      end associate !-- FV, etc.

    end do !-- iC

  end subroutine ComputeFromTemperature


  subroutine ComputeFromPrimitive ( FS_CS, CS )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_CS
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeFromPrimitive', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, CS % Atlas % nCharts

      associate &
        (   CSV  =>  FS_CS % Storage ( iC ) % Value, &
          M_Ref  =>  CS % BaryonMass, &
          N_Min  =>  CS % BaryonDensityMin, &
          E_Min  =>  CS % EnergyDensityMin, &
          Gamma  =>  CS % AdiabaticIndex, &
          C_V    =>  CS % SpecificHeatVolume, &
          N_0    =>  CS % FiducialBaryonDensity, &
          P_0    =>  CS % FiducialPressure )
      associate &
        ( M    =>  CSV ( :, CS % BARYON_MASS ), &
          N    =>  CSV ( :, CS % BARYON_DENSITY_C ), &
          V_1  =>  CSV ( :, CS % VELOCITY_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_U_3 ), &
          D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
          S_1  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
          E    =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
          G    =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          P    =>  CSV ( :, CS % PRESSURE ), &
          T    =>  CSV ( :, CS % TEMPERATURE ), &
          SB   =>  CSV ( :, CS % ENTROPY_PER_BARYON ), &
          DS   =>  CSV ( :, CS % ENTROPY_DENSITY_B ), &
          SS   =>  CSV ( :, CS % SOUND_SPEED ), &
          Gmm  =>  CSV ( :, CS % ADIABATIC_INDEX ) )

      call Apply_EOS_I_E_A_Kernel &
             ( M, N, E, P, T, SB, SS, Gmm, M_Ref, N_Min, E_Min, Gamma, C_V, &
               N_0, P_0, UseDeviceOption = CS % DeviceMemory )

      select type ( Gn  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Gn % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, Gn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Gn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Gn % METRIC_F_DD_33 ) )

        call CS % Compute_D_S_G_DS_G_Kernel & 	 	 
               ( N, V_1, V_2, V_3, E, SB, M, SS, M_DD_11, M_DD_22, M_DD_33, &
                 N_Min, E_Min, D, S_1, S_2, S_3, G, DS, &
                 UseDeviceOption = CS % DeviceMemory )

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_P_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromPrimitive', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Gn

      end associate !-- M, etc.
      end associate !-- CSV, etc.

    end do !-- iC

  end subroutine ComputeFromPrimitive


  subroutine ComputeFromBalancedAll ( CS, T_Option, DetectFeaturesOption )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    logical ( KDL ), intent ( in ), optional :: &
      DetectFeaturesOption

    integer ( KDI ) :: &
      iC
    logical ( KDL ) :: &
      DetectFeatures
    type ( TimerForm ), pointer :: &
      T_G, &
      T_K

    call Show ( 'ComputeFromBalancedAll', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    if ( present ( T_Option ) ) then
      T_K  =>  PROGRAM_HEADER % Timer &
                 ( Handle = CS % iTimer_CFB, &
                   Name = trim ( T_Option % Name ) // '_Krnl', &
                   Level = T_Option % Level + 1 )
    else
      T_K  =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )
    do iC  =  1, CS % Atlas % nCharts

      associate &
        (   CSV  =>  CS % Storage ( iC ) % Value, &
          M_Ref  =>  CS % BaryonMass, &
          N_Min  =>  CS % BaryonDensityMin, &
          E_Min  =>  CS % EnergyDensityMin, &
          Gamma  =>  CS % AdiabaticIndex, &
          C_V    =>  CS % SpecificHeatVolume, &
          N_0    =>  CS % FiducialBaryonDensity, &
          P_0    =>  CS % FiducialPressure, &
          Threshold  =>  CS % EntropyEnergyThreshold )
      associate &
        ( M    =>  CSV ( :, CS % BARYON_MASS ), &
          N    =>  CSV ( :, CS % BARYON_DENSITY_C ), &
          V_1  =>  CSV ( :, CS % VELOCITY_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_U_3 ), &
          D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
          S_1  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
          E    =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
          G    =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          P    =>  CSV ( :, CS % PRESSURE ), &
          T    =>  CSV ( :, CS % TEMPERATURE ), &
          SB   =>  CSV ( :, CS % ENTROPY_PER_BARYON ), &
          DS   =>  CSV ( :, CS % ENTROPY_DENSITY_B ), &
          SS   =>  CSV ( :, CS % SOUND_SPEED ), &
          Gmm  =>  CSV ( :, CS % ADIABATIC_INDEX ) )
   
      select type ( Gn  =>  CS % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  Gn % Storage ( iC ) % Value )
        associate &
          ( M_UU_11  =>  GSV ( :, Gn % METRIC_F_UU_11 ), &
            M_UU_22  =>  GSV ( :, Gn % METRIC_F_UU_22 ), &
            M_UU_33  =>  GSV ( :, Gn % METRIC_F_UU_33 ), &
            M_DD_11  =>  GSV ( :, Gn % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, Gn % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, Gn % METRIC_F_DD_33 ) )

        call CS % Compute_N_V_E_SB_G_A_Kernel &
               ( D, S_1, S_2, S_3, G, DS, M, M_UU_11, M_UU_22, M_UU_33, &
                 N_Min, E_Min, N, V_1, V_2, V_3, E, SB, &
                 UseDeviceOption = CS % DeviceMemory )

        if ( CS % UseEntropy ) then
          select type ( F  =>  CS % Features )
            class is ( Features_F_P_Form )
          associate ( FV  =>  CS % Features % Storage ( iC ) % Value )
          associate ( Shock  =>  FV ( :, F % SHOCK ) )
          call Apply_EOS_I_E_SB_A_Kernel &
                 ( M, N, E, P, T, SB, SS, Gmm, Shock, Threshold, &
                   M_Ref, N_Min, E_Min, Gamma, C_V, N_0, P_0, &
                   UseDeviceOption = CS % DeviceMemory )
          call CS % Compute_D_S_G_DS_G_Kernel & 	 	 
                 ( N, V_1, V_2, V_3, E, SB, M, SS, M_DD_11, M_DD_22, M_DD_33, &
                   N_Min, E_Min, D, S_1, S_2, S_3, G, DS, &
                   UseDeviceOption = CS % DeviceMemory )
          end associate !-- Shock
          end associate !-- FV
          end select !-- F
        else !-- use only energy
          call Apply_EOS_I_E_A_Kernel &
                 ( M, N, E, P, T, SB, SS, Gmm, &
                   M_Ref, N_Min, E_Min, Gamma, C_V, N_0, P_0, &
                   UseDeviceOption = CS % DeviceMemory )
        end if !-- UseEntropy

        end associate !-- M_UU_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Fluid_P_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeFromBalancedAll', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Gn

      end associate !-- M, etc.
      end associate !-- CSV, etc.

    end do !-- iC
    if ( associated ( T_K ) ) call T_K % Stop ( )

    DetectFeatures = .false.
    if ( present ( DetectFeaturesOption ) ) &
      DetectFeatures = DetectFeaturesOption
    if ( DetectFeatures ) &
      call CS % Features % Detect ( )

  end subroutine ComputeFromBalancedAll


  subroutine ComputeFromBalancedSingle ( CS, iC, iV, T_Option )

    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      CS
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

    type ( TimerForm ), pointer :: &
      T_K

    call Show ( 'ComputeFromBalancedAll', CONSOLE % INFO_6 )
    call Show ( CS % Name, 'Fluid', CONSOLE % INFO_6 )

    if ( present ( T_Option ) ) then
      T_K  =>  PROGRAM_HEADER % Timer &
                 ( Handle = CS % iTimer_CFB, &
                   Name = trim ( T_Option % Name ) // '_Krnl', &
                   Level = T_Option % Level + 1 )
    else
      T_K  =>  null ( )
    end if

    if ( associated ( T_K ) ) call T_K % Start ( )

    associate &
      (   CSV  =>  CS % Storage ( iC ) % Value, &
        M_Ref  =>  CS % BaryonMass, &
        N_Min  =>  CS % BaryonDensityMin, &
        E_Min  =>  CS % EnergyDensityMin, &
        Gamma  =>  CS % AdiabaticIndex, &
        C_V    =>  CS % SpecificHeatVolume, &
        N_0    =>  CS % FiducialBaryonDensity, &
        P_0    =>  CS % FiducialPressure )
    associate &
      ( M    =>  CSV ( :, CS % BARYON_MASS ), &
        N    =>  CSV ( :, CS % BARYON_DENSITY_C ), &
        V_1  =>  CSV ( :, CS % VELOCITY_U_1 ), &
        V_2  =>  CSV ( :, CS % VELOCITY_U_2 ), &
        V_3  =>  CSV ( :, CS % VELOCITY_U_3 ), &
        D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
        S_1  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
        S_2  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
        S_3  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
        E    =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
        G    =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
        P    =>  CSV ( :, CS % PRESSURE ), &
        T    =>  CSV ( :, CS % TEMPERATURE ), &
        SB   =>  CSV ( :, CS % ENTROPY_PER_BARYON ), &
        DS   =>  CSV ( :, CS % ENTROPY_DENSITY_B ), &
        SS   =>  CSV ( :, CS % SOUND_SPEED ) )
 
    select type ( Gn  =>  CS % Geometry )
    class is ( Gravitation_G_Form )

      associate &
        ( GSV  =>  Gn % Storage ( iC ) % Value )
      associate &
        ( M_UU_11  =>  GSV ( :, Gn % METRIC_F_UU_11 ), &
          M_UU_22  =>  GSV ( :, Gn % METRIC_F_UU_22 ), &
          M_UU_33  =>  GSV ( :, Gn % METRIC_F_UU_33 ) )

      call CS % Compute_N_V_E_SB_G_S_Kernel &
             ( D, S_1, S_2, S_3, G, DS, M, M_UU_11, M_UU_22, M_UU_33, &
               N_Min, E_Min, iV, N, V_1, V_2, V_3, E, SB )

      end associate !-- M_UU_11, etc.
      end associate !-- GSV

    class default
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_P_I__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFromBalancedSingle', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Gn

    call Apply_EOS_I_E_S_Kernel &
           ( M, N, E, P, T, SB, SS, M_Ref, N_Min, E_Min, Gamma, C_V, &
             N_0, P_0, iV )

    end associate !-- M, etc.
    end associate !-- CSV, etc.

    if ( associated ( T_K ) ) call T_K % Stop ( )

  end subroutine ComputeFromBalancedSingle


  impure elemental subroutine Finalize ( F )

    type ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

  end subroutine Finalize


end module Fluid_P_I__Form
