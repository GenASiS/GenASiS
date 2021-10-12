module Fluid_P_I__Form

  !-- Fluid_Perfect_Ideal__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_P__Form

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
    procedure, private, pass :: &
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
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromPrimitiveCommon
  !   procedure, public, pass ( C ) :: &
  !     ComputeFromConservedCommon
  !   procedure, public, pass ( C ) :: &
  !     ComputeRawFluxes
    procedure, public, nopass :: &
      Apply_EOS_I_T_Kernel
  !   procedure, public, nopass :: &
  !     Apply_EOS_I_SB_E_Kernel
  !   procedure, public, nopass :: &
  !     Apply_EOS_I_E_Kernel
    final :: &
      Finalize
  end type Fluid_P_I_Form


  interface
  
    module subroutine Apply_EOS_I_T_Kernel &
                 ( P, E, SB, CS, M, N, T, Gamma, C_V, N0, P0, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        P, &
        E, &
        SB, &
        CS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        T
      real ( KDR ), intent ( in ) :: &
        Gamma, &
        C_V, &
        N0, &
        P0
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Apply_EOS_I_T_Kernel

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
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaPrimitiveOption, &
      iaBalancedOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

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

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( Fluid_P_I_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % BARYON_DENSITY_C, CS % VELOCITY_U, &
                    CS % ENERGY_DENSITY_C, CS % PRESSURE, CS % TEMPERATURE, &
                    CS % MACH_NUMBER, CS % ENTROPY_PER_BARYON ] )

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( Fluid_P_I_Form ), intent ( in ) :: &
      FS

    call FS % Fluid_D_Form % Show ( )

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

    call Show ( 'ComputeFromBalanced', CONSOLE % INFO_6 )
    call Show ( F % Name, 'Fluid', CONSOLE % INFO_6 )

    do iC  =  1, F % Atlas % nCharts

      associate &
        (    FV  =>  F % Storage ( iC ) % Value, &
          M_Ref  =>  F % BaryonMass, &
          N_Min  =>  F % BaryonDensityMin )
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
          CS   =>  FV ( :, F % SOUND_SPEED ), &
          MN   =>  FV ( :, F % MACH_NUMBER ) )
   
      call F % Compute_M_Kernel &
             ( M_Ref, M, UseDeviceOption = F % DeviceMemory )
      call F % Apply_EOS_I_T_Kernel &
             ( P, E, SB, CS, M, N, T, F % AdiabaticIndex, &
               F % SpecificHeatVolume, F % FiducialBaryonDensity, &
               F % FiducialPressure, UseDeviceOption = F % DeviceMemory )

      end associate !-- M, etc.
      end associate !-- FV, etc.

    end do !-- iC

  end subroutine ComputeFromTemperature


  impure elemental subroutine Finalize ( F )

    type ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

  end subroutine Finalize


end module Fluid_P_I__Form
