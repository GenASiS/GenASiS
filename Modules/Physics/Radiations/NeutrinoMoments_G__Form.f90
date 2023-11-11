module NeutrinoMoments_G__Form

  !-- NeutrinoMoments_Grey__Form

  use Basics
  use Mathematics
  use Fluids
  use Units_R__Form
  use Interactions_BM__Form
  use PhotonMoments_G__Form

  implicit none
  private

  integer ( KDI ), private, parameter :: &
    N_FIELDS_NM    = 6, &
    N_PRIMITIVE_NM = 1, &
    N_BALANCED_NM  = 1

  type, public, extends ( PhotonMoments_G_Form ) :: NeutrinoMoments_G_Form
    integer ( KDI ) :: &
      N_FIELDS_NM    = N_FIELDS_NM, &
      N_PRIMITIVE_NM = N_PRIMITIVE_NM, &
      N_BALANCED_NM  = N_BALANCED_NM
    integer ( KDI ) :: &
      NUMBER_DENSITY_C    = 0, &  !-- Comoving
      NUMBER_DENSITY_C_EQ = 0, &
      NUMBER_DENSITY_B    = 0     !-- Balanced
    integer ( KDI ) :: &
      DEGENERACY_GREY, &
      ENERGY_AVERAGE, &
      OCCUPANCY_AVERAGE
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    final :: &
      Finalize
    procedure, public, pass ( CS ) :: &
      SetStream
    procedure, public, pass :: &
      ComputeSpectralParameters
    procedure, public, pass :: &
      ComputeEquilibrium
  end type NeutrinoMoments_G_Form

    private :: &
      Compute_SP_Kernel, &
      Compute_Eq_Kernel
    
    interface

      module subroutine Compute_SP_Kernel &
               ( T_R, Eta_R, E_Ave, F_Ave, J, N, UseDeviceOption )
        !-- Compute_SpectralParameters_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          T_R, Eta_R, &
          E_Ave, F_Ave
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J, N
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_SP_Kernel
 
      module subroutine Compute_Eq_Kernel &
               ( J_Eq, N_Eq, T, Mu_E, Mu_NP, Sign, UseDeviceOption )
        !-- Compute_Equilibrium_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          J_Eq, N_Eq
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          T, &
          Mu_E, Mu_NP
        real ( KDR ), intent ( in ) :: &
          Sign
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_Eq_Kernel

    end interface


contains


  subroutine InitializeAllocate_RM &
               ( RM, G, Units_R, RadiationType, FieldOption, VectorOption, &
                 NameOption, UnitOption, VectorIndicesOption, &
                 iaPrimitiveOption, iaBalancedOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    character ( * ), intent ( in ) :: &
      RadiationType
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

    integer ( KDI ) :: &
      iC, &  !-- iChart
      oF, &  !-- oField
      oP, &  !-- oPrimitive
      oB, &  !-- oBalanced
      nFields, &
      nPrimitive, &
      nBalanced
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaPrimitive, &
      iaBalanced
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RM % Type  ==  '' ) &
      RM % Type  =  'a NeutrinoMoments_G' 
    
    Name  =  'NeutrinoMoments'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  RM % N_FIELDS_CS  +  RM % N_FIELDS_RM  +  RM % N_FIELDS_PM

    RM % NUMBER_DENSITY_C     =  oF + 1
    RM % NUMBER_DENSITY_C_EQ  =  oF + 2
    RM % NUMBER_DENSITY_B     =  oF + 3
    RM % DEGENERACY_GREY      =  oF + 4
    RM % ENERGY_AVERAGE       =  oF + 5
    RM % OCCUPANCY_AVERAGE    =  oF + 6

    nFields  =  oF  +  RM % N_FIELDS_NM
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RM % N_FIELDS_NM ) &
      = [ 'NumberDensity_C   ', &
          'NumberDensity_C_Eq', &
          'NumberDensity_B   ', &
          'DegeneracyGrey    ', &
          'EnergyAverage     ', &
          'OccupancyAverage  ' ]

    !-- Units

    associate ( nC  =>  G % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC
      FieldUnit ( RM % NUMBER_DENSITY_C, iC ) &
        =  Units_R ( iC ) % NumberDensity
      FieldUnit ( RM % NUMBER_DENSITY_C_EQ, iC ) &
        =  Units_R ( iC ) % NumberDensity
      FieldUnit ( RM % NUMBER_DENSITY_B, iC ) &
        =  Units_R ( iC ) % NumberDensity
      FieldUnit ( RM % ENERGY_AVERAGE, iC ) &
        =  Units_R ( iC ) % Coordinate_MS ( 1 )
    end do !-- iC

    end associate !-- nC

    !-- Primitive fields

    oP  =  RM % N_PRIMITIVE_CS  +  RM % N_PRIMITIVE_RM

    if ( present ( iaPrimitiveOption ) ) then
      nPrimitive  =  size ( iaPrimitiveOption )
      allocate ( iaPrimitive, source = iaPrimitiveOption )
    else
      nPrimitive  =  oP  +  RM % N_PRIMITIVE_NM
      allocate ( iaPrimitive ( nPrimitive ) )
    end if !-- iaPrimitiveOption

    iaPrimitive ( oP  +  1 : oP  +  RM % N_PRIMITIVE_NM )  &
      =  [ RM % NUMBER_DENSITY_C ]

    !-- Balanced fields

    oB  =  RM % N_BALANCED_CS  +  RM % N_BALANCED_RM

    if ( present ( iaBalancedOption ) ) then
      nBalanced  =  size ( iaBalancedOption )
      allocate ( iaBalanced, source = iaBalancedOption )
    else
      nBalanced  =  oB  +  RM % N_BALANCED_NM
      allocate ( iaBalanced ( nBalanced ) )
    end if !-- iaBalancedOption

    iaBalanced ( oB  +  1 : oB  +  RM % N_BALANCED_NM )  &
      =  [ RM % NUMBER_DENSITY_B ]

    !-- PhotonMoments_G

    call RM % PhotonMoments_G_Form % Initialize &
           ( G, Units_R, RadiationType, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             UnitOption = FieldUnit, &
             VectorIndicesOption = VectorIndicesOption, &
             iaPrimitiveOption = iaPrimitive, &
             iaBalancedOption = iaBalanced, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_RM


  impure elemental subroutine Finalize ( PM )

    type ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      PM

  end subroutine Finalize


  subroutine SetStream ( S, CS )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      CS

    call S % AddFieldSet &
           ( CS, &
             iaSelectedOption &
               =  [ CS % ENERGY_DENSITY_C, &
                    CS % MOMENTUM_DENSITY_C_U, &
                    CS % FLUX_FACTOR, &
                    CS % STRESS_FACTOR, &
                    CS % TEMPERATURE_GREY, &
                    CS % NUMBER_DENSITY_C, &
                    CS % DEGENERACY_GREY, &
                    CS % ENERGY_AVERAGE, &
                    CS % OCCUPANCY_AVERAGE ] )

  end subroutine SetStream


  subroutine ComputeSpectralParameters ( RM )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeSpectralParameters', CONSOLE % INFO_6 )
    call Show ( RM % Name, 'NeutrinoMoments', CONSOLE % INFO_6 )

    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RMV  =>  RM % Storage ( iC ) % Value )
      associate &
        (   T_R    =>  RMV ( :, RM % TEMPERATURE_GREY ), &
          Eta_R    =>  RMV ( :, RM % DEGENERACY_GREY ), &
            E_Ave  =>  RMV ( :, RM % ENERGY_AVERAGE ), &
            F_Ave  =>  RMV ( :, RM % OCCUPANCY_AVERAGE ), &
            J      =>  RMV ( :, RM % ENERGY_DENSITY_C ), &
            N      =>  RMV ( :, RM % NUMBER_DENSITY_C ) )
               
      call Compute_SP_Kernel &
             ( T_R, Eta_R, E_Ave, F_Ave, J, N, &
               UseDeviceOption  =  RM % DeviceMemory )

      end associate !-- T_R, etc.
      end associate !-- RV, etc.
    end do !-- iC

  end subroutine ComputeSpectralParameters


  subroutine ComputeEquilibrium ( RM )

    class ( NeutrinoMoments_G_Form ), intent ( inout ) :: &
      RM

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeEquilibrium', CONSOLE % INFO_6 )
    call Show ( RM % Name, 'NeutrinoMoments', CONSOLE % INFO_6 )

    select type ( I  =>  RM % Interactions )
      class is ( Interactions_BM_Form )
    select type ( F  =>  I % Fluid )
      class is ( Fluid_P_HN_Form )

    do iC  =  1, RM % Atlas % nCharts
      associate &
        ( RMV  =>  RM % Storage ( iC ) % Value, &
           FV  =>   F % Storage ( iC ) % Value )
      associate &
        (  J_Eq  =>  RMV ( :, RM % ENERGY_DENSITY_C_EQ ), &
           N_Eq  =>  RMV ( :, RM % NUMBER_DENSITY_C_EQ ), &
           T     =>   FV ( :,  F % TEMPERATURE ), &
          Mu_E   =>   FV ( :,  F % CHEMICAL_POTENTIAL_E ), &
          Mu_NP  =>   FV ( :,  F % CHEMICAL_POTENTIAL_N_P ) )

      select case ( trim ( RM % RadiationType ) )
      case ( 'NEUTRINOS_E' )
        call Compute_Eq_Kernel &
               ( J_Eq, N_Eq, T, Mu_E, Mu_NP, Sign = +1.0_KDR, &
                 UseDeviceOption = RM % DeviceMemory )
      case ( 'NEUTRINOS_E_BAR' )
        call Compute_Eq_Kernel &
               ( J_Eq, N_Eq, T, Mu_E, Mu_NP, Sign = -1.0_KDR, &
                 UseDeviceOption = RM % DeviceMemory )
      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( RM % RadiationType, 'RadiationType', CONSOLE % ERROR )
        call Show ( 'NeutrinoMoments_G__Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeEquilibrium', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- Name

      end associate !-- T_R, etc.
      end associate !-- RV, etc.
    end do !-- iC

    end select !-- F
    end select !-- I

  end subroutine ComputeEquilibrium


end module NeutrinoMoments_G__Form
