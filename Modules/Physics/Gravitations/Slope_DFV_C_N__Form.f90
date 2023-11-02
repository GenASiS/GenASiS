module Slope_DFV_C_N__Form

  !-- Slope_DivergenceFiniteVolume_Curvili_Newton__Form

  use Basics
  use Mathematics
  use Gravitation_N_H__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_C_N_Form
    integer ( KDI ) :: &
      iBaryonMass_F, &
      iBaryonDensity_F, &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iVelocity_F, &
      iMomentum_B
    class ( CurrentSetForm ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_C_N
    generic, public :: &
      Initialize => InitializeAllocate_C_N
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_DFV_C_N_Form

    private :: &
      ComputeMomentumKernel, &
      ComputeEnergyKernel

    interface
    
      module subroutine ComputeMomentumKernel &
               ( ProperCell, M, N, dPhi_1, dPhi_2, dPhi_3, &
                 S_S_1, S_S_2, S_S_3, UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M, &
          N, &
          dPhi_1, dPhi_2, dPhi_3
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_S_1, S_S_2, S_S_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeMomentumKernel

      module subroutine ComputeEnergyKernel &
               ( ProperCell, V_1, V_2, V_3, S_S_1, S_S_2, S_S_3, S_G, &
                 UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          V_1, V_2, V_3, &
          S_S_1, S_S_2, S_S_3
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_G
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeEnergyKernel

    end interface


contains


  subroutine InitializeAllocate_C_N &
               ( S, Fluid, iVelocity_F, iMomentum_B, iBaryonMass_F, &
                 iBaryonDensity_F, iEnergy_B )

    class ( Slope_DFV_C_N_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      Fluid
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iVelocity_F, &
      iMomentum_B
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass_F, &
      iBaryonDensity_F, &
      iEnergy_B

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_C_N' 
    
    Name  =  trim ( Fluid % Name ) // '_Slp_DFV_C_N'

    S % Fluid  =>  Fluid

    S % iBaryonMass_F     =  iBaryonMass_F
    S % iBaryonDensity_F  =  iBaryonDensity_F
    S % iVelocity_F       =  iVelocity_F
    S % iMomentum_B       =  iMomentum_B
    S % iEnergy_B         =  iEnergy_B

    associate ( F  =>  S % Fluid )

    call S % Slope_H_Form % Initialize &
           ( F % Atlas, &
             FieldOption = F % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             nFieldsOption = F % nBalanced, &
             IgnorabilityOption = F % IGNORABILITY + 1 )

    end associate !-- F

  end subroutine InitializeAllocate_C_N


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_DFV_C_N_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate &
      ( F  =>  S % Fluid )
    select type ( G  =>  F % Geometry )
      class is ( Gravitation_N_H_Form )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate &
        ( FV  =>  F % Storage ( iC ) % Value, &
          GV  =>  G % Storage ( iC ) % Value, &
          SV  =>  S % Storage ( iC ) % Value )
      associate &
        ( M       =>  FV ( :, S % iBaryonMass_F ), &
          N       =>  FV ( :, S % iBaryonDensity_F ), &
          dPhi_1  =>  GV ( :, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
          dPhi_2  =>  GV ( :, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
          dPhi_3  =>  GV ( :, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
          S_S_1   =>  SV ( :, S % iMomentum_B ( 1 ) ), &
          S_S_2   =>  SV ( :, S % iMomentum_B ( 2 ) ), &
          S_S_3   =>  SV ( :, S % iMomentum_B ( 3 ) ) )

      call ComputeMomentumKernel &
             ( C % ProperCell, M, N, dPhi_1, dPhi_2, dPhi_3, &
               S_S_1, S_S_2, S_S_3, UseDeviceOption = S % DeviceMemory )

      if ( S % iEnergy_B  >  0 ) then
        associate &
          ( V_1  =>  FV ( :, S % iVelocity_F ( 1 ) ), &
            V_2  =>  FV ( :, S % iVelocity_F ( 2 ) ), &
            V_3  =>  FV ( :, S % iVelocity_F ( 3 ) ), &
            S_G  =>  SV ( :, S % iEnergy_B ) )

        call ComputeEnergyKernel &
               ( C % ProperCell, V_1, V_2, V_3, S_S_1, S_S_2, S_S_3, S_G, &
                 UseDeviceOption = S % DeviceMemory )

        end associate !-- V_1, etc.
      end if

      end associate !-- M, etc.
      end associate !-- FV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_DFV_C_N__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    class default
      call Show ( 'Newtonian gravitation expected', CONSOLE % ERROR )
      call Show ( 'Slope_DFV_C_N__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    end associate !-- F

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_C_N_Form ), intent ( inout ) :: &
      S

    nullify ( S % Fluid )
    
  end subroutine Finalize


end module Slope_DFV_C_N__Form
