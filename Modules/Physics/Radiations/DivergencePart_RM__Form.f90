module DivergencePart_RM__Form

  !-- DivergencePart_RadiationMoments__Form

  use Basics
  use Mathematics
  use Gravitations
  use RadiationMoments_BM__Form

  implicit none
  private

  type, public, extends ( DivergencePart_CS_Form ) :: DivergencePart_RM_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass ( DP ) :: &
      ComputeFluxes
    ! procedure, public, pass ( DP ) :: &
    !   ComputeStresses
    final :: &
      Finalize
  end type DivergencePart_RM_Form

    private :: &
      Compute_FS_G_Kernel, &
      Compute_S_UD_Kernel

    interface

      module subroutine Compute_FS_G_Kernel &
               ( J, H_1, H_2, H_3, H_Dim, S_Dim, SF, &
                 M_DD_11, M_DD_22, M_DD_33, M_UU_Dim, iDim, &
                 F_E, F_S_1, F_S_2, F_S_3, UseDeviceOption )
        !-- Compute_FluxSet_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J, &
          H_1, H_2, H_3, &
          H_Dim, &
          S_Dim, &
          SF, &
          M_DD_11, M_DD_22, M_DD_33, &
          M_UU_Dim
        integer ( KDI ), intent ( in ) :: &
          iDim
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          F_E, &
          F_S_1, F_S_2, F_S_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_FS_G_Kernel

      module subroutine Compute_S_UD_Kernel &
               ( )
        !-- Compute_Stress_UD_Kernel
      end subroutine Compute_S_UD_Kernel

    end interface


contains


  subroutine Initialize ( DP, CS, NameOption, IgnorabilityOption )

    class ( DivergencePart_RM_Form ), intent ( inout ) :: &
      DP
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( DP % Type  ==  '' ) &
      DP % Type  =  'a DivergencePart_RM' 

    Name  =  'S'  !-- Streaming; is this even used? Compare 'V' in Fluid
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call DP % DivergencePart_CS_Form % Initialize &
           ( CS, NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine Initialize


  subroutine ComputeFluxes ( FS_F, DP, FS_CS, iC, iD )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_F  !-- Fluxes
    class ( DivergencePart_RM_Form ), intent ( in ) :: &
      DP
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    integer ( KDI ) :: &
      iEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    select type ( CS  =>  DP % CurrentSet )
      class is ( RadiationMoments_BM_Form )

    call Search &
           ( CS % iaBalanced, CS % ENERGY_DENSITY_B, iEnergy )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_B_D_1, iMomentum ( 1 ) )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_B_D_2, iMomentum ( 2 ) )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_B_D_3, iMomentum ( 3 ) )

    associate &
      ( FSV  =>  FS_F  % Storage ( iC ) % Value, &
        CSV  =>  FS_CS % Storage ( iC ) % Value )
    associate &
      ( F_E      =>  FSV ( :, iEnergy ), &
        F_S_1    =>  FSV ( :, iMomentum ( 1 ) ), &
        F_S_2    =>  FSV ( :, iMomentum ( 2 ) ), &
        F_S_3    =>  FSV ( :, iMomentum ( 3 ) ), &
          J      =>  CSV ( :, CS % ENERGY_DENSITY_C ), &
          H_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_1 ), &
          H_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_2 ), &
          H_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U_3 ), &
          H_Dim  =>  CSV ( :, CS % MOMENTUM_DENSITY_C_U ( iD ) ), &
          S_Dim  =>  CSV ( :, CS % MOMENTUM_DENSITY_B_D ( iD ) ), &
         SF      =>  CSV ( :, CS % STRESS_FACTOR ) )
 
    select type ( G  =>  CS % Geometry )
    class is ( Gravitation_G_Form )

      associate &
        ( GSV  =>  G % Storage ( iC ) % Value )
      associate &
        ( M_DD_11   =>  GSV ( :, G % METRIC_F_DD_11 ), &
          M_DD_22   =>  GSV ( :, G % METRIC_F_DD_22 ), &
          M_DD_33   =>  GSV ( :, G % METRIC_F_DD_33 ), &
          M_UU_Dim  =>  GSV ( :, G % METRIC_F_UU ( iD ) ) )

      call Compute_FS_G_Kernel &
             ( J, H_1, H_2, H_3, H_Dim, S_Dim, SF, &
               M_DD_11, M_DD_22, M_DD_33, M_UU_Dim, iD, &
               F_E, F_S_1, F_S_2, F_S_3, &
               UseDeviceOption = CS % DeviceMemory )
  
      end associate !-- M_UU_Dim
      end associate !-- GSV

    class default
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'DivergencePart_RM__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeFluxes', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

    end associate !-- F_E, etc.
    end associate !-- FSV, etc.


    end select !-- CS

  end subroutine ComputeFluxes


  impure elemental subroutine Finalize ( DP )

    type ( DivergencePart_RM_Form ), intent ( inout ) :: &
      DP

  end subroutine Finalize


end module DivergencePart_RM__Form
