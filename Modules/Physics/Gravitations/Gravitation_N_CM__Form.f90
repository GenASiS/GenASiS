module Gravitation_N_CM__Form

  !-- Gravitation_Newton_CentralMass__Form

  use Basics
  use Mathematics
  use Gravitation_N_H__Form

  implicit none
  private

  type, public, extends ( Gravitation_N_H_Form ) :: Gravitation_N_CM_Form
    real ( KDR ) :: &
      GravitationalConstant, &
      Mass
  contains
    procedure, private, pass :: &
      InitializeAllocate_N_CM
    generic, public :: &
      Initialize => InitializeAllocate_N_CM
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
  end type Gravitation_N_CM_Form

    private :: &
      SolveKernel

    interface
    
      module subroutine SolveKernel &
               ( Phi, GradPhi_1, GradPhi_2, GradPhi_3, R, G, M, nD, &
                 UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Phi, &
          GradPhi_1, GradPhi_2, GradPhi_3
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          R
        real ( KDR ), intent ( in ) :: &
          G, &  !-- Gravitational constant
          M     !-- Central mass
        integer ( KDI ), intent ( in ) :: &
          nD
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine SolveKernel

    end interface


contains


  subroutine InitializeAllocate_N_CM &
               ( G, A, GravitationalConstant, Mass, FieldOption, VectorOption, &
                 NameOption, DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, UnitOption, &
                 VectorIndicesOption, nFieldsOption, IgnorabilityOption )

    class ( Gravitation_N_CM_Form ), intent ( inout ), target :: &
      G
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    real ( KDR ), intent ( in ) :: &
      GravitationalConstant, &
      Mass
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption, &
      AssociateFieldsOption
    type ( MeasuredValueForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    if ( G % Type  ==  '' ) &
      G % Type  =  'a Gravitation_N_CM'

    call G % Gravitation_N_H_Form % Initialize &
           ( A, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, AssociateFieldsOption, UnitOption, &
             VectorIndicesOption, nFieldsOption, IgnorabilityOption )

    G % GravitationalConstant  =  GravitationalConstant
    G % Mass                   =  Mass

  end subroutine InitializeAllocate_N_CM


  subroutine Show_FS ( FS )

    class ( Gravitation_N_CM_Form ), intent ( in ) :: &
      FS

    integer ( KDI ) :: &
      iD

    call FS % Gravitation_N_H_Form % Show ( )

    call Show ( FS % Mass, 'Mass' )

  end subroutine Show_FS


  subroutine Solve ( G, F, iBaryonMass, iBaryonDensity, T_Option )

    class ( Gravitation_N_CM_Form ), intent ( inout ) :: &
      G
    class ( FieldSetForm ), intent ( in ) :: &
      F  !-- Fluid
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass, &
      iBaryonDensity
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iChart

    associate ( nC  =>  G % Atlas % nCharts )
    do iC  =  1,  nC

      associate &
        ( C    =>  G % Atlas % Chart ( iC ) % Element, &
          GSV  =>  G % Storage ( iC ) % Value )
      associate &
        (     Phi    =>  GSV ( :, G % POTENTIAL ), &
          GradPhi_1  =>  GSV ( :, G % POTENTIAL_GRADIENT_D_1 ), &
          GradPhi_2  =>  GSV ( :, G % POTENTIAL_GRADIENT_D_2 ), &
          GradPhi_3  =>  GSV ( :, G % POTENTIAL_GRADIENT_D_3 ), &
              R      =>  GSV ( :, G % CENTER_U_1 ) )

      select type ( C )
      class is ( Chart_GS_C_Form )

      call SolveKernel &
             ( Phi, GradPhi_1, GradPhi_2, GradPhi_3, R, &
               G = G % GravitationalConstant, M = G % Mass, &
               nD = C % nDimensions, UseDeviceOption = G % DeviceMemory )

      class default
        call Show ( 'Chart type not implemented', CONSOLE % ERROR )
        call Show ( C % Name, 'Name', CONSOLE % ERROR )
        call Show ( C % Type, 'Type', CONSOLE % ERROR )
        call Show ( 'Gravitation_N_CM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Solve', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select    

      end associate !-- Phi, etc.
      end associate !-- C, etc.

    end do !-- iC
    end associate !-- nC

  end subroutine Solve


  impure elemental subroutine Finalize ( G )

    type ( Gravitation_N_CM_Form ), intent ( inout ) :: &
      G

  end subroutine Finalize


end module Gravitation_N_CM__Form
