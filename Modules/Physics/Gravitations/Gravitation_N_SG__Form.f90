module Gravitation_N_SG__Form

  !-- Gravitation_Newton_SelfGravity__Form

  use Basics
  use Mathematics
  use Gravitation_N_H__Form

  implicit none
  private

  type, public, extends ( Gravitation_N_H_Form ) :: Gravitation_N_SG_Form
    integer ( KDI ) :: &
      iTimer_S = 0, &  !-- Source
      iTimer_G = 0     !-- Gradient
    real ( KDR ) :: &
      GravitationalConstant
    type ( FieldSetForm ), allocatable :: &
      Source, &
      Solution
    type ( FieldSetForm ), dimension ( : ), allocatable :: &
      SolutionGradient
    type ( GradientForm ), allocatable :: &
      Gradient
    type ( Poisson_ASCG_Form ), allocatable :: &
      Poisson
  contains
    procedure, private, pass :: &
      InitializeAllocate_N_SG
    generic, public :: &
      Initialize => InitializeAllocate_N_SG
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
  end type Gravitation_N_SG_Form

    private :: &
      ComputeSourceKernel

    interface

      module subroutine ComputeSourceKernel &
               ( M, N, G, S, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M, &
          N
        real ( KDR ), intent ( in ) :: &
          G
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeSourceKernel

   end interface

   
contains


  subroutine InitializeAllocate_N_SG &
               ( G, A, GravitationalConstant, FieldOption, VectorOption, &
                 NameOption, DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, UnitOption, &
                 VectorIndicesOption, nFieldsOption, IgnorabilityOption )

    class ( Gravitation_N_SG_Form ), intent ( inout ), target :: &
      G
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    real ( KDR ), intent ( in ) :: &
      GravitationalConstant
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

    integer ( KDI ) :: &
      iD, &
      MaxDegree
    character ( 1 ) :: &
      DimensionNumber

    if ( G % Type  ==  '' ) &
      G % Type  =  'a Gravitation_N_SG'

    call G % Gravitation_N_H_Form % Initialize &
           ( A, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, AssociateFieldsOption, UnitOption, &
             VectorIndicesOption, nFieldsOption, IgnorabilityOption )

    G % GravitationalConstant  =  GravitationalConstant

    !-- Source

    allocate ( G % Source )
    associate ( S  =>  G % Source )
    call S % Initialize &
           ( A, &
             FieldOption = [ 'Source' ], &
             NameOption = 'GravitationSource', &
             DeviceMemoryOption = DeviceMemoryOption, &
             nFieldsOption = 1, &
             IgnorabilityOption = A % IGNORABILITY + 1 )
    end associate !-- S

    !-- Solution and SolutionGradient

    allocate ( G % Solution )
    associate &
      ( S  =>  G % Solution )

    call S % Initialize &
           ( G, &
             iaSelected = [ G % POTENTIAL ], &
             NameOption = 'GravitationSolution', &
             IgnorabilityOption = A % IGNORABILITY + 1 )

    allocate ( G % SolutionGradient ( 3 ) )
    do iD  =  1, 3
      write ( DimensionNumber, fmt = '(i1.1)' ) iD
      associate &
        ( SG  =>  G % SolutionGradient ( iD ) )
      call SG % Initialize &
             ( G, &
               iaSelected = [ G % POTENTIAL_GRADIENT_D ( iD ) ], &
               NameOption = 'GravitationGradient_' // DimensionNumber, &
               IgnorabilityOption = A % IGNORABILITY + 1 )
      end associate !-- iD
    end do !-- iD

    select type ( A  =>  S % Atlas )
    class is ( Atlas_SCG_CC_Form )
      call S % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
      call S % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
      call S % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Gravitation_N_SG__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeAllocate_N_SG', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- S

    !-- Gradient

    allocate ( G % Gradient )
    associate ( Gt  =>  G % Gradient )
    call Gt % Initialize ( G, G % Solution )
    end associate !-- Gt

    !-- Poisson

    MaxDegree  =  12
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( G % Poisson )
    associate ( PA  =>  G % Poisson )
    call PA % Initialize ( G, 'MULTIPOLE', MaxDegree )
    end associate !-- PA

  end subroutine InitializeAllocate_N_SG


  subroutine Show_FS ( FS )

    class ( Gravitation_N_SG_Form ), intent ( in ) :: &
      FS

    integer ( KDI ) :: &
      iD

    call FS % Gravitation_N_H_Form % Show ( )

    call Show ( FS % GravitationalConstant, 'GravitationalConstant' )

    call FS % Source % Show ( )
    call FS % Solution % Show ( )
    do iD = 1, 3
      call FS % SolutionGradient ( iD ) % Show ( )
    end do !-- iD

    call FS % Poisson % Show ( )

  end subroutine Show_FS


  subroutine Solve ( G, F, iBaryonMass, iBaryonDensity, T_Option )

    class ( Gravitation_N_SG_Form ), intent ( inout ) :: &
      G
    class ( FieldSetForm ), intent ( in ) :: &
      F  !-- Fluid
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass, &
      iBaryonDensity
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC, &
      iD
    type ( TimerForm ), pointer :: &
      T_S, &  !-- Source
      T_P, &  !-- Poisson
      T_G     !-- Gradient

    !-- Source

    associate &
      (  S  =>  G % Source, &
        nC  =>  G % Atlas % nCharts )

    if ( present ( T_Option ) ) then
      T_S  =>  PROGRAM_HEADER % Timer &
                 ( Handle = G % iTimer_S, &
                   Name = trim ( T_Option % Name ) // '_Src', &
                   Level = T_Option % Level + 1 )
    else
      T_S  =>  null ( )
    end if

    if ( associated ( T_S ) ) call T_S % Start ( )
    do iC  =  1, nC
      associate &
        ( FV  =>  F % Storage ( iC ) % Value, &
          SV  =>  S % Storage ( iC ) % Value )

      call ComputeSourceKernel &
             ( M = FV ( :, iBaryonMass ), &
               N = FV ( :, iBaryonDensity ), &
               G = G % GravitationalConstant, &
               S = SV ( :, S % iaSelected ( 1 ) ), &
               UseDeviceOption = S % DeviceMemory )

      end associate !-- FV, etc.
    end do !-- iC
    if ( associated ( T_S ) ) call T_S % Stop ( )

    end associate !-- S, etc.

    !-- Poisson

    associate ( P  =>  G % Poisson )

    if ( present ( T_Option ) ) then
      T_P  =>  P % Timer ( Level = T_Option % Level + 1 )
      call T_P % Start ( )
      call P % Solve ( G % Solution, G % Source, T_Option = T_P )
      call T_P % Stop ( )
    else
      call P % Solve ( G % Solution, G % Source )
    end if

    end associate !-- P

    !-- Gradient

    associate &
      (  Gt  =>  G % Gradient, &
        nC   =>  G % Atlas % nCharts )

    if ( present ( T_Option ) ) then
      T_G  =>  PROGRAM_HEADER % Timer &
                 ( Handle = G % iTimer_G, &
                   Name = trim ( T_Option % Name ) // '_Grdnt', &
                   Level = T_Option % Level + 1 )
    else
      T_G  =>  null ( )
    end if

    if ( associated ( T_G ) ) call T_G % Start ( )
    do iC  =  1, nC
      associate ( nD  =>  G % Atlas % Chart ( iC ) % Element % nDimensions )
      do iD  =  1, nD
        associate ( SG  =>  G % SolutionGradient ( iD ) )
        call Gt % Compute ( iD )
        call Gt % Copy ( SG )
        end associate !-- SG
      end do !-- iD
      end associate !-- nD
    end do !-- iC
    if ( associated ( T_G ) ) call T_G % Stop ( )

    end associate !-- Gt, etc.

  end subroutine Solve


  impure elemental subroutine Finalize ( G )

    type ( Gravitation_N_SG_Form ), intent ( inout ) :: &
      G

    if ( allocated ( G % Poisson ) ) &
      deallocate ( G % Poisson )
    if ( allocated ( G % Gradient ) ) &
      deallocate ( G % Gradient )
    if ( allocated ( G % SolutionGradient ) ) &
      deallocate ( G % SolutionGradient )
    if ( allocated ( G % Solution ) ) &
      deallocate ( G % Solution )
    if ( allocated ( G % Source ) ) &
      deallocate ( G % Source )

  end subroutine Finalize


end module Gravitation_N_SG__Form
