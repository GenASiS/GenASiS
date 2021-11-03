module Gravitation_N_SG__Form

  !-- Gravitation_Newton_SelfGravity__Form

  use Basics
  use Mathematics
  use Gravitation_N_H__Form

  implicit none
  private

  type, public, extends ( Gravitation_N_H_Form ) :: Gravitation_N_SG_Form
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
      InitializeAllocate_FS
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


  subroutine InitializeAllocate_FS &
               ( FS, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( Gravitation_N_SG_Form ), intent ( inout ), target :: &
      FS
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
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

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a Gravitation_N_SG'

    call FS % Gravitation_N_H_Form % Initialize &
           ( A, FieldOption, VectorOption, NameOption, &
             DeviceMemoryOption, PinnedMemoryOption, &
             DevicesCommunicateOption, AssociateFieldsOption, UnitOption, &
             VectorIndicesOption, nFieldsOption, IgnorabilityOption )

    !-- Source

    allocate ( FS % Source )
    associate ( S  =>  FS % Source )
    call S % Initialize &
           ( A, &
             FieldOption = [ 'Source' ], &
             NameOption = 'GravitationSource', &
             DeviceMemoryOption = DeviceMemoryOption, &
             nFieldsOption = 1, &
             IgnorabilityOption = A % IGNORABILITY + 1 )
    end associate !-- S

    !-- Solution and SolutionGradient

    allocate ( FS % Solution )
    associate &
      ( S  =>  FS % Solution )

    call S % Initialize &
           ( FS, &
             iaSelected = [ FS % POTENTIAL ], &
             NameOption = 'GravitationSolution', &
             IgnorabilityOption = A % IGNORABILITY + 1 )

    allocate ( FS % SolutionGradient ( 3 ) )
    do iD  =  1, 3
      write ( DimensionNumber, fmt = '(i1.1)' ) iD
      associate &
        ( SG  =>  FS % SolutionGradient ( iD ) )
      call SG % Initialize &
             ( FS, &
               iaSelected = [ FS % POTENTIAL_GRADIENT_D ( iD ) ], &
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
      call Show ( 'InitializeAllocate_FS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- S

    !-- Gradient

    allocate ( FS % Gradient )
    associate ( Gt  =>  FS % Gradient )
    call Gt % Initialize ( FS, FS % Solution )
    end associate !-- Gt

    !-- Poisson

    MaxDegree  =  12
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( FS % Poisson )
    associate ( PA  =>  FS % Poisson )
    call PA % Initialize ( FS, 'MULTIPOLE', MaxDegree )
    end associate !-- PA

  end subroutine InitializeAllocate_FS


  subroutine Show_FS ( FS )

    class ( Gravitation_N_SG_Form ), intent ( in ) :: &
      FS

    integer ( KDI ) :: &
      iD

    call FS % Gravitation_N_H_Form % Show ( )

    call FS % Source % Show ( )
    call FS % Solution % Show ( )
    do iD = 1, 3
      call FS % SolutionGradient ( iD ) % Show ( )
    end do !-- iD

    call FS % Poisson % Show ( )

  end subroutine Show_FS


  subroutine Solve ( G, F, Constant_G, iBaryonMass, iBaryonDensity )

    class ( Gravitation_N_SG_Form ), intent ( inout ) :: &
      G
    class ( FieldSetForm ), intent ( in ) :: &
      F  !-- Fluid
    real ( KDR ), intent ( in ) :: &
      Constant_G  !-- Gravitational
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass, &
      iBaryonDensity

    integer ( KDI ) :: &
      iC, &
      iD

    associate &
      (  S  =>  G % Source, &
        nC  =>  G % Atlas % nCharts )
    do iC  =  1, nC
      associate &
        ( FV  =>  F % Storage ( iC ) % Value, &
          SV  =>  S % Storage ( iC ) % Value )

      call ComputeSourceKernel &
             ( M = FV ( :, iBaryonMass ), &
               N = FV ( :, iBaryonDensity ), &
               G = Constant_G, &
               S = SV ( :, S % iaSelected ( 1 ) ), &
               UseDeviceOption = S % DeviceMemory )

      end associate !-- FV, etc.
    end do !-- iC
    end associate !-- S, etc.

    associate ( P  =>  G % Poisson )
    call P % Solve ( G % Solution, G % Source )
    end associate !-- P

    associate &
      (  Gt  =>  G % Gradient, &
        nC   =>  G % Atlas % nCharts )
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
