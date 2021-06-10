module Gravitation_NSG_A__Form

  !-- Gravitation_NewtonSelfGravity_Atlas_Form

  use Basics
  use Mathematics
  use Gravitation_NH_C__Form
  use Gravitation_NH_A__Form

  implicit none
  private

  type, public, extends ( Gravitation_NH_A_Form ) :: Gravitation_NSG_A_Form
    type ( FieldSet_A_Form ), allocatable :: &
      Source_A, &
      Solution_A
    type ( Poisson_A_Form ), allocatable :: &
      Poisson_A
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass ( GA ) :: &
      SetStream
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Gravitation_NSG_A_Form


contains


  subroutine InitializeAllocate_FS &
               ( FSA, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Gravitation_NSG_A_Form ), intent ( inout ), target :: &
      FSA
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
      DevicesCommunicateOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      MaxDegree

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a Gravitation_NSG_A'

    call FSA % Gravitation_NH_A_Form % Initialize &
               ( A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    !-- Source_A

    allocate ( FSA % Source_A )
    associate ( SA  =>  FSA % Source_A )
    call SA % Initialize &
           ( A, &
             FieldOption = [ 'Source' ], &
             NameOption = 'GravitationSource', &
             DeviceMemoryOption = DeviceMemoryOption, &
             nFieldsOption = 1 )
    end associate !-- SA

    !-- Solution_A

    allocate ( FSA % Solution_A )
    associate &
      ( SA  =>  FSA % Solution_A )
    select type ( GC  =>  FSA % FieldSet_C ( 1 ) % Element )
      class is ( Gravitation_NH_C_Form )

    call SA % Initialize &
           ( FSA, &
             iaSelected = [ GC % POTENTIAL ], &
             NameOption = 'GravitationSolution' )

    select type ( A  =>  SA % Atlas )
    class is ( Atlas_SCG_CC_Form )
      call SA % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
      call SA % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Gravitation_NSG_A__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeAllocate_FS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end select !-- GC
    end associate !-- SA

    !-- Poisson_A

    MaxDegree  =  8
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( FSA % Poisson_A )
    associate ( PA  =>  FSA % Poisson_A )
    call PA % Initialize ( FSA, 'MULTIPOLE', MaxDegree )
    end associate !-- PA

  end subroutine InitializeAllocate_FS


  subroutine Show_FS ( FSA )

    class ( Gravitation_NSG_A_Form ), intent ( in ) :: &
      FSA

    call FSA % Gravitation_NH_A_Form % Show ( )
    call FSA % Source_A % Show ( )
    call FSA % Solution_A % Show ( )
    call FSA % Poisson_A % Show ( )

  end subroutine Show_FS


  subroutine SetStream ( SA, GA, iaAdditionalOption )

    class ( Stream_A_Form ), intent ( inout ) :: &
      SA
    class ( Gravitation_NSG_A_Form ), intent ( in ) :: &
      GA
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaAdditionalOption

    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
      class is ( Gravitation_NH_C_Form )

    call GA % Gravitation_NH_A_Form % SetStream &
           ( SA, iaAdditionalOption = [ GC % POTENTIAL ] )

    end select !-- GC

  end subroutine SetStream


  subroutine Compute ( GA )

    class ( Gravitation_NSG_A_Form ), intent ( inout ) :: &
      GA

  end subroutine Compute


  impure elemental subroutine Finalize ( GA )

    type ( Gravitation_NSG_A_Form ), intent ( inout ) :: &
      GA

    if ( allocated ( GA % Poisson_A ) ) &
      deallocate ( GA % Poisson_A )
    if ( allocated ( GA % Solution_A ) ) &
      deallocate ( GA % Solution_A )
    if ( allocated ( GA % Source_A ) ) &
      deallocate ( GA % Source_A )

  end subroutine Finalize


end module Gravitation_NSG_A__Form
