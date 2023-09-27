module Gravitation_N_H__Form

  !-- Gravitation_Newton_Header__Form

  use Basics
  use Mathematics
  use Gravitation_G__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_N  = 4, &
      N_VECTORS_N = 1

  type, public, extends ( Gravitation_G_Form ) :: Gravitation_N_H_Form
    integer ( KDI ) :: &
      N_FIELDS_N = N_FIELDS_N, &
      N_VECTORS_N = N_VECTORS_N
    integer ( KDI ) :: &
      POTENTIAL = 0, &
      POTENTIAL_GRADIENT_D_1 = 0, &
      POTENTIAL_GRADIENT_D_2 = 0, &
      POTENTIAL_GRADIENT_D_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      POTENTIAL_GRADIENT_D = 0
    integer ( KDI ) :: &
      iTimer = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, public, pass ( G ) :: &
      SetStream
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      Solve
    final :: &
      Finalize
  end type Gravitation_N_H_Form


contains


  subroutine InitializeAllocate_FS &
               ( FS, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, AssociateFieldsOption, &
                 UnitOption, VectorIndicesOption, nFieldsOption, &
                 IgnorabilityOption )

    class ( Gravitation_N_H_Form ), intent ( inout ), target :: &
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
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      iV, &
      oF, &  !-- oField
      oV, &  !-- oVector
      nFields, &
      nVectors
    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field, &
      Vector

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a Gravitation_N_H' 
    
    Name  =  'Gravitation'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  FS % N_FIELDS_F

    FS % POTENTIAL               =  oF + 1
    FS % POTENTIAL_GRADIENT_D_1  =  oF + 2
    FS % POTENTIAL_GRADIENT_D_2  =  oF + 3
    FS % POTENTIAL_GRADIENT_D_3  =  oF + 4

    nFields  =  oF  +  FS % N_FIELDS_N
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    FS % POTENTIAL_GRADIENT_D  =  [ FS % POTENTIAL_GRADIENT_D_1, &
                                    FS % POTENTIAL_GRADIENT_D_2, &
                                    FS % POTENTIAL_GRADIENT_D_3 ]

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + FS % N_FIELDS_N ) &
      = [ 'Potential            ', &
          'PotentialGradient_D_1', &
          'PotentialGradient_D_2', &
          'PotentialGradient_D_3' ]

    !-- Units

    !-- Vector indices

    oV  =  FS % N_VECTORS_F

    if ( present ( VectorIndicesOption ) ) then
      nVectors  =  size ( VectorIndicesOption )
      allocate ( VectorIndices ( nVectors ) )
      do iV  =  oV  +  FS % N_VECTORS_N  +  1,  nVectors 
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do !-- iV
    else
      nVectors  =  oV  +  FS % N_VECTORS_N
      allocate ( VectorIndices ( nVectors ) )
    end if

    call VectorIndices ( oV + 1 ) % Initialize ( FS % POTENTIAL_GRADIENT_D )

    !-- Vector names

    if ( present ( VectorOption ) ) then
      allocate ( Vector, source = VectorOption )
    else
      allocate ( Vector ( nVectors ) )
    end if !-- FieldOption

    Vector ( oV  +  1 : oV  +  FS % N_VECTORS_N ) &
      = [ 'PotentialGradient_D' ]

    !-- Geometry_F

    call FS % Gravitation_G_Form % Initialize &
           ( A, &
             FieldOption = Field, &
             VectorOption = VectorOption, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             AssociateFieldsOption = AssociateFieldsOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_FS


  subroutine SetStream ( S, G, iaAdditionalOption )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Gravitation_N_H_Form ), intent ( in ) :: &
      G
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaAdditionalOption

    call G % Geometry_F_Form % SetStream &
           ( S, iaAdditionalOption = [ G % POTENTIAL, &
                                       G % POTENTIAL_GRADIENT_D ] )

  end subroutine SetStream


  function Timer ( G, Level ) result ( T )

    class ( Gravitation_N_H_Form ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = G % iTimer, &
               Name = trim ( G % Name ) // '_Slv', &
               Level = Level )

  end function Timer


  subroutine Solve ( G, F, iBaryonMass, iBaryonDensity, T_Option )

    class ( Gravitation_N_H_Form ), intent ( inout ) :: &
      G
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      F  !-- Fluid
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass, &
      iBaryonDensity
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    ! call Show ( 'Solve should be overridden', CONSOLE % WARNING )
    ! call Show ( G % Name, 'Name', CONSOLE % WARNING )
    ! call Show ( 'Gravitation_N_H__Form', 'module', CONSOLE % WARNING )
    ! call Show ( 'Solve', 'subroutine', CONSOLE % WARNING )

  end subroutine Solve


  impure elemental subroutine Finalize ( G )

    type ( Gravitation_N_H_Form ), intent ( inout ) :: &
      G

  end subroutine Finalize


end module Gravitation_N_H__Form
