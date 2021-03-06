module Sources_RM__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_RM  = 8, &
      N_VECTORS_RM = 3

  type, public, extends ( Sources_C_Form ) :: Sources_RM_Form
    integer ( KDI ) :: &
      N_FIELDS_RM    = N_FIELDS_RM, &
      N_VECTORS_RM   = N_VECTORS_RM, &
      INTERACTIONS_J = 0, &
      INTERACTIONS_N = 0
    integer ( KDI ), dimension ( 3 ) :: &
      INTERACTIONS_H_D = 0, &
      CURVILINEAR_S_D  = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_SRM
    generic, public :: &
      Initialize => InitializeAllocate_SRM
    final :: &
      Finalize
  end type Sources_RM_Form

    private :: &
      InitializeBasics, &
      SetUnits

contains


  subroutine InitializeAllocate_SRM &
               ( SRM, RM, TimeUnit, EnergyUnit, VariableOption, VectorOption, &
                 NameOption, ClearOption, PinnedOption, UnitOption, &
                 VectorIndicesOption )

    class ( Sources_RM_Form ), intent ( inout ) :: &
      SRM
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit, &
      EnergyUnit   
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption, &
      PinnedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
        VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector

    SRM % N_FIELDS_C = RM % N_CONSERVED

    call InitializeBasics &
           ( SRM, RM, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call SetUnits ( VariableUnit, SRM, RM, TimeUnit, EnergyUnit )

    call SRM % Sources_C_Form % Initialize &
           ( RM, TimeUnit, RM % iaConserved, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, PinnedOption = PinnedOption, &
             UnitOption = VariableUnit, VectorIndicesOption = VectorIndices )

  end subroutine InitializeAllocate_SRM


  impure elemental subroutine Finalize ( SRM )

    type ( Sources_RM_Form ), intent ( inout ) :: &
      SRM

    !-- Trigger finalization of parent

  end subroutine Finalize


  subroutine InitializeBasics &
               ( SRM, RM, Variable, Vector, Name, VariableUnit, &
                 VectorIndices, VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( Sources_RM_Form ), intent ( inout ) :: &
      SRM
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    character ( LDF ), intent ( out ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    !-- FIXME: intent(out) here caused ICE with Intel Compiler 15
    !          Temporarily set to intent(inout)
    !type ( Integer_1D_Form ), dimension ( : ), allocatable, &
    !  intent ( out ) :: &
    type ( Integer_1D_Form ), dimension ( : ), allocatable, &
      intent ( inout ) :: &
        VectorIndices
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      iV, &  !-- iVector
      oF, &  !-- oField
      oV     !-- oVector
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    if ( SRM % Type == '' ) &
      SRM % Type = 'a Sources_RM'

    Name = 'Sources_RM'
    if ( present ( NameOption ) ) &
      Name = NameOption

    !-- variable indices

    oF = SRM % N_FIELDS_C
    if ( SRM % N_FIELDS == 0 ) &
      SRM % N_FIELDS = oF + SRM % N_FIELDS_RM

    SRM % INTERACTIONS_J    =  oF + 1
    SRM % INTERACTIONS_N    =  oF + 2
    SRM % INTERACTIONS_H_D  =  oF + [ 3, 4, 5 ]
    SRM % CURVILINEAR_S_D   =  oF + [ 6, 7, 8 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( SRM % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( oF + 1 : oF + SRM % N_FIELDS_RM ) &
      = [ 'Interactions_J    ', &
          'Interactions_N    ', &
          'Interactions_H_D_1', &
          'Interactions_H_D_2', &
          'Interactions_H_D_3', &
          'Curvilinear_S_D_1 ', &
          'Curvilinear_S_D_2 ', &
          'Curvilinear_S_D_3 ' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( SRM % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    oV = SRM % N_VECTORS_C
    if ( SRM % N_VECTORS == 0 ) &
      SRM % N_VECTORS = oV + SRM % N_VECTORS_RM

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( SRM % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( oV + 1 : oV + SRM % N_VECTORS_RM ) &
      = [ 'Div_F_S_D       ', &
          'Curvilinear_S_D ', &
          'Interactions_H_D' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = oV + SRM % N_VECTORS_RM + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( SRM % N_VECTORS ) )
    end if

    do iD = 1, 3
      call Search ( RM % iaConserved, &
                    RM % CONSERVED_MOMENTUM_D ( iD ), &
                    iMomentum ( iD ) )
    end do

    call VectorIndices ( oV + 1 ) % Initialize ( iMomentum )
    call VectorIndices ( oV + 2 ) % Initialize ( SRM % CURVILINEAR_S_D )
    call VectorIndices ( oV + 3 ) % Initialize ( SRM % INTERACTIONS_H_D )

  end subroutine InitializeBasics


  subroutine SetUnits ( VariableUnit, SRM, RM, TimeUnit, EnergyUnit )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      VariableUnit
    class ( Sources_RM_Form ), intent ( in ) :: &
      SRM
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit, &
      EnergyUnit

    integer ( KDI ) :: &
      iD

    VariableUnit ( SRM % INTERACTIONS_J )  &
      =  RM % Unit ( RM % COMOVING_ENERGY )  /  TimeUnit
    VariableUnit ( SRM % INTERACTIONS_N )  &
      =  RM % Unit ( RM % COMOVING_ENERGY )  /  EnergyUnit  /  TimeUnit

    do iD = 1, 3
      VariableUnit ( SRM % CURVILINEAR_S_D ( iD ) )  &
        =  RM % Unit ( RM % CONSERVED_MOMENTUM_D ( iD ) )  /  TimeUnit
      VariableUnit ( SRM % INTERACTIONS_H_D ( iD ) )  &
        =  RM % Unit ( RM % CONSERVED_MOMENTUM_D ( iD ) )  /  TimeUnit
    end do

  end subroutine SetUnits


end module Sources_RM__Form
