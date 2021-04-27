module CurrentSet_CH__Form

  !-- CurrentSet_ChartHeader_Form

  use Basics
  use Manifolds

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_PRIMITIVE_CS = 0, &
      N_CONSERVED_CS = 0, &
      N_FIELDS_CS    = 6, &
      N_VECTORS_CS   = 2

  type, public :: CurrentSet_CH_Form
    integer ( KDI ) :: &
      N_PRIMITIVE_CS = N_PRIMITIVE_CS, &
      N_CONSERVED_CS = N_CONSERVED_CS, &
      N_FIELDS_CS    = N_FIELDS_CS, &
      N_VECTORS_CS   = N_VECTORS_CS
    integer ( KDI ) :: &
      FAST_EIGENSPEED_PLUS_U_1  = 0, &
      FAST_EIGENSPEED_PLUS_U_2  = 0, &
      FAST_EIGENSPEED_PLUS_U_3  = 0, &
      FAST_EIGENSPEED_MINUS_U_1 = 0, &
      FAST_EIGENSPEED_MINUS_U_2 = 0, &
      FAST_EIGENSPEED_MINUS_U_3 = 0
    integer ( KDI ), dimension ( 3 ) :: &
      FAST_EIGENSPEED_PLUS_U, &
      FAST_EIGENSPEED_MINUS_U
  contains
    procedure, public, pass :: &
      Initialize_H
    final :: &
      Finalize
  end type CurrentSet_CH_Form


contains


  subroutine Initialize_H &
               ( CSC, C, NameOption, nFieldsOption, FieldOption, UnitOption )

    class ( CurrentSet_CH_Form ), intent ( inout ) :: &
      CSC
    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    character ( * ), intent ( inout ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
    character ( * ), dimension ( : ), intent ( out ), allocatable, optional :: &
      FieldOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( out ), allocatable, &
      optional :: &
        UnitOption

    integer ( KDI ) :: &
      nFields
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    Name  =  'Current'
    if ( present ( NameOption ) ) then
      if ( NameOption  ==  '' ) then
        NameOption  =  Name
      else
        Name  =  NameOption
      end if
    end if

    !-- Field indices

    CSC % FAST_EIGENSPEED_PLUS_U_1   =  1
    CSC % FAST_EIGENSPEED_PLUS_U_2   =  2
    CSC % FAST_EIGENSPEED_PLUS_U_3   =  3
    CSC % FAST_EIGENSPEED_MINUS_U_1  =  4
    CSC % FAST_EIGENSPEED_MINUS_U_2  =  5
    CSC % FAST_EIGENSPEED_MINUS_U_3  =  6

    nFields  =  CSC % N_FIELDS_CS
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    CSC % FAST_EIGENSPEED_PLUS_U  &
      =  [ CSC % FAST_EIGENSPEED_PLUS_U_1, &
           CSC % FAST_EIGENSPEED_PLUS_U_2, &
           CSC % FAST_EIGENSPEED_PLUS_U_3 ]
    CSC % FAST_EIGENSPEED_MINUS_U  &
      =  [ CSC % FAST_EIGENSPEED_MINUS_U_1, &
           CSC % FAST_EIGENSPEED_MINUS_U_2, &
           CSC % FAST_EIGENSPEED_MINUS_U_3 ]

    !-- Field names

    allocate ( Field ( nFields ) )

    Field ( 1 : CSC % N_FIELDS_CS ) &
      = [ 'FastEigenspeedPlus_U_1 ', &
          'FastEigenspeedPlus_U_2 ', &
          'FastEigenspeedPlus_U_3 ', &
          'FastEigenspeedMinus_U_1', &
          'FastEigenspeedMinus_U_2', &
          'FastEigenspeedMinus_U_3' ]
          
    if ( present ( FieldOption ) ) &
      allocate ( FieldOption, source = Field )

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( CSC )

    type ( CurrentSet_CH_Form ), intent ( inout ) :: &
      CSC

!    if ( allocated ( GC % FieldSet ) ) &
!      deallocate ( GC % FieldSet )

  end subroutine Finalize

  
end module CurrentSet_CH__Form
