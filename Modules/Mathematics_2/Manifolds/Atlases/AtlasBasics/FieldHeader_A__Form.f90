!-- FieldHeader_A__Form handles metadata for a set of related fields on 
!   an Atlas.

module FieldHeader_A__Form

  !-- FieldHeader_Atlas__Form

  use Basics
  use AtlasHeader_Form

  implicit none
  private

  type, public :: FieldHeader_A_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    logical ( KDL ) :: &
      UsePinnedMemory
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
    class ( AtlasHeaderForm ), pointer :: &
      Atlas => null ( )
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    final :: &
      Finalize
  end type FieldHeader_A_Form


contains


  subroutine Initialize_H &
               ( FA, A, NameShort, UsePinnedMemoryOption, IgnorabilityOption )

    class ( FieldHeader_A_Form ), intent ( inout ) :: &
      FA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    FA % IGNORABILITY = A % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FA % IGNORABILITY = IgnorabilityOption

    if ( FA % Type == '' ) &
      FA % Type = 'a Field_A' 
    
    FA % UsePinnedMemory = .false.
    if ( present ( UsePinnedMemoryOption ) ) &
      FA % UsePinnedMemory = UsePinnedMemoryOption
    
    FA % Name = trim ( NameShort ) // '_' // trim ( A % Name ) 

    call Show ( 'Initializing ' // trim ( FA % Type ), FA % IGNORABILITY )
    call Show ( FA % Name, 'Name', FA % IGNORABILITY )
   
    FA % NameShort = NameShort
    call Show ( FA % NameShort, 'NameShort', FA % IGNORABILITY )

    FA % Atlas => A

    ! call FA % SetField ( )

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FA )

    type ( FieldHeader_A_Form ), intent ( inout ) :: &
      FA

    nullify ( FA % Atlas )

    call Show ( 'Finalizing ' // trim ( FA % Type ), FA % IGNORABILITY )
    call Show ( FA % Name, 'Name', FA % IGNORABILITY )
   
  end subroutine Finalize


end module FieldHeader_A__Form
