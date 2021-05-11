module FieldSet_A__Form

  !-- FieldSet_Atlas__Form

  use Basics
  use Manifolds
  use FieldSet_C__Form

  implicit none
  private

  type, public :: FieldSet_A_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( FieldSet_C_Element ), dimension ( : ), allocatable :: &
      FieldSet_C
  contains
    procedure, private, pass :: &
      InitializeAllocate_FS
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeAllocate_FS, InitializeClone
    procedure, public, pass :: &
      Clear => Clear_FSA
    procedure, public, pass :: &
      Show => Show_FSA
    final :: &
      Finalize
  end type FieldSet_A_Form

  type, public :: FieldSet_A_Element
    !-- FieldSet_Atlas_Element
    class ( FieldSet_A_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type FieldSet_A_Element


contains


  subroutine InitializeAllocate_FS &
               ( FSA, A, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( FieldSet_A_Form ), intent ( inout ), target :: &
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
      iC  !-- iChart

    FSA % IGNORABILITY  =  A % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FSA % IGNORABILITY  =  IgnorabilityOption

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a FieldSet_A'

    FSA % Name  =  'Fields'
    if ( present ( NameOption ) ) &
      FSA % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( FSA % Type ), A % IGNORABILITY )
    call Show ( FSA % Name, 'Name', A % IGNORABILITY )

    FSA % Atlas  =>  A

    if ( .not. allocated ( FSA % FieldSet_C ) ) then
      associate ( nC  =>  A % nCharts )
      allocate ( FSA % FieldSet_C ( nC ) )
      do iC  =  1, nC
        allocate ( FSA % FieldSet_C ( iC ) % Element )
        associate &
          ( FSC  =>  FSA % FieldSet_C ( iC ) % Element, &
              C  =>    A %      Chart ( iC ) % Element )
        call FSC % Initialize &
               ( C, FieldOption, VectorOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, UnitOption, VectorIndicesOption, &
                 nFieldsOption, IgnorabilityOption )
        end associate !-- FSC
      end do !-- iC
      end associate !-- nC
    end if !-- allocated FSA % FieldSet_C

  end subroutine InitializeAllocate_FS


  subroutine InitializeClone &
               ( FSA_T, FSA_S, iaSelected, NameOption, IgnorabilityOption )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      FSA_T  !-- FSA_Target
    class ( FieldSet_A_Form ), intent ( in ) :: &
      FSA_S  !-- FSA_Source
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaSelected
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC  !-- iChart

    associate ( nC  =>  FSA_S % Atlas % nCharts )
    allocate ( FSA_T % FieldSet_C ( nC ) )
 
    call FSA_T % Initialize &
           ( FSA_S % Atlas, &
             NameOption = NameOption, &
             IgnorabilityOption = IgnorabilityOption )

    do iC  =  1, nC
      allocate ( FSA_T % FieldSet_C ( iC ) % Element )
      associate &
        ( FSC_T  =>  FSA_T % FieldSet_C ( iC ) % Element, &
          FSC_S  =>  FSA_S % FieldSet_C ( iC ) % Element )
      call FSC_T % Initialize &
             ( FSC_S, iaSelected, NameOption, IgnorabilityOption )
      end associate !-- FSC_T, etc.
    end do !-- iC

    end associate !-- nC

  end subroutine InitializeClone


  subroutine Clear_FSA ( FSA )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      FSA

   integer ( KDI ) :: &
     iC  !-- iC

    associate ( A  =>  FSA % Atlas )
    do iC  =  1, A % nCharts
      if ( allocated ( FSA % FieldSet_C ( iC ) % Element ) ) then
        associate ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
        call FSC % Clear ( )
        end associate !-- FSC
      end if  
    end do !-- iC
    end associate  !-- A

  end subroutine Clear_FSA


  subroutine Show_FSA ( FSA )

    class ( FieldSet_A_Form ), intent ( in ) :: &
      FSA

   integer ( KDI ) :: &
     iC  !-- iC
   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( FSA % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', FSA % IGNORABILITY )

    associate ( A  =>  FSA % Atlas )

    call Show ( FSA % Name, 'Name',  FSA % IGNORABILITY )
    call Show (   A % Name, 'Atlas', FSA % IGNORABILITY )

    do iC  =  1, A % nCharts
      if ( allocated ( FSA % FieldSet_C ( iC ) % Element ) ) then
        associate ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
        call FSC % Show ( )
        end associate !-- FSC
      end if  
    end do !-- iC

    end associate  !-- A

  end subroutine Show_FSA


  impure elemental subroutine Finalize ( FSA )

    type ( FieldSet_A_Form ), intent ( inout ) :: &
      FSA

    if ( allocated ( FSA % FieldSet_C ) ) &
      deallocate ( FSA % FieldSet_C )  

    nullify ( FSA % Atlas )

    call Show ( 'Finalizing ' // trim ( FSA % Type ), FSA % IGNORABILITY )
    call Show ( FSA % Name, 'Name', FSA % IGNORABILITY )

  end subroutine Finalize


  impure elemental subroutine Finalize_E ( FSE )
    
    type ( FieldSet_A_Element ), intent ( inout ) :: &
      FSE

    if ( allocated ( FSE % Element ) ) &
      deallocate ( FSE % Element )

  end subroutine Finalize_E


end module FieldSet_A__Form
