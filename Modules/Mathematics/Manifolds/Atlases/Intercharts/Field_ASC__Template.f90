!-- Field_ASC is a template for a set of fields on an Atlas_SC.

module Field_ASC__Template

  !-- Field_AtlasSingleChart_Template

  use Basics
  use AtlasBasics
  use Charts

  implicit none
  private

  type, public, extends ( FieldAtlasTemplate ), abstract :: &
    Field_ASC_Template
      class ( FieldChartTemplate ), allocatable :: &
        Chart
  contains
    procedure, public, pass :: &
      InitializeTemplate_ASC
    procedure, public, pass :: &
      AllocateDevice_ASC_Template
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_ASC_Template
    procedure, private, pass :: &
      Clear_FA
    procedure, public, pass :: &
      FinalizeTemplate_ASC
    !-- FIXME: This should be automatically inherited from FieldAtlasTemplate
    !          but XL compiler got confused
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type Field_ASC_Template
  
    abstract interface 
      subroutine SF ( FA )
        use Basics
        import Field_ASC_Template
        class ( Field_ASC_Template ), intent ( inout ) :: &
          FA
      end subroutine
    end interface

contains


  subroutine InitializeTemplate_ASC &
               ( FA, A, NameShort, UsePinnedMemoryOption, IgnorabilityOption )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FA % Type == '' ) &
      FA % Type = 'a Field_ASC' 

    call FA % InitializeTemplate &
           ( A, NameShort, UsePinnedMemoryOption, IgnorabilityOption )

  end subroutine InitializeTemplate_ASC
  
  
  subroutine AllocateDevice_ASC_Template ( FA, AssociateVariablesOption )
    
    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA
    logical ( KDL ), intent ( in ), optional :: &
      AssociateVariablesOption
    
    select type ( FC => FA % Chart ) 
    class is ( Field_CSL_Template )
      call FC % AllocateDevice ( AssociateVariablesOption )
    class default
      call Show ( 'Field type not implemented', CONSOLE % ERROR )
      call Show ( 'AllocateDevice_ASC', 'subroutine',  CONSOLE % ERROR )
      call Show ( 'Field_ASC_Template', 'module',  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select
  
  end subroutine AllocateDevice_ASC_Template


  subroutine Clear_FA ( FA )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA

    call FA % Chart % Clear ( )

  end subroutine Clear_FA


  impure elemental subroutine FinalizeTemplate_ASC ( FA )

    class ( Field_ASC_Template ), intent ( inout ) :: &
      FA

    if ( allocated ( FA % Chart ) ) &
      deallocate ( FA % Chart )

    call FA % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_ASC


end module Field_ASC__Template
