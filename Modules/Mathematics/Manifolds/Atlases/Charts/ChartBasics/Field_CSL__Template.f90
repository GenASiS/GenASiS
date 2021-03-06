!-- Field_CSL is a template for a set of fields on a Chart_SL.

module Field_CSL__Template

  !-- Field_ChartSingleLevel_Template

  use Basics
  use Chart_Template
  use FieldChart_Template

  implicit none
  private

  type, public, extends ( FieldChartTemplate ), abstract :: Field_CSL_Template
    integer ( KDI ) :: &
      nValues
    class ( StorageForm ), allocatable :: &
      Field, &
      FieldOutput
    type ( MessageIncoming_1D_R_Form ), allocatable :: &
      IncomingFace_L_R, &
      IncomingFace_R_L, &
      IncomingEdge_LL_RR, &
      IncomingEdge_RR_LL, &
      IncomingEdge_LR_RL, &
      IncomingEdge_RL_LR
    type ( MessageOutgoing_1D_R_Form ), allocatable :: &
      OutgoingFace_L_R, &
      OutgoingFace_R_L, &
      OutgoingEdge_LL_RR, &
      OutgoingEdge_RR_LL, &
      OutgoingEdge_LR_RL, &
      OutgoingEdge_RL_LR
  contains
    procedure, public, pass :: &
      InitializeTemplate_CSL
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_CSL
    procedure, public, pass :: &
      UpdateHost => UpdateHost_CSL
    procedure, private, pass :: &
      Clear_FC
    procedure, public, pass :: &
      FinalizeTemplate_CSL
    !-- FIXME: This should be automatically inherited from FieldChartTemplate
    !          but XL compiler got confused
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type Field_CSL_Template
  
    !-- FIXME: see above
    abstract interface 
      subroutine SF ( FC )
        import Field_CSL_Template
        class ( Field_CSL_Template ), intent ( inout ) :: &
          FC
      end subroutine
    end interface

contains


  subroutine InitializeTemplate_CSL &
               ( FC, C, NameShort, UsePinnedMemory, nValues, &
                 IgnorabilityOption )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ) :: &
      UsePinnedMemory
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FC % Type == '' ) &
      FC % Type = 'a Field_CSL' 

    FC % nValues = nValues

    call FC % InitializeTemplate &
           ( C, NameShort, UsePinnedMemory, IgnorabilityOption )

  end subroutine InitializeTemplate_CSL
  
  
  subroutine AllocateDevice_CSL ( FC, AssociateVariablesOption )
  
    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC
    logical ( KDL ), intent ( in ), optional :: &
      AssociateVariablesOption
    
    call FC % Field % AllocateDevice ( AssociateVariablesOption )
  
  end subroutine AllocateDevice_CSL
  

  subroutine UpdateHost_CSL ( FC )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC

    call FC % Field % UpdateHost ( )

  end subroutine UpdateHost_CSL


  subroutine Clear_FC ( FC )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC

    call FC % Field % Clear ( )

  end subroutine Clear_FC


  impure elemental subroutine FinalizeTemplate_CSL ( FC )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC

    if ( allocated ( FC % FieldOutput ) ) &
      deallocate ( FC % FieldOutput )
    if ( allocated ( FC % Field ) ) &
      deallocate ( FC % Field )

    call FC % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_CSL


end module Field_CSL__Template
