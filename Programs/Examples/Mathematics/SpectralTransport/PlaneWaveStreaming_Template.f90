module PlaneWaveStreaming_Template

  use Basics
  use Mathematics
  use RadiationMoments_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( Integrator_C_Template ), abstract :: &
    PlaneWaveStreamingTemplate
      type ( RadiationMoments_BSLL_ASC_CSLD_Form ), allocatable :: &
        RadiationMoments_BSLL_ASC_CSLD, &
        Reference_BSLL_ASC_CSLD, &
        Difference_BSLL_ASC_CSLD
  contains 
    procedure, public, pass :: &
      InitializeTemplate_PWS
    procedure, public, pass :: &
      FinalizeTemplate_PWS
  end type PlaneWaveStreamingTemplate

contains


  subroutine InitializeTemplate_PWS ( PWS, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nEnergyCells

    if ( PWS % Type == '' ) &
      PWS % Type = 'a PlaneWaveStreaming'

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: PWS % PositionSpace )
    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( PWS % Geometry_ASC )
    associate ( GA => PWS % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: PWS % MomentumSpace )
    select type ( MS => PWS % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, Name )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    nEnergyCells = 4
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( CoordinateSystemOption = 'SPHERICAL', &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Radiation

    allocate ( PWS % RadiationMoments_BSLL_ASC_CSLD )
    associate ( RMB => PWS % RadiationMoments_BSLL_ASC_CSLD )
    call RMB % Initialize ( MS, 'GENERIC' )

    !-- Cleanup

    end associate !-- RMB
    end select !-- MS
    end select !-- PS

  end subroutine InitializeTemplate_PWS


  subroutine FinalizeTemplate_PWS ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS

    if ( allocated ( PWS % Difference_BSLL_ASC_CSLD ) ) &
      deallocate ( PWS % Difference_BSLL_ASC_CSLD )
    if ( allocated ( PWS % Reference_BSLL_ASC_CSLD ) ) &
      deallocate ( PWS % Reference_BSLL_ASC_CSLD )
    if ( allocated ( PWS % RadiationMoments_BSLL_ASC_CSLD ) ) &
      deallocate ( PWS % RadiationMoments_BSLL_ASC_CSLD )

    call PWS % FinalizeTemplate_C ( )

  end subroutine FinalizeTemplate_PWS


end module PlaneWaveStreaming_Template
