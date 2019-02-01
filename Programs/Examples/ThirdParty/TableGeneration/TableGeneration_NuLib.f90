program TableGeneration_NuLib

  use nulib
  use GenASiS

  implicit none

  type ( Atlas_SC_Form ), allocatable :: &
    MomentumSpace

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'TableGeneration_NuLib' )

  allocate ( MomentumSpace )
  call MomentumSpace % Initialize &
         ( 'MomentumSpace', nDimensionsOption = 1 )

  call SetEnergyGrid ( MomentumSpace )
  call PrepareNuLib ( MomentumSpace )
  call CreateTable ( )

  deallocate ( MomentumSpace )
  deallocate ( PROGRAM_HEADER )

end program TableGeneration_NuLib


subroutine SetEnergyGrid ( MS )

  use GenASiS

  type ( Atlas_SC_Form ), intent ( inout ) :: &
    MS

  integer ( KDI ) :: &
    nEnergyCells
  real ( KDR ) :: &
    EnergyScale
  real ( KDR ), dimension ( 1 ) :: &
    Scale
  type ( MeasuredValueForm ), dimension ( 1 ) :: &
    CoordinateUnit
  character ( LDL ), dimension ( 1 ) :: &
    Spacing, &
    CoordinateLabel

  call Show ( 'Creating GenASiS energy grid' )

  nEnergyCells = 20
  call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

  EnergyScale  =  10.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

  Scale            =  EnergyScale
  CoordinateUnit   =  UNIT % MEGA_ELECTRON_VOLT
  Spacing          =  'COMPACTIFIED'
  CoordinateLabel  =  'Energy'
  
  call MS % CreateChart &
         ( SpacingOption          = Spacing, &
           CoordinateLabelOption  = CoordinateLabel, &
           CoordinateSystemOption = 'SPHERICAL', &
           CoordinateUnitOption   = CoordinateUnit, &
           ScaleOption            = Scale, &
           nCellsOption           = [ nEnergyCells ], &
           nGhostLayersOption     = [ 0, 0, 0 ] )
  call MS % SetGeometry ( )

end subroutine SetEnergyGrid


subroutine PrepareNuLib ( MS )

  use nulib, only: &
    energies, &
    bin_widths, &
    bin_bottom, &
    bin_top, &
    initialize_nulib
  use GenASiS

  type ( Atlas_SC_Form ), intent ( inout ) :: &
    MS

  integer ( KDI ) :: &
    NeutrinoScheme, &
    NumberSpecies, &
    NumberGroups
  class ( GeometryFlatForm ), pointer :: &
    G

  call Show ( 'Preparing NuLib' )

  G => MS % Geometry ( )

  associate &
    ( Energy      =>  G % Value ( :, G % CENTER_U ( 1 ) ), &
      WidthLeft   =>  G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
      WidthRight  =>  G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
      EnergyUnit  =>  MS % Chart % CoordinateUnit ( 1 ) )

  call Show ( Energy, EnergyUnit, 'Energy' )
  call Show ( WidthLeft, EnergyUnit, 'WidthLeft' )
  call Show ( WidthRight, EnergyUnit, 'WidthRight' )

  !-- Initialize nulib

  NeutrinoScheme  =  1  !-- Nu_e, Nu_e_bar, Nu_x
  call PROGRAM_HEADER % GetParameter ( NeutrinoScheme, 'NeutrinoScheme' )

  NumberSpecies  =  6  !-- Required by nulib

  select type ( C => MS % Chart )
  class is ( Chart_SL_Template )
  NumberGroups  =  C % nCells ( 1 )
  end select !-- C

  call initialize_nulib &
         ( neutrino_scheme_in = NeutrinoScheme, &
           number_species_in  = NumberSpecies, &
           number_groups_in   = NumberGroups )

  !-- Set nulib module variables

  energies    =  Energy
  bin_widths  =  WidthLeft + WidthRight
  bin_bottom  =  Energy - WidthLeft
  bin_top     =  Energy + WidthRight  

  end associate !-- Energy, etc.
  nullify ( G )

end subroutine PrepareNuLib


subroutine CreateTable ( )

  use nulib, only: &
    number_groups, &
    energies, &
    bin_widths, &
    bin_bottom, &
    bin_top
  use GenASiS 

  call Show ( 'Creating NuLib table' )
  call Show ( number_groups, 'number_groups' )
  call Show ( energies, 'energies' )
  call Show ( bin_widths, 'bin_widths' )
  call Show ( bin_bottom, 'bin_bottom' )
  call Show ( bin_top, 'bin_top' )

end subroutine CreateTable
