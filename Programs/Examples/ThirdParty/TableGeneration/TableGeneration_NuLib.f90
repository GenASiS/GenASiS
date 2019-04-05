program TableGeneration_NuLib

  use GenASiS

  implicit none

  character ( LDL ), parameter :: outdir = "./"
  character ( LDL ), parameter :: base = "NuLib"
  character ( LDL ), parameter :: vnum = "1.0"
  character ( LDF ), parameter :: &
    eos_table_filename = "./LS220_234r_136t_50y_analmu_20091212_SVNr26.h5"
  integer ( KDI ),   parameter :: table_size_D = 82 ! Mass Density
  integer ( KDI ),   parameter :: table_size_T = 65 ! Temperature
  integer ( KDI ),   parameter :: table_size_Y = 51 ! Electron Fraction
  integer ( KDI ),   parameter :: table_size_X = 61 ! Electron Degeneracy
  real ( KDR ),      parameter :: min_LogD =  6.0_KDR
  real ( KDR ),      parameter :: max_LogD = 15.5_KDR
  real ( KDR ),      parameter :: min_LogT = log10( 0.050_KDR )
  real ( KDR ),      parameter :: max_LogT = log10( 150.0_KDR )
  real ( KDR ),      parameter :: min_Y = 0.035_KDR
  real ( KDR ),      parameter :: max_Y = 0.550_KDR
  real ( KDR ),      parameter :: min_LogX = log10( 0.100_KDR )
  real ( KDR ),      parameter :: max_LogX = log10( 100.0_KDR )

  logical ( KDL ) :: &
    include_NES, &
    include_Pair
  integer ( KDI ) :: &
    number_local_species
  real ( KDR ), dimension( : ), allocatable :: &
    table_D, &
    table_T, &
    table_Y, &
    table_X, &
    eos_variables
  real ( KDR ), dimension( :, :, :, :, : ), allocatable :: &
    table_emission, &
    table_absopacity, &
    table_scatopacity, &
    table_nes_phi_0, &
    table_nes_phi_1
  real ( KDR ), dimension( :, :, :, :, :, : ), allocatable :: &
    table_pair_phi_0, &
    table_pair_phi_1
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


contains


subroutine SetEnergyGrid ( MS )

  type ( Atlas_SC_Form ), intent ( inout ) :: &
    MS

  integer ( KDI ) :: &
    nEnergyCells
  type ( MeasuredValueForm ), dimension ( 1 ) :: &
    CoordinateUnit
  character ( LDL ) :: &
    EnergyGrid
  character ( LDL ), dimension ( 1 ) :: &
    CoordinateLabel

  nEnergyCells = 18
  call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

  CoordinateUnit   =  UNIT % MEGA_ELECTRON_VOLT
  CoordinateLabel  =  'Energy'

  EnergyGrid = 'GEOMETRIC'
  call PROGRAM_HEADER % GetParameter ( EnergyGrid, 'EnergyGrid' )

  call Show ( 'Selecting EnergyGrid' )
  call Show ( EnergyGrid, 'EnergyGrid' )

  select case ( trim ( EnergyGrid ) )
  case ( 'NU_LIB' )
    call SetEnergyGridNuLib &
           ( MS, CoordinateLabel, CoordinateUnit, nEnergyCells )
  case ( 'GEOMETRIC' )
    call SetEnergyGridGeometric &
           ( MS, CoordinateLabel, CoordinateUnit, nEnergyCells )
  case ( 'COMPACTIFIED' )
    call SetEnergyGridCompactified &
           ( MS, CoordinateLabel, CoordinateUnit, nEnergyCells )
  case default
    call Show ( 'EnergyGrid not recognized', CONSOLE % ERROR )
    call Show ( EnergyGrid, 'EnergyGrid', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )
  end select !-- EnergyGrid

end subroutine SetEnergyGrid


subroutine SetEnergyGridNuLib &
             ( MS, CoordinateLabel, CoordinateUnit, nEnergyCells )

  type ( Atlas_SC_Form ), intent ( inout ) :: &
    MS
  character ( LDL ), dimension ( 1 ), intent ( in ) :: &
    CoordinateLabel
  type ( MeasuredValueForm ), dimension ( 1 ), intent ( in ) :: &
    CoordinateUnit
  integer ( KDI ), intent ( in ) :: &
    nEnergyCells

  integer ( KDI ) :: &
    iE
  real ( KDR ) :: &
    MinWidth, &
    ZoomFactor
  real ( KDR ), dimension ( nEnergyCells ) :: &
    Bottom, &
    Top, &
    Width, &
    Center
  type ( Real_1D_Form ), dimension ( 1 ) :: &
    Edge
  class ( GeometryFlatForm ), pointer :: &
    G

  call Show ( 'Creating NuLib energy grid' )

  call MS % CreateChart &
         ( CoordinateLabelOption  = CoordinateLabel, &
           CoordinateSystemOption = 'SPHERICAL', &
           CoordinateUnitOption   = CoordinateUnit, &
           nCellsOption           = [ nEnergyCells ], &
           nGhostLayersOption     = [ 0, 0, 0 ] )
  call MS % SetGeometry ( )

  MinWidth                 =  2.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  Bottom ( 1 )             =  0.0_KDR
  Bottom ( 2 )             =  MinWidth
  Bottom ( 3 )             =  Bottom ( 2 )  +  MinWidth
  Bottom ( nEnergyCells )  =  250.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT

  !-- Compute zoom factor

  call nulib_series2 &
         ( nEnergyCells - 1, Bottom ( 2 ), Bottom ( nEnergyCells ), MinWidth, &
           ZoomFactor )

  !-- Compute bin bottom, top, width, center

  do iE = 4, nEnergyCells
    Bottom ( iE )  =  Bottom ( iE - 1 )  +  &
                      ( Bottom ( iE - 1 ) - Bottom ( iE - 2 ) ) * ZoomFactor
  end do

  do iE = 1, nEnergyCells - 1
    Top    ( iE )        =  Bottom ( iE + 1 )
    Width  ( iE )        =  Top ( iE )  -  Bottom ( iE )
    Center ( iE )        =  ( Bottom ( iE ) + Top ( iE ) )  /  2.0_KDR
  end do
  Width ( nEnergyCells )  &
    =  Width ( nEnergyCells - 1 )  *  ZoomFactor
  Top ( nEnergyCells )  &
    =  Bottom ( nEnergyCells )  +  Width ( nEnergyCells )
  Center ( nEnergyCells )  &
    =  ( Bottom ( nEnergyCells )  +  Top ( nEnergyCells ) )  /  2.0_KDR

  !-- Hack geometry

  call Edge ( 1 ) % Initialize ( nEnergyCells + 1 )
  Edge ( 1 ) % Value  =  [ Bottom, Top ( nEnergyCells ) ]
  select type ( MSC => MS % Chart )
  class is ( Chart_SL_Template )
    call MSC % ResetGeometry ( Edge )
  end select !-- MSC

  G => MS % Geometry ( )

  G % Value ( :, G % CENTER_U ( 1 ) )       =  Center
  G % Value ( :, G % WIDTH_LEFT_U ( 1 ) )   =  Center - Bottom
  G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) )  =  Top - Center

  nullify ( G )

end subroutine SetEnergyGridNuLib


subroutine SetEnergyGridGeometric &
             ( MS, CoordinateLabel, CoordinateUnit, nEnergyCells )

  type ( Atlas_SC_Form ), intent ( inout ) :: &
    MS
  character ( LDL ), dimension ( 1 ), intent ( in ) :: &
    CoordinateLabel
  type ( MeasuredValueForm ), dimension ( 1 ), intent ( in ) :: &
    CoordinateUnit
  integer ( KDI ), intent ( in ) :: &
    nEnergyCells

  real ( KDR ) :: &
    MinWidth, &
    MaxEnergy
  real ( KDR ), dimension ( 1 ) :: &
    MaxCoordinate, &
    Scale
  character ( LDL ), dimension ( 1 ) :: &
    Spacing

  call Show ( 'Creating GenASiS GEOMETRIC energy grid' )

  MinWidth   =  2.0_KDR    *  UNIT % MEGA_ELECTRON_VOLT
  MaxEnergy  =  310.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  call PROGRAM_HEADER % GetParameter ( MinWidth, 'MinWidth' )
  call PROGRAM_HEADER % GetParameter ( MaxEnergy, 'MaxEnergy' )

  MaxCoordinate  =  MaxEnergy
  Scale          =  MinWidth
  Spacing        =  'GEOMETRIC'
  
  call MS % CreateChart &
         ( SpacingOption          = Spacing, &
           CoordinateLabelOption  = CoordinateLabel, &
           CoordinateSystemOption = 'SPHERICAL', &
           CoordinateUnitOption   = CoordinateUnit, &
           MaxCoordinateOption    = MaxCoordinate, &
           ScaleOption            = Scale, &
           nCellsOption           = [ nEnergyCells ], &
           nGhostLayersOption     = [ 0, 0, 0 ] )
  call MS % SetGeometry ( )

end subroutine SetEnergyGridGeometric


subroutine SetEnergyGridCompactified &
             ( MS, CoordinateLabel, CoordinateUnit, nEnergyCells )

  type ( Atlas_SC_Form ), intent ( inout ) :: &
    MS
  character ( LDL ), dimension ( 1 ), intent ( in ) :: &
    CoordinateLabel
  type ( MeasuredValueForm ), dimension ( 1 ), intent ( in ) :: &
    CoordinateUnit
  integer ( KDI ), intent ( in ) :: &
    nEnergyCells

  real ( KDR ) :: &
    EnergyScale
  real ( KDR ), dimension ( 1 ) :: &
    Scale
  character ( LDL ), dimension ( 1 ) :: &
    Spacing

  call Show ( 'Creating GenASiS COMPACTIFIED energy grid' )

  EnergyScale  =  10.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT
  call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

  Scale    =  EnergyScale
  Spacing  =  'COMPACTIFIED'
  
  call MS % CreateChart &
         ( SpacingOption          = Spacing, &
           CoordinateLabelOption  = CoordinateLabel, &
           CoordinateSystemOption = 'SPHERICAL', &
           CoordinateUnitOption   = CoordinateUnit, &
           ScaleOption            = Scale, &
           nCellsOption           = [ nEnergyCells ], &
           nGhostLayersOption     = [ 0, 0, 0 ] )
  call MS % SetGeometry ( )

end subroutine SetEnergyGridCompactified


subroutine PrepareNuLib ( MS )

  use nulib, only: &
    energies, &
    bin_widths, &
    bin_bottom, &
    bin_top, &
    initialize_nulib

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
    neutrino_scheme, &
    number_groups, &
    energies, &
    bin_widths, &
    bin_bottom, &
    bin_top, &
    total_eos_variables, &
    rhoindex, &
    tempindex, &
    yeindex, &
    single_point_return_all, &
    single_ipoint_return_all, &
    single_epannihil_kernel_point_return_all, &
    add_nue_Iscattering_electrons, &
    add_anue_Iscattering_electrons, &
    add_numu_Iscattering_electrons, &
    add_anumu_Iscattering_electrons, &
    add_nutau_Iscattering_electrons, &
    add_anutau_Iscattering_electrons, &
    add_nue_kernel_epannihil, &
    add_anue_kernel_epannihil, &
    add_numu_kernel_epannihil, &
    add_anumu_kernel_epannihil, &
    add_nutau_kernel_epannihil, &
    add_anutau_kernel_epannihil
  integer ( KDI ) :: &
    iD, iT, iY, iX, iS, iG, iGp
  real( KDR ) :: &
    dLogD, dLogT, dY, dLogX
  real ( KDR ), dimension( :, : ), allocatable :: &
    emission, &
    absopacity, &
    scatopacity, &
    nes_phi_0, &
    nes_phi_1
  real ( KDR ), dimension( :, :, : ), allocatable :: &
    pair_phi_0, &
    pair_phi_1

  select case ( neutrino_scheme )

    case ( 1 )

      number_local_species = 3

    case ( 2 )

      number_local_species = 4

    case ( 3 )

      number_local_species = 6

    case default

      call Show ( 'incorrect neutrino scheme' )

  end select

  include_NES &
  ! --- Set in requested_interactions.inc ---
    = any( [ add_nue_Iscattering_electrons,   &
             add_anue_Iscattering_electrons,  &
             add_numu_Iscattering_electrons,  &
             add_anumu_Iscattering_electrons, &
             add_nutau_Iscattering_electrons, &
             add_anutau_Iscattering_electrons ] )

  include_Pair &
    = any( [ add_nue_kernel_epannihil,   &
             add_anue_kernel_epannihil,  &
             add_numu_kernel_epannihil,  &
             add_anumu_kernel_epannihil, &
             add_nutau_kernel_epannihil, &
             add_anutau_kernel_epannihil ] )

  call Show ( 'Creating NuLib table' )
  call Show ( neutrino_scheme, 'neutrino_scheme' )
  call Show ( number_local_species, 'number_local_species' )
  call Show ( number_groups, 'number_groups' )
  call Show ( energies, 'energies' )
  call Show ( bin_widths, 'bin_widths' )
  call Show ( bin_bottom, 'bin_bottom' )
  call Show ( bin_top, 'bin_top' )
  call Show ( include_NES, 'include_NES' )
  call Show ( include_Pair, 'include_Pair' )

  call read_eos_table( eos_table_filename )

  ! --- Density Grid ---

  allocate( table_D(table_size_D) )

  dLogD = ( max_LogD - min_LogD ) / dble( table_size_D - 1 )
  do iD = 1, table_size_D

    table_D(iD) = 10.0d0 ** ( min_LogD + dble( iD - 1 ) * dLogD )

  end do

  call Show ( table_D, 'Table Densities' )

  ! --- Temperature Grid ---

  allocate( table_T(table_size_T) )

  dLogT = ( max_LogT - min_LogT ) / dble( table_size_T - 1 )
  do iT = 1, table_size_T

    table_T(iT) = 10.0d0 ** ( min_LogT + dble( iT - 1 ) * dLogT )

  end do

  call Show ( table_T, 'Table Temperatures' )

  ! --- Electron Fraction Grid ---

  allocate( table_Y(table_size_Y) )

  dY = ( max_Y - min_Y ) / dble( table_size_Y - 1 )
  do iY = 1, table_size_Y

    table_Y(iY) = min_Y + dble( iY - 1 ) * dY

  end do

  call Show ( table_Y, 'Table Electron Fractions' )

  ! --- Electron Degeneracy Parameter Grid ---
  
  allocate( table_X(table_size_X) )

  dLogX = ( max_LogX - min_LogX ) / dble( table_size_X - 1 )
  do iX = 1, table_size_X

    table_X(iX) = 10.0d0 ** ( min_LogX + dble( iX - 1 ) * dLogX )

  enddo

  call Show ( table_X, 'Table Electron Degen. Parameters' )

  allocate &
    ( table_emission &
        (table_size_D,table_size_T,table_size_Y,number_local_species, &
         number_groups) )
  allocate &
    ( table_absopacity &
        (table_size_D,table_size_T,table_size_Y,number_local_species, &
         number_groups) )
  allocate &
    ( table_scatopacity &
        (table_size_D,table_size_T,table_size_Y,number_local_species, &
         number_groups) )

  allocate( eos_variables(total_eos_variables) )

  allocate( emission   (number_local_species,number_groups) )
  allocate( absopacity (number_local_species,number_groups) )
  allocate( scatopacity(number_local_species,number_groups) )

  do iD = 1, table_size_D
  do iT = 1, table_size_T
  do iY = 1, table_size_Y

    call Show( [ iD, iT, iY ], '[ iD, iT, iY ]' )

    eos_variables = 0.0_KDR

    eos_variables(rhoindex)  = table_D(iD)
    eos_variables(tempindex) = table_T(iT)
    eos_variables(yeindex)   = table_Y(iY)

    call Show( [ table_D(iD), table_T(iT), table_Y(iY) ], '[ D, T, Y ]' )

    call set_eos_variables( eos_variables )

    call single_point_return_all &
           ( eos_variables, emission, absopacity, scatopacity, neutrino_scheme )

    do iS = 1, number_local_species
    do iG = 1, number_groups

      table_emission   (iD,iT,iY,iS,iG) = emission   (iS,iG)
      table_absopacity (iD,iT,iY,iS,iG) = absopacity (iS,iG)
      table_scatopacity(iD,iT,iY,iS,iG) = scatopacity(iS,iG)

    end do
    end do

  end do
  end do
  end do

  deallocate( emission, absopacity, scatopacity )
  deallocate( eos_variables )

  if( include_NES )then

    call Show ( 'Computing NES Kernels' )

    allocate &
      ( table_nes_phi_0 &
          (table_size_T,table_size_X,number_groups,number_local_species, &
           number_groups) )

    allocate &
      ( table_nes_phi_1 &
          (table_size_T,table_size_X,number_groups,number_local_species, &
           number_groups) )

    allocate( nes_phi_0(number_local_species,number_groups) )
    allocate( nes_phi_1(number_local_species,number_groups) )

    do iT = 1, table_size_T
    do iX = 1, table_size_X

      call Show( [ iT, iX ], '[ iT, iX ]' )

      do iGp = number_groups, 1, - 1

        call single_Ipoint_return_all &
               ( iGp, table_X(iX), table_T(iT), nes_phi_0, nes_phi_1, &
                 neutrino_scheme )

        do iS = 1, number_local_species
        do iG = iGp+1, number_groups

          nes_phi_0(iS,iG) &
            = exp( - ( energies(iG) - energies(iGp) ) / table_T(iT) ) &
                * table_nes_phi_0(iT,iX,iG,iS,iGp)

          nes_phi_1(iS,iG) &
            = exp( - ( energies(iG) - energies(iGp) ) / table_T(iT) ) &
                * table_nes_phi_1(iT,iX,iG,iS,iGp)

        end do
        end do

        do iS = 1, number_local_species
        do iG = 1, number_groups

          table_nes_phi_0(iT,iX,iGp,iS,iG) = nes_phi_0(iS,iG)
          table_nes_phi_1(iT,iX,iGp,iS,iG) = nes_phi_1(iS,iG)

        end do
        end do

      end do

    end do
    end do

    deallocate( nes_phi_0, nes_phi_1 )

  end if

  if( include_Pair )then

    call Show ( 'Computing Pair Kernels' )

    allocate &
      ( table_pair_phi_0 &
          (table_size_T,table_size_X,number_groups,number_local_species, &
           number_groups,2) )

    allocate &
      ( table_pair_phi_1 &
          (table_size_T,table_size_X,number_groups,number_local_species, &
           number_groups,2) )

    allocate( pair_phi_0(number_local_species,number_groups,2) )
    allocate( pair_phi_1(number_local_species,number_groups,2) )

    do iT = 1, table_size_T
    do iX = 1, table_size_X

      call Show( [ iT, iX ], '[ iT, iX ]' )

      do iGp = number_groups, 1, - 1

        call single_epannihil_kernel_point_return_all &
               ( iGp, table_X(iX), table_T(iT), pair_phi_0, pair_phi_1, &
                 neutrino_scheme )

        do iS = 1, number_local_species
        do iG = 1, number_groups

          table_pair_phi_0(iT,iX,iGp,iS,iG,1) = pair_phi_0(iS,iG,1)
          table_pair_phi_0(iT,iX,iGp,iS,iG,2) = pair_phi_0(iS,iG,2)
          table_pair_phi_1(iT,iX,iGp,iS,iG,1) = pair_phi_1(iS,iG,1)
          table_pair_phi_1(iT,iX,iGp,iS,iG,2) = pair_phi_1(iS,iG,2)

        end do
        end do

      end do

    end do
    end do

    deallocate( pair_phi_0, pair_phi_1 )

  end if

  call WriteTable

  deallocate( table_D )
  deallocate( table_T )
  deallocate( table_Y )

  deallocate( table_emission )
  deallocate( table_absopacity )
  deallocate( table_scatopacity )

  if( include_NES )then

    deallocate( table_nes_phi_0, table_nes_phi_1 )

  end if

  if( include_Pair )then

    deallocate( table_pair_phi_0, table_pair_phi_1 )

  end if

end subroutine CreateTable


subroutine WriteTable

  use nulib, only: &
    number_groups, &
    bin_bottom, &
    bin_top, &
    bin_widths, &
    energies
    
  use hdf5

  character ( 8 )     :: date
  character ( LDL )   :: string_D
  character ( LDL )   :: string_T
  character ( LDL )   :: string_Y
  character ( LDL )   :: string_X
  character ( LDL )   :: string_ng
  character ( LDL )   :: string_ns
  character ( LDF )   :: table_filename
  integer ( KDI )     :: error, cerror, rank, values(8)
  integer ( HID_T )   :: file_id, dset_id, dspace_id
  integer ( HSIZE_T ) :: dims1(1), dims5(5), dims6(6)
  real ( KDR )        :: timestamp

  call date_and_time( DATE = date, VALUES = values )

  write( string_D,  * ) table_size_D
  write( string_T,  * ) table_size_T
  write( string_Y,  * ) table_size_Y
  write( string_ng, * ) number_groups
  write( string_ns, * ) number_local_species
  write( string_X,  * ) table_size_X

  timestamp = dble( values(1) ) * 10000.0d0 &
                + dble( values(2) ) * 100.0 &
                + dble( values(3) ) &
                + ( dble( values(5) ) + dble( values(6) ) / 60.0d0 &
                      + dble( values(7) ) / 3600.0d0 ) / 24.0

  if( include_NES .or. include_Pair )then

    table_filename &
      = trim( adjustl( outdir ) ) &
        // trim( adjustl( base ) ) &
        // "_rho"   // trim( adjustl( string_D ) ) &
        // "_temp"  // trim( adjustl( string_T ) ) &
        // "_ye"    // trim( adjustl( string_Y ) ) &
        // "_ng"    // trim( adjustl( string_ng ) ) &
        // "_ns"    // trim( adjustl( string_ns ) ) &
        // "_Itemp" // trim( adjustl( string_T ) ) &
        // "_Ieta"  // trim( adjustl( string_X ) ) &
        // "_version" // trim( adjustl( vnum ) ) &
        // "_" // trim( adjustl( date ) ) // ".h5"

  else

    table_filename &
      = trim( adjustl( outdir ) ) &
        // trim( adjustl( base ) ) &
        // "_rho"  // trim( adjustl( string_D ) ) &
        // "_temp" // trim( adjustl( string_T ) ) &
        // "_ye"   // trim( adjustl( string_Y ) ) &
        // "_ng"   // trim( adjustl( string_ng ) ) &
        // "_ns"   // trim( adjustl( string_ns ) ) &
        // "_version" // trim( adjustl( vnum ) ) &
        // "_" // trim( adjustl( date ) ) // ".h5"

  end if

  call Show ( 'Writing Table' )
  call Show ( table_filename, 'table_filename' )
  call Show ( timestamp, 'timestamp' )
  call Show ( number_local_species, 'number_local_species' )

  cerror = 0

  call h5open_f( error )
  cerror = cerror + error 

  call h5fcreate_f( table_filename, H5F_ACC_TRUNC_F, file_id, error )
  cerror = cerror + error

  ! --- write scalars (rank = 1, dims1(1) = 1) ---
  rank = 1; dims1 = 1

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "number_species", H5T_NATIVE_INTEGER, dspace_id, &
           dset_id, error)
  call h5dwrite_f &
         ( dset_id, H5T_NATIVE_INTEGER, number_local_species, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "timestamp", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, timestamp, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "number_groups", H5T_NATIVE_INTEGER, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, number_groups, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "nrho", H5T_NATIVE_INTEGER, dspace_id,dset_id, error )
  call h5dwrite_f &
         ( dset_id, H5T_NATIVE_INTEGER, table_size_D, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "ntemp", H5T_NATIVE_INTEGER, dspace_id,dset_id, error )
  call h5dwrite_f &
         ( dset_id, H5T_NATIVE_INTEGER, table_size_T, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "nye", H5T_NATIVE_INTEGER, dspace_id,dset_id, error )
  call h5dwrite_f &
         ( dset_id, H5T_NATIVE_INTEGER, table_size_Y, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  ! --- write energy grid ---
  rank = 1; dims1 = number_groups

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "neutrino_energies", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, energies, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "bin_widths", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, bin_widths, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error   
    
  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "bin_bottom", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE,bin_bottom, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "bin_top", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE,bin_top, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  ! --- write density grid ---
  rank = 1; dims1 = table_size_D

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "rho_points", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_D, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error     

  ! --- write temperature grid ---
  rank = 1; dims1 = table_size_T

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "temp_points", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_T, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  ! --- write electron fraction grid ---
  rank = 1; dims1 = table_size_Y

  call h5screate_simple_f( rank, dims1, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "ye_points", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_Y, dims1, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  ! --- write emissivities, absorption scattering opacities ---
  rank = 5
  dims5(1) = table_size_D
  dims5(2) = table_size_T
  dims5(3) = table_size_Y
  dims5(4) = number_local_species
  dims5(5) = number_groups
    
  call h5screate_simple_f( rank, dims5, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "emissivities", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_emission, dims5, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error   

  call h5screate_simple_f( rank, dims5, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "absorption_opacity", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_absopacity, dims5, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error   

  call h5screate_simple_f( rank, dims5, dspace_id, error )
  call h5dcreate_f &
         ( file_id, "scattering_opacity", H5T_NATIVE_DOUBLE, dspace_id, &
           dset_id, error )
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE,table_scatopacity, dims5, error )
  call h5dclose_f( dset_id, error )
  call h5sclose_f( dspace_id, error )
  cerror = cerror + error

  if( include_NES .or. include_Pair )then

    rank = 1; dims1 = 1

    call h5screate_simple_f( rank, dims1, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "Itemp", H5T_NATIVE_INTEGER, dspace_id, dset_id, error )
    call h5dwrite_f &
           ( dset_id, H5T_NATIVE_INTEGER, table_size_T, dims1, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error

    rank = 1; dims1 = 1
    call h5screate_simple_f( rank, dims1, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "Ieta", H5T_NATIVE_INTEGER, dspace_id, dset_id, error )
    call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, table_size_X, dims1, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error

    rank = 1; dims1 = table_size_T
    call h5screate_simple_f( rank, dims1, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "temp_Ipoints", H5T_NATIVE_DOUBLE, dspace_id, &
             dset_id, error )
    call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_T, dims1, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error

    rank = 1; dims1 = table_size_X
    call h5screate_simple_f( rank, dims1, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "eta_Ipoints", H5T_NATIVE_DOUBLE, dspace_id, &
             dset_id, error )
    call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, table_X, dims1, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error   

  end if

  if( include_NES )then

    rank = 5
    dims5(1) = table_size_T
    dims5(2) = table_size_X
    dims5(3) = number_groups
    dims5(4) = number_local_species  
    dims5(5) = number_groups  

    call h5screate_simple_f( rank, dims5, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "inelastic_phi0", H5T_NATIVE_DOUBLE, dspace_id, &
             dset_id, error )
    call h5dwrite_f &
           ( dset_id, H5T_NATIVE_DOUBLE, table_nes_phi_0, dims5, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error   

    call h5screate_simple_f( rank, dims5, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "inelastic_phi1", H5T_NATIVE_DOUBLE, dspace_id, &
             dset_id, error )
    call h5dwrite_f &
           ( dset_id, H5T_NATIVE_DOUBLE, table_nes_phi_1, dims5, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error

  end if

  if( include_Pair )then

    rank = 6
    dims6(1) = table_size_T
    dims6(2) = table_size_X
    dims6(3) = number_groups
    dims6(4) = number_local_species  
    dims6(5) = number_groups
    dims6(6) = 2
       
    call h5screate_simple_f( rank, dims6, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "epannihil_phi0", H5T_NATIVE_DOUBLE, dspace_id, &
             dset_id, error )
    call h5dwrite_f &
           ( dset_id, H5T_NATIVE_DOUBLE, table_pair_phi_0, dims6, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error   

    call h5screate_simple_f( rank, dims6, dspace_id, error )
    call h5dcreate_f &
           ( file_id, "epannihil_phi1", H5T_NATIVE_DOUBLE, dspace_id, &
             dset_id, error )
    call h5dwrite_f &
           ( dset_id, H5T_NATIVE_DOUBLE, table_pair_phi_1, dims6, error )
    call h5dclose_f( dset_id, error )
    call h5sclose_f( dspace_id, error )
    cerror = cerror + error

  end if

  ! --- close h5 files, check for error ---

  if( cerror .ne. 0 )then
    call Show ( 'Errors Writing Table' )
    call Show ( cerror, 'cerror' )
    stop
  else
    call Show ( 'Finished Writing Table' )
  endif

  call h5fclose_f( file_id, error )

  call h5close_f( error )

end subroutine WriteTable


end program TableGeneration_NuLib
