module Tally_F_P_MHN__Form

  !-- Tally_Perfect_MeanHeavyNucleus_Form

  use Basics
  use Mathematics
  use Fluid_P_MHN__Form
  use Tally_F_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_MEAN_HEAVY_NUCLEUS = 1
  
  type, public, extends ( Tally_F_P_Form ) :: Tally_F_P_MHN_Form
    integer ( KDI ) :: &
      N_INTEGRALS_MEAN_HEAVY_NUCLEUS = N_INTEGRALS_MEAN_HEAVY_NUCLEUS, &
      PROTON_NUMBER = 0
    real ( KDR ), dimension ( : ), pointer :: &
      GravitationalPotential => null ( )
    procedure ( ), pointer, nopass :: &
      ComputeGravitationalPotential => null ( )
  contains
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      SetGravitationalPotential
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrandGalilean
    procedure, public, pass :: &
      ComputeBoundaryIntegrandGalilean_CSL
  end type Tally_F_P_MHN_Form


contains


  subroutine InitializeFluid &
               ( T, A, MassUnitOption, EnergyUnitOption, MomentumUnitOption, &
                 AngularMomentumUnitOption )
    
    class ( Tally_F_P_MHN_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      MassUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

    integer ( KDI ) :: &
      oI     !-- oIntegral

    oI  =  T % N_INTEGRALS_DUST  +  T % N_INTEGRALS_PERFECT
    if ( T % N_INTEGRALS == 0 ) &
      T % N_INTEGRALS = oI + T % N_INTEGRALS_MEAN_HEAVY_NUCLEUS
    
    call T % Tally_F_P_Form % Initialize &
           ( A, MassUnitOption, EnergyUnitOption, MomentumUnitOption, &
             AngularMomentumUnitOption )
    
    T % PROTON_NUMBER  =  oI + 1
    
    T % Variable ( oI + 1 : oI + T % N_INTEGRALS_MEAN_HEAVY_NUCLEUS ) &
      = [ 'ProtonNumber' ]

    if ( present ( MassUnitOption ) ) then
      T % Unit ( oI + 1 : oI + T % N_INTEGRALS_MEAN_HEAVY_NUCLEUS ) &
        = [ MassUnitOption ]
    else
      T % Unit ( oI + 1 : oI + T % N_INTEGRALS_MEAN_HEAVY_NUCLEUS ) &
        = [ UNIT % IDENTITY ]
    end if

    call T % SelectVariables ( A )

    T % Atlas => A

  end subroutine InitializeFluid

  
  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_F_P_MHN_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    class ( GeometryFlatForm ), pointer :: &
      G

    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )

    select type ( A )
    class is ( Atlas_SC_Form )

    G => A % Geometry ( )
    
    select type ( G )
    type is ( GeometryFlatForm )
      T % nSelected = 13
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, T % PROTON_NUMBER, T % MOMENTUM, &
            T % FLUID_ENERGY, T % INTERNAL_ENERGY, T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, T % GRAVITATIONAL_ENERGY, T % TOTAL_ENERGY ]
    
    ! type is ( GeometryNewtonianForm )
    !   allocate ( T % iaSelected ( 8 ) )
    !   T % iaSelected &
    !     = [ T % BARYON_NUMBER, T % MOMENTUM, T % KINETIC_ENERGY, &
    !         T % FLUID_ENERGY, T % GRAVITATIONAL_ENERGY, &
    !         T % TOTAL_ENERGY ]
    
    class default 
      call Show ( 'This type is not implemented yet', CONSOLE % WARNING )
      call Show ( 'Tally_F_P_MHN__Form % SelectVariables', &
                  'Subroutine', CONSOLE % WARNING )
            
    end select !-- G

    class default 
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_P_MHN__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    nullify ( G )

  end subroutine SelectVariables
  

  subroutine SetGravitationalPotential ( T, Phi )

    class ( Tally_F_P_MHN_Form ), intent ( inout ) :: &
      T
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Phi

    T % GravitationalPotential  =>  Phi

  end subroutine SetGravitationalPotential


  impure elemental subroutine Finalize ( T )
  
    type ( Tally_F_P_MHN_Form ), intent ( inout ) :: &
      T

    !-- Trigger finalization of parent

  end subroutine Finalize
  

  subroutine ComputeInteriorIntegrandGalilean &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_F_P_MHN_Form ), intent ( inout ) :: &
      T
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    type ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI     !-- iIntegral

    select type ( C )
    class is ( Fluid_P_MHN_Form )

    call T % Tally_F_P_Form &
           % ComputeInteriorIntegrandGalilean ( Integrand, C, G, nDimensions )

    call T % ComputeGravitationalPotential ( C, T % Atlas, G )

    associate &
      ( D   => C % Value ( :, C % CONSERVED_DENSITY ), &
        DP  => C % Value ( :, C % CONSERVED_PROTON_DENSITY ), &
        CE  => C % Value ( :, C % CONSERVED_ENERGY ), &
        Phi => T % GravitationalPotential )
    
    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % PROTON_NUMBER ) then
        call Copy ( DP, Integrand ( iS ) % Value )
      else if ( iI == T % GRAVITATIONAL_ENERGY ) then
        Integrand ( iS ) % Value  =  0.5_KDR * D * Phi
      else if ( iI == T % TOTAL_ENERGY ) then
        Integrand ( iS ) % Value  =  CE  +  0.5_KDR * D * Phi
      end if !-- iI
    end do !-- iS

    end associate !-- DP
    end select !-- F

  end subroutine ComputeInteriorIntegrandGalilean


  subroutine ComputeBoundaryIntegrandGalilean_CSL &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_F_P_MHN_Form ), intent ( inout ) :: &
      T
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: & 
      C
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    type ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence

    integer ( KDI ) :: &
      iF, &   !-- iFace
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
      iC, &   !-- iConnectivity
      iD, &   !-- iDimension
      iFluence, &
      iProton, &
      iDensity, &
      iEnergy

    select type ( C )
    class is ( Fluid_P_MHN_Form )

    call T % Tally_F_P_Form % ComputeBoundaryIntegrandGalilean_CSL &
           ( Integrand, C, CSL, G, BoundaryFluence )

    do iFluence = 1, C % N_CONSERVED
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_PROTON_DENSITY ) &
        iProton = iFluence
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_DENSITY ) &
        iDensity = iFluence
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_ENERGY ) &
        iEnergy = iFluence
    end do !-- iFluence

    associate ( Cnnct => CSL % Atlas % Connectivity )

    do iF = 1, Cnnct % nFaces
      associate ( CE => BoundaryFluence ( iProton, iF ) % Value )
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % PROTON_NUMBER ) then
          call CopyCollapse ( CE, Integrand ( iS, iF ) % Value )
        end if !-- iI
      end do !-- iS
      end associate !-- CE
    end do !-- iF

    !-- outer radial boundary, assume 1D
    iD  =  1
    iC  =  Cnnct % iaOuter ( iD )
    associate &
      ( D   => BoundaryFluence ( iDensity, iC ) % Value ( 1, 1, 1 ), &
        CE  => BoundaryFluence ( iEnergy, iC ) % Value ( 1, 1, 1 ), &
        Phi => T % GravitationalPotential ( C % nValues - 1 ) )
    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % GRAVITATIONAL_ENERGY ) then
        Integrand ( iS, iC ) % Value ( 1, 1, 1 ) =  D * Phi
      else if ( iI == T % TOTAL_ENERGY ) then
        Integrand ( iS, iC ) % Value ( 1, 1, 1 ) =  CE  +  D * Phi 
      end if !-- iI
    end do !-- iS

    end associate !-- D, etc.

    end associate !-- Cnnct

    end select !-- C

  end subroutine ComputeBoundaryIntegrandGalilean_CSL


end module Tally_F_P_MHN__Form
