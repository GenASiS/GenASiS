module Tally_F_P_HN__Form

  !-- Tally_Perfect_MeanHeavyNucleus_Form

  use Basics
  use Mathematics
  use Spaces
  use Fluid_P_HN__Form
  use Tally_F_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_HEAVY_NUCLEUS = 1
  
  type, public, extends ( Tally_F_P_Form ) :: Tally_F_P_HN_Form
    integer ( KDI ) :: &
      N_INTEGRALS_HEAVY_NUCLEUS = N_INTEGRALS_HEAVY_NUCLEUS, &
      PROTON_NUMBER = 0
  contains
    procedure, public, pass :: &
      InitializeFluid
    procedure, public, pass :: &
      SelectVariables
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_CSL
  end type Tally_F_P_HN_Form


contains


  subroutine InitializeFluid &
               ( T, A, NumberUnitOption, EnergyUnitOption, &
                 MomentumUnitOption, AngularMomentumUnitOption )
    
    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      NumberUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

    integer ( KDI ) :: &
      oI     !-- oIntegral

    oI  =  T % N_INTEGRALS_DUST  +  T % N_INTEGRALS_PERFECT
    if ( T % N_INTEGRALS == 0 ) &
      T % N_INTEGRALS = oI + T % N_INTEGRALS_HEAVY_NUCLEUS
    
    call T % Tally_F_P_Form % Initialize &
           ( A, NumberUnitOption, EnergyUnitOption, MomentumUnitOption, &
             AngularMomentumUnitOption )
    
    T % PROTON_NUMBER  =  oI + 1
    
    T % Variable ( oI + 1 : oI + T % N_INTEGRALS_HEAVY_NUCLEUS ) &
      = [ 'ProtonNumber' ]

    if ( present ( NumberUnitOption ) ) then
      T % Unit ( oI + 1 : oI + T % N_INTEGRALS_HEAVY_NUCLEUS ) &
        = [ NumberUnitOption ]
    else
      T % Unit ( oI + 1 : oI + T % N_INTEGRALS_HEAVY_NUCLEUS ) &
        = [ UNIT % IDENTITY ]
    end if

    call T % SelectVariables ( A )

    T % Atlas => A

  end subroutine InitializeFluid

  
  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
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
    type is ( Geometry_G_Form )
      T % nSelected = 11
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, T % PROTON_NUMBER, T % MOMENTUM, &
            T % FLUID_ENERGY, T % INTERNAL_ENERGY, T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM ]
    
    type is ( Geometry_N_Form )
      T % nSelected = 13
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, T % PROTON_NUMBER, T % MOMENTUM, &
            T % FLUID_ENERGY, T % INTERNAL_ENERGY, T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, T % GRAVITATIONAL_ENERGY, &
            T % TOTAL_ENERGY ]
    
    class default 
      call Show ( 'This type is not implemented yet', CONSOLE % WARNING )
      call Show ( 'Tally_F_P_HN__Form % SelectVariables', &
                  'Subroutine', CONSOLE % WARNING )
            
    end select !-- G

    class default 
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_P_HN__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    nullify ( G )

  end subroutine SelectVariables
  

  impure elemental subroutine Finalize ( T )
  
    type ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T

    !-- Trigger finalization of parent

  end subroutine Finalize
  

  subroutine ComputeInteriorIntegrand ( T, Integrand, C, G, nDimensions )

    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI     !-- iIntegral

    select type ( C )
    class is ( Fluid_P_HN_Form )

    call T % Tally_F_P_Form &
           % ComputeInteriorIntegrand ( Integrand, C, G, nDimensions )

    associate ( DP => C % Value ( :, C % CONSERVED_PROTON_DENSITY ) )
    
    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % PROTON_NUMBER ) then
        call Copy ( DP, Integrand ( iS ) % Value )
      end if !-- iI
    end do !-- iS

    end associate !-- DP
    end select !-- C

  end subroutine ComputeInteriorIntegrand


  subroutine ComputeBoundaryIntegrand_CSL &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_F_P_HN_Form ), intent ( inout ) :: &
      T
    type ( Real_3D_Form ), dimension ( :, : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: & 
      C
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BoundaryFluence

    integer ( KDI ) :: &
      iF, &   !-- iFace
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
!       iC, &   !-- iConnectivity
!       iD, &   !-- iDimension
      iFluence, &
      iProton

    select type ( C )
    class is ( Fluid_P_HN_Form )

    call T % Tally_F_P_Form % ComputeBoundaryIntegrand_CSL &
           ( Integrand, C, CSL, G, BoundaryFluence )

    do iFluence = 1, C % N_CONSERVED
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_PROTON_DENSITY ) &
        iProton = iFluence
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

!     !-- outer radial boundary, assume 1D
!     iD  =  1
!     iC  =  Cnnct % iaOuter ( iD )
!     associate &
!       ( D   => BoundaryFluence ( iDensity, iC ) % Value ( 1, 1, 1 ), &
!         CE  => BoundaryFluence ( iEnergy, iC ) % Value ( 1, 1, 1 ), &
!         Phi => 0.5_KDR * ( T % GravitationalPotential ( C % nValues - 1 ) &
!                            + T % GravitationalPotential ( C % nValues - 2 ) ) )
!     do iS = 1, T % nSelected
!       iI = T % iaSelected ( iS )
!       if ( iI == T % GRAVITATIONAL_ENERGY ) then
!         Integrand ( iS, iC ) % Value ( 1, 1, 1 ) =  D * Phi
!       else if ( iI == T % TOTAL_ENERGY ) then
!         Integrand ( iS, iC ) % Value ( 1, 1, 1 ) =  CE  +  D * Phi 
!       end if !-- iI
!     end do !-- iS

!     end associate !-- D, etc.

    end associate !-- Cnnct

    end select !-- C

  end subroutine ComputeBoundaryIntegrand_CSL


end module Tally_F_P_HN__Form
