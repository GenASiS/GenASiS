module Tally_FM__Form

  !-- Tally_FishboneMoncrief_Form

  use Basics
  use Mathematics
  use Fluid_P_P__Form
  use Tally_F_P__Form
  
  implicit none
  private

  type, public, extends ( Tally_F_P_Form ) :: Tally_FM_Form
    real ( KDR ) :: &
      CentralMass
  contains
    procedure, public, pass :: &
      SetCentralMass
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      ComputeInteriorIntegrandGalilean
    procedure, public, pass :: &
      ComputeBoundaryIntegrandGalilean_CSL
  end type Tally_FM_Form

contains


  subroutine SetCentralMass ( T, CentralMass )

    class ( Tally_FM_Form ), intent ( inout ) :: &
      T
    real ( KDR ), intent ( in ) :: &
      CentralMass

    T % CentralMass = CentralMass
    
  end subroutine SetCentralMass


  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_FM_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    if ( allocated ( T % iaSelected ) ) &
      deallocate ( T % iaSelected )

    T % nSelected = 12

    allocate ( T % iaSelected ( T % nSelected ) )

    T % iaSelected &
      = [ T % BARYON_NUMBER, T % MOMENTUM, T % FLUID_ENERGY, &
          T % INTERNAL_ENERGY, T % KINETIC_ENERGY, T % ANGULAR_MOMENTUM, &
          T % GRAVITATIONAL_ENERGY, T % TOTAL_ENERGY ]
    
  end subroutine SelectVariables
  

  subroutine ComputeInteriorIntegrandGalilean &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_FM_Form ), intent ( inout ) :: &
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
    class is ( Fluid_P_P_Form )

    call T % Tally_F_P_Form &
           % ComputeInteriorIntegrandGalilean ( Integrand, C, G, nDimensions )

    associate &
      ( CE => C % Value ( :, C % CONSERVED_ENERGY ), &
        N  => C % Value ( :, C % COMOVING_DENSITY ), &
        R  => G % Value ( :, G % CENTER ( 1 ) ), &
        M  => T % CentralMass, &
        GC => CONSTANT % GRAVITATIONAL )

    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % GRAVITATIONAL_ENERGY ) then
        Integrand ( iS ) % Value  =  - GC * M * N / R
      else if ( iI == T % TOTAL_ENERGY ) then
        Integrand ( iS ) % Value  =  CE  -  GC * M * N / R
      end if !-- iI
    end do !-- iS

    end associate !-- CE, etc.
    end select !-- F

  end subroutine ComputeInteriorIntegrandGalilean


  subroutine ComputeBoundaryIntegrandGalilean_CSL &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_FM_Form ), intent ( inout ) :: &
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
      iD, &   !-- iDimension
      iF, &   !-- iFace
      iC, &   !-- iConnectivity
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
      iFluence, &
      iDensity, &
      iEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      nB    !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      X_1, X_2, X_3

    call T % Tally_F_P_Form % ComputeBoundaryIntegrandGalilean_CSL &
           ( Integrand, C, CSL, G, BoundaryFluence )

    select type ( C )
    class is ( Fluid_P_P_Form )
       
    do iFluence = 1, C % N_CONSERVED
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_DENSITY ) &
        iDensity = iFluence
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_ENERGY ) &
        iEnergy = iFluence
    end do !-- iFluence

    associate ( Cnnct => CSL % Atlas % Connectivity )
    do iD = 1, CSL % nDimensions

      nB = shape ( BoundaryFluence ( 1, Cnnct % iaInner ( iD ) ) % Value )

      !-- Geometry
      allocate ( X_1 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_2 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_3 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

      do iF = 1, 2

        if ( iF == 1 ) then
          iC = Cnnct % iaInner ( iD )
        else if ( iF == 2 ) then
          iC = Cnnct % iaOuter ( iD )
        end if

        associate &
          ( D  => BoundaryFluence ( iDensity, iC ) % Value, &
            CE => BoundaryFluence ( iEnergy, iC ) % Value, &
            M  => T % CentralMass, &
            GC => CONSTANT % GRAVITATIONAL )

        call T % ComputeFacePositions ( CSL, G, iD, iF, X_1, X_2, X_3 )

        do iS = 1, T % nSelected
          iI = T % iaSelected ( iS )
          if ( iI == T % GRAVITATIONAL_ENERGY ) then
            Integrand ( iS, iC ) % Value  =  - GC * M * D / X_1
          else if ( iI == T % TOTAL_ENERGY ) then
            Integrand ( iS, iC ) % Value  =  CE  -  GC * M * D / X_1
          end if !-- iI
        end do !-- iS

        end associate !-- D, etc.
      end do !-- iF
      deallocate ( X_1, X_2, X_3 )
    end do !-- iD
    end associate !-- Cnnct

    end select !-- C

  end subroutine ComputeBoundaryIntegrandGalilean_CSL


end module Tally_FM__Form
