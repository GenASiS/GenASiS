module Tally_F_D__Form

  !-- Tally_Fluid_Dust__Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergyBasics
  use Fluid_D__Form

  implicit none
  private
  
    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_DUST = 10

  type, public, extends ( Tally_C_Form ) :: Tally_F_D_Form
    integer ( KDI ) :: &
      N_INTEGRALS_DUST     = N_INTEGRALS_DUST, &
      BARYON_NUMBER        = 0, &
      KINETIC_ENERGY       = 0, &
      GRAVITATIONAL_ENERGY = 0, &
      TOTAL_ENERGY         = 0
    integer ( KDI ), dimension ( 3 ) :: &
      MOMENTUM = 0, &
      ANGULAR_MOMENTUM = 0
  contains
    procedure, private, pass :: &
      InitializeFluid
    generic, public :: &
      Initialize => InitializeFluid
    procedure, public, pass :: &
      SelectVariables
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeInteriorIntegrand
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_CSL
    procedure, public, pass :: &
      ComputeInteriorIntegrand_G
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_CSL_G
    procedure, public, pass :: &
      ComputeInteriorIntegrand_N
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_CSL_N
  end type Tally_F_D_Form

    private :: &
      ComputeDensity_KE, &
      ComputeDensity_AM_Rectangular, &
      ComputeDensity_LM_CylindricalHorizontal, &
      ComputeDensity_AM_CylindricalHorizontal, &
      ComputeDensity_LM_SphericalHorizontal, &
      ComputeDensity_LM_SphericalVertical, &
      ComputeDensity_AM_SphericalHorizontal, &
      ComputeFluence_AM_Rectangular, &
      ComputeFluence_LM_CylindricalHorizontal, &
      ComputeFluence_AM_CylindricalHorizontal, &
      ComputeFluence_LM_SphericalVertical, &
      ComputeFluence_AM_SphericalHorizontal, &
      ComputeFluence_GE_Spherical

contains


  subroutine InitializeFluid ( T, A, Units )
    
    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units

    integer ( KDI ) :: &
      iU  !-- iUnit

    if ( T % N_INTEGRALS == 0 ) &
      T % N_INTEGRALS = T % N_INTEGRALS_DUST
    
    T % BARYON_NUMBER        = 1
    T % MOMENTUM             = [ 2, 3, 4 ]
    T % KINETIC_ENERGY       = 5
    T % ANGULAR_MOMENTUM     = [ 6, 7, 8 ]
    T % GRAVITATIONAL_ENERGY = 9
    T % TOTAL_ENERGY         = 10
    
    if ( .not. allocated ( T % Value ) ) then
      allocate ( T % Value ( T % N_INTEGRALS ) )
      call Clear ( T % Value )
    end if
    
    if ( .not. allocated ( T % Variable ) ) &
      allocate ( T % Variable ( T % N_INTEGRALS ) )
    
    T % Variable ( 1 : T % N_INTEGRALS_DUST ) &
      = [ 'BaryonNumber       ', &
          'Momentum_1         ', &
          'Momentum_2         ', &
          'Momentum_3         ', &
          'KineticEnergy      ', &
          'AngularMomentum_1  ', &
          'AngularMomentum_2  ', &
          'AngularMomentum_3  ', &
          'GravitationalEnergy', &
          'TotalEnergy        ' ]                    
          
    if ( .not. allocated ( T % Unit ) ) &
      allocate ( T % Unit ( T % N_INTEGRALS ) )
    
    T % Unit ( 1 : T % N_INTEGRALS_DUST ) &
      = [ Units % Number, spread ( Units % Momentum, 1, 3 ), Units % Energy, &
          spread ( Units % AngularMomentum, 1, 3 ), Units % Energy, &
          Units % Energy ]
    
    call T % SelectVariables ( A )

    T % Atlas => A

  end subroutine InitializeFluid
  
  
  subroutine SelectVariables ( T, A ) 
    
    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ) :: &
      A

    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( A )
    class is ( Atlas_SC_Form )

    G => A % Geometry ( )

    select type ( G )
    type is ( Geometry_G_Form )
      T % nSelected = 8
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, T % MOMENTUM, T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM ]
    
    class is ( Geometry_N_Form )
      T % nSelected = 10
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, T % MOMENTUM, T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, T % GRAVITATIONAL_ENERGY, &
            T % TOTAL_ENERGY ]
    
    class default 
      call Show ( 'Geometry type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
  
    class default 
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    nullify ( G )

  end subroutine SelectVariables
  

  impure elemental subroutine Finalize ( T )

    type ( Tally_F_D_Form ), intent ( inout ) :: &
      T

    !-- Trigger finalization of parent

  end subroutine Finalize


  subroutine ComputeInteriorIntegrand ( T, Integrand, C, G, nDimensions )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Integrand
    class ( CurrentTemplate ), intent ( in ) :: &
      C
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    select type ( G )
    type is ( Geometry_G_Form )
      call T % ComputeInteriorIntegrand_G &
             ( Integrand, C, G, nDimensions )
    class is ( Geometry_N_Form )
      call T % ComputeInteriorIntegrand_N &
             ( Integrand, C, G, nDimensions )
    class default 
      call Show ( 'Geometry type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeInteriorIntegrand', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

  end subroutine ComputeInteriorIntegrand


  subroutine ComputeBoundaryIntegrand_CSL &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
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

    select type ( G )
    type is ( Geometry_G_Form )
      call T % ComputeBoundaryIntegrand_CSL_G &
             ( Integrand, C, CSL, G, BoundaryFluence )
    class is ( Geometry_N_Form )
      call T % ComputeBoundaryIntegrand_CSL_N &
             ( Integrand, C, CSL, G, BoundaryFluence )
    class default 
      call Show ( 'Geometry type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeInteriorIntegrand', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    
  end subroutine ComputeBoundaryIntegrand_CSL


  subroutine ComputeInteriorIntegrand_G &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
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
    class is ( Fluid_D_Form )
       
    associate &
      ( D   => C % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
        V_1 => C % Value ( :, C % VELOCITY_U ( 1 ) ), &
        V_2 => C % Value ( :, C % VELOCITY_U ( 2 ) ), &
        V_3 => C % Value ( :, C % VELOCITY_U ( 3 ) ), &
        S_1 => C % Value ( :, C % MOMENTUM_DENSITY_D ( 1 ) ), & 
        S_2 => C % Value ( :, C % MOMENTUM_DENSITY_D ( 2 ) ), & 
        S_3 => C % Value ( :, C % MOMENTUM_DENSITY_D ( 3 ) ), &
        X_1 => G % Value ( :, G % CENTER_U ( 1 ) ), &
        X_2 => G % Value ( :, G % CENTER_U ( 2 ) ), &
        X_3 => G % Value ( :, G % CENTER_U ( 3 ) ) )

    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % BARYON_NUMBER ) then
        call Copy ( D, Integrand ( iS ) % Value )
      else if ( iI == T % KINETIC_ENERGY ) then
        call ComputeDensity_KE &
               ( S_1, S_2, S_3, V_1, V_2, V_3, Integrand ( iS ) % Value )
      end if !-- iI
    end do !-- iS

    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % TOTAL_ENERGY ) then
        call Copy ( Integrand ( T % KINETIC_ENERGY ) % Value, &
                    Integrand ( iS ) % Value )
      end if !-- iI
    end do !-- iS

    select case ( trim ( G % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % MOMENTUM ( 1 ) ) then
          call Copy ( S_1, Integrand ( iS ) % Value )
        else if ( iI == T % MOMENTUM ( 2 ) ) then
          call Copy ( S_2, Integrand ( iS ) % Value )
        else if ( iI == T % MOMENTUM ( 3 ) ) then
          call Copy ( S_3, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
          call ComputeDensity_AM_Rectangular &
                 ( X_2, X_3, S_2, S_3, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
          call ComputeDensity_AM_Rectangular &
                 ( X_3, X_1, S_3, S_1, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
          call ComputeDensity_AM_Rectangular &
                 ( X_1, X_2, S_1, S_2, Integrand ( iS ) % Value )
        end if !-- iI
      end do !-- iS
    case ( 'CYLINDRICAL' )
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % MOMENTUM ( 1 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_LM_CylindricalHorizontal &
                   ( X_1, X_3, S_1, S_3, 1, Integrand ( iS ) % Value )
        else if ( iI == T % MOMENTUM ( 2 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_LM_CylindricalHorizontal &
                   ( X_1, X_3, S_1, S_3, 2, Integrand ( iS ) % Value )
        else if ( iI == T % MOMENTUM ( 3 ) ) then
          if ( nDimensions > 1 ) &
            call Copy ( S_2, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_AM_CylindricalHorizontal &
                   ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
                     Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_AM_CylindricalHorizontal &
                   ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
                     Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
          if ( nDimensions > 1 ) &
            call Copy ( S_3, Integrand ( iS ) % Value )
        end if !-- iI
      end do !-- iS
    case ( 'SPHERICAL' )
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % MOMENTUM ( 1 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_LM_SphericalHorizontal &
                   ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
                     Integrand ( iS ) % Value )
        else if ( iI == T % MOMENTUM ( 2 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_LM_SphericalHorizontal &
                   ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
                     Integrand ( iS ) % Value )
        else if ( iI == T % MOMENTUM ( 3 ) ) then
          if ( nDimensions > 1 ) &
            call ComputeDensity_LM_SphericalVertical &
                   ( X_1, X_2, S_1, S_2, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_AM_SphericalHorizontal &
                   ( X_2, X_3, S_2, S_3, 1, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
          if ( nDimensions > 2 ) &
            call ComputeDensity_AM_SphericalHorizontal &
                   ( X_2, X_3, S_2, S_3, 2, Integrand ( iS ) % Value )
        else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
          if ( nDimensions > 1 ) &
            call Copy ( S_3, Integrand ( iS ) % Value )
        end if !-- iI
      end do !-- iS
    end select !-- CoordinateSystem

    end associate !-- N, etc.
    end select !-- C
 
  end subroutine ComputeInteriorIntegrand_G


  subroutine ComputeBoundaryIntegrand_CSL_G &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
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
      iD, &   !-- iDimension
      iF, &   !-- iFace
      iC, &   !-- iConnectivity
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
      iFluence, &
      iDensity, &
      iMomentum_1, &
      iMomentum_2, &
      iMomentum_3
    integer ( KDI ), dimension ( 3 ) :: &
      nB    !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      X_1, X_2, X_3

    select type ( C )
    class is ( Fluid_D_Form )
       
    do iFluence = 1, C % N_CONSERVED
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_BARYON_DENSITY ) &
        iDensity = iFluence
      if ( C % iaConserved ( iFluence ) == C % MOMENTUM_DENSITY_D ( 1 ) ) &
        iMomentum_1 = iFluence
      if ( C % iaConserved ( iFluence ) == C % MOMENTUM_DENSITY_D ( 2 ) ) &
        iMomentum_2 = iFluence
      if ( C % iaConserved ( iFluence ) == C % MOMENTUM_DENSITY_D ( 3 ) ) &
        iMomentum_3 = iFluence
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
          ( D   => BoundaryFluence ( iDensity,    iC ) % Value, &
            S_1 => BoundaryFluence ( iMomentum_1, iC ) % Value, &
            S_2 => BoundaryFluence ( iMomentum_2, iC ) % Value, &
            S_3 => BoundaryFluence ( iMomentum_3, iC ) % Value )
        call T % ComputeFacePositions ( CSL, G, iD, iF, X_1, X_2, X_3 )
        select case ( trim ( G % CoordinateSystem ) )
        case ( 'RECTANGULAR' )
          do iS = 1, T % nSelected
            iI = T % iaSelected ( iS )
            if ( iI == T % BARYON_NUMBER ) then
              call CopyCollapse ( D, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 1 ) ) then
              call CopyCollapse ( S_1, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 2 ) ) then
              call CopyCollapse ( S_2, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 3 ) ) then
              call CopyCollapse ( S_3, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
              call ComputeFluence_AM_Rectangular &
                     ( X_2, X_3, S_2, S_3, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
              call ComputeFluence_AM_Rectangular &
                     ( X_3, X_1, S_3, S_1, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
              call ComputeFluence_AM_Rectangular &
                     ( X_1, X_2, S_1, S_2, Integrand ( iS, iC ) % Value )
            end if !-- iI
          end do !-- iS
        case ( 'CYLINDRICAL' )
          do iS = 1, T % nSelected
            iI = T % iaSelected ( iS )
            if ( iI == T % BARYON_NUMBER ) then
              call CopyCollapse ( D, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 1 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_LM_CylindricalHorizontal &
                       ( X_1, X_3, S_1, S_3, 1, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 2 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_LM_CylindricalHorizontal &
                       ( X_1, X_3, S_1, S_3, 2, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 3 ) ) then
              if ( CSL % nDimensions > 1 ) &
                call Copy ( S_2, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_AM_CylindricalHorizontal &
                       ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
                         Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_AM_CylindricalHorizontal &
                       ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
                         Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
              if ( CSL % nDimensions > 1 ) &
                call Copy ( S_3, Integrand ( iS, iC ) % Value )
            end if !-- iI
          end do !-- iS
        case ( 'SPHERICAL' )
          do iS = 1, T % nSelected
            iI = T % iaSelected ( iS )
            if ( iI == T % BARYON_NUMBER ) then
              call CopyCollapse ( D, Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 1 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_LM_SphericalHorizontal &
                       ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
                         Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 2 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_LM_SphericalHorizontal &
                       ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
                         Integrand ( iS, iC ) % Value )
            else if ( iI == T % MOMENTUM ( 3 ) ) then
              if ( CSL % nDimensions > 1 ) &
                call ComputeFluence_LM_SphericalVertical &
                       ( X_1, X_2, S_1, S_2, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_AM_SphericalHorizontal &
                       ( X_2, X_3, S_2, S_3, 1, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
              if ( CSL % nDimensions > 2 ) &
                call ComputeFluence_AM_SphericalHorizontal &
                       ( X_2, X_3, S_2, S_3, 2, Integrand ( iS, iC ) % Value )
            else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
              if ( CSL % nDimensions > 1 ) &
                call CopyCollapse ( S_3, Integrand ( iS, iC ) % Value )
            end if !-- iI
          end do !-- iS
        end select !-- CoordinateSystem
        end associate !-- D, etc.

      end do !-- iF
      deallocate ( X_1, X_2, X_3 )
    end do !-- iD
    end associate !-- Cnnct

    end select !-- C
       
  end subroutine ComputeBoundaryIntegrand_CSL_G


  subroutine ComputeInteriorIntegrand_N &
               ( T, Integrand, C, G, nDimensions )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
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

    call T % ComputeInteriorIntegrand_G ( Integrand, C, G, nDimensions )

    select type ( A => T % Atlas )
    class is ( Atlas_SC_Form )    
    select type ( GA => A % Geometry_ASC )
    class is ( Geometry_ASC_Form )

    select type ( C )
    class is ( Fluid_D_Form )
    select type ( G )
    class is ( Geometry_N_Form )

    associate &
      ( M   => C % Value ( :, C % BARYON_MASS ), &
        D   => C % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
        Phi => G % Value ( :, G % POTENTIAL ) )

    select case ( trim ( GA % GravitySolverType ) )
    case ( 'UNIFORM', 'CENTRAL_MASS' )  !-- External potential
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % GRAVITATIONAL_ENERGY ) then
          Integrand ( iS ) % Value  =  M * D * Phi
        else if ( iI == T % TOTAL_ENERGY ) then
          Integrand ( iS ) % Value  &
            =  Integrand ( iS ) % Value  +  M * D * Phi
        end if !-- iI
      end do !-- iS     
    case default
      do iS = 1, T % nSelected
        iI = T % iaSelected ( iS )
        if ( iI == T % GRAVITATIONAL_ENERGY ) then
          Integrand ( iS ) % Value  =  0.5_KDR * M * D * Phi
        else if ( iI == T % TOTAL_ENERGY ) then
          Integrand ( iS ) % Value  &
            =  Integrand ( iS ) % Value  +  0.5_KDR * M * D * Phi
        end if !-- iI
      end do !-- iS     
    end select !-- GravitySolverType

    end associate !-- M, etc.

    class default 
      call Show ( 'Geometry type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeInteriorIntegrand_N', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    end select !-- C
    end select !-- GA
    end select !-- A

  end subroutine ComputeInteriorIntegrand_N


  subroutine ComputeBoundaryIntegrand_CSL_N &
               ( T, Integrand, C, CSL, G, BoundaryFluence )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
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
      iD, &   !-- iDimension
      iF, &   !-- iFace
      iC, &   !-- iConnectivity
      iS, &   !-- iSelected
      iI, &   !-- iIntegral
      iFluence, &
      iDensity, &
      iGravity
    integer ( KDI ), dimension ( 3 ) :: &
      nB    !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      X_1, X_2, X_3
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Phi, &
      R, &
      M

    call T % ComputeBoundaryIntegrand_CSL_G &
           ( Integrand, C, CSL, G, BoundaryFluence )

    select type ( C )
    class is ( Fluid_D_Form )
    select type ( G )
    type is ( Geometry_N_Form )

    do iFluence = 1, C % N_CONSERVED
      if ( C % iaConserved ( iFluence ) == C % CONSERVED_BARYON_DENSITY ) &
        iDensity = iFluence
    end do !-- iFluence

    call CSL % SetVariablePointer ( G % Value ( :, G % POTENTIAL ), Phi )
    call CSL % SetVariablePointer ( C % Value ( :, C % BARYON_MASS ), M )

    associate ( Cnnct => CSL % Atlas % Connectivity )
    DimensionLoop: do iD = 1, CSL % nDimensions
      nB = shape ( BoundaryFluence ( 1, Cnnct % iaInner ( iD ) ) % Value )

      !-- Geometry
      allocate ( X_1 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_2 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_3 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

      FaceLoop: do iF = 1, 2
        if ( iF == 1 ) then
          iC = Cnnct % iaInner ( iD )
        else if ( iF == 2 ) then
          iC = Cnnct % iaOuter ( iD )
        end if

        associate ( D => BoundaryFluence ( iDensity, iC ) % Value )
        call T % ComputeFacePositions ( CSL, G, iD, iF, X_1, X_2, X_3 )
        select case ( trim ( G % CoordinateSystem ) )
        case ( 'SPHERICAL' )

          if ( iD /= 1 ) &
            exit DimensionLoop
!          if ( iF /= 2 ) &
!            cycle FaceLoop
          do iS = 1, T % nSelected
            iI = T % iaSelected ( iS )
            call CSL % SetVariablePointer &
                   ( G % Value ( :, G % CENTER_U ( 1 ) ), R )
            if ( iI == T % GRAVITATIONAL_ENERGY ) then
              call ComputeFluence_GE_Spherical &
                     ( CSL, M, Phi, R, D, X_1, Integrand ( iS, iC ) % Value )
              iGravity = iS
            end if !-- iI
          end do !-- iS

          do iS = 1, T % nSelected
            iI = T % iaSelected ( iS )
            if ( iI == T % TOTAL_ENERGY ) then
              Integrand ( iS, iF ) % Value  &
                =  Integrand ( iS, iF ) % Value  &
                   +  Integrand ( iGravity, iF ) % Value
            end if !-- iI
          end do !-- iS

        end select !-- CoordinateSystem

        end associate !-- D, etc.
      end do FaceLoop !-- iF
      deallocate ( X_1, X_2, X_3 )
    end do DimensionLoop !-- iD
    end associate !-- Cnnct

    end select !-- G
    end select !-- C
    nullify ( Phi )

  end subroutine ComputeBoundaryIntegrand_CSL_N


  subroutine ComputeDensity_KE ( S_1, S_2, S_3, V_1, V_2, V_3, I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      S_1, S_2, S_3, &
      V_1, V_2, V_3
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( I )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      I ( iV ) = 0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                             +  S_2 ( iV ) * V_2 ( iV )  &
                             +  S_3 ( iV ) * V_3 ( iV )  )
    end do
    !$OMP end parallel do

  end subroutine ComputeDensity_KE


  subroutine ComputeDensity_AM_Rectangular ( X_J, X_K, S_J, S_K, I_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_J, X_K, &
      S_J, S_K
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( I_I )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      I_I ( iV ) =  X_J ( iV ) * S_K ( iV ) -  X_K ( iV ) * S_J ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeDensity_AM_Rectangular


  subroutine ComputeDensity_LM_CylindricalHorizontal &
               ( X_1, X_3, S_1, S_3, iD, I_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1, X_3, &
      S_1, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = size ( I_I )
    
    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     cos ( X_3 ( iV ) )  *  S_1 ( iV )  &
                       -  sin ( X_3 ( iV ) )  &
                          /  max ( X_1 ( iV ), SqrtTiny )  &
                          *  S_3 ( iV ) 
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     sin ( X_3 ( iV ) )  *  S_1 ( iV )  &
                       +  cos ( X_3 ( iV ) )  &
                          /  max ( X_1 ( iV ), SqrtTiny )  &
                          *  S_3 ( iV ) 
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeDensity_LM_CylindricalHorizontal


  subroutine ComputeDensity_AM_CylindricalHorizontal &
               ( X_1, X_2, X_3, S_1, S_2, S_3, iD, I_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1, X_2, X_3, &
      S_1, S_2, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = size ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =  sin ( X_3 ( iV ) )  *  (   X_1 ( iV )  *  S_2 ( iV )  &
                                                - X_2 ( iV )  *  S_1 ( iV ) ) &
                       - X_2 ( iV )  *  cos ( X_3 ( iV ) )  &
                         /  max ( X_1 ( iV ), SqrtTiny )  &
                         *  S_3 ( iV )                         
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =  cos ( X_3 ( iV ) )  *  (   X_2 ( iV )  *  S_1 ( iV )  &
                                                - X_1 ( iV )  *  S_2 ( iV ) ) &
                       - X_2 ( iV )  *  sin ( X_3 ( iV ) )  &
                         /  max ( X_1 ( iV ), SqrtTiny )  &
                         *  S_3 ( iV )                         
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeDensity_AM_CylindricalHorizontal


  subroutine ComputeDensity_LM_SphericalHorizontal &
               ( X_1, X_2, X_3, S_1, S_2, S_3, iD, I_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1, X_2, X_3, &
      S_1, S_2, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = size ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     sin ( X_2 ( iV ) )  *  cos ( X_3 ( iV ) ) &
                          * S_1 ( iV ) &
                       +  cos ( X_2 ( iV ) )  *  cos ( X_3 ( iV ) ) &
                          /  max ( X_1 ( iV ), SqrtTiny )  &
                          *  S_2 ( iV ) &
                       -  sin ( X_3 ( iV ) )  &
                          /  max ( X_1 ( iV ) * sin ( X_2 ( iV ) ), &
                                   SqrtTiny )  &
                          *  S_3 ( iV ) 
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     sin ( X_2 ( iV ) )  *  sin ( X_3 ( iV ) ) &
                          *  S_1 ( iV ) &
                       +  cos ( X_2 ( iV ) )  *  sin ( X_3 ( iV ) ) &
                          /  max ( X_1 ( iV ), SqrtTiny )  &
                          *  S_2 ( iV ) &
                       +  cos ( X_3 ( iV ) )  &
                          /  max ( X_1 ( iV ) * sin ( X_2 ( iV ) ), &
                                   SqrtTiny )  &
                          *  S_3 ( iV ) 
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeDensity_LM_SphericalHorizontal


  subroutine ComputeDensity_LM_SphericalVertical ( X_1, X_2, S_1, S_2, I_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_1, X_2, &
      S_1, S_2
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = size ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      I_I ( iV )  =  cos ( X_2 ( iV ) )  *  S_1 ( iV ) &
                     -  ( sin ( X_2 ( iV ) )  &
                          /  max ( X_1 ( iV ), SqrtTiny ) )  &
                        *  S_2 ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeDensity_LM_SphericalVertical


  subroutine ComputeDensity_AM_SphericalHorizontal &
               ( X_2, X_3, S_2, S_3, iD, I_I )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X_2, X_3, &
      S_2, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = size ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =  - sin ( X_3 ( iV ) )  *  S_2 ( iV )  &
                       - ( cos ( X_2 ( iV ) )  *  cos ( X_3 ( iV ) )  &
                           /  sin ( max ( X_2 ( iV ), SqrtTiny ) ) )  &
                         *  S_3 ( iV )  
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =    cos ( X_3 ( iV ) )  *  S_2 ( iV )  &
                       - ( cos ( X_2 ( iV ) )  *  sin ( X_3 ( iV ) )  &
                           /  sin ( max ( X_2 ( iV ), SqrtTiny ) ) )  &
                         *  S_3 ( iV )  
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeDensity_AM_SphericalHorizontal


  subroutine ComputeFluence_AM_Rectangular ( X_J, X_K, S_J, S_K, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      X_J, X_K, &
      S_J, S_K
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( I_I )

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          I_I ( iV, jV, kV ) &
            =  X_J ( iV, jV, kV ) * S_K ( iV, jV, kV ) &
               -  X_K ( iV, jV, kV ) * S_J ( iV, jV, kV )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeFluence_AM_Rectangular


  subroutine ComputeFluence_LM_CylindricalHorizontal &
               ( X_1, X_3, S_1, S_3, iD, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      X_1, X_3, &
      S_1, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = shape ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =  cos ( X_3 ( iV, jV, kV ) )  *  S_1 ( iV, jV, kV )  &
                 -  sin ( X_3 ( iV, jV, kV ) )  &
                    /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
                    *  S_3 ( iV, jV, kV ) 
          end do
        end do
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =  sin ( X_3 ( iV, jV, kV ) )  *  S_1 ( iV, jV, kV )  &
                 +  cos ( X_3 ( iV, jV, kV ) )  &
                    /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
                    *  S_3 ( iV, jV, kV ) 
          end do
        end do
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeFluence_LM_CylindricalHorizontal


  subroutine ComputeFluence_AM_CylindricalHorizontal &
               ( X_1, X_2, X_3, S_1, S_2, S_3, iD, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      X_1, X_2, X_3, &
      S_1, S_2, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = shape ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =  sin ( X_3 ( iV, jV, kV ) )  &
                   *  (   X_1 ( iV, jV, kV )  *  S_2 ( iV, jV, kV )  &
                        - X_2 ( iV, jV, kV )  *  S_1 ( iV, jV, kV ) ) &
                 -  X_2 ( iV, jV, kV )  *  cos ( X_3 ( iV, jV, kV ) )  &
                    /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
                    *  S_3 ( iV, jV, kV )                         
          end do
        end do
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =  cos ( X_3 ( iV, jV, kV ) )  &
                   *  (   X_2 ( iV, jV, kV )  *  S_1 ( iV, jV, kV )  &
                        - X_1 ( iV, jV, kV )  *  S_2 ( iV, jV, kV ) ) &
                 -  X_2 ( iV, jV, kV )  *  sin ( X_3 ( iV, jV, kV ) )  &
                    /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
                    *  S_3 ( iV, jV, kV )                         
          end do
        end do
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeFluence_AM_CylindricalHorizontal


  subroutine ComputeFluence_LM_SphericalHorizontal &
               ( X_1, X_2, X_3, S_1, S_2, S_3, iD, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      X_1, X_2, X_3, &
      S_1, S_2, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = shape ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  =     sin ( X_2 ( iV, jV, kV ) )  &
                                      *  cos ( X_3 ( iV, jV, kV ) ) &
                                      *  S_1 ( iV, jV, kV ) &
                                   +  cos ( X_2 ( iV, jV, kV ) )  &
                                      *  cos ( X_3 ( iV, jV, kV ) ) &
                                      /  max ( X_1 ( iV, jV, kV ), SqrtTiny ) &
                                      *  S_2 ( iV, jV, kV ) &
                                   -  sin ( X_3 ( iV, jV, kV ) )  &
                                      /  ( max ( X_1 ( iV, jV, kV )  &
                                                 * sin ( X_2 ( iV, jV, kV ) ),&
                                                 SqrtTiny ) )  &
                                      *  S_3 ( iV, jV, kV ) 
          end do
        end do
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  =     sin ( X_2 ( iV, jV, kV ) )  &
                                      *  sin ( X_3 ( iV, jV, kV ) ) &
                                      *  S_1 ( iV, jV, kV ) &
                                   +  cos ( X_2 ( iV, jV, kV ) )  &
                                      *  sin ( X_3 ( iV, jV, kV ) ) &
                                      /  max ( X_1 ( iV, jV, kV ), &
                                               SqrtTiny )  &
                                      *  S_2 ( iV, jV, kV ) &
                                   +  cos ( X_3 ( iV, jV, kV ) )  &
                                      /  ( max ( X_1 ( iV, jV, kV )  &
                                                 * sin ( X_2 ( iV, jV, kV ) ),&
                                                 SqrtTiny ) )  &
                                      *  S_3 ( iV, jV, kV ) 
          end do
        end do
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeFluence_LM_SphericalHorizontal


  subroutine ComputeFluence_LM_SphericalVertical ( X_1, X_2, S_1, S_2, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      X_1, X_2, &
      S_1, S_2
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = shape ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          I_I ( iV, jV, kV ) &
            =  cos ( X_2 ( iV, jV, kV ) )  *   S_1 ( iV, jV, kV ) &
               -  sin ( X_2 ( iV, jV, kV ) ) &
                    /  max ( X_1 ( iV, jV, kV ), SqrtTiny ) &
                  *  S_2 ( iV, jV, kV )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeFluence_LM_SphericalVertical

  
  subroutine ComputeFluence_AM_SphericalHorizontal &
               ( X_2, X_3, S_2, S_3, iD, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      X_2, X_3, &
      S_2, S_3
    integer ( KDI ), intent ( in ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny

    nV = shape ( I_I )

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =  - sin ( X_3 ( iV, jV, kV ) )  *  S_2 ( iV, jV, kV )  &
                 - ( cos ( X_2 ( iV, jV, kV ) )  &
                     *  cos ( X_3 ( iV, jV, kV ) )  &
                     /  max ( sin ( X_2 ( iV, jV, kV ) ), SqrtTiny ) )  &
                   *  S_3 ( iV, jV, kV )  
          end do
        end do
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =    cos ( X_3 ( iV, jV, kV ) )  *  S_2 ( iV, jV, kV )  &
                 - ( cos ( X_2 ( iV, jV, kV ) )  &
                     *  sin ( X_3 ( iV, jV, kV ) )  &
                     /  max ( sin ( X_2 ( iV, jV, kV ) ), SqrtTiny ) )  &
                   *  S_3 ( iV, jV, kV )  
          end do
        end do
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeFluence_AM_SphericalHorizontal


  subroutine ComputeFluence_GE_Spherical ( CSL, M, Phi, R, D, R_I, I_I )

    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL    
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      M, &
      Phi, &
      R, &
      D, &
      R_I
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      oI, oJ, oK, &
      jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    nV  =  shape ( I_I )

    oI  =  CSL % nGhostLayers ( 1 )  +  CSL % nCellsBrick ( 1 )  -  1
    oJ  =  CSL % nGhostLayers ( 2 )
    oK  =  CSL % nGhostLayers ( 3 )

    !$OMP parallel do private ( jV, kV ) collapse ( 2 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
          I_I ( 1, jV, kV )  &
            =  Phi ( oI + 1, oJ + jV, oK + kV )  &
               *  R ( oI + 1, oJ + jV, oK + kV )  &
                  /  max ( R_I ( 1, jV, kV ), SqrtTiny )  &
               *  M ( oI + 1, oJ + jV, oK + kV )  *  D ( 1, jV, kV )
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeFluence_GE_Spherical


end module Tally_F_D__Form
