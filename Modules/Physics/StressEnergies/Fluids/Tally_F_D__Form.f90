module Tally_F_D__Form

  !-- Tally_Fluid_Dust__Form

  use Basics
  use Mathematics
  use Spaces
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
      ComputeFluence_AM_SphericalHorizontal

contains


  subroutine InitializeFluid &
               ( T, A, NumberUnitOption, EnergyUnitOption, &
                 MomentumUnitOption, AngularMomentumUnitOption )
    
    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      NumberUnitOption, &
      EnergyUnitOption, &
      MomentumUnitOption, &
      AngularMomentumUnitOption

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
    
    if ( present ( NumberUnitOption ) .and. present ( MomentumUnitOption ) &
         .and. present ( AngularMomentumUnitOption ) ) &
    then
      T % Unit ( 1 : T % N_INTEGRALS_DUST ) &
        = [ NumberUnitOption, spread ( MomentumUnitOption, 1, 3 ), &
            EnergyUnitOption, spread ( AngularMomentumUnitOption, 1, 3 ), &
            EnergyUnitOption, EnergyUnitOption ]
    else
!-- FIXME: NAG chokes on this, so use explicit loop
!      T % Unit ( 1 : T % N_FIELDS_DUST ) = UNIT % IDENTITY
      do iU = 1, T % N_INTEGRALS_DUST
        T % Unit ( iU ) = UNIT % IDENTITY
      end do
    end if
    
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
    
    type is ( Geometry_N_Form )
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
    type is ( Geometry_N_Form )
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
    type is ( Geometry_N_Form )
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

    select type ( C )
    class is ( Fluid_D_Form )
    select type ( G )
    type is ( Geometry_N_Form )

    associate &
      ( M   => C % Value ( :, C % BARYON_MASS ), &
        D   => C % Value ( :, C % CONSERVED_BARYON_DENSITY ), &
        Phi => G % Value ( :, G % GRAVITATIONAL_POTENTIAL ) )

    do iS = 1, T % nSelected
      iI = T % iaSelected ( iS )
      if ( iI == T % GRAVITATIONAL_ENERGY ) then
        Integrand ( iS ) % Value  =  0.5_KDR * M * D * Phi
      else if ( iI == T % TOTAL_ENERGY ) then
        Integrand ( iS ) % Value  &
          =  Integrand ( T % KINETIC_ENERGY ) % Value  +  0.5_KDR * M * D * Phi
      end if !-- iI
    end do !-- iS     

    end associate !-- D, etc.
    end select !-- G
    end select !-- C

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

    call T % ComputeBoundaryIntegrand_CSL_G &
           ( Integrand, C, CSL, G, BoundaryFluence )

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

    nV = size ( I_I )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     cos ( X_3 ( iV ) )  *  S_1 ( iV )  &
                       -  sin ( X_3 ( iV ) )  /  X_1 ( iV )  &
                          *  S_3 ( iV ) 
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     sin ( X_3 ( iV ) )  *  S_1 ( iV )  &
                       +  cos ( X_3 ( iV ) )  /  X_1 ( iV )  &
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

    nV = size ( I_I )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =  sin ( X_3 ( iV ) )  *  (   X_1 ( iV )  *  S_2 ( iV )  &
                                                - X_2 ( iV )  *  S_1 ( iV ) ) &
                       - X_2 ( iV )  *  cos ( X_3 ( iV ) )  /  X_1 ( iV )  &
                         *  S_3 ( iV )                         
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =  cos ( X_3 ( iV ) )  *  (   X_2 ( iV )  *  S_1 ( iV )  &
                                                - X_1 ( iV )  *  S_2 ( iV ) ) &
                       - X_2 ( iV )  *  sin ( X_3 ( iV ) )  /  X_1 ( iV )  &
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

    nV = size ( I_I )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     sin ( X_2 ( iV ) )  *  cos ( X_3 ( iV ) ) &
                            * S_1 ( iV ) &
                       +  cos ( X_2 ( iV ) )  *  cos ( X_3 ( iV ) ) &
                            /  X_1 ( iV )  *  S_2 ( iV ) &
                       -  sin ( X_3 ( iV ) )  &
                            /  ( X_1 ( iV ) * sin ( X_2 ( iV ) ) )  &
                            *  S_3 ( iV ) 
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =     sin ( X_2 ( iV ) )  *  sin ( X_3 ( iV ) ) &
                            * S_1 ( iV ) &
                       +  cos ( X_2 ( iV ) )  *  sin ( X_3 ( iV ) ) &
                            /  X_1 ( iV )  *  S_2 ( iV ) &
                       +  cos ( X_3 ( iV ) )  &
                            /  ( X_1 ( iV ) * sin ( X_2 ( iV ) ) )  &
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

    nV = size ( I_I )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      I_I ( iV )  =  cos ( X_2 ( iV ) )  *  S_1 ( iV ) &
                     -  ( sin ( X_2 ( iV ) )  /  X_1 ( iV ) )  *  S_2 ( iV )
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

    nV = size ( I_I )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =  - sin ( X_3 ( iV ) )  *  S_2 ( iV )  &
                       - ( cos ( X_2 ( iV ) )  *  cos ( X_3 ( iV ) )  &
                           /  sin ( X_2 ( iV ) ) )  *  S_3 ( iV )  
      end do
      !$OMP end parallel do
    case ( 2 )
      !$OMP parallel do private ( iV )
      do iV = 1, nV
        I_I ( iV )  =    cos ( X_3 ( iV ) )  *  S_2 ( iV )  &
                       - ( cos ( X_2 ( iV ) )  *  sin ( X_3 ( iV ) )  &
                           /  sin ( X_2 ( iV ) ) )  *  S_3 ( iV )  
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

    nV = shape ( I_I )

    select case ( iD )
    case ( 1 )
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            I_I ( iV, jV, kV )  &
              =  cos ( X_3 ( iV, jV, kV ) )  *  S_1 ( iV, jV, kV )  &
                 -  sin ( X_3 ( iV, jV, kV ) )  /  X_1 ( iV, jV, kV )  &
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
                 +  cos ( X_3 ( iV, jV, kV ) )  /  X_1 ( iV, jV, kV )  &
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

    nV = shape ( I_I )

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
                    /  X_1 ( iV, jV, kV )  *  S_3 ( iV, jV, kV )                         
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
                    /  X_1 ( iV, jV, kV )  *  S_3 ( iV, jV, kV )                         
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

    nV = shape ( I_I )

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
                                      /  X_1 ( iV, jV, kV )  &
                                      *  S_2 ( iV, jV, kV ) &
                                   -  sin ( X_3 ( iV, jV, kV ) )  &
                                      /  ( X_1 ( iV, jV, kV ) &
                                           * sin ( X_2 ( iV, jV, kV ) ) )  &
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
                                      /  X_1 ( iV, jV, kV )  &
                                      *  S_2 ( iV, jV, kV ) &
                                   +  cos ( X_3 ( iV, jV, kV ) )  &
                                      /  ( X_1 ( iV, jV, kV ) &
                                           * sin ( X_2 ( iV, jV, kV ) ) )  &
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

    nV = shape ( I_I )

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          I_I ( iV, jV, kV ) &
            =  cos ( X_2 ( iV, jV, kV ) )  *   S_1 ( iV, jV, kV ) &
               -  sin ( X_2 ( iV, jV, kV ) )  /  X_1 ( iV, jV, kV ) &
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

    nV = shape ( I_I )

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
                   /  sin ( X_2 ( iV, jV, kV ) ) )  *  S_3 ( iV, jV, kV )  
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
                   /  sin ( X_2 ( iV, jV, kV ) ) )  *  S_3 ( iV, jV, kV )  
          end do
        end do
      end do
      !$OMP end parallel do
    end select

  end subroutine ComputeFluence_AM_SphericalHorizontal


end module Tally_F_D__Form
