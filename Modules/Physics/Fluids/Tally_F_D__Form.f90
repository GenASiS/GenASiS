module Tally_F_D__Form

  !-- Tally_Fluid_Dust__Form

  use Basics
  use Mathematics
  use Gravitations
  use Units_F__Form
  use Fluid_D__Form

  implicit none
  private
  
    integer ( KDI ), private, parameter :: &
      N_INTEGRALS_D = 10

  type, public, extends ( Tally_CS_Form ) :: Tally_F_D_Form
    integer ( KDI ) :: &
      N_INTEGRALS_D        = N_INTEGRALS_D, &
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
      ComputeBoundaryIntegrand
    procedure, public, pass :: &
      ComputeInteriorIntegrand_G
    procedure, public, pass :: &
      ComputeInteriorIntegrand_N
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_G
    procedure, public, pass :: &
      ComputeBoundaryIntegrand_N
  end type Tally_F_D_Form

    private :: &
      ComputeDensity_KE, &
      ComputeDensity_AM_Rectangular, &
!       ComputeDensity_LM_CylindricalHorizontal, &
!       ComputeDensity_AM_CylindricalHorizontal, &
      ComputeDensity_LM_SphericalHorizontal, &
      ComputeDensity_LM_SphericalVertical, &
      ComputeDensity_AM_SphericalHorizontal, &
      ComputeFluence_KE_Rectangular, &
      ComputeFluence_AM_Rectangular, &
!       ComputeFluence_LM_CylindricalHorizontal, &
!       ComputeFluence_AM_CylindricalHorizontal, &
      ComputeFluence_LM_SphericalVertical, &
      ComputeFluence_KE_Spherical, &
      ComputeFluence_AM_SphericalHorizontal, &
      ComputeFluence_GE_Spherical


contains


  subroutine InitializeFluid ( T, G, Units )
    
    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G
    class ( Units_F_Form ), intent ( in ) :: &
      Units

    if ( T % nIntegrals  ==  0 ) &
      T % nIntegrals  =  T % N_INTEGRALS_D
    
    T % BARYON_NUMBER        = 1
    T % MOMENTUM             = [ 2, 3, 4 ]
    T % KINETIC_ENERGY       = 5
    T % ANGULAR_MOMENTUM     = [ 6, 7, 8 ]
    T % GRAVITATIONAL_ENERGY = 9
    T % TOTAL_ENERGY         = 10
    
    if ( .not. allocated ( T % Value ) ) then
      allocate ( T % Value ( T % nIntegrals ) )
      call Clear ( T % Value )
    end if
    
    if ( .not. allocated ( T % Variable ) ) &
      allocate ( T % Variable ( T % nIntegrals ) )
    
    T % Variable ( 1 : T % N_INTEGRALS_D ) &
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
      allocate ( T % Unit ( T % nIntegrals ) )
    
    T % Unit ( 1 : T % N_INTEGRALS_D ) &
      = [ Units % Number, &
          spread ( Units % Momentum, 1, 3 ), &
          Units % Energy, &
          spread ( Units % AngularMomentum, 1, 3 ), &
          Units % Energy, &
          Units % Energy ]
    
    T % Geometry  =>  G

    call T % SelectVariables ( )

  end subroutine InitializeFluid
  
  
  subroutine SelectVariables ( T ) 
    
    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T

    select type ( G  =>  T % Geometry )
    type is ( Gravitation_G_Form )
      T % nSelected = 8
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, &
            T % MOMENTUM, &
            T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM ]    
    class is ( Gravitation_N_H_Form )
      T % nSelected = 10
      allocate ( T % iaSelected ( T % nSelected ) )
      T % iaSelected &
        = [ T % BARYON_NUMBER, &
            T % MOMENTUM, &
            T % KINETIC_ENERGY, &
            T % ANGULAR_MOMENTUM, &
            T % GRAVITATIONAL_ENERGY, &
            T % TOTAL_ENERGY ]    
    class default 
      call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SelectVariables', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
  
  end subroutine SelectVariables
  

  impure elemental subroutine Finalize ( T )

    type ( Tally_F_D_Form ), intent ( inout ) :: &
      T

  end subroutine Finalize


  subroutine ComputeInteriorIntegrand ( T, CS )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS

    select type ( G  =>  T % Geometry )
    type is ( Gravitation_G_Form )
      call T % ComputeInteriorIntegrand_G ( CS )
    class is ( Gravitation_N_H_Form )
      call T % ComputeInteriorIntegrand_N ( CS )
    class default 
      call Show ( 'Geometry type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeInteriorIntegrand', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G

  end subroutine ComputeInteriorIntegrand


  subroutine ComputeBoundaryIntegrand ( T, CS, C, BF )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BF

    select type ( G  =>  T % Geometry )
    type is ( Gravitation_G_Form )
      call T % ComputeBoundaryIntegrand_G ( CS, C, BF )
    class is ( Gravitation_N_H_Form )
      call T % ComputeBoundaryIntegrand_N ( CS, C, BF )
    class default 
      call Show ( 'Geometry type not recognized', CONSOLE % ERROR )
      call Show ( 'Tally_F_D__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeBoundaryIntegrand', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- G
    
  end subroutine ComputeBoundaryIntegrand


  subroutine ComputeInteriorIntegrand_G ( T, CS )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI, &  !-- iIntegral
      iKineticEnergy

    select type ( CS )
      class is ( Fluid_D_Form )
    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( G    =>  T % Geometry, &
        C    =>  A % Chart_GS, &
        CSV  =>  CS % Storage_GS % Value, &
         GV  =>  T % Geometry % Storage_GS % Value, &
         IV  =>  T % InteriorIntegral % Integrand % Storage_GS % Value )
    associate &
      ( D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
        V_1  =>  CSV ( :, CS % VELOCITY_U_1 ), &
        V_2  =>  CSV ( :, CS % VELOCITY_U_2 ), &
        V_3  =>  CSV ( :, CS % VELOCITY_U_3 ), &
        S_1  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), & 
        S_2  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), & 
        S_3  =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
        X_1  =>   GV ( :, G % CENTER_U_1 ), &
        X_2  =>   GV ( :, G % CENTER_U_2 ), &
        X_3  =>   GV ( :, G % CENTER_U_3 ) )

    do iS  =  1,  T % nSelected
      iI  =  T % iaSelected ( iS )
      if ( iI  ==  T % BARYON_NUMBER ) then
        call Copy ( D, IV ( :, iS ) )
      else if ( iI  ==  T % KINETIC_ENERGY ) then
        iKineticEnergy  =  iS
        call ComputeDensity_KE &
               ( S_1, S_2, S_3, V_1, V_2, V_3, IV ( :, iS ) )
      end if !-- iI
    end do !-- iS

    do iS  =  1, T % nSelected
      iI  =  T % iaSelected ( iS )
      if ( iI  ==  T % TOTAL_ENERGY ) then
        call Copy ( IV ( :, iKineticEnergy ), IV ( :, iS ) )
      end if !-- iI
    end do !-- iS

    select case ( trim ( C % CoordinateSystem ) )
    case ( 'RECTANGULAR' )
      do iS  =  1, T % nSelected
        iI  =  T % iaSelected ( iS )
        if ( iI  ==  T % MOMENTUM ( 1 ) ) then
          call Copy ( S_1, IV ( :, iS ) )
        else if ( iI  ==  T % MOMENTUM ( 2 ) ) then
          call Copy ( S_2, IV ( :, iS ) )
        else if ( iI  ==  T % MOMENTUM ( 3 ) ) then
          call Copy ( S_3, IV ( :, iS ) )
        else if ( iI  ==  T % ANGULAR_MOMENTUM ( 1 ) ) then
          call ComputeDensity_AM_Rectangular &
                 ( X_2, X_3, S_2, S_3, IV ( :, iS ) )
        else if ( iI  ==  T % ANGULAR_MOMENTUM ( 2 ) ) then
          call ComputeDensity_AM_Rectangular &
                 ( X_3, X_1, S_3, S_1, IV ( :, iS ) )
        else if ( iI  ==  T % ANGULAR_MOMENTUM ( 3 ) ) then
          call ComputeDensity_AM_Rectangular &
                 ( X_1, X_2, S_1, S_2, IV ( :, iS ) )
        end if !-- iI
      end do !-- iS
!     case ( 'CYLINDRICAL' )
!       do iS = 1, T % nSelected
!         iI = T % iaSelected ( iS )
!         if ( iI == T % MOMENTUM ( 1 ) ) then
!           if ( nDimensions > 2 ) &
!             call ComputeDensity_LM_CylindricalHorizontal &
!                    ( X_1, X_3, S_1, S_3, 1, Integrand ( iS ) % Value )
!         else if ( iI == T % MOMENTUM ( 2 ) ) then
!           if ( nDimensions > 2 ) &
!             call ComputeDensity_LM_CylindricalHorizontal &
!                    ( X_1, X_3, S_1, S_3, 2, Integrand ( iS ) % Value )
!         else if ( iI == T % MOMENTUM ( 3 ) ) then
!           if ( nDimensions > 1 ) &
!             call Copy ( S_2, Integrand ( iS ) % Value )
!         else if ( iI == T % ANGULAR_MOMENTUM ( 1 ) ) then
!           if ( nDimensions > 2 ) &
!             call ComputeDensity_AM_CylindricalHorizontal &
!                    ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
!                      Integrand ( iS ) % Value )
!         else if ( iI == T % ANGULAR_MOMENTUM ( 2 ) ) then
!           if ( nDimensions > 2 ) &
!             call ComputeDensity_AM_CylindricalHorizontal &
!                    ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
!                      Integrand ( iS ) % Value )
!         else if ( iI == T % ANGULAR_MOMENTUM ( 3 ) ) then
!           if ( nDimensions > 1 ) &
!             call Copy ( S_3, Integrand ( iS ) % Value )
!         end if !-- iI
!       end do !-- iS
    case ( 'SPHERICAL' )
      do iS  =  1, T % nSelected
        iI  =  T % iaSelected ( iS )
        if ( iI  ==  T % MOMENTUM ( 1 ) ) then
          if ( C % nDimensions > 2 ) &
            call ComputeDensity_LM_SphericalHorizontal &
                   ( X_1, X_2, X_3, S_1, S_2, S_3, 1, IV ( :, iS ) )
        else if ( iI  ==  T % MOMENTUM ( 2 ) ) then
          if ( C % nDimensions > 2 ) &
            call ComputeDensity_LM_SphericalHorizontal &
                   ( X_1, X_2, X_3, S_1, S_2, S_3, 2, IV ( :, iS ) )
        else if ( iI  ==  T % MOMENTUM ( 3 ) ) then
          if ( C % nDimensions > 1 ) &
            call ComputeDensity_LM_SphericalVertical &
                   ( X_1, X_2, S_1, S_2, IV ( :, iS ) )
        else if ( iI  ==  T % ANGULAR_MOMENTUM ( 1 ) ) then
          if ( C % nDimensions > 2 ) &
            call ComputeDensity_AM_SphericalHorizontal &
                   ( X_2, X_3, S_2, S_3, 1, IV ( :, iS ) )
        else if ( iI  ==  T % ANGULAR_MOMENTUM ( 2 ) ) then
          if ( C % nDimensions > 2 ) &
            call ComputeDensity_AM_SphericalHorizontal &
                   ( X_2, X_3, S_2, S_3, 2, IV ( :, iS ) )
        else if ( iI  ==  T % ANGULAR_MOMENTUM ( 3 ) ) then
          if ( C % nDimensions > 1 ) &
            call Copy ( S_3, IV ( :, iS ) )
        end if !-- iI
      end do !-- iS
    end select !-- CoordinateSystem

    end associate !-- D, etc.
    end associate !-- CSV, etc.
    end select !-- A
    end select !-- CS
 
  end subroutine ComputeInteriorIntegrand_G


  subroutine ComputeInteriorIntegrand_N ( T, CS )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iI     !-- iIntegral

    call T % ComputeInteriorIntegrand_G ( CS )

    select type ( CS )
      class is ( Fluid_D_Form )
    select type ( G  =>  T % Geometry )
      class is ( Gravitation_N_H_Form )
    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( CSV  =>  CS % Storage_GS % Value, &
         GV  =>  G % Storage_GS % Value, &
         IV  =>  T % InteriorIntegral % Integrand % Storage_GS % Value )
    associate &
      ( M    =>  CSV ( :, CS % BARYON_MASS ), &
        D    =>  CSV ( :, CS % BARYON_DENSITY_B ), &
        Phi  =>   GV ( :, G % POTENTIAL ) )

!    select case ( trim ( GA % GravitySolverType ) )
!    case ( 'UNIFORM', 'CENTRAL_MASS' )  !-- External potential
!       do iS = 1, T % nSelected
!         iI = T % iaSelected ( iS )
!         if ( iI == T % GRAVITATIONAL_ENERGY ) then
!           Integrand ( iS ) % Value  =  M * D * Phi
!         else if ( iI == T % TOTAL_ENERGY ) then
!           Integrand ( iS ) % Value  &
!             =  Integrand ( iS ) % Value  +  M * D * Phi
!         end if !-- iI
!       end do !-- iS     
!     case default

      do iS  =  1,  T % nSelected
        iI  =  T % iaSelected ( iS )
        if ( iI  ==  T % GRAVITATIONAL_ENERGY ) then
          IV ( :, iS )  =  0.5_KDR * M * D * Phi
        else if ( iI  ==  T % TOTAL_ENERGY ) then
          IV ( :, iS )  =  IV ( :, iS )  +  0.5_KDR * M * D * Phi
        end if !-- iI
      end do !-- iS     

!     end select !-- GravitySolverType

    end associate !-- M, etc.
    end associate !-- CSV, etc.
    end select !-- A
    end select !-- G
    end select !-- CS

  end subroutine ComputeInteriorIntegrand_N


  subroutine ComputeBoundaryIntegrand_G ( T, CS, C, BF )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BF

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
      iMomentum_3, &
      iKineticEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      nB    !-- nBoundary
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      X_1, X_2, X_3
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      M

    select type ( CS )
      class is ( Fluid_D_Form )
    associate &
      ( G   =>  T % Geometry, &
        I   =>  T % BoundaryIntegral % Integrand, &
        Cy  =>  C % Connectivity )
       
    do iFluence  =  1,  CS % nBalanced
      if ( CS % iaBalanced ( iFluence )  ==  CS % BARYON_DENSITY_B ) &
        iDensity = iFluence
      if ( CS % iaBalanced ( iFluence )  ==  CS % MOMENTUM_DENSITY_D ( 1 ) ) &
        iMomentum_1 = iFluence
      if ( CS % iaBalanced ( iFluence )  ==  CS % MOMENTUM_DENSITY_D ( 2 ) ) &
        iMomentum_2 = iFluence
      if ( CS % iaBalanced ( iFluence )  ==  CS % MOMENTUM_DENSITY_D ( 3 ) ) &
        iMomentum_3 = iFluence
    end do !-- iFluence

    associate ( CSV  =>  CS % Storage_GS % Value )
    call C % SetFieldPointer ( CSV ( :, CS % BARYON_MASS ), M )
    end associate !-- CSV

    do iD  =  1,  C % nDimensions

      nB = shape ( BF ( 1, Cy % iaInner ( iD ) ) % Value )

      !-- Geometry
      allocate ( X_1 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_2 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_3 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

      do iF = 1, 2

        if ( iF == 1 ) then
          iC  =  Cy % iaInner ( iD )
        else if ( iF == 2 ) then
          iC  =  Cy % iaOuter ( iD )
        end if

        associate &
          ( D    =>  BF ( iDensity,    iC ) % Value, &
            S_1  =>  BF ( iMomentum_1, iC ) % Value, &
            S_2  =>  BF ( iMomentum_2, iC ) % Value, &
            S_3  =>  BF ( iMomentum_3, iC ) % Value )

        call T % ComputeFacePositions ( G, C, iD, iF, X_1, X_2, X_3 )

        select case ( trim ( C % CoordinateSystem ) )
        case ( 'RECTANGULAR' )
          do iS  =  1,  T % nSelected
            iI = T % iaSelected ( iS )
            if ( iI  ==  T % BARYON_NUMBER ) then
              call Copy ( D, I ( iS, iC ) % Value )
            else if ( iI  ==  T % MOMENTUM ( 1 ) ) then
              call Copy ( S_1, I ( iS, iC ) % Value )
            else if ( iI  ==  T % MOMENTUM ( 2 ) ) then
              call Copy ( S_2, I ( iS, iC ) % Value )
            else if ( iI  ==  T % MOMENTUM ( 3 ) ) then
              call Copy ( S_3, I ( iS, iC ) % Value )
            else if ( iI  ==  T % KINETIC_ENERGY ) then
              iKineticEnergy  =  iS
              call ComputeFluence_KE_Rectangular &
                     ( S_1, S_2, S_3, M, D, I ( iS, iC ) % Value )
            else if ( iI  ==  T % ANGULAR_MOMENTUM ( 1 ) ) then
              call ComputeFluence_AM_Rectangular &
                     ( X_2, X_3, S_2, S_3, I ( iS, iC ) % Value )
            else if ( iI  ==  T % ANGULAR_MOMENTUM ( 2 ) ) then
              call ComputeFluence_AM_Rectangular &
                     ( X_3, X_1, S_3, S_1, I ( iS, iC ) % Value )
            else if ( iI  ==  T % ANGULAR_MOMENTUM ( 3 ) ) then
              call ComputeFluence_AM_Rectangular &
                     ( X_1, X_2, S_1, S_2, I ( iS, iC ) % Value )
            end if !-- iI
          end do !-- iS
!         case ( 'CYLINDRICAL' )
!           do iS = 1, T % nSelected
!             iI = T % iaSelected ( iS )
!             if ( iI  ==  T % BARYON_NUMBER ) then
!               call Copy ( D, I ( iS, iC ) % Value )
!             else if ( iI  ==  T % MOMENTUM ( 1 ) ) then
!               if ( CSL % nDimensions > 2 ) &
!                 call ComputeFluence_LM_CylindricalHorizontal &
!                        ( X_1, X_3, S_1, S_3, 1, I ( iS, iC ) % Value )
!             else if ( iI  ==  T % MOMENTUM ( 2 ) ) then
!               if ( CSL % nDimensions > 2 ) &
!                 call ComputeFluence_LM_CylindricalHorizontal &
!                        ( X_1, X_3, S_1, S_3, 2, I ( iS, iC ) % Value )
!             else if ( iI  ==  T % MOMENTUM ( 3 ) ) then
!               if ( CSL % nDimensions > 1 ) &
!                 call Copy ( S_2, I ( iS, iC ) % Value )
!             else if ( iI  ==  T % ANGULAR_MOMENTUM ( 1 ) ) then
!               if ( CSL % nDimensions > 2 ) &
!                 call ComputeFluence_AM_CylindricalHorizontal &
!                        ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
!                          I ( iS, iC ) % Value )
!             else if ( iI  ==  T % ANGULAR_MOMENTUM ( 2 ) ) then
!               if ( CSL % nDimensions > 2 ) &
!                 call ComputeFluence_AM_CylindricalHorizontal &
!                        ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
!                          I ( iS, iC ) % Value )
!             else if ( iI  ==  T % ANGULAR_MOMENTUM ( 3 ) ) then
!               if ( CSL % nDimensions > 1 ) &
!                 call Copy ( S_3, I ( iS, iC ) % Value )
!             end if !-- iI
!           end do !-- iS
        case ( 'SPHERICAL' )
          do iS = 1, T % nSelected
            iI  =  T % iaSelected ( iS )
            if ( iI  ==  T % BARYON_NUMBER ) then
              call Copy ( D, I ( iS, iC ) % Value )
            else if ( iI  ==  T % MOMENTUM ( 1 ) ) then
              if ( C % nDimensions > 2 ) &
                call ComputeFluence_LM_SphericalHorizontal &
                       ( X_1, X_2, X_3, S_1, S_2, S_3, 1, &
                         I ( iS, iC ) % Value )
            else if ( iI  ==  T % MOMENTUM ( 2 ) ) then
              if ( C % nDimensions > 2 ) &
                call ComputeFluence_LM_SphericalHorizontal &
                       ( X_1, X_2, X_3, S_1, S_2, S_3, 2, &
                         I ( iS, iC ) % Value )
            else if ( iI  ==  T % MOMENTUM ( 3 ) ) then
              if ( C % nDimensions > 1 ) &
                call ComputeFluence_LM_SphericalVertical &
                       ( X_1, X_2, S_1, S_2, I ( iS, iC ) % Value )
           else if ( iI  ==  T % KINETIC_ENERGY ) then
             iKineticEnergy  =  iS
             call ComputeFluence_KE_Spherical &
                    ( S_1, M, D, I ( iS, iC ) % Value )
            else if ( iI  ==  T % ANGULAR_MOMENTUM ( 1 ) ) then
              if ( C % nDimensions > 2 ) &
                call ComputeFluence_AM_SphericalHorizontal &
                       ( X_2, X_3, S_2, S_3, 1, I ( iS, iC ) % Value )
            else if ( iI  ==  T % ANGULAR_MOMENTUM ( 2 ) ) then
              if ( C % nDimensions > 2 ) &
                call ComputeFluence_AM_SphericalHorizontal &
                       ( X_2, X_3, S_2, S_3, 2, I ( iS, iC ) % Value )
            else if ( iI  ==  T % ANGULAR_MOMENTUM ( 3 ) ) then
              if ( C % nDimensions > 1 ) &
                call Copy ( S_3, I ( iS, iC ) % Value )
            end if !-- iI
          end do !-- iS
        end select !-- CoordinateSystem
        end associate !-- D, etc.

        do iS  =  1, T % nSelected
          iI  =  T % iaSelected ( iS )
          if ( iI  ==  T % TOTAL_ENERGY ) then
            call Copy ( I ( iKineticEnergy, iC ) % Value, &
                        I ( iS, iC ) % Value )
          end if !-- iI
        end do !-- iS

      end do !-- iF

      deallocate ( X_1, X_2, X_3 )

    end do !-- iD

    end associate !-- G, etc.
    end select !-- CS
       
  end subroutine ComputeBoundaryIntegrand_G


  subroutine ComputeBoundaryIntegrand_N ( T, CS, C, BF )

    class ( Tally_F_D_Form ), intent ( inout ) :: &
      T
    class ( FieldSetForm ), intent ( in ) :: &
      CS
    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    type ( Real_3D_Form ), dimension ( :, : ), intent ( in ) :: &
      BF

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

    call T % ComputeBoundaryIntegrand_G ( CS, C, BF )

    select type ( CS )
      class is ( Fluid_D_Form )
    select type ( G  =>  T % Geometry )
      class is ( Gravitation_N_H_Form )
    associate &
      ( CSV  =>  CS % Storage_GS % Value, &
         GV  =>  G  % Storage_GS % Value, &
         I   =>  T % BoundaryIntegral % Integrand, &
         Cy  =>  C % Connectivity )

    do iFluence  =  1,  CS % nBalanced
      if ( CS % iaBalanced ( iFluence )  ==  CS % BARYON_DENSITY_B ) &
        iDensity = iFluence
    end do !-- iFluence

    call C % SetFieldPointer (  GV ( :, G % POTENTIAL ), Phi )
    call C % SetFieldPointer ( CSV ( :, CS % BARYON_MASS ), M )

    DimensionLoop: do iD  =  1,  C % nDimensions

      nB = shape ( BF ( 1, Cy % iaInner ( iD ) ) % Value )

      !-- Geometry
      allocate ( X_1 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_2 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )
      allocate ( X_3 ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

      FaceLoop: do iF = 1, 2

        if ( iF == 1 ) then
          iC  =  Cy % iaInner ( iD )
        else if ( iF == 2 ) then
          iC  =  Cy % iaOuter ( iD )
        end if

        associate ( D  =>  BF ( iDensity, iC ) % Value )

        call T % ComputeFacePositions ( G, C, iD, iF, X_1, X_2, X_3 )

        select case ( trim ( C % CoordinateSystem ) )
        case ( 'SPHERICAL' )

          if ( iD /= 1 ) &
            exit DimensionLoop
!          if ( iF /= 2 ) &
!            cycle FaceLoop

          call C % SetFieldPointer ( GV ( :, G % CENTER_U ( 1 ) ), R )

          do iS  =  1,  T % nSelected
            iI  =  T % iaSelected ( iS )
            if ( iI  ==  T % GRAVITATIONAL_ENERGY ) then
              call ComputeFluence_GE_Spherical &
                     ( C, M, Phi, R, D, X_1, I ( iS, iC ) % Value )
              iGravity = iS
            end if !-- iI
          end do !-- iS

          do iS  =  1,  T % nSelected
            iI  =  T % iaSelected ( iS )
            if ( iI  ==  T % TOTAL_ENERGY ) then
              I ( iS, iF ) % Value  &
                =  I ( iS, iF ) % Value  +  I ( iGravity, iF ) % Value
            end if !-- iI
          end do !-- iS

        end select !-- CoordinateSystem

        end associate !-- D

      end do FaceLoop !-- iF

      deallocate ( X_1, X_2, X_3 )

    end do DimensionLoop !-- iD

    end associate !-- CSV, etc.
    end select !-- G
    end select !-- CS

  end subroutine ComputeBoundaryIntegrand_N


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

    !$OMP parallel do
    do iV  =  1,  nV
      I ( iV )  =  0.5_KDR  *  (    S_1 ( iV ) * V_1 ( iV )  &
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

    !$OMP parallel do
    do iV = 1, nV
      I_I ( iV )  =  X_J ( iV ) * S_K ( iV )  -  X_K ( iV ) * S_J ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeDensity_AM_Rectangular


!   subroutine ComputeDensity_LM_CylindricalHorizontal &
!                ( X_1, X_3, S_1, S_3, iD, I_I )

!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       X_1, X_3, &
!       S_1, S_3
!     integer ( KDI ), intent ( in ) :: &
!       iD
!     real ( KDR ), dimension ( : ), intent ( out ) :: &
!       I_I

!     integer ( KDI ) :: &
!       iV, &
!       nV
!     real ( KDR ) :: &
!       SqrtTiny

!     nV = size ( I_I )
    
!     SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

!     select case ( iD )
!     case ( 1 )
!       !$OMP parallel do private ( iV )
!       do iV = 1, nV
!         I_I ( iV )  =     cos ( X_3 ( iV ) )  *  S_1 ( iV )  &
!                        -  sin ( X_3 ( iV ) )  &
!                           /  max ( X_1 ( iV ), SqrtTiny )  &
!                           *  S_3 ( iV ) 
!       end do
!       !$OMP end parallel do
!     case ( 2 )
!       !$OMP parallel do private ( iV )
!       do iV = 1, nV
!         I_I ( iV )  =     sin ( X_3 ( iV ) )  *  S_1 ( iV )  &
!                        +  cos ( X_3 ( iV ) )  &
!                           /  max ( X_1 ( iV ), SqrtTiny )  &
!                           *  S_3 ( iV ) 
!       end do
!       !$OMP end parallel do
!     end select

!   end subroutine ComputeDensity_LM_CylindricalHorizontal


!   subroutine ComputeDensity_AM_CylindricalHorizontal &
!                ( X_1, X_2, X_3, S_1, S_2, S_3, iD, I_I )

!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       X_1, X_2, X_3, &
!       S_1, S_2, S_3
!     integer ( KDI ), intent ( in ) :: &
!       iD
!     real ( KDR ), dimension ( : ), intent ( out ) :: &
!       I_I

!     integer ( KDI ) :: &
!       iV, &
!       nV
!     real ( KDR ) :: &
!       SqrtTiny

!     nV = size ( I_I )

!     SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

!     select case ( iD )
!     case ( 1 )
!       !$OMP parallel do private ( iV )
!       do iV = 1, nV
!         I_I ( iV )  =  sin ( X_3 ( iV ) )  *  (   X_1 ( iV )  *  S_2 ( iV )  &
!                                                 - X_2 ( iV )  *  S_1 ( iV ) ) &
!                        - X_2 ( iV )  *  cos ( X_3 ( iV ) )  &
!                          /  max ( X_1 ( iV ), SqrtTiny )  &
!                          *  S_3 ( iV )                         
!       end do
!       !$OMP end parallel do
!     case ( 2 )
!       !$OMP parallel do private ( iV )
!       do iV = 1, nV
!         I_I ( iV )  =  cos ( X_3 ( iV ) )  *  (   X_2 ( iV )  *  S_1 ( iV )  &
!                                                 - X_1 ( iV )  *  S_2 ( iV ) ) &
!                        - X_2 ( iV )  *  sin ( X_3 ( iV ) )  &
!                          /  max ( X_1 ( iV ), SqrtTiny )  &
!                          *  S_3 ( iV )                         
!       end do
!       !$OMP end parallel do
!     end select

!   end subroutine ComputeDensity_AM_CylindricalHorizontal


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
      !$OMP parallel do
      do iV  =  1,  nV
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


  subroutine ComputeFluence_KE_Rectangular ( S_1, S_2, S_3, M, D, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      S_1, S_2, S_3, &
      M, &
      D
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( I_I )

    !$OMP parallel do collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          if ( D ( iV, jV, kV )  >  0.0_KDR ) &
            I_I ( iV, jV, kV ) &
              =  0.5_KDR  *  (    S_1 ( iV, jV, kV ) * S_1 ( iV, jV, kV ) &
                               +  S_2 ( iV, jV, kV ) * S_2 ( iV, jV, kV ) &
                               +  S_2 ( iV, jV, kV ) * S_2 ( iV, jV, kV ) ) &
                          /  ( D ( iV, jV, kV ) )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeFluence_KE_Rectangular


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

    !$OMP parallel do collapse ( 3 )
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


!   subroutine ComputeFluence_LM_CylindricalHorizontal &
!                ( X_1, X_3, S_1, S_3, iD, I_I )

!     real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
!       X_1, X_3, &
!       S_1, S_3
!     integer ( KDI ), intent ( in ) :: &
!       iD
!     real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
!       I_I

!     integer ( KDI ) :: &
!       iV, jV, kV
!     integer ( KDI ), dimension ( 3 ) :: &
!       nV
!     real ( KDR ) :: &
!       SqrtTiny

!     nV = shape ( I_I )

!     SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

!     select case ( iD )
!     case ( 1 )
!       !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
!       do kV = 1, nV ( 3 )
!         do jV = 1, nV ( 2 )
!           do iV = 1, nV ( 1 )
!             I_I ( iV, jV, kV )  &
!               =  cos ( X_3 ( iV, jV, kV ) )  *  S_1 ( iV, jV, kV )  &
!                  -  sin ( X_3 ( iV, jV, kV ) )  &
!                     /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
!                     *  S_3 ( iV, jV, kV ) 
!           end do
!         end do
!       end do
!       !$OMP end parallel do
!     case ( 2 )
!       !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
!       do kV = 1, nV ( 3 )
!         do jV = 1, nV ( 2 )
!           do iV = 1, nV ( 1 )
!             I_I ( iV, jV, kV )  &
!               =  sin ( X_3 ( iV, jV, kV ) )  *  S_1 ( iV, jV, kV )  &
!                  +  cos ( X_3 ( iV, jV, kV ) )  &
!                     /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
!                     *  S_3 ( iV, jV, kV ) 
!           end do
!         end do
!       end do
!       !$OMP end parallel do
!     end select

!   end subroutine ComputeFluence_LM_CylindricalHorizontal


!   subroutine ComputeFluence_AM_CylindricalHorizontal &
!                ( X_1, X_2, X_3, S_1, S_2, S_3, iD, I_I )

!     real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
!       X_1, X_2, X_3, &
!       S_1, S_2, S_3
!     integer ( KDI ), intent ( in ) :: &
!       iD
!     real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
!       I_I

!     integer ( KDI ) :: &
!       iV, jV, kV
!     integer ( KDI ), dimension ( 3 ) :: &
!       nV
!     real ( KDR ) :: &
!       SqrtTiny

!     nV = shape ( I_I )

!     SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

!     select case ( iD )
!     case ( 1 )
!       !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
!       do kV = 1, nV ( 3 )
!         do jV = 1, nV ( 2 )
!           do iV = 1, nV ( 1 )
!             I_I ( iV, jV, kV )  &
!               =  sin ( X_3 ( iV, jV, kV ) )  &
!                    *  (   X_1 ( iV, jV, kV )  *  S_2 ( iV, jV, kV )  &
!                         - X_2 ( iV, jV, kV )  *  S_1 ( iV, jV, kV ) ) &
!                  -  X_2 ( iV, jV, kV )  *  cos ( X_3 ( iV, jV, kV ) )  &
!                     /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
!                     *  S_3 ( iV, jV, kV )                         
!           end do
!         end do
!       end do
!       !$OMP end parallel do
!     case ( 2 )
!       !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 )
!       do kV = 1, nV ( 3 )
!         do jV = 1, nV ( 2 )
!           do iV = 1, nV ( 1 )
!             I_I ( iV, jV, kV )  &
!               =  cos ( X_3 ( iV, jV, kV ) )  &
!                    *  (   X_2 ( iV, jV, kV )  *  S_1 ( iV, jV, kV )  &
!                         - X_1 ( iV, jV, kV )  *  S_2 ( iV, jV, kV ) ) &
!                  -  X_2 ( iV, jV, kV )  *  sin ( X_3 ( iV, jV, kV ) )  &
!                     /  max ( X_1 ( iV, jV, kV ), SqrtTiny )  &
!                     *  S_3 ( iV, jV, kV )                         
!           end do
!         end do
!       end do
!       !$OMP end parallel do
!     end select

!   end subroutine ComputeFluence_AM_CylindricalHorizontal


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

  
  subroutine ComputeFluence_KE_Spherical ( S_1, M, D, I_I )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      S_1, &
      M, &
      D
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      I_I

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( I_I )

    !$OMP parallel do collapse ( 3 )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          if ( D ( iV, jV, kV )  >  0.0_KDR ) &
            I_I ( iV, jV, kV ) &
              =  0.5_KDR  *  ( S_1 ( iV, jV, kV ) * S_1 ( iV, jV, kV ) ) &
                          /  ( D ( iV, jV, kV ) )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ComputeFluence_KE_Spherical


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


  subroutine ComputeFluence_GE_Spherical ( C, M, Phi, R, D, R_I, I_I )

    class ( Chart_GS_Form ), intent ( in ) :: &
      C    
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

    oI  =  C % nGhostLayers ( 1 )  +  C % nCellsBrick ( 1 )  -  1
    oJ  =  C % nGhostLayers ( 2 )
    oK  =  C % nGhostLayers ( 3 )

    !$OMP parallel do collapse ( 2 )
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
