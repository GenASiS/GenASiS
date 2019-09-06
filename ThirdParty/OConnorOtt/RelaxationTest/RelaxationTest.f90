PROGRAM RelaxationTest

  USE nulibtable, ONLY: &
    nulibtable_nrho, &
    nulibtable_ntemp, &
    nulibtable_nye, &
    nulibtable_number_groups, &
    nulibtable_number_species, &
    nulibtable_logrho, &
    nulibtable_logtemp, &
    nulibtable_ye, &
    nulibtable_number_easvariables, &
    nulibtable_energies, &
    nulibtable_ewidths, &
    nulibtable_ebottom, &
    nulibtable_etop
  USE nulibtable_interface, ONLY: &
    nulibtable_reader

  USE HDF5

  IMPLICIT NONE

  ! --- Program Parameters ---

  LOGICAL,  PARAMETER :: UpdateMatter = .TRUE.
  LOGICAL,  PARAMETER :: Include_EC   = .TRUE.
  LOGICAL,  PARAMETER :: Include_NES  = .TRUE.
  LOGICAL,  PARAMETER :: Include_Pair = .FALSE.
  INTEGER,  PARAMETER :: dp = KIND( 1.d0 )
  INTEGER,  PARAMETER :: nx = 213
  INTEGER,  PARAMETER :: keyeint = 0
  INTEGER,  PARAMETER :: keytemp = 1
  REAL(dp), PARAMETER :: Log1d100 = LOG( 1.0d100 )
  REAL(dp), PARAMETER :: FourPi  = 4.0_DP * ACOS( - 1.0_DP )
  REAL(dp), PARAMETER :: mBcgs   = 1.6605402d-24 ! [ g ]
  REAL(dp), PARAMETER :: kMeV    = 8.61733d-11   ! [ MeV K^{-1} ]
  REAL(dp), PARAMETER :: ccgs    = 2.997924d+10  ! [ cm s^{-1} ]
  REAL(dp), PARAMETER :: hMeV    = 4.135667d-21  ! [ MeV s ]
  REAL(dp), PARAMETER :: MeV2Erg = 1.602d-6      ! [ Erg MeV^{-1} ]

  ! --- Timers ---

  REAL(dp) :: Timer_Total
  REAL(dp) :: Timer_IO
  REAL(dp) :: Timer_Eos
  REAL(dp) :: Timer_Opacity
  REAL(dp) :: Timer_FP_AA

  ! --- NuLib Arrays ---

  REAL(dp), ALLOCATABLE :: eas_energy(:,:)
  REAL(dp), ALLOCATABLE :: ies_energy_energy(:,:,:)
  REAL(dp), ALLOCATABLE :: epannihil_energy_energy(:,:,:)

  ! --- Program Variables ---

  CHARACTER(16) :: InitialProfile
  INTEGER       :: ix, iS, iE, nE
  INTEGER       :: iCycle, iCycleD, iCycleW
  INTEGER       :: keyerr
  INTEGER       :: FileNumber = 0
  INTEGER       :: nIterations(nx) = 0
  REAL(dp)      :: Exponent
  REAL(dp)      :: G_Shift, FD_Shift
  REAL(dp)      :: dt, time, t_end, temp_in
  REAL(dp)      :: xR(nx), xD(nx), xT(nx), xY(nx), xEta(nx)
  REAL(dp)      :: xEnr(nx), xPrs(nx), xEnt(nx), xCs2(nx), xdEdT(nx)
  REAL(dp)      :: xdPdE_D(nx), xdPdD_E(nx), xXa(nx), xXh(nx), xXn(nx)
  REAL(dp)      :: xXp(nx), xAbar(nx), xZbar(nx), xMu_e(nx), xMu_n(nx)
  REAL(dp)      :: xMu_p(nx), xMu_hat(nx), xMu_nu(2,nx), xQ_Y(nx), xQ_Enr(nx)
  REAL(dp), ALLOCATABLE, DIMENSION(:)     :: E, dE, dV_E
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: Chi_EC
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: Chi
  REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: Eta
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Phi_0_In__NES
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Phi_0_Out_NES
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Phi_0_In__Pair
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Phi_0_Out_Pair
  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: xJ, xJnew

  CALL InitializeTimers

  CALL Timer_Start( Timer_Total )

  InitialProfile = 'FermiDirac'

  G_Shift  = 1.0_DP ! --- Shift Factor for Gaussian    Initial Profile
  FD_Shift = 0.8_DP ! --- Shift Factor for Fermi-Dirac Initial Profile

  dt      = 1.0d-6
  time    = 0.0d0
  t_end   = 1.0d3 * dt
  iCycleD = 1
  iCycleW = 10

  WRITE(*,*)
  WRITE(*,'(A4,A)') '', 'RelaxationTest'
  WRITE(*,*)

  CALL Timer_Start( Timer_IO )

  CALL ReadMatterProfiles( nx, xR, xD, xT, xY )

  CALL Timer_Stop( Timer_IO )

  ! --- Read EOS Table ---

  CALL Timer_Start( Timer_IO )

  CALL read_eos_table( "NuLibEosTable.h5" )

  CALL Timer_Stop( Timer_IO )

  ! --- Read Opacity Table ---

  CALL Timer_Start( Timer_IO )

  CALL nulibtable_reader &
         ( filename &
             = 'NuLibOpacityTable.h5', &
           include_Ielectron &
             = Include_NES, &
           include_epannihil_kernels &
             = Include_Pair )

  CALL Timer_Stop( Timer_IO )

  CALL AllocateRelaxation

  ! --- Initialize Radiation Field ---

  DO ix = 1, nx

    ! --- Evaluate Eos (with D,T,Y) ---

    CALL Timer_Start( Timer_Eos )

    temp_in = kMeV * xT(ix); xEnr(ix) = 0.0_DP

    CALL nuc_eos_full &
           ( xD(ix), temp_in, xY(ix), xEnr(ix),      &
             xPrs(ix), xEnt(ix), xCs2(ix), xdEdT(ix),      &
             xdPdE_D(ix), xdPdD_E(ix), xXa(ix), xXh(ix),   &
             xXn(ix), xXp(ix), xAbar(ix), xZbar(ix),       &
             xMu_e(ix), xMu_n(ix), xMu_p(ix), xMu_hat(ix), &
             keytemp, keyerr, 1.0d-9 )

    xEta(ix)     = xMu_e(ix) / ( kMeV * xT(ix) )
    xMu_nu(1,ix) = xMu_e(ix) + xMu_p(ix) - xMu_n(ix)
    xMu_nu(2,ix) = - xMu_nu(1,ix)

    CALL Timer_Stop( Timer_Eos )

    DO iS = 1, 2

      IF( TRIM( InitialProfile ) == 'Gaussian' )THEN

        DO iE = 1, nE

          Exponent = MIN( MAX( - ( E(iE) - G_Shift * kMeV * xT(ix) )**2 / ( 2.0d2 ), - Log1d100 ), Log1d100 )

          xJ(iE,iS,ix) = 0.99_dp * E(iE) * EXP( Exponent )

        END DO

      ELSEIF( TRIM( InitialProfile ) == 'FermiDirac' )THEN

        xJ(:,iS,ix) = E / ( EXP( ( E - FD_Shift * xMu_nu(iS,ix) ) / ( kMeV * xT(ix) ) ) + 1.0_DP )

      ELSE

        WRITE(*,*)
        WRITE(*,'(A4,A,A)') '', 'Invalid Initial Profile: ', TRIM( InitialProfile )
        WRITE(*,*)
        STOP

      END IF

    END DO

  END DO

  xQ_Y   = 0.0_dp
  xQ_Enr = 0.0_dp

  CALL Timer_Start( Timer_IO )

  CALL WriteOutput ! --- Initial Condition

  CALL Timer_Stop( Timer_IO )

  iCycle = 0
  DO WHILE( time < t_end )

    iCycle = iCycle + 1

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A4,A8,I8.8,A8,ES9.3E2,A6,ES9.3E2)') &
        '', 'cycle = ', iCycle, 'time = ', time, 'dt = ', dt

    END IF

    DO ix = 1, nx

      ! --- Compute Opacities ---

      CALL Timer_Start( Timer_Opacity )

      DO iS = 1, 2

        ! --- Emission/Absorption/Scattering ---

        CALL nulibtable_single_species_range_energy &
               ( xD(ix), kMeV * xT(ix), xY(ix), iS, eas_energy, &
                 nulibtable_number_groups, nulibtable_number_easvariables )

        Chi_EC(:,iS) = eas_energy(:,2)

        IF( Include_NES )THEN

          ! --- Neutrino-Electron Scattering ---

          CALL nulibtable_inelastic_single_species_range_energy &
                 ( kMeV * xT(ix), xEta(ix), iS, ies_energy_energy, &
                   nulibtable_number_groups, nulibtable_number_groups, 2 )

          Phi_0_In__NES(:,:,iS) = TRANSPOSE( ies_energy_energy(:,:,1) )
          Phi_0_Out_NES(:,:,iS) =            ies_energy_energy(:,:,1)

        ELSE

          Phi_0_In__NES(:,:,iS) = 0.0_DP
          Phi_0_Out_NES(:,:,iS) = 0.0_DP

        END IF

        IF( Include_Pair )THEN

          ! --- Pair Creation and Annihilation ---

          CALL nulibtable_epannihil_single_species_range_energy &
                 ( kMeV * xT(ix), xEta(ix), iS, epannihil_energy_energy, &
                   nulibtable_number_groups, nulibtable_number_groups, 4 )

          Phi_0_In__Pair(:,:,iS) = epannihil_energy_energy(:,:,1)
          Phi_0_Out_Pair(:,:,iS) = epannihil_energy_energy(:,:,2)

        ELSE

          Phi_0_In__Pair(:,:,iS) = 0.0_DP
          Phi_0_Out_Pair(:,:,iS) = 0.0_DP

        END IF

      END DO

      CALL Timer_Stop( Timer_Opacity )

      CALL Timer_Start( Timer_FP_AA )

      CALL Solve_FP_AA &
             ( dt, xD(ix), kMeV * xT(ix), xY(ix), xMu_nu(1:2,ix), &
               E(1:nE), Chi_EC, xJ(1:nE,1:2,ix), xJnew(1:nE,1:2,ix), &
               xQ_Y(ix), xQ_Enr(ix), nIterations(ix) )

!!$      CALL SolveEmAb_FP_AA & ! --- Fully Implicit Em-Ab Solver
!!$             ( dt, xD(ix), kMeV * xT(ix), xY(ix), xEnr(ix), &
!!$               xMu_nu(1:2,ix), E(1:nE), Chi_EC, xJ(1:nE,1:2,ix), &
!!$               xJnew(1:nE,1:2,ix), xQ_Y(ix), xQ_Enr(ix), nIterations(ix) )

      xJ(1:nE,1:2,ix) = xJnew(1:nE,1:2,ix)

      IF( UpdateMatter )THEN

        xY  (ix) = xY  (ix) - xQ_Y  (ix)
        xEnr(ix) = xEnr(ix) - xQ_Enr(ix)

      END IF

      CALL Timer_Stop( Timer_FP_AA )

      ! --- Evaluate Eos (with D,E,Y) ---

      CALL Timer_Start( Timer_Eos )

      temp_in = kMeV * xT(ix) ! K -> MeV

      CALL nuc_eos_full &
             ( xD(ix), temp_in, xY(ix), xEnr(ix),      &
               xPrs(ix), xEnt(ix), xCs2(ix), xdEdT(ix),      &
               xdPdE_D(ix), xdPdD_E(ix), xXa(ix), xXh(ix),   &
               xXn(ix), xXp(ix), xAbar(ix), xZbar(ix),       &
               xMu_e(ix), xMu_n(ix), xMu_p(ix), xMu_hat(ix), &
               keyeint, keyerr, 1.0d-9 )

      xT(ix) = temp_in / kMeV ! MeV -> K

      xEta(ix)     = xMu_e(ix) / ( kMeV * xT(ix) )
      xMu_nu(1,ix) = xMu_e(ix) + xMu_p(ix) - xMu_n(ix)
      xMu_nu(2,ix) = - xMu_nu(1,ix)

      CALL Timer_Stop( Timer_Eos )

    END DO

    time = time + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL Timer_Start( Timer_IO )

      CALL WriteOutput

      CALL Timer_Stop( Timer_IO )

    END IF

  END DO

  CALL Timer_Start( Timer_IO )

  CALL WriteOutput

  CALL Timer_Stop( Timer_IO )

  CALL DeallocateRelaxation

  CALL Timer_Stop( Timer_Total )

  CALL FinalizeTimers

CONTAINS


  SUBROUTINE ReadMatterProfiles( nx, xR, xD, xT, xY )

    INTEGER,  INTENT(in)  :: nx
    REAL(dp), INTENT(out) :: xR(:), xD(:), xT(:), xY(:)

    CHARACTER(04) :: Format1 = "(I5)"
    CHARACTER(09) :: Format2 = "(4ES12.3)"
    INTEGER       :: ix
    REAL(dp)      :: buffer(4*nx)

    OPEN( 1, FILE = "MatterProfiles.d", FORM = "formatted", ACTION = 'read' )

    READ( 1, Format2 ) buffer

    DO ix = 1, nx

      xR(nx-ix+1) = buffer((ix-1)*4+1)
      xD(nx-ix+1) = buffer((ix-1)*4+2)
      xT(nx-ix+1) = buffer((ix-1)*4+3)
      xY(nx-ix+1) = buffer((ix-1)*4+4)

    END DO

    CLOSE( 1, STATUS = 'keep' )

  END SUBROUTINE ReadMatterProfiles


  SUBROUTINE AllocateRelaxation

    nE = nulibtable_number_groups

    ALLOCATE( E(nE), dE(nE), dV_E(nE) )

    E    = nulibtable_energies
    dE   = nulibtable_ewidths
    dV_E = ( nulibtable_etop**3 - nulibtable_ebottom**3 ) / 3.0_DP

    ALLOCATE( Chi_EC(nE,2) )
    ALLOCATE( Chi   (nE,2) )
    ALLOCATE( Eta   (nE,2) )

    ALLOCATE( Phi_0_In__NES(nE,nE,2) )
    ALLOCATE( Phi_0_Out_NES(nE,nE,2) )

    ALLOCATE( Phi_0_In__Pair(nE,nE,2) )
    ALLOCATE( Phi_0_Out_Pair(nE,nE,2) )

    ALLOCATE( xJ   (nE,2,nx) )
    ALLOCATE( xJnew(nE,2,nx) )

    ! --- NuLib Arrays ---

    ALLOCATE &
      ( eas_energy &
          (nulibtable_number_groups,nulibtable_number_easvariables) )
    ALLOCATE &
      ( ies_energy_energy &
          (nulibtable_number_groups,nulibtable_number_groups,2) )
    ALLOCATE &
      ( epannihil_energy_energy &
          (nulibtable_number_groups,nulibtable_number_groups,4) )

  END SUBROUTINE AllocateRelaxation


  SUBROUTINE DeallocateRelaxation

    DEALLOCATE( E, dE, dV_E )

    DEALLOCATE( Chi_EC, Chi, Eta )

    DEALLOCATE( Phi_0_In__NES )
    DEALLOCATE( Phi_0_Out_NES )

    DEALLOCATE( Phi_0_In__Pair )
    DEALLOCATE( Phi_0_Out_Pair )

    DEALLOCATE( xJ, xJnew )

    ! --- NuLib Arrays ---

    DEALLOCATE( eas_energy )
    DEALLOCATE( ies_energy_energy )
    DEALLOCATE( epannihil_energy_energy )

  END SUBROUTINE DeallocateRelaxation


  SUBROUTINE Solve_FP_AA( dt, D, T, Y, Mnu, E, Chi_EC, Jold, Jnew, Q_Y, Q_Enr, nIter )

    REAL(DP), INTENT(in)  :: dt, D, T, Y, Mnu(2)
    REAL(DP), INTENT(in)  :: E     (1:nE)
    REAL(DP), INTENT(in)  :: Chi_EC(1:nE,1:2)
    REAL(DP), INTENT(in)  :: Jold  (1:nE,1:2)
    REAL(DP), INTENT(out) :: Jnew  (1:nE,1:2)
    REAL(DP), INTENT(out) :: Q_Y, Q_Enr
    INTEGER,  INTENT(out) :: nIter

    INTEGER,  PARAMETER :: M = 3
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIter = 100
    REAL(DP), PARAMETER :: Rtol = 1.0d-8

    LOGICAL  :: CONVERGED
    INTEGER  :: OS_1, OS_2, iS, jS, k, mk, i, INFO
    REAL(DP) :: C_EC, C_NES, C_Pair
    REAL(DP) :: W(1:nE), Alpha(M), WORK(LWORK)
    REAL(DP) :: J0      (1:nE,1:2), Eta_EC  (1:nE,1:2)
    REAL(DP) :: Chi_NES (1:nE,1:2), Eta_NES (1:nE,1:2) ! --- NES   Opacities
    REAL(DP) :: Chi_Pair(1:nE,1:2), Eta_Pair(1:nE,1:2) ! --- Pair  Opacities
    REAL(DP) :: Chi     (1:nE,1:2), Eta     (1:nE,1:2) ! --- Total Opacities
    REAL(DP) :: GVEC(2*nE,M), GVECm(2*nE)
    REAL(DP) :: FVEC(2*nE,M), FVECm(2*nE)
    REAL(DP) :: BVEC(2*nE), AMAT(2*nE,M)

    IF( Include_EC )THEN
      C_EC = 1.0_DP
    ELSE
      C_EC = 0.0_DP
    END IF

    IF( Include_NES )THEN
      C_NES = 1.0_DP
    ELSE
      C_NES = 0.0_DP
    END IF

    IF( Include_Pair )THEN
      C_Pair = 1.0_DP
    ELSE
      C_Pair = 0.0_DP
    END IF

    OS_1 = 0
    OS_2 = nE

    ! --- Energy Integration Weights ---

    W = FourPi / ( ccgs * ( ccgs * hMeV )**3 ) * dV_E / E

    DO iS = 1, 2

      ! --- Equilibrium Spectral Distribution ---

      J0(:,iS) = E / ( EXP( ( E - Mnu(iS) ) / T ) + 1.0_DP )

      ! --- Electron Capture Emissivity ---

      Eta_EC(:,iS) = Chi_EC(:,iS) * J0(:,iS)

    END DO

    ! --- Initial Guess ---

    Jnew = Jold

    ! --- Start Fixed-Point Iteration ---

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k = k + 1

      mk = MIN( k, M )

      DO iS = 1, 2

        jS = (1+iS)-2*(iS-1) ! --- Anitineutrino Index

        ! --- NES Emissivity and Opacity ---

        CALL DGEMV( 'N', nE, nE, 1.0_DP, Phi_0_In__NES(:,:,iS), nE, &
                    W * Jnew(:,iS),     1, 0.0_DP, Eta_NES(:,iS), 1 )

        Chi_NES(:,iS) = Eta_NES(:,iS) ! --- First Part of Chi_NES

        CALL DGEMV( 'N', nE, nE, 1.0_DP, Phi_0_Out_NES(:,:,iS), nE, &
                    W * (E-Jnew(:,iS)), 1, 1.0_DP, Chi_NES(:,iS), 1 )

        ! --- Pair Emissivity and Opacity ---

        CALL DGEMV( 'N', nE, nE, 1.0_DP, Phi_0_In__Pair(:,:,iS), nE, &
                    W * (E-Jnew(:,jS)), 1, 0.0_DP, Eta_Pair(:,iS), 1 )

        Chi_Pair(:,iS) = Eta_Pair(:,iS) ! --- First Part of Chi_Pair

        CALL DGEMV( 'N', nE, nE, 1.0_DP, Phi_0_Out_Pair(:,:,iS), nE, &
                    W * Jnew(:,jS),     1, 1.0_DP, Chi_Pair(:,iS), 1 )

        ! --- Total Emissivity and Opacity ---

        Eta(:,iS) = C_EC * Eta_EC(:,iS) + E * ( C_NES * Eta_NES(:,iS) + C_Pair * Eta_Pair(:,iS) )
        Chi(:,iS) = C_EC * Chi_EC(:,iS) +     ( C_NES * Chi_NES(:,iS) + C_Pair * Chi_Pair(:,iS) )

      END DO

      ! --- Right-Hand Side Vectors ---

      ! --- Neutrinos ---

      GVEC(OS_1+1:OS_1+nE,mk) &
        = ( Jold(:,1) + ccgs * dt * Eta(:,1) ) &
           / ( 1.0_DP + ccgs * dt * Chi(:,1) )

      ! --- Antineutrinos ---

      GVEC(OS_2+1:OS_2+nE,mk) &
        = ( Jold(:,2) + ccgs * dt * Eta(:,2) ) &
           / ( 1.0_DP + ccgs * dt * Chi(:,2) )

      ! --- Residual Vectors ---

      ! --- Neutrinos ---

      FVEC(OS_1+1:OS_1+nE,mk) = GVEC(OS_1+1:OS_1+nE,mk) - Jnew(:,1)

      ! --- Antineutrinos ---

      FVEC(OS_2+1:OS_2+nE,mk) = GVEC(OS_2+1:OS_2+nE,mk) - Jnew(:,2)

      IF( mk == 1 )THEN

        GVECm = GVEC(:,mk)

      ELSE

        BVEC(:) &
          = - FVEC(:,mk)
        AMAT(:,1:mk-1) &
          = FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

        CALL DGELS( 'N', 2*nE, mk-1, 1, AMAT(:,1:mk-1), 2*nE, BVEC, 2*nE, WORK, LWORK, INFO )

        Alpha(1:mk-1) = BVEC(1:mk-1)
        Alpha(mk) = 1.0_DP - SUM( Alpha(1:mk-1) )

        GVECm = 0.0_DP
        DO i = 1, mk

          GVECm = GVECm + Alpha(i) * GVEC(:,i)

        END DO

      END IF

      FVECm(OS_1+1:OS_1+nE) = GVECm(OS_1+1:OS_1+nE) - Jnew(:,1)
      FVECm(OS_2+1:OS_2+nE) = GVECm(OS_2+1:OS_2+nE) - Jnew(:,2)

      IF( ENORM( FVECm(OS_1+1:OS_1+nE) ) <= Rtol * ENORM( Jold(:,1) ) .AND. &
          ENORM( FVECm(OS_2+1:OS_2+nE) ) <= Rtol * ENORM( Jold(:,2) ) ) &
      THEN

        CONVERGED = .TRUE.

        nIter = k

      END IF

      Jnew(:,1) = GVECm(OS_1+1:OS_1+nE)
      Jnew(:,2) = GVECm(OS_2+1:OS_2+nE)

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

      END IF

    END DO

    ! --- Compute Matter Sources ---

    Q_Y   = ( FourPi * mBcgs )   / ( D * ( ccgs * hMeV )**3 ) &
              * ( SUM  ( ( Jnew(:,1) - Jold(:,1) ) * dV_E / E ) &
                  - SUM( ( Jnew(:,2) - Jold(:,2) ) * dV_E / E ) )

    Q_Enr = ( FourPi * MeV2Erg ) / ( D * ( ccgs * hMeV )**3 ) &
              * ( SUM  ( ( Jnew(:,1) - Jold(:,1) ) * dV_E     ) &
                  + SUM( ( Jnew(:,2) - Jold(:,2) ) * dV_E     ) )

  END SUBROUTINE Solve_FP_AA


  SUBROUTINE SolveEmAb_FP_AA( dt, D, T, Y, U, Mnu, E, Chi_EC, Jold, Jnew, Q_Y, Q_Enr, nIter )

    REAL(DP), INTENT(in)  :: dt, D, T, Y, U, Mnu(2)
    REAL(DP), INTENT(in)  :: E     (1:nE)
    REAL(DP), INTENT(in)  :: Chi_EC(1:nE,1:2)
    REAL(DP), INTENT(in)  :: Jold  (1:nE,1:2)
    REAL(DP), INTENT(out) :: Jnew  (1:nE,1:2)
    REAL(DP), INTENT(out) :: Q_Y, Q_Enr
    INTEGER,  INTENT(out) :: nIter

    INTEGER,  PARAMETER :: M = 2
    INTEGER,  PARAMETER :: LWORK = 2 * M
    INTEGER,  PARAMETER :: MaxIter = 100
    INTEGER,  PARAMETER :: iY = 1
    INTEGER,  PARAMETER :: iU = 2
    REAL(DP), PARAMETER :: Rtol = 1.0d-8

    LOGICAL  :: CONVERGED
    INTEGER  :: i, k, mk, iS, INFO
    REAL(DP) :: Alpha(M), WORK(LWORK)
    REAL(DP) :: Dold, Told, Yold, Uold, Mold(2)
    REAL(DP) :: Dnew, Tnew, Ynew, Unew, Mnew(2)
    REAL(DP) :: zPrs, zEnt, zCs2, zdEdT, zdPdE_D, zdPdD_E, zXa, zXh, zXn, zXp
    REAL(DP) :: zAbar, zZbar, zMu_e, zMu_n, zMu_p, zMu_hat
    REAL(DP) :: GVEC(2,M), GVECm(2)
    REAL(DP) :: FVEC(2,M), FVECm(2)
    REAL(DP) :: J0(1:nE,1:2)
    REAL(DP) :: BVEC(2), AMAT(2,M)

    Dold = D;   Dnew = Dold
    Told = T;   Tnew = Told
    Yold = Y;   Ynew = Yold
    Uold = U;   Unew = Uold
    Mold = Mnu; Mnew = Mold

!!$    CALL WriteVector( nE, E(:), 'E.dat' )
!!$    CALL WriteVector( nE, Chi_EC(:,1), 'Chi_1.dat' )
!!$    CALL WriteVector( nE, Chi_EC(:,2), 'Chi_2.dat' )

    k = 0
    CONVERGED = .FALSE.
    DO WHILE( .NOT. CONVERGED .AND. k < MaxIter )

      k = k + 1

      mk = MIN( k, M )

      DO iS = 1, 2

        J0(:,iS) = E / ( EXP( ( E - Mnew(iS) ) / Tnew ) + 1.0_DP )

      END DO

      GVEC(iY,mk) &
        = Yold - ( FourPi * mBcgs )   / ( D * ( ccgs * hMeV )**3 ) &
            * ( SUM  ( ccgs*dt*Chi_EC(:,1) * ( J0(:,1) - Jold(:,1) ) / ( 1.0_DP + ccgs*dt*Chi_EC(:,1) ) * dV_E/E ) &
                - SUM( ccgs*dt*Chi_EC(:,2) * ( J0(:,2) - Jold(:,2) ) / ( 1.0_DP + ccgs*dt*Chi_EC(:,2) ) * dV_E/E ) )

      FVEC(iY,mk) = GVEC(iY,mk) - Ynew

      GVEC(iU,mk) &
        = Uold - ( FourPi * MeV2Erg ) / ( D * ( ccgs * hMeV )**3 ) &
            * ( SUM  ( ccgs*dt*Chi_EC(:,1) * ( J0(:,1) - Jold(:,1) ) / ( 1.0_DP + ccgs*dt*Chi_EC(:,1) ) * dV_E ) &
                + SUM( ccgs*dt*Chi_EC(:,2) * ( J0(:,2) - Jold(:,2) ) / ( 1.0_DP + ccgs*dt*Chi_EC(:,2) ) * dV_E ) )

      FVEC(iU,mk) = GVEC(iU,mk) - Unew

      IF( mk == 1 )THEN

        GVECm = GVEC(:,mk)

      ELSE

        BVEC           = - FVEC(:,mk)
        AMAT(:,1:mk-1) = + FVEC(:,1:mk-1) - SPREAD( FVEC(:,mk), DIM = 2, NCOPIES = mk-1 )

        CALL DGELS( 'N', 2, mk-1, 1, AMAT(:,1:mk-1), 2, BVEC, 2, WORK, LWORK, INFO )

        IF( INFO .NE. 0 )THEN
          PRINT*, "INFO  = ", INFO
        END IF

        Alpha(1:mk-1) = BVEC(1:mk-1)
        Alpha(mk) = 1.0_DP - SUM( Alpha(1:mk-1) )

        GVECm = 0.0_DP
        DO i = 1, mk
          GVECm = GVECm + Alpha(i) * GVEC(:,i)
        END DO

      END IF

      FVECm(iY) = GVECm(iY) - Ynew
      FVECm(iU) = GVECm(iU) - Unew

      IF( ENORM( [ FVECm(iY) ] ) <= Rtol * ENORM( [ Yold ] ) .AND. &
          ENORM( [ FVECm(iU) ] ) <= Rtol * ENORM( [ Uold ] ) )THEN

        CONVERGED = .TRUE.

        nIter = k

      END IF

!!$      PRINT*, "  ", k, ABS( FVECm(iY) )/(Rtol*ABS(Yold)), ABS( FVECm(iU) )/(Rtol*ABS(Uold))

      Ynew = GVECm(iY)
      Unew = GVECm(iU)

      IF( mk == M .AND. .NOT. CONVERGED )THEN

        GVEC = CSHIFT( GVEC, SHIFT = + 1, DIM = 2 )
        FVEC = CSHIFT( FVEC, SHIFT = + 1, DIM = 2 )

      END IF

      CALL nuc_eos_full &
             ( Dnew, Tnew, Ynew, Unew, zPrs, zEnt, zCs2, zdEdT, zdPdE_D, &
               zdPdD_E, zXa, zXh, zXn, zXp, zAbar, zZbar, zMu_e, zMu_n,  &
               zMu_p, zMu_hat, keyeint, keyerr, 1.0d-9 )

      Mnew(1) = zMu_e + zMu_p - zMu_n
      Mnew(2) = - Mnew(1)

    END DO

    DO iS = 1, 2

      J0(:,iS) = E / ( EXP( ( E - Mnew(iS) ) / Tnew ) + 1.0_DP )

      Jnew(:,iS) = ( Jold(:,iS) + ccgs * dt * Chi_EC(:,iS) * J0(:,iS) ) / ( 1.0_DP + ccgs * dt * Chi_EC(:,iS) )

    END DO

    ! --- Compute Matter Sources ---

    Q_Y   = ( FourPi * mBcgs )   / ( D * ( ccgs * hMeV )**3 ) &
              * ( SUM  ( ( Jnew(:,1) - Jold(:,1) ) * dV_E / E ) &
                  - SUM( ( Jnew(:,2) - Jold(:,2) ) * dV_E / E ) )

    Q_Enr = ( FourPi * MeV2Erg ) / ( D * ( ccgs * hMeV )**3 ) &
              * ( SUM  ( ( Jnew(:,1) - Jold(:,1) ) * dV_E     ) &
                  + SUM( ( Jnew(:,2) - Jold(:,2) ) * dV_E     ) )

  END SUBROUTINE SolveEmAb_FP_AA


  PURE REAL(DP) FUNCTION ENORM( X )

    REAL(DP), DIMENSION(:), INTENT(in) :: X

    ENORM = SQRT( DOT_PRODUCT( X, X ) )

    RETURN
  END FUNCTION ENORM


  SUBROUTINE WriteOutput

    USE HDF5

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName
    CHARACTER(256) :: GroupName
    INTEGER        :: HDFERR
    INTEGER(HID_T) :: FILE_ID

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName = './Output/Relaxation_' // FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1D_REAL( [ time ], DatasetName, FILE_ID )

    ! --- Write Energy Grid ---

    GroupName = 'Energy Grid'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/E'

    CALL WriteDataset1D_REAL( E, DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/R'

    CALL WriteDataset1D_REAL( xR / 1.0d5, DatasetName, FILE_ID )

    ! --- Write Matter ---

    GroupName = 'Matter'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/D'

    CALL WriteDataset1D_REAL( xD, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/T'

    CALL WriteDataset1D_REAL( xT * kMeV, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Y'

    CALL WriteDataset1D_REAL( xY, DatasetName, FILE_ID )

    ! --- Write Neutrino-Matter Sources ---

    GroupName = 'Neutrino-Matter Sources'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Q_Y'

    CALL WriteDataset1D_REAL( xQ_Y, DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Q_E'

    CALL WriteDataset1D_REAL( xQ_Enr, DatasetName, FILE_ID )

    ! --- Write Radiation ---

    GroupName = 'Radiation'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Mu_1'

    CALL WriteDataset1D_REAL( xMu_nu(1,:), DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/J_1'

    CALL WriteDataset2D_REAL( xJ(:,1,:), DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/Mu_2'

    CALL WriteDataset1D_REAL( xMu_nu(2,:), DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/J_2'

    CALL WriteDataset2D_REAL( xJ(:,2,:), DatasetName, FILE_ID )

    ! --- Solver Data ---

    GroupName = 'Solver'

    CALL CreateGroupHDF( FileName, GroupName , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/nIterations'

    CALL WriteDataset1D_REAL &
           ( DBLE( nIterations ), DatasetName, FILE_ID )

    ! --- Close Files ---

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteOutput


  SUBROUTINE CreateGroupHDF( FileName, GroupName, FILE_ID )

    USE HDF5

    CHARACTER(len=*), INTENT(in) :: FileName
    CHARACTER(len=*), INTENT(in) :: GroupName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER        :: HDFERR
    INTEGER(HID_T) :: GROUP_ID

    CALL H5GCREATE_F( FILE_ID, TRIM( GroupName ), GROUP_ID, HDFERR )

    CALL H5GCLOSE_F( GROUP_ID, HDFERR )

  END SUBROUTINE CreateGroupHDF


  SUBROUTINE WriteDataset1D_REAL( Dataset, DatasetName, FILE_ID )

    USE HDF5

    REAL(DP),         INTENT(in) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER          :: HDFERR
    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SIZE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 1, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset1D_REAL


  SUBROUTINE WriteDataset2D_REAL( Dataset, DatasetName, FILE_ID )

    USE HDF5

    REAL(DP),         INTENT(in) :: Dataset(:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER          :: HDFERR
    INTEGER(HSIZE_T) :: DATASIZE(2)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 2, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset2D_REAL


  SUBROUTINE InitializeTimers

    Timer_Total   = 0.0_DP
    Timer_IO      = 0.0_DP
    Timer_Eos     = 0.0_DP
    Timer_Opacity = 0.0_DP
    Timer_FP_AA   = 0.0_DP

  END SUBROUTINE InitializeTimers


  SUBROUTINE FinalizeTimers

    WRITE(*,*)
    WRITE(*,'(A4,A16,ES10.4E2)') &
      '', 'Timer_Total = ',   Timer_Total
    WRITE(*,*)
    WRITE(*,'(A4,A16,ES10.4E2)') &
      '', 'Timer_IO = ',      Timer_IO
    WRITE(*,'(A4,A16,ES10.4E2)') &
      '', 'Timer_Eos = ',     Timer_Eos
    WRITE(*,'(A4,A16,ES10.4E2)') &
      '', 'Timer_Opacity = ', Timer_Opacity
    WRITE(*,'(A4,A16,ES10.4E2)') &
      '', 'Timer_FP_AA = ',   Timer_FP_AA
    WRITE(*,'(A4,A16,ES10.4E2)') &
      '', 'Sum = ',           Timer_IO + Timer_Eos + Timer_Opacity + Timer_FP_AA
    WRITE(*,*)

  END SUBROUTINE FinalizeTimers


  SUBROUTINE Timer_Start( Timer )

    REAL(DP), INTENT(inout) :: Timer

    Timer = Timer - TimerWTIME()

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP), INTENT(inout) :: Timer

    Timer = Timer + TimerWTIME()

  END SUBROUTINE Timer_Stop


  REAL(DP) FUNCTION TimerWTIME()

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64

    INTEGER(INT64) :: clock_read
    INTEGER(INT64) :: clock_rate
    INTEGER(INT64) :: clock_max

    CALL SYSTEM_CLOCK( clock_read, clock_rate, clock_max )
    TimerWTIME = REAL( clock_read, DP ) / REAL( clock_rate, DP )

    RETURN
  END FUNCTION TimerWTIME


  SUBROUTINE WriteVector( N, Vec, FileName )

    INTEGER,                INTENT(in) :: N
    REAL(DP), DIMENSION(N), INTENT(in) :: Vec
    CHARACTER(LEN=*),       INTENT(in) :: FileName

    INTEGER :: FUNIT

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) N, Vec(1:N)

    CLOSE( FUNIT )

  END SUBROUTINE WriteVector


  SUBROUTINE WriteMatrix( M, N, Mat, FileName )

    INTEGER,                  INTENT(in) :: M, N
    REAL(DP), DIMENSION(M,N), INTENT(in) :: Mat
    CHARACTER(LEN=*),         INTENT(in) :: FileName

    INTEGER :: FUNIT

    OPEN( NEWUNIT = FUNIT, FILE = TRIM( FileName ) )

    WRITE( FUNIT, * ) M, N, Mat(1:M,1:N)

    CLOSE( FUNIT )

  END SUBROUTINE WriteMatrix


END PROGRAM RelaxationTest
