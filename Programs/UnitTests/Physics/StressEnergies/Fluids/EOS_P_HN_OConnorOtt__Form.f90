module EOS_P_HN_OConnorOtt__Form

  use HDF5
  use Basics
  
  implicit none
  private
  
  type, public :: EOS_P_HN_OConnorOtt_Form
    integer ( KDI ) :: &
      LOG_PRESSURE              = 1, &
      LOG_ENERGY                = 2, &
      ENTROPY                   = 3, &
      CHEMICAL_POTENTIAL_NU     = 4, &
      SOUND_SPEED_SQUARE        = 5, &
      DE_DT                     = 6, &
      DP_DRHO_E                 = 7, &
      DP_DE_RHO                 = 8, &
      CHEMICAL_POTENTIAL_HAT    = 9, &
      CHEMICAL_POTENTIAL_E      = 10, &
      CHEMICAL_POTENTIAL_P      = 11, &
      CHEMICAL_POTENTIAL_N      = 12, &
      MASS_FRACTION_A           = 13, &
      MASS_FRACTION_H           = 14, &
      MASS_FRACTION_N           = 15, &
      MASS_FRACTION_P           = 16, &
      MASS_NUMBER_BAR           = 17, &
      ATOMIC_NUMBER_BAR         = 18, &
      GAMMA                     = 19, &
      N_VARIABLES               = 19
    integer ( KDI ) :: &
      nDensity, nTemperature, nElectronFraction, &
      MinDensity, MaxDensity, &
      MinElectronFraction, MaxElectronFraction, &
      MinTemperature, MaxTemperature
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaFluidOutput, &
      iaSelected, &
      Error
    real ( KDR ) :: & 
      EnergyShift
    real ( KDR ), dimension ( : ), allocatable :: &
      LogDensity, &
      LogTemperature, &
      ElectronFraction
    real ( KDR ), dimension ( :, :, :, : ), allocatable :: &
      Table
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      ComputeFromTemperature
    final :: &
      Finalize
  end type EOS_P_HN_OConnorOtt_Form
  
  real ( KDR ), parameter :: &
    T_MAX_HACK          = 240.0_KDR, &
    mev_to_erg          = 1.60217733e-6_KDR, &
    amu_cgs             = 1.66053873e-24_KDR, &
    amu_mev             = 931.49432e0_KDR, &
    pi                  = 3.14159265358979e0_KDR, &
    ggrav               = 6.672e-8_KDR, &
    temp_mev_to_kelvin  = 1.1604447522806e10_KDR, &
    clight              = 2.99792458e10_KDR, &
    kb_erg              = 1.380658e-16_KDR, &
    kb_mev              = 8.61738568e-1_KDR


contains

  subroutine Initialize ( E )
  
    class ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
    
    integer ( KDI ) :: &
      error, &
      accerr
    integer ( HID_T ) :: &
      file_id, &
      dset_id, dspace_id
    integer ( HSIZE_T ), dimension ( 1 ) :: &
      dims1
    integer ( HSIZE_T ), dimension ( 3 ) :: &
      dims3
    character ( LDF ) :: &
      Filename
    
    !-- FIXME: Hardcoded filename for now
    Filename = '../Parameters/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5'
        
    call h5open_f(error)

    call h5fopen_f (trim(adjustl(Filename)), H5F_ACC_RDONLY_F, file_id, error)

    associate &
      ( nrho          => E % nDensity, &
        ntemp         => E % nTemperature, &
        nye           => E % nElectronFraction, &
        nvars         => E % N_VARIABLES, &
        energy_shift  => E % EnergyShift )
        
    ! read scalars
    dims1(1)=1
    call h5dopen_f(file_id, "pointsrho", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
    call h5dclose_f(dset_id,error)

    dims1(1)=1
    call h5dopen_f(file_id, "pointstemp", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
    call h5dclose_f(dset_id,error)

    dims1(1)=1
    call h5dopen_f(file_id, "pointsye", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
    call h5dclose_f(dset_id,error)

    if(error.ne.0) then
       stop "Could not read EOS table file"
    endif

    call Show ( Filename, 'Reading Equation of State' )
    call Show ( [ E % nDensity, E % nTemperature, E % nElectronFraction ], &
                'Table sizes' )
    call Show ( E % N_VARIABLES, 'nVariables' )

    
    allocate ( E % Table ( nrho, ntemp, nye, nvars ) )
    
    associate ( alltables => E % Table )
    
    call Show ( shape ( E % Table ), 'Table Shape' ) 
    
    dims3(1)=nrho
    dims3(2)=ntemp
    dims3(3)=nye

    call h5dopen_f(file_id, "logpress", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,1), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "logenergy", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,2), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "entropy", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,3), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "munu", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,4), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "cs2", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,5), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "dedt", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,6), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "dpdrhoe", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,7), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error
    call h5dopen_f(file_id, "dpderho", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,8), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

  ! chemical potentials
    call h5dopen_f(file_id, "muhat", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,9), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "mu_e", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,10), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "mu_p", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,11), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "mu_n", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,12), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

  ! compositions
    call h5dopen_f(file_id, "Xa", dset_id, error) 
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,13), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "Xh", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,14), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "Xn", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,15), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "Xp", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,16), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    ! Gamma
    call h5dopen_f(file_id, "gamma", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,19), dims3, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    allocate ( E % LogDensity ( E % nDensity ) )
    dims1(1)=nrho
    call h5dopen_f(file_id, "logrho", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, E % LogDensity, dims1, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    allocate ( E % LogTemperature ( E % nTemperature ) )
    dims1(1)=ntemp
    call h5dopen_f(file_id, "logtemp", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, E % LogTemperature, dims1, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    allocate ( E % ElectronFraction ( E % nElectronFraction ) )
    dims1(1)=nye
    call h5dopen_f(file_id, "ye", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, E % ElectronFraction, dims1, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5dopen_f(file_id, "energy_shift", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, energy_shift, dims1, error)
    call h5dclose_f(dset_id,error)
    accerr=accerr+error

    call h5fclose_f (file_id,error)

    call h5close_f (error)

    end associate !-- alltables
    end associate !-- logrho, etc
      
  end subroutine Initialize
  
  
  subroutine SelectVariables ( E, iaFluidOutput, iaSelected )
    
    class ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaFluidOutput, &
      iaSelected
    
    allocate ( E % iaFluidOutput, source = ( iaFluidOutput ) )
    allocate ( E % iaSelected, source = ( iaSelected ) )
  
  end subroutine SelectVariables

  
  subroutine ComputeFromTemperature ( E, Fluid, iaFluidInput )
  
    class ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
    class ( StorageForm ), intent ( inout ) :: &
      Fluid
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaFluidInput
    
    call Interpolate_3D_Kernel &
           ( Fluid % Value, E % Table, E % LogDensity, E % LogTemperature, &
             E % ElectronFraction, E % EnergyShift, iaFluidInput, &
             E % iaFluidOutput, E % iaSelected )
  
  end subroutine ComputeFromTemperature
  
  
  subroutine Finalize ( E )
  
    type ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
  
  
  end subroutine Finalize
  
  
  
  subroutine Interpolate_3D_Kernel &
               ( F, T, XT, YT, ZT, E_Shift, ia_F_I, ia_F_O, ia_E )
    
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      F
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      T
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      XT, &      !-- LogDensity (typically)
      YT, &      !-- LogTemperature (typically)
      ZT         !-- ElectronFraction  (typically)
    real ( KDR ), intent ( in ) :: &
      E_Shift
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      ia_F_I, &  !-- iaFluidInput
      ia_F_O, &  !-- iaFluidOutput
      ia_E       !-- iaEOS
    
    integer ( KDI ) :: &
      iValue, &
      iV, &  !-- iVariable
      iS, &  !-- iSelected
      iF, &  !-- iFluid
      nx,ny,nz,nvars, &
      nValues, &
      ix,iy,iz
    real ( KDR ) :: &
      L_x, &
      L_y, &
      delx, dely, delz, &
      dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi,&
      a1, a2, a3, a4, a5, a6, a7, a8 
    real ( KDR ), dimension ( 8 ) :: &
      fh

    call Show ( 'Interpolate_3D' )
    
    nx = size ( XT )
    ny = size ( YT )
    nz = size ( ZT )
    nvars = size ( T, dim = 4 )
    nValues = size ( F, dim = 1 )

    !--  determine spacing parameters of (equidistant!!!) table
    dx    = (xt(nx) - xt(1)) / real(nx-1, kind=KDR)
    dy    = (yt(ny) - yt(1)) / real(ny-1, kind=KDR)
    dz    = (zt(nz) - zt(1)) / real(nz-1, kind=KDR)

    dxi   = 1. / dx
    dyi   = 1. / dy
    dzi   = 1. / dz

    dxyi  = dxi * dyi
    dxzi  = dxi * dzi
    dyzi  = dyi * dzi

    dxyzi = dxi * dyi * dzi
    
    associate ( ft => T )
    
    do  iValue = 1, nValues
      
      !-- Convert to log space for the x and y
      L_x = log10 ( F ( iValue, ia_F_I ( 1 ) ) )
      L_y = log10 ( F ( iValue, ia_F_I ( 2 ) ) )
      
      !-- determine location in (equidistant!!!) table 
      ix = 2 + INT( (L_x - xt(1) - 1.e-10_KDR) * dxi )
      iy = 2 + INT( (L_y - yt(1) - 1.e-10_KDR) * dyi )
      iz = 2 + INT( ( F ( iValue, ia_F_I ( 3 ) ) - zt(1) - 1.e-10_KDR) * dzi )
                                                 
      ix = MAX( 2, MIN( ix, nx ) )
      iy = MAX( 2, MIN( iy, ny ) )   
      iz = MAX( 2, MIN( iz, nz ) )

      !-- set-up auxiliary arrays for Lagrange interpolation
                                                             
      delx = xt(ix) - L_x
      dely = yt(iy) - L_y
      delz = zt(iz) - F ( iValue, ia_F_I ( 3 ) )
  
      do iS = 1, size ( ia_E )
        
        iV = ia_E ( iS )
        iF = ia_F_O ( iS )
        
        fh(1) = ft(ix  , iy  , iz,   iv)
        fh(2) = ft(ix-1, iy  , iz,   iv)   
        fh(3) = ft(ix  , iy-1, iz,   iv)   
        fh(4) = ft(ix  , iy  , iz-1, iv)
        fh(5) = ft(ix-1, iy-1, iz,   iv)
        fh(6) = ft(ix-1, iy  , iz-1, iv)
        fh(7) = ft(ix  , iy-1, iz-1, iv)
        fh(8) = ft(ix-1, iy-1, iz-1, iv)
        
        !-- set up coefficients of the interpolation polynomial and 
        !   evaluate function values 
        a1 = fh(1)
        a2 = dxi   * ( fh(2) - fh(1) )
        a3 = dyi   * ( fh(3) - fh(1) )
        a4 = dzi   * ( fh(4) - fh(1) )
        a5 = dxyi  * ( fh(5) - fh(2) - fh(3) + fh(1) )
        a6 = dxzi  * ( fh(6) - fh(2) - fh(4) + fh(1) )
        a7 = dyzi  * ( fh(7) - fh(3) - fh(4) + fh(1) )
        a8 = dxyzi * ( fh(8) - fh(1) + fh(2) + fh(3) + &
             fh(4) - fh(5) - fh(6) - fh(7) )

        f(iValue, iF)  &
          = a1 +  a2 * delx               &
             +  a3 * dely                      &  
             +  a4 * delz                      &  
             +  a5 * delx * dely               &  
             +  a6 * delx * delz               &
             +  a7 * dely * delz               &
             +  a8 * delx * dely * delz
        
        !-- T ( :, :, :, 1 ) and T ( :, :, :, 2 ) is in log space
        if ( iV == 1 .or. iV == 2 ) then
          f ( iValue, iF ) = 10.0e0_KDR ** f ( iValue, iF )
        end if
        if ( iV == 2 ) then
          f ( iValue, iF ) = f ( iValue, iF ) - E_Shift
        end if
                  
      enddo
      
    enddo
                      
    end associate
  
  end subroutine Interpolate_3D_Kernel


end module EOS_P_HN_OConnorOtt__Form
