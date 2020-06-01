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
    logical ( KDL ) :: &
      AllocatedDevice
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_EOS_P_HN
    procedure, public, pass :: &
      SelectVariables
    procedure, public, pass :: &
      ComputeFromTemperature, &
      ComputeFromEnergy
    final :: &
      Finalize
  end type EOS_P_HN_OConnorOtt_Form
  
    private :: &
      Interpolate_3D_Kernel, &
      InterpolateTableKernel
  
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


  interface
  
    module subroutine Interpolate_3D_Kernel &
                 ( F, T, XT, YT, ZT, E_Shift, ia_F_I, ia_F_O, ia_E, &
                   UseDeviceOption )
      use Basics      
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
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine Interpolate_3D_Kernel
  
  end interface

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
    
    allocate ( E % iaFluidOutput, source = iaFluidOutput )
    allocate ( E % iaSelected, source = iaSelected )
    
  end subroutine SelectVariables
  
  
  subroutine AllocateDevice_EOS_P_HN( E )
  
    class ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
      
    type ( c_ptr ) :: &
      D_P  !-- Device Pointer

    call AllocateDevice ( E % Error, D_P )
    call AssociateHost  ( D_P, E % Error )

    call AllocateDevice ( E % Table, D_P )
    call AssociateHost  ( D_P, E % Table )
    call UpdateDevice   ( E % Table, D_P )

    call AllocateDevice ( E % LogDensity, D_P )
    call AssociateHost  ( D_P, E % LogDensity )
    call UpdateDevice   ( E % LogDensity, D_P )

    call AllocateDevice ( E % LogTemperature, D_P )
    call AssociateHost  ( D_P, E % LogTemperature )
    call UpdateDevice   ( E % LogTemperature, D_P )

    call AllocateDevice ( E % ElectronFraction, D_P )
    call AssociateHost  ( D_P, E % ElectronFraction )
    call UpdateDevice   ( E % ElectronFraction, D_P )
    
    E % AllocatedDevice = .true.

  end subroutine AllocateDevice_EOS_P_HN

  
  subroutine ComputeFromTemperature ( E, Fluid, iaFluidInput )
  
    class ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
    class ( StorageForm ), intent ( inout ) :: &
      Fluid
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaFluidInput
    
    logical ( KDL ) :: &
      UseDevice 
      
    call Interpolate_3D_Kernel &
           ( Fluid % Value, E % Table, E % LogDensity, E % LogTemperature, &
             E % ElectronFraction, E % EnergyShift, iaFluidInput, &
             E % iaFluidOutput, E % iaSelected, &
             UseDeviceOption = UseDevice )
  
  end subroutine ComputeFromTemperature
  
  
  subroutine ComputeFromEnergy ( E, Fluid, iaFluidInput, iSolve )
  
    class ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
    class ( StorageForm ), intent ( inout ) :: &
      Fluid
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaFluidInput
    integer ( KDI ), intent ( in ) :: &
      iSolve
    
    integer ( KDI ) :: &
      iSelected
    
    call Search ( E % iaFluidOutput, iSolve, iSelected )
    
    call FindTemperatureKernel &
           ( Fluid % Value, E % Table, E % LogDensity, E % LogTemperature, &
             E % ElectronFraction, iaFluidInput, iSolve, &
             E % iaSelected ( iSelected ), ShiftOption = E % EnergyShift, &
             LogScaleOption = .true. )
    
    call E % ComputeFromTemperature ( Fluid, iaFluidInput )
    
  end subroutine ComputeFromEnergy
    
  
  subroutine Finalize ( E )
  
    type ( EOS_P_HN_OConnorOtt_Form ), intent ( inout ) :: &
      E
    
    if ( allocated ( E % Table ) ) &
      deallocate ( E % Table )
   
    if ( allocated ( E % ElectronFraction ) ) &
      deallocate ( E % ElectronFraction )
    
    if ( allocated ( E % LogTemperature ) ) &
      deallocate ( E % LogTemperature ) 
    
    if ( allocated ( E % LogDensity ) ) &
      deallocate ( E % LogDensity )
    
    if ( allocated ( E % Error ) ) &
      deallocate ( E % Error )
      
    if ( allocated ( E % iaSelected ) ) &
      deallocate ( E % iaSelected )
    
    if ( allocated ( E % iaFluidOutput ) ) &
      deallocate ( E % iaFluidOutput )
  
  end subroutine Finalize
  
  
  subroutine FindTemperatureKernel &
               ( F, T, T_L_N, T_L_T, T_Ye, ia_F_I, i_SF, i_ST, &
                 ToleranceOption, nIterationsOption, ShiftOption, &
                 LogScaleOption )
  
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      F
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      T
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      T_L_N, &      !-- TableLogDensity
      T_L_T, &      !-- TableLogTemperature
      T_Ye          !-- TableElectronFraction
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      ia_F_I        !-- iaFluidInput
    integer ( KDI ), intent ( in ) :: & 
      i_SF, &       !-- index of Fluid to solve from
      i_ST          !-- index of the corresponding quantity in EOS table 
    logical ( KDL ), intent ( in ), optional :: &
      LogScaleOption
    real ( KDR ), intent ( in ), optional :: &
      ShiftOption, &
      ToleranceOption
    integer ( KDI ), intent ( in ), optional :: &
      nIterationsOption
    
    integer ( KDI ) :: &
      iV, &
      iI, &   !-- iIteration
      nValues, &
      nIterations
    real ( KDR ) :: &
      Shift, &
      Tolerance, &
      T_L_T_Max, T_L_T_Min, &
      L_N, Ye, &
      SV_R, D2, &
      L_T, L_T_1, &   !-- LogTemperature (from Fluid input)
      SV_0, SV_1, &
      L_DT
    !-- FIXME: temporary
    real ( KDR ) :: &
      D1, D3
    logical ( KDL ) :: &
      LogScale
    
    
    Shift = 0.0_KDR
    if ( present ( ShiftOption ) ) &
      Shift = ShiftOption
    
    nIterations = 20
    if ( present ( nIterationsOption ) ) &
      nIterations = nIterationsOption
    
    Tolerance = 1e-10_KDR
    if ( present ( ToleranceOption ) ) &
      Tolerance = ToleranceOption
      
    LogScale = .false.
    if ( present ( LogScaleOption ) ) &
      LogScale = .true.
    
    nValues = size ( F, dim = 1 )
    
    T_L_T_Max = T_L_T ( size ( T_L_T ) )
    T_L_T_Min = T_L_T ( 1 )
    
    !call Show ( [ T_L_T_Max, T_L_T_Min ], 'T_L_T_Max, T_L_T_Min' )
    !call Show ( i_ST , 'i_ST' )
    
    Value: do iV = 1, nValues
      
      L_N   = log10 ( F ( iV, ia_F_I ( 1 ) ) )
      L_T   = log10 ( F ( iV, ia_F_I ( 2 ) ) )
      Ye    = F ( iV, ia_F_I ( 3 )  )
      
      L_T_1 = L_T
      
      if ( LogScale ) then
        SV_0 = log10 ( max ( ( F ( iV, i_SF ) + Shift ), 1.0_KDR ) )
      else
        SV_0 = F ( iV, i_SF ) + Shift
      end if 
      SV_1 = SV_0
      
      call InterpolateTableKernel &
             ( L_N, L_T, Ye, T, T_L_N, T_L_T, T_Ye, i_ST, SV_R, D2 )
      
      if ( abs ( SV_R - SV_0 ) < Tolerance * abs ( SV_0 ) ) then
        cycle Value
      end if
      
      do iI = 1, nIterations
        !call Show ( iI, 'iIteration' )
        
        L_DT  = - ( SV_R - SV_0 ) / D2
        L_T_1 = L_T
        L_T   = max ( min ( ( L_T + L_DT ), T_L_T_Max ), T_L_T_Min )
        SV_1  = SV_R
        
        !call Show ( [ D2, L_DT, L_T ], 'D2, L_DT, L_T' )
        !call Show ( [ D2, L_DT, L_T ], 'D2, L_DT, L_T' )
        
        call InterpolateTableKernel &
               ( L_N, L_T, Ye, T, T_L_N, T_L_T, T_Ye, i_ST, SV_R, D2 )
      
        if ( abs ( SV_R - SV_0 )  <  Tolerance * abs ( SV_0 ) ) then
          !call Show ( 'SUCCESS: Tolerance satisfied' )
          F ( iV, ia_F_I ( 2 ) ) = 10.0_KDR ** L_T
          cycle Value
        endif 
        
        ! if we are closer than 10^-2  to the 
        ! root (eps-eps0)=0, we are switching to 
        ! the secant method, since the table is rather coarse and the
        ! derivatives may be garbage.

        if ( abs ( SV_R - SV_0 )  <  1e-3_KDR * abs ( SV_0 ) ) then
          !call Show ( 'Switching to Secand Method' )          
          D2 = ( SV_R - SV_1 ) / ( L_T - L_T_1 )
        end if
      
        if ( iI == nIterations ) &
          call Show ( 'Error, Max iteration reached' )
      end do 

    end do Value
    
  end subroutine FindTemperatureKernel
  
  
  subroutine InterpolateTableKernel &
               ( X, Y, Z, T, XT, YT, ZT, i_ST, &
                 SV_R, D2 )
    
    real ( KDR ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
      T
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      XT, &      !-- TableLogDensity
      YT, &      !-- TableLogTemperature
      ZT         !-- TableElectronFraction
    integer ( KDI ), intent ( in ) :: &
      i_ST          !-- index of table variable
    real ( KDR ), intent ( out ) :: &
      SV_R, &
      D2
    
    integer ( KDI ) :: &
      iValue, &
      nx,ny,nz,nvars, &
      ix,iy,iz
    real ( KDR ) :: &
      delx, dely, delz, &
      dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi, &
      a1, a2, a3, a4, a5, a6, a7, a8  
    real ( KDR ), dimension ( 8 ) :: &
      fh

    nx = size ( XT )
    ny = size ( YT )
    nz = size ( ZT )
    nvars = size ( T, dim = 4 )  

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
      
    !-- determine location in (equidistant!!!) table 
    ix = 2 + INT( (X - xt(1) - 1.e-10_KDR) * dxi ) 
    iy = 2 + INT( (Y - yt(1) - 1.e-10_KDR) * dyi ) 
    iz = 2 + INT( (Z - zt(1) - 1.e-10_KDR) * dzi )
                                               
    ix = MAX( 2, MIN( ix, nx ) )               
    iy = MAX( 2, MIN( iy, ny ) )               
    iz = MAX( 2, MIN( iz, nz ) )               

    !-- set-up auxiliary arrays for Lagrange interpolation
                                                          
    delx = xt(ix) - X                                   
    dely = yt(iy) - Y
    delz = zt(iz) - Z

    fh(1) = T(ix  , iy  , iz,   i_ST)
    fh(2) = T(ix-1, iy  , iz,   i_ST)
    fh(3) = T(ix  , iy-1, iz,   i_ST)
    fh(4) = T(ix  , iy  , iz-1, i_ST)
    fh(5) = T(ix-1, iy-1, iz,   i_ST)
    fh(6) = T(ix-1, iy  , iz-1, i_ST)
    fh(7) = T(ix  , iy-1, iz-1, i_ST)
    fh(8) = T(ix-1, iy-1, iz-1, i_ST)

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

    D2 = -a3
    SV_R  &
      = a1 +  a2 * delx               &
         +  a3 * dely                      &
         +  a4 * delz                      &
         +  a5 * delx * dely               &
         +  a6 * delx * delz               &
         +  a7 * dely * delz               &
         +  a8 * delx * dely * delz

  end subroutine InterpolateTableKernel
  
  
end module EOS_P_HN_OConnorOtt__Form
