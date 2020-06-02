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
      ComputeFromEnergy, &
      ComputeFromEntropy
    final :: &
      Finalize
  end type EOS_P_HN_OConnorOtt_Form
  
    private :: &
      Interpolate_3D_Kernel, &
      FindTemperatureKernel, &
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
    
    module subroutine FindTemperatureKernel &
               ( F, T, T_L_N, T_L_T, T_Ye, ia_F_I, i_SF, i_ST, &
                 LogScaleOption, UseDeviceOption, ShiftOption, &
                 ToleranceOption, nIterationsOption )
      use Basics
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
        LogScaleOption, &
        UseDeviceOption
      real ( KDR ), intent ( in ), optional :: &
        ShiftOption, &
        ToleranceOption
      integer ( KDI ), intent ( in ), optional :: &
        nIterationsOption
    end subroutine FindTemperatureKernel
  
    module subroutine InterpolateTableKernel &
                 ( X, Y, Z, T, XT, YT, ZT, i_ST, &
                   SV_R, D2 )
      use Basics
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
    end subroutine InterpolateTableKernel
  
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
    
    UseDevice = ( E % AllocatedDevice .and. Fluid % AllocatedDevice )
    
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
    logical ( KDL ) :: &
      UseDevice 
    
    UseDevice = ( E % AllocatedDevice .and. Fluid % AllocatedDevice )
    
    call Search ( E % iaFluidOutput, iSolve, iSelected )
    
    call FindTemperatureKernel &
           ( Fluid % Value, E % Table, E % LogDensity, E % LogTemperature, &
             E % ElectronFraction, iaFluidInput, iSolve, &
             E % iaSelected ( iSelected ), ShiftOption = E % EnergyShift, &
             LogScaleOption = .true., UseDeviceOption = UseDevice )
    
    call E % ComputeFromTemperature ( Fluid, iaFluidInput )
    
  end subroutine ComputeFromEnergy
    
  
  subroutine ComputeFromEntropy ( E, Fluid, iaFluidInput, iSolve )
  
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
    logical ( KDL ) :: &
      UseDevice 
    
    UseDevice = ( E % AllocatedDevice .and. Fluid % AllocatedDevice )
    
    call Search ( E % iaFluidOutput, iSolve, iSelected )
    
    call FindTemperatureKernel &
           ( Fluid % Value, E % Table, E % LogDensity, E % LogTemperature, &
             E % ElectronFraction, iaFluidInput, iSolve, &
             E % iaSelected ( iSelected ), UseDeviceOption = UseDevice )
    
    call E % ComputeFromTemperature ( Fluid, iaFluidInput )
    
  end subroutine ComputeFromEntropy
    
  
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
  
  
end module EOS_P_HN_OConnorOtt__Form
