#include "Preprocessor"

submodule ( Fluid_P_HN__Form ) Fluid_P_HN__Kernel

  use Basics
  
  implicit none
  
contains

  
  module procedure Apply_EOS_PrologueKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
           
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( P )

    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nValues

        M ( iV )  =  M_Ref
        N ( iV )  =  max ( N ( iV ), N_Min )

        E ( iV )  =  ( E ( iV )  /  ( M ( iV )  *  N ( iV ) )  -  OR_Shift ) &
                     /  SpecificEnergy_CGS
        N ( iV )  =  M ( iV )  *  N ( iV )  /  MassDensity_CGS
        T ( iV )  =  T ( iV )  /  MeV
        P ( iV )  =  P ( iV )  /  Pressure_CGS

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nValues

        M ( iV )  =  M_Ref
        N ( iV )  =  max ( N ( iV ), N_Min )

        E ( iV )  =  ( E ( iV )  /  ( M ( iV )  *  N ( iV ) )  -  OR_Shift ) &
                     /  SpecificEnergy_CGS
        N ( iV )  =  M ( iV )  *  N ( iV )  /  MassDensity_CGS
        T ( iV )  =  T ( iV )  /  MeV
        P ( iV )  =  P ( iV )  /  Pressure_CGS

      end do
      !$OMP end parallel do
    
    end if
    
  end procedure Apply_EOS_PrologueKernel
  
  
end submodule Fluid_P_HN__Kernel
