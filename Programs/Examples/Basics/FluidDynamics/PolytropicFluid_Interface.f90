  interface
    
    subroutine ComputeConservedPolytropic_C &
                 ( G, E, N, V_1, V_2, V_3, nValues ) &
                 bind ( c, name = 'ComputeConservedPolytropic_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        G, &
        E, &
        N, &
        V_1, V_2, V_3
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeConservedPolytropic_C
 
    subroutine ComputePrimitivePolytropic_C &
                 ( E, G, N, V_1, V_2, V_3, nValues ) &
                 bind ( c, name = 'ComputePrimitivePolytropic_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        E, &
        G, &
        N, &
        V_1, V_2, V_3
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputePrimitivePolytropic_C
 
    subroutine ComputeAuxiliaryPolytropic_C &
                 ( P, K, N, E, Gamma, nValues ) &
                 bind ( c, name = 'ComputeAuxiliaryPolytropic_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        P, &
        K, &
        N, &
        E, &
        Gamma
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeAuxiliaryPolytropic_C

    subroutine ComputeEigenspeedsPolytropic_C &
                 ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, CS, N, &
                   V_1, V_2, V_3, P, Gamma, nValues ) &
                 bind ( c, name = 'ComputeEigenspeedsPolytropic_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        FEP_1, FEP_2, FEP_3, &
        FEM_1, FEM_2, FEM_3, &
        CS, &
        N, &
        V_1, V_2, V_3, &
        P, &
        Gamma
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeEigenspeedsPolytropic_C
!
!    subroutine ApplyBoundaryConditionsReflectingPolytropic_C &
!                 ( E_E, Gamma_E, E_I, Gamma_I, nB, oBE, oBI, nSizes ) &
!                 bind ( c, name = 'ApplyBoundaryConditionsReflectingPolytropic_C' )
!      use iso_c_binding
!      implicit none
!      type ( c_ptr ), value :: &
!        E_E, &
!        Gamma_E, &
!        E_I, &
!        Gamma_I, &
!        nB, &
!        oBE, oBI
!      integer ( c_int ), dimension ( 3 ) :: &
!        nSizes
!    end subroutine ApplyBoundaryConditionsReflectingPolytropic_C 
! 
    subroutine ComputeRawFluxesPolytropic_C &
                 ( F_D, F_S_1, F_S_2, F_S_3, F_S_Dim, F_G, D, &
                   S_1, S_2, S_3, G, P, V_Dim, nValues ) &
                 bind ( c, name = 'ComputeRawFluxesPolytropic_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        F_D, &
        F_S_1, F_S_2, F_S_3, &
        F_S_Dim, &
        F_G, &
        D, &
        S_1, S_2, S_3, &
        G, &
        P, &
        V_Dim
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeRawFluxesPolytropic_C  

  end interface
