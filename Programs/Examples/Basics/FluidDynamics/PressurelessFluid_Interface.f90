  interface
    
    subroutine ComputeConservedPressureless_C &
                ( D, S_1, S_2, S_3, N, V_1, V_2, V_3, nValues ) &    
                 bind ( c, name = 'ComputeConservedPressureless_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        D, &
        S_1, S_2, S_3, &
        N, &
        V_1, V_2, V_3
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeConservedPressureless_C
 
    subroutine ComputePrimitivePressureless_C  &
                ( N, V_1, V_2, V_3, D, S_1, S_2, S_3, nValues ) &
                 bind ( c, name = 'ComputePrimitivePressureless_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        N, &
        V_1, V_2, V_3, &
        D, &
        S_1, S_2, S_3
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputePrimitivePressureless_C
 
    subroutine ComputeEigenspeedsPressureless_C &
                 ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                   V_1, V_2, V_3, nValues ) &
                 bind ( c, name = 'ComputeEigenspeedsPressureless_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        FEP_1, FEP_2, FEP_3, &
        FEM_1, FEM_2, FEM_3, &
        V_1, V_2, V_3
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeEigenspeedsPressureless_C

    subroutine ApplyBoundaryConditionsReflectingPressureless_C &
                 ( N_E, VI_E, VJ_E, VK_E, N_I, VI_I, VJ_I, VK_I, &
                   nB, oBE, oBI, nSizes ) &
                 bind ( c, name = 'ApplyBoundaryConditionsReflectingPressureless_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        N_E, &
        VI_E, VJ_E, VK_E, &
        N_I, &
        VI_I, VJ_I, VK_I, &
        nB, &
        oBE, oBI
      integer ( c_int ), dimension ( 3 ) :: &
        nSizes
    end subroutine ApplyBoundaryConditionsReflectingPressureless_C 
 
    subroutine ComputeRiemannSolverInputPressureless_C &
                ( AP_I, AP_O, AM_I, AM_O, LP_I, LP_O, LM_I, LM_O, &
                  lV, uV, iaS_M, iaS_P, nSizes ) &
                 bind ( c, name = 'ComputeRiemannSolverInputPressureless_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        AP_I, AP_O, &
        AM_I, AM_O, &
        LP_I, LP_O, &
        LM_I, LM_O
      integer ( c_int ), dimension ( 3 ) :: &
        lV, uV, &
        iaS_M, iaS_P, &
        nSizes
    end subroutine ComputeRiemannSolverInputPressureless_C 
  
  end interface
