  interface 

    subroutine ComputeDifferences_C &
                 ( V, lV, uV, iaS, iD, dV_Left, dV_Right, nSizes ) &
                 bind ( c, name = 'ComputeDifferences_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &   
        V
      integer ( c_int ), dimension ( 3 ) :: &
        lV, uV, &
        iaS
      integer ( c_int ), value :: &  
        iD
      type ( c_ptr ), value :: &
        dV_Left, dV_Right
      integer ( c_int ), dimension ( 3 ) :: &
        nSizes
    end subroutine ComputeDifferences_C
    
    subroutine ComputeReconstruction_C &
                 ( V, dV_Left, dV_Right, Theta, V_Inner, V_Outer, &
                   nValues ) &
                 bind ( c, name = 'ComputeReconstruction_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        V, &
        dV_Left, dV_Right
      real ( c_double ), value :: &
        Theta
      type ( c_ptr ), value :: &
        V_Inner, V_Outer
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeReconstruction_C
    
    subroutine ComputeFluxes_C & 
                 ( AP_I, AP_O, AM_I, AM_O, RF_I, RF_O, &
                   U_I, U_O, lV, uV, iaS, iD, F_I, F_O, nSizes ) &
                 bind ( c, name = 'ComputeFluxes_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        AP_I, AP_O, &
        AM_I, AM_O, &
        RF_I, RF_O, &
        U_I, U_O
      integer ( c_int ), dimension ( 3 ) :: &
        lV, uV, &
        iaS
      integer ( c_int ), value :: &
        iD
      type ( c_ptr ), value :: &
        F_I, F_O
      integer ( c_int ), dimension ( 3 ) :: &
        nSizes
    end subroutine ComputeFluxes_C

    subroutine ComputeUpdate_C &
                 ( dU, F_I, F_O, V, A, dT, nValues ) &
                 bind ( c, name = 'ComputeUpdate_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &  
        dU, &
        F_I, F_O
      real ( c_double ), value :: &
        V, &
        A, &
        dT
      integer ( c_int ), value :: &
        nValues
    end subroutine ComputeUpdate_C
     
    subroutine AddUpdate_C ( O, U, C, nValues ) &
                 bind ( c, name = 'AddUpdate_C' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        O, &
        U
      type ( c_ptr ), value :: &
        C 
      integer ( c_int ), value :: &
        nValues
    end subroutine AddUpdate_C

    subroutine CombineUpdates_C ( C, O, U, nValues ) &
                 bind ( c, name = 'CombineUpdates_C' )
      use iso_c_binding
      implicit none 
      type ( c_ptr ), value :: &
         C 
      type ( c_ptr ), value :: &
         O, &
         U
      integer ( c_int ), value :: &
        nValues
    end subroutine CombineUpdates_C

  end interface
