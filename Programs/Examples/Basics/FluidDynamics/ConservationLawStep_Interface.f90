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
