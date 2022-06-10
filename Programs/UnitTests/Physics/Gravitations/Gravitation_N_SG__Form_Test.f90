program Gravitation_N_SG__Form_Test

  !-- Gravitation_NewtonSelfGravity_Atlas__Form_Test

  use Basics
  use Mathematics
  use Gravitations

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( Gravitation_N_SG_Form ), allocatable :: &
    G

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Gravitation_N_SG__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A, GravitationalConstant = 1.0_KDR )
  call G % SetStream ( S )
  call S % AddFieldSet ( G % Source )

  call A % Show ( )
  call G % Show ( )

  call TestHomogeneousSpheres ( )

  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestHomogeneousSpheres ( )

    integer ( KDI ) :: &
      iHS
    real ( KDR ), dimension ( 3 ) :: &
      Radius, &
      Density
    type ( FieldSetForm ), allocatable :: &
      Fluid, &
      Reference, &
      Difference
    type ( TimerForm ), pointer :: &
      T_G, &
      T_W

    call Show ( 'Testing homogeneous spheres' )

    associate ( C  =>  A % Chart_GS )
    
    allocate ( Fluid )
    call Fluid % Initialize &
           ( A, &
             FieldOption = [ 'BaryonMass   ', 'BaryonDensity' ], &
             NameOption = 'Fluid', &
             nFieldsOption = 2 )
    
    allocate ( Reference )
    call Reference % Initialize &
           ( A, &
             FieldOption = [ 'Potential' ], &
             NameOption = 'Reference', &
             nFieldsOption = 1 )
    
    allocate ( Difference )
    call Difference % Initialize &
           ( A, &
             FieldOption = [ 'Potential' ], &
             NameOption = 'Difference', &
             nFieldsOption = 1 )
    
    call S % AddFieldSet ( Fluid )
    call S % AddFieldSet ( Reference )
    call S % AddFieldSet ( Difference )
    call S % Show ( )

    Radius  =  C % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ]
    call PROGRAM_HEADER % GetParameter ( Radius, 'Radius' )

    Density  =  1.0_KDR  /  Radius ** 3
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iHS  =  1, 3

      call SetHomogeneousSphere &
             ( Fluid, Reference, G, Density ( iHS ), Radius ( iHS )  )
    
      call Fluid % UpdateDevice ( )

      T_G  =>  G % Timer ( Level = 1 )
      call T_G % Start ( )
      call G % Solve ( Fluid, iBaryonMass = 1, iBaryonDensity = 2 )
      call T_G % Stop ( )

      call Show ( Radius ( iHS ), 'Radius', nLeadingLinesOption = 2 )
      call Show ( Density ( iHS ), 'Density' )
      call ComputeError ( Difference, G % Solution, Reference )

      T_W  =>  S % TimerWrite ( Level = 1 )
      call T_W % Start ( )
      call GIS % Open ( GIS % ACCESS_CREATE )
      call S % Write ( )
      call GIS % Close ( )
      call T_W % Stop ( )

    end do !-- iHS

    end associate !-- C

  end subroutine TestHomogeneousSpheres


  subroutine ComputeError ( Difference, Solution, Reference )

    class ( FieldSetForm ), intent ( inout ) :: &
      Difference, &         
      Solution, &
      Reference

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    associate &
      ( SV  =>  Solution % Storage_GS % Value &
                  ( :, Solution % iaSelected ( 1 ) ), &
        RV  =>  Reference  % Storage_GS % Value ( :, 1 ), &
        DV  =>  Difference % Storage_GS % Value ( :, 1 ) )

    call MultiplyAdd ( SV, RV, -1.0_KDR, DV )

    select type ( A  =>  S % Atlas )
      class is ( Atlas_SCG_Form ) 
    associate &
      ( C  =>  A % Chart_GS )

    call CO % Initialize &
           ( C % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 )  &
      =  sum ( abs ( DV ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 )  &
      =  sum ( abs ( RV ), mask = C % ProperCell )
    call CO % Reduce ( REDUCTION % SUM )

    associate &
      ( Norm_D  =>  CO % Incoming % Value ( 1 ), &
        Norm_R  =>  CO % Incoming % Value ( 2 ) )

    L1  =  Norm_D / Norm_R

    end associate !-- Norm_D_1, etc.

    call Show ( L1, '*** L1 error', nTrailingLinesOption = 2 )

    end associate !-- C
    end select !-- A
    end associate !-- SV, etc.

  end subroutine ComputeError


  subroutine SetHomogeneousSphere &
               ( Fluid, Reference, Geometry, Density, Radius )

    class ( FieldSetForm ), intent ( inout ) :: &
      Fluid, &
      Reference
    class ( Geometry_F_Form ), intent ( in ) :: &
      Geometry
    real ( KDR ), intent ( in ) :: &
      Density, &
      Radius

    !-- Geometry

    associate &
      ( GV  =>  Geometry % Storage_GS % Value )
    associate &
      ( R_E  =>  GV ( :, G % EDGE_I_U ( 1 ) ), &
        R_W  =>  GV ( :, G % WIDTH_U  ( 1 ) ), &
        R_C  =>  GV ( :, G % CENTER_U ( 1 ) ) )

    !-- Fluid

    associate &
      ( FV  =>  Fluid % Storage_GS % Value )
    associate &
      ( M  =>  FV ( :, 1 ), &
        D  =>  FV ( :, 2 ) )

    call SetDensityKernel ( R_E, R_W, Radius, Density, M, D )

    end associate !-- M, etc.
    end associate !-- FV

    !-- Reference

    associate &
      ( RV  =>  Reference % Storage_GS % Value )
    associate &
      ( Phi      =>  RV ( :, 1 ), &
        FourPi   =>  4.0_KDR * CONSTANT % PI )

    where ( R_C  <  Radius )
      Phi  =  1.0_KDR / 6.0_KDR  *  FourPi  *  Density  *  R_C ** 2  &
              -  1.0_KDR / 2.0_KDR  *  FourPi  *  Density  &
                                    *  Radius ** 2
    elsewhere
      Phi  =  - 1.0_KDR / 3.0_KDR  *  FourPi  *  Density  &
                                   *  Radius ** 3  /  R_C
    end where

    end associate !-- Phi, etc.
    end associate !-- RV

    !-- Cleanup

    end associate !-- R_E, etc.
    end associate !-- GV

  end subroutine SetHomogeneousSphere


  subroutine SetDensityKernel ( R_E, R_W, RD, Density, M, D )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_E, &
      R_W
    real ( KDR ), intent ( in ) :: &
      RD, &
      Density
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      M, &
      D

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      R_I, R_O

    do iV  =  1, size ( R_E )

      M ( iV )  =  1.0_KDR

      R_I  =  R_E ( iV )
      R_O  =  R_E ( iV )  +  R_W ( iV )
      if ( R_O  <=  RD ) then
        D ( iV )  =  Density
      else if ( R_I  <  RD .and. R_O  >  RD ) then
        D ( iV )  =  Density * ( RD ** 3  -  R_I ** 3 ) &
                     / ( R_O ** 3  -  R_I ** 3 )
      else
        D ( iV )  =  0.0_KDR
      end if

    end do !-- iV

  end subroutine SetDensityKernel


end program Gravitation_N_SG__Form_Test
