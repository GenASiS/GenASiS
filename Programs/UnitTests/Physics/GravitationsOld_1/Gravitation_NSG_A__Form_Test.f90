program Gravitation_NSG_A__Form_Test

  !-- Gravitation_NewtonSelfGravity_Atlas__Form_Test

  use Basics
  use Mathematics
  use Gravitations

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Gravitation_NSG_A_Form ), allocatable :: &
    GA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Gravitation_NSG_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )
  call GA % SetStream ( SA )
  call SA % AddFieldSet ( GA % Source_A )

  call  A % Show ( )
  call GA % Show ( )

  call TestHomogeneousSpheres ( )

  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestHomogeneousSpheres ( )

    integer ( KDI ) :: &
      iHS
    real ( KDR ), dimension ( 3 ) :: &
      RadiusDensity, &
      Density
    type ( FieldSet_A_Form ), allocatable :: &
      Fluid_A, &
      Reference_A, &
      Difference_A

    call Show ( 'Testing homogeneous spheres' )

    associate ( C  =>  A % Chart_GS )
    
    allocate ( Fluid_A )
    call Fluid_A % Initialize &
           ( A, &
             FieldOption = [ 'BaryonMass   ', 'BaryonDensity' ], &
             NameOption = 'Fluid', &
             nFieldsOption = 2 )
    
    allocate ( Reference_A )
    call Reference_A % Initialize &
           ( A, &
             FieldOption = [ 'Potential' ], &
             NameOption = 'Reference', &
             nFieldsOption = 1 )
    
    allocate ( Difference_A )
    call Difference_A % Initialize &
           ( A, &
             FieldOption = [ 'Potential' ], &
             NameOption = 'Difference', &
             nFieldsOption = 1 )
    
    call SA % AddFieldSet ( Fluid_A )
    call SA % AddFieldSet ( Reference_A )
    call SA % AddFieldSet ( Difference_A )
    call SA % Show ( )

    RadiusDensity  =  C % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ]
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density  =  1.0_KDR  /  RadiusDensity ** 3
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iHS  =  1, 3

      call SetHomogeneousSphere &
             ( Fluid_A, Reference_A, GA, &
               Density ( iHS ), RadiusDensity ( iHS )  )
    
      call Fluid_A % UpdateDevice ( )

      call GA % Solve &
             ( Fluid_A, Constant_G = 1.0_KDR, iBaryonMass = 1, &
               iBaryonDensity = 2 )

      call Show ( RadiusDensity ( iHS ), 'Radius', nLeadingLinesOption = 2 )
      call Show ( Density ( iHS ), 'Density' )
      call ComputeError ( Difference_A, GA % Solution_A, Reference_A )

      call GIS % Open ( GIS % ACCESS_CREATE )
      call SA % Write ( )
      call GIS % Close ( )

    end do !-- iHS

    end associate !-- C

  end subroutine TestHomogeneousSpheres


  subroutine ComputeError ( Difference_A, Solution_A, Reference_A )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Difference_A, &         
      Solution_A, &
      Reference_A

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    associate &
      ( SC  =>  Solution_A   % FieldSet_C ( 1 ) % Element, &
        RC  =>  Reference_A  % FieldSet_C ( 1 ) % Element, &
        DC  =>  Difference_A % FieldSet_C ( 1 ) % Element )
    associate &
      ( SV  =>  SC % Storage_FSC % Storage &
                   % Value ( :, SC % iaSelected ( 1 ) ), &
        RV  =>  RC % Storage_FSC % Storage % Value ( :, 1 ), &
        DV  =>  DC % Storage_FSC % Storage % Value ( :, 1 ) )

    call MultiplyAdd ( SV, RV, -1.0_KDR, DV )

    select type ( C  =>  SC % Chart )
    class is ( Chart_GS_Form ) 

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

    end select !-- C
    end associate !-- SV, etc.
    end associate !-- SC, etc.

  end subroutine ComputeError


  subroutine SetHomogeneousSphere &
               ( Fluid_A, Reference_A, Geometry_A, &
                 Density, RadiusDensity )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Fluid_A, &
      Reference_A
    class ( Geometry_F_A_Form ), intent ( in ) :: &
      Geometry_A
    real ( KDR ), intent ( in ) :: &
      Density, &
      RadiusDensity

    !-- Geometry

    select type ( GC  =>  Geometry_A % FieldSet_C ( 1 ) % Element )
      class is ( Geometry_F_C_Form )
    associate &
      ( GV  =>  GC % Storage_FSC % Storage % Value )
    associate &
      ( R_E  =>  GV ( :, GC % EDGE_I_U ( 1 ) ), &
        R_W  =>  GV ( :, GC % WIDTH_U  ( 1 ) ), &
        R_C  =>  GV ( :, GC % CENTER_U ( 1 ) ) )

    !-- Fluid

    associate &
      ( FC  =>  Fluid_A % FieldSet_C ( 1 ) % Element )
    associate &
      ( FV  =>  FC % Storage_FSC % Storage % Value )
    associate &
      ( M  =>  FV ( :, 1 ), &
        D  =>  FV ( :, 2 ) )

    call SetDensityKernel ( R_E, R_W, RadiusDensity, Density, M, D )

    end associate !-- D
    end associate !-- SV
    end associate !-- SC

    !-- Reference

    associate &
      ( RC  =>  Reference_A % FieldSet_C ( 1 ) % Element )
    associate &
      ( RV  =>  RC % Storage_FSC % Storage % Value )
    associate &
      ( Phi      =>  RV ( :, 1 ), &
        FourPi   =>  4.0_KDR * CONSTANT % PI )

    where ( R_C  <  RadiusDensity )
      Phi  =  1.0_KDR / 6.0_KDR  *  FourPi  *  Density  *  R_C ** 2  &
              -  1.0_KDR / 2.0_KDR  *  FourPi  *  Density  &
                                    *  RadiusDensity ** 2
    elsewhere
      Phi  =  - 1.0_KDR / 3.0_KDR  *  FourPi  *  Density  &
                                   *  RadiusDensity ** 3  /  R_C
    end where

    end associate !-- Phi, etc.
    end associate !-- RV
    end associate !-- RC

    !-- Cleanup

    end associate !-- R_E, etc.
    end associate !-- GV
    end select !-- GC

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


end program Gravitation_NSG_A__Form_Test
