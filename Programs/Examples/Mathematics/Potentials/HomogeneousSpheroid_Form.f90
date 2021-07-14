module HomogeneousSpheroid_Form

  use Basics
  use Mathematics

  implicit none
  private

  type, public :: HomogeneousSpheroidForm
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( Atlas_SCG_CC_Form ), allocatable :: &
      Atlas
    type ( FieldSetForm ), allocatable :: &
      Source, &
      Solution, &
      Reference, &
      Difference
    type ( StreamForm ), allocatable :: &
      Stream
    type ( Geometry_F_Form ), allocatable :: &
      Geometry
    type ( GradientForm ), allocatable :: &
      GradSolution
    type ( Poisson_ASCG_Form ), allocatable :: &
      Poisson
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, private, pass :: &
      SetHomogeneousSpheroids
  end type HomogeneousSpheroidForm

    private :: &
      SetHomogeneousSpheroidKernel, &
      ComputeError


contains


  subroutine Initialize ( HS )

    class ( HomogeneousSpheroidForm ), intent ( inout ) :: &
      HS

  integer ( KDI ) :: &
    nEquations, &
    MaxDegree
  real ( KDR ) :: &
    Density
  real ( KDR ), dimension ( : ), allocatable :: &
    SemiMajor, &
    Eccentricity
  character ( LDL ), dimension ( : ), allocatable :: &
    Field

    call Show ( 'Initializing a HomogeneousSpheroid' )

    allocate ( HS % GridImageStream )
    associate ( GIS  =>  HS % GridImageStream )
    call GIS % Initialize &
           ( PROGRAM_HEADER % Name, &
             CommunicatorOption = PROGRAM_HEADER % Communicator )

    allocate ( HS % Atlas )
    associate ( A  =>  HS % Atlas )
    call A % Initialize &
           ( RadiusMax = 10.0_KDR, &
             RadiusCore = 10.0_KDR / 8.0_KDR, &
             CommunicatorOption = PROGRAM_HEADER % Communicator )

    allocate ( HS % Stream )
    associate ( S  =>  HS % Stream )
    call S % Initialize ( A, GIS )

    allocate ( HS % Geometry )
    associate ( G  =>  HS % Geometry )
    call G % Initialize ( A )
    call G % SetStream ( S )

    nEquations = 3

    MaxDegree = 12
    call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

    allocate ( HS % Poisson )
    associate ( P  =>  HS % Poisson )
    call P % Initialize ( G, 'MULTIPOLE', MaxDegree, nEquations )

    allocate ( Field ( nEquations ) )
    Field  =  [ 'OblateSpheroid_1', &
                'OblateSpheroid_2', &
                'ProlateSpheroid ' ]

    allocate ( HS % Source )
    call HS % Source % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Source', &
             DeviceMemoryOption = P % Laplacian_M % DeviceMemory, &
             nFieldsOption = nEquations )
    
    allocate ( HS % Solution )
    call HS % Solution % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Solution', &
             DeviceMemoryOption = P % Laplacian_M % DeviceMemory, &
             nFieldsOption = nEquations )
    call HS % Solution % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
    call HS % Solution % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    call HS % Solution % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    
    allocate ( HS % Reference )
    call HS % Reference % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Reference', &
             nFieldsOption = nEquations )
    
    allocate ( HS % Difference )
    call HS % Difference % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Difference', &
             nFieldsOption = nEquations )

    call S % AddFieldSet ( HS % Source )
    call S % AddFieldSet ( HS % Solution )
    call S % AddFieldSet ( HS % Reference )
    call S % AddFieldSet ( HS % Difference )

    allocate ( HS % GradSolution )
    associate ( GS  =>  HS % GradSolution )
    call GS % Initialize ( G, HS % Solution )
    call GS % SetStream ( S )

    call A % Show ( )
    call G % Show ( )
    call P % Show ( )
    call S % Show ( )

    call Show ( 'Spheroid Parameters' )

    associate ( C  =>  A % Chart_GS )
    allocate ( SemiMajor ( nEquations ) )
    SemiMajor  =  C % MaxCoordinate ( 1 ) / 2
    call PROGRAM_HEADER % GetParameter ( SemiMajor, 'SemiMajor' )
    call Show ( SemiMajor, 'SemiMajor' )
    end associate !-- C

    allocate ( Eccentricity ( nEquations ) )
    Eccentricity  =  sqrt ( 1.0_KDR &
                            - [ 0.7_KDR ** 2, 0.2_KDR ** 2, 0.2_KDR ** 2 ] )
    call PROGRAM_HEADER % GetParameter ( Eccentricity, 'Eccentricity' )
    call Show ( Eccentricity, 'Eccentricity' )

    Density  =  1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  )
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )
    call Show ( Density, 'Density' )

    call HS % SetHomogeneousSpheroids &
           ( SemiMajor, Eccentricity, Density, nEquations )

    end associate !-- GS
    end associate !-- P
    end associate !-- G
    end associate !-- S
    end associate !-- A
    end associate !-- GIS

  end subroutine Initialize


  subroutine Compute ( HS )

    class ( HomogeneousSpheroidForm ), intent ( inout ) :: &
      HS

    integer ( KDI ) :: &
      iD

    associate ( P  =>  HS % Poisson )
    call P % Solve ( HS % Solution, HS % Source )
    end associate !-- P

    call ComputeError ( HS % Difference, HS % Solution, HS % Reference )

    associate ( nD  =>  HS % Atlas % Chart ( 1 ) % Element % nDimensions )
    do iD  =  1, nD
      associate ( GS  =>  HS % GradSolution )
      call GS % Compute ( iD )
      end associate !-- GS
    end do !-- iD
    end associate !-- nD

    associate &
      ( GIS  =>  HS % GridImageStream, &
          S  =>  HS % Stream )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call S % Write ( )
    call GIS % Close ( )
    end associate !-- GIS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( HS )

    type ( HomogeneousSpheroidForm ), intent ( inout ) :: &
      HS

    if ( allocated ( HS % Poisson ) ) &
      deallocate ( HS % Poisson )
    if ( allocated ( HS % GradSolution ) ) &
      deallocate ( HS % GradSolution )
    if ( allocated ( HS % Geometry ) ) &
      deallocate ( HS % Geometry )
    if ( allocated ( HS % Stream ) ) &
      deallocate ( HS % Stream )
    if ( allocated ( HS % Difference ) ) &
      deallocate ( HS % Difference )
    if ( allocated ( HS % Reference ) ) &
      deallocate ( HS % Reference )
    if ( allocated ( HS % Solution ) ) &
      deallocate ( HS % Solution )
    if ( allocated ( HS % Source ) ) &
      deallocate ( HS % Source )
    if ( allocated ( HS % Atlas ) ) &
      deallocate ( HS % Atlas )
    if ( allocated ( HS % GridImageStream ) ) &
      deallocate ( HS % GridImageStream )

    call Show ( 'Finalizing a HomogeneousSpheroid' )

  end subroutine Finalize


  subroutine SetHomogeneousSpheroids &
               ( HS, SemiMajor, Eccentricity, Density, nEquations )

    class ( HomogeneousSpheroidForm ), intent ( inout ) :: &
      HS
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      SemiMajor, &
      Eccentricity
    real ( KDR ), intent ( in ) :: &
      Density
    integer ( KDI ), intent ( in ) :: &
      nEquations

    integer ( KDI ) :: &
      iE !-- iEquation
    real ( KDR ) :: &
      ProlateMajor
    real ( KDR ), dimension ( nEquations ) :: &
      SemiMinor

    SemiMinor  =  sqrt ( 1.0_KDR - Eccentricity ** 2 ) * SemiMajor

    !-- Prolate Spheroid 
    ProlateMajor              =  SemiMajor ( nEquations )
    SemiMajor ( nEquations )  =  SemiMinor ( nEquations )
    SemiMinor ( nEquations )  =  ProlateMajor

    do iE  =  1, nEquations
      call SetHomogeneousSpheroidKernel &
             ( HS % Source, HS % Reference, HS % Geometry, &
               Density, SemiMajor ( iE), SemiMinor ( iE ), iE )
    end do

    call HS % Source % UpdateDevice ( )

  end subroutine SetHomogeneousSpheroids


  subroutine SetHomogeneousSpheroidKernel &
               ( Source, Reference, G, Density, a_1, a_3, &
                 iField )

    class ( FieldSetForm ), intent ( inout ) :: &
      Source, &
      Reference
    class ( Geometry_F_Form ), intent ( in ) :: &
      G
    real ( KDR ), intent ( in ) :: &
      Density, &
      a_1, &
      a_3
    integer ( KDI ), intent ( in ) :: &
      iField

    integer ( KDI ) :: &
      iC, &       !-- iCell
      iS, jS, kS  !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ) :: &
      dVS, &  !-- dVolumeSubcell
       VS, &   !-- VolumeSubcell
      e, &
      a_1_sq, &
      a_3_sq, &
      C_I, &
      C_A, &
      C_B, &
      rho_sq_II, &
      rho_sq_IO, &
      rho_sq_OI, &
      rho_sq_OO, &
      Z_sq_II, &
      Z_sq_IO, &
      Z_sq_OI, &
      Z_sq_OO, &
      rho_S_sq, &
      Z_s_sq, &
      Pi, &
      FourPi
    real ( KDR ), dimension ( 3 ) :: &
      X_I, &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS,  &  !--  X_Subcell
      X, &
      dX
    real ( KDR ), dimension ( : ), allocatable :: &
      VF, &  !-- VolumeFraction
      rho_sq, & !-- X^2 + Y^2 
      Z_sq, &
      l, &
      C_I_vec

        Pi  =  CONSTANT % PI
    FourPi  =  4.0_KDR * Pi

    !-- Geometry

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        GV  =>  G % Storage_GS % Value )
    associate &
      ( dV  =>  GV ( :, G % VOLUME ) )

    allocate ( VF ( size ( dV ) ) )
    call Clear ( VF )

    nSubCells  =  1
    nSubcells ( : C % nDimensions ) = 20

    a_1_sq  =  a_1 ** 2
    a_3_sq  =  a_3 ** 2

    do iC  =  1, size ( dV )

      if ( .not. C % ProperCell ( iC ) ) cycle

       X  =  GV ( iC, G % CENTER_U_1 : G % CENTER_U_3 )
      dX  =  GV ( iC, G %  WIDTH_U_1 : G %  WIDTH_U_3 )

      X_I  =  X  -  0.5_KDR * dX
      X_O  =  X  +  0.5_KDR * dX

      rho_sq_II  =  ( X_I ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
      rho_sq_IO  =  ( X_O ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
      rho_sq_OI  =  ( X_I ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2
      rho_sq_OO  =  ( X_O ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2

      Z_sq_II  =  ( X_I ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
      Z_sq_IO  =  ( X_O ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
      Z_sq_OI  =  ( X_I ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2
      Z_sq_OO  =  ( X_O ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2

      if (       rho_sq_II / a_1_sq  +  Z_sq_II / a_3_sq  <=  1.0_KDR &
           .and. rho_sq_IO / a_1_sq  +  Z_sq_IO / a_3_sq  <=  1.0_KDR &
           .and. rho_sq_OI / a_1_sq  +  Z_sq_OI / a_3_sq  <=  1.0_KDR &
           .and. rho_sq_OO / a_1_sq  +  Z_sq_OO / a_3_sq  <=  1.0_KDR ) &
      then 
        VF ( iC ) = 1.0_KDR
        cycle
      end if

      if (       rho_sq_II / a_1_sq  +  Z_sq_II / a_3_sq  >  1.0_KDR &
           .and. rho_sq_IO / a_1_sq  +  Z_sq_IO / a_3_sq  >  1.0_KDR &
           .and. rho_sq_OI / a_1_sq  +  Z_sq_OI / a_3_sq  >  1.0_KDR &
           .and. rho_sq_OO / a_1_sq  +  Z_sq_OO / a_3_sq  >  1.0_KDR ) &
      then 
        VF ( iC ) = 0.0_KDR
        cycle
      end if

      dXS  =  ( X_O  -  X_I ) / nSubcells

      VS = 0.0_KDR
      do kS = 1, nSubcells ( 3 )
        do jS = 1, nSubcells ( 2 )
          do iS = 1, nSubcells ( 1 )

            XS  =  X_I  +  ( [ iS, jS, kS ] - 0.5_KDR ) * dXS
            rho_S_sq  =  ( XS ( 1 ) * sin ( XS ( 2 ) ) ) ** 2
              Z_S_sq  =  ( XS ( 1 ) * cos ( XS ( 2 ) ) ) ** 2

            select case ( C % nDimensions )
            case ( 2 )
              dVS  =  2 * Pi * XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                      * dXS ( 1 ) * dXS ( 2 )
            case ( 3 )
              dVS  =  XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                      * dXS ( 1 ) * dXS ( 2 ) * dXS ( 3 )
            end select !-- nDimensions

            VS  =  VS + dVS
            if ( rho_S_sq / a_1_sq  +  Z_S_sq / a_3_sq  <=  1.0_KDR ) &
              VF ( iC )  =  VF ( iC )  +  dVS

          end do !-- iS
        end do !-- jS
      end do !-- kS
      VF ( iC )  =  VF ( iC )  /  VS

    end do !-- iC

    associate &
      (     R  =>  GV ( :, G % CENTER_U ( 1 ) ), &
        Theta  =>  GV ( :, G % CENTER_U ( 2 ) ) )

    allocate &
      ( rho_sq  ( size ( R ) ) , &
        Z_sq    ( size ( R ) ), &
        l       ( size ( R ) ), &
        C_I_vec ( size ( R ) ) )
         
    rho_sq  =  ( R * sin ( Theta ) ) ** 2
      Z_sq  =  ( R * cos ( Theta ) ) ** 2

    end associate !-- R, etc.
    end associate !-- dV
    end associate !-- C, etc.
    end select !-- A

    !-- Source

    associate &
      ( SV  =>  Source % Storage_GS % Value )
    associate &
      ( D  =>  SV ( :, iField ) )

    D  =  Density * VF * FourPi

    end associate !-- D
    end associate !-- SV

    !-- Reference

    e  =  sqrt ( 1.0_KDR  -  min ( ( a_3 / a_1 ), ( a_1 / a_3 ) )  ** 2 )

    if ( a_3 < a_1 ) then
      C_I = 2 * sqrt ( 1.0_KDR - e ** 2 ) / e * asin ( e )
      C_A = sqrt ( 1.0_KDR - e ** 2 ) / e ** 3 * asin ( e ) &
            - ( 1.0_KDR - e ** 2 ) / e ** 2
      C_B = 2.0_KDR / e ** 2 &
            - 2 * sqrt ( 1.0_KDR - e ** 2 ) / e ** 3 * asin ( e )
    else
      C_I  = 1.0_KDR / e * log ( ( 1 + e ) / ( 1 - e ) )
      C_A = 1.0_KDR / e ** 2 &
            - ( 1.0_KDR - e ** 2 ) / ( 2 * e ** 3 ) &
                * log ( ( 1 + e ) / ( 1 - e ) )
      C_B = ( 1.0_KDR - e ** 2 ) / e ** 3  * log ( ( 1 + e ) / ( 1 - e ) ) &
            - 2 * ( 1.0_KDR - e ** 2 ) / e ** 2
    end if

    associate &
      ( RV  =>  Reference % Storage_GS % Value )
    associate &
      ( Phi  =>  RV ( :, iField ) )

    where ( rho_sq / a_1_sq + Z_sq / a_3_sq <= 1.0_KDR )
      Phi  =   - Density * Pi &
                 * ( C_I * a_1_sq  - C_A * rho_sq - C_B * Z_sq )
    end where

    l = 0.0_KDR

    where ( rho_sq / a_1_sq + Z_sq / a_3_sq > 1.0_KDR )
      l = 0.5_KDR * ( ( rho_sq + Z_sq - a_1_sq - a_3 **2  ) &
            + sqrt ( ( a_1_sq + a_3_sq - Z_sq - rho_sq) ** 2 &
                     - 4.0_KDR * ( a_1_sq * a_3_sq - rho_sq * a_3_sq &
                                   - Z_sq * a_1_sq ) ) )
    end where

    where ( l < 0.0_KDR )
      l = 0.5_KDR * ( ( rho_sq + Z_sq - a_1_sq - a_3 **2  ) &
          - sqrt ( ( a_1_sq + a_3_sq - Z_sq - rho_sq) ** 2 &
                   - 4.0_KDR * ( a_1_sq * a_3_sq - rho_sq * a_3_sq &
                                 - Z_sq * a_1_sq ) ) )
    end where

    if ( a_3 < a_1 ) then
      C_I_vec = Pi / sqrt ( a_1_sq - a_3_sq ) &
            - 2.0_KDR / sqrt ( a_1_sq - a_3_sq ) &
            * atan ( sqrt ( ( a_3_sq + l ) / ( a_1_sq - a_3_sq ) ) )
    else
      C_I_vec = -1.0_KDR / sqrt ( a_3_sq - a_1_sq ) &
                * log ( ( sqrt ( a_3_sq + l ) &
                           - sqrt ( a_3_sq - a_1_sq ) ) ** 2 &
                             / ( a_1_sq + l ) )
    end if

    where (  rho_sq / a_1_sq + Z_sq / a_3_sq > 1.0_KDR )
      Phi  = - Density * Pi * a_1_sq * a_3 &
               * ( ( 1.0_KDR + rho_sq / ( 2 * ( a_3_sq - a_1_sq ) ) &
                   - Z_sq / ( a_3_sq - a_1_sq ) ) * C_I_vec &
             - rho_sq * sqrt ( a_3_sq + l ) &
               / ( ( a_3_sq - a_1_sq ) * ( a_1_sq + l ) ) &
             - Z_sq * ( 2.0_KDR &
                        / ( ( a_1_sq + l ) * sqrt ( a_3_sq + l ) ) &
                        - 2.0_KDR  * sqrt ( a_3_sq + l ) &
                          / ( ( a_3_sq - a_1_sq ) * ( a_1_sq + l ) ) ))
    end where

    end associate !-- Phi
    end associate !-- RV

  end subroutine SetHomogeneousSpheroidKernel


  subroutine ComputeError ( Difference, Solution, Reference )

    class ( FieldSetForm ), intent ( inout ) :: &
      Difference, &         
      Solution, &
      Reference

    integer ( KDI ) :: &
      nEquations
    real ( KDR ) :: &
      L1_1, &
      L1_2, &
      L1_3
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    nEquations  =  Solution % nFields

    call Difference % MultiplyAdd ( Solution, Reference, -1.0_KDR )

    select type ( A  =>  Difference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        RV  =>  Reference % Storage_GS % Value, &
        DV  =>  Difference % Storage_GS % Value )

    call CO % Initialize &
           ( C % Communicator, [ 2 * nEquations ], [ 2 * nEquations ] )
    CO % Outgoing % Value ( 1 )  &
      =  sum ( abs ( DV ( :, 1 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 )  &
      =  sum ( abs ( DV ( :, 2 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 3 )  &
      =  sum ( abs ( DV ( :, 3 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 4 )  &
      =  sum ( abs ( RV ( :, 1 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 5 )  &
      =  sum ( abs ( RV ( :, 2 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 6 )  &
      =  sum ( abs ( RV ( :, 3 ) ), mask = C % ProperCell )
    call CO % Reduce ( REDUCTION % SUM )

    associate &
      ( Norm_D_1  =>  CO % Incoming % Value ( 1 ), &
        Norm_D_2  =>  CO % Incoming % Value ( 2 ), &
        Norm_D_3  =>  CO % Incoming % Value ( 3 ), &
        Norm_R_1  =>  CO % Incoming % Value ( 4 ), &
        Norm_R_2  =>  CO % Incoming % Value ( 5 ), &
        Norm_R_3  =>  CO % Incoming % Value ( 6 ) )

    L1_1  =  Norm_D_1 / Norm_R_1
    L1_2  =  Norm_D_2 / Norm_R_2
    L1_3  =  Norm_D_3 / Norm_R_3

    end associate !-- Norm_D_1, etc.

    call Show ( L1_1, '*** L1_1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_2, '*** L1_2 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_3, '*** L1_3 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    ! Difference % Value = abs ( Difference % Value / Reference % Value )

    end associate !-- C, etc.
    end select !-- A

  end subroutine ComputeError


end module HomogeneousSpheroid_Form
