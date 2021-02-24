!-- Geometry_ASC represents the geometry on a position space represented
!   by a single chart.

module Geometry_ASC__Form

  !-- Geometry_AtlasSingleChart_Form

  use Basics
  use Mathematics
  use Geometry_G__Form
  use Geometry_N__Form
  use Geometry_CSL__Form

  implicit none
  private

  type, public, extends ( GeometryFlat_ASC_Form ) :: Geometry_ASC_Form
    real ( KDR ) :: &
      GravitationalConstant, &
      UniformAcceleration, &
      CentralMass
    type ( MeasuredValueForm ) :: &
      CentralMassUnit
    character ( LDL ) :: &
      GravitySolverType = ''
    type ( Storage_ASC_Form ), allocatable :: &
      Storage_ASC
    type ( GradientForm ), allocatable :: &
      Gradient
    type ( Poisson_ASC_Form ), allocatable :: &
      Poisson_ASC
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AllocateDevice => AllocateDevice_G_ASC 
    procedure, public, pass :: &
      UpdateHost => UpdateHost_G_ASC 
    procedure, private, pass :: &
      Geometry_G_CSL
    generic, public :: &
      Geometry_G => Geometry_G_CSL
    procedure, private, pass :: &
      Geometry_N_CSL
    generic, public :: &
      Geometry_N => Geometry_N_CSL
    procedure, public, pass :: &
      ComputeGravity
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Geometry_ASC_Form

    private :: &
      ComputeGravityUniform, &
      ComputeGravityCentralMass, &
      ComputeGravityMultipole

      private :: &
        ComputeGravityUniformKernel, &
        ComputeGravityCentralMassKernel, &
        ComputeGravitySource
  
  interface
  
    module subroutine ComputeGravityUniformKernel &
             ( Phi, GradPhi_1, GradPhi_2, GradPhi_3, X, Y, Z, &
               Acceleration, nDimensions, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Phi, &
        GradPhi_1, GradPhi_2, GradPhi_3
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        X, Y, Z
      real ( KDR ), intent ( in ) :: &
        Acceleration
      integer ( KDI ), intent ( in ) :: &
        nDimensions
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeGravityUniformKernel


    module subroutine ComputeGravityCentralMassKernel &
             ( Phi, GradPhi_1, GradPhi_2, GradPhi_3, X_1, X_2, X_3, G, M, &
               UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Phi, &
        GradPhi_1, GradPhi_2, GradPhi_3
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        X_1, X_2, X_3
      real ( KDR ), intent ( in ) :: &
        G, &
        M
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeGravityCentralMassKernel


    module subroutine ComputeGravitySource &
             ( M, N, G, S, UseDeviceOption )
      use Basics
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N
      real ( KDR ), intent ( in ) :: &
        G
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        S
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeGravitySource

  end interface

contains


  subroutine Initialize &
               ( GA, A, GeometryType, NameShortOption, &
                 GravitySolverTypeOption, UsePinnedMemoryOption, &
                 CentralMassUnitOption, GravitationalConstantOption, &
                 UniformAccelerationOption, CentralMassOption, &
                 IgnorabilityOption )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      GeometryType
    character ( * ), intent ( in ), optional :: &
      NameShortOption, &
      GravitySolverTypeOption
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      CentralMassUnitOption
    real ( KDR ), intent ( in ), optional :: &
      GravitationalConstantOption, &
      UniformAccelerationOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      NameShort
    class ( Geometry_N_Form ), pointer :: &
      G

    if ( GA % Type == '' ) &
      GA % Type = 'a Geometry_ASC'

    GA % GeometryType = GeometryType    

    call GA % InitializeFlat &
          ( A, NameShortOption = NameShortOption, &
            UsePinnedMemoryOption = UsePinnedMemoryOption, &
            IgnorabilityOption = IgnorabilityOption )

    select case ( trim ( GA % GeometryType ) )
    case ( 'NEWTONIAN', 'NEWTONIAN_STRESS' )

      GA % GravitationalConstant  =  CONSTANT % GRAVITATIONAL
      if ( present ( GravitationalConstantOption ) ) &
        GA % GravitationalConstant  =  GravitationalConstantOption

      if ( present ( GravitySolverTypeOption ) ) then
        GA % GravitySolverType = GravitySolverTypeOption
        call Show ( GA % GravitySolverType, 'GravitySolverType', &
                    GA % IGNORABILITY )
      else
        call Show ( 'NEWTONIAN geometry requires GravitySolverTypeOption', &
                    CONSOLE % ERROR )
        call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      select case ( trim ( GA % GravitySolverType ) )
      case ( 'UNIFORM' )

        if ( present ( UniformAccelerationOption ) ) then
          GA % UniformAcceleration = UniformAccelerationOption
          call Show ( GA % UniformAcceleration, 'UniformAcceleration', &
                      GA % IGNORABILITY )
        else
          call Show ( 'GravitySolverType UNIFORM requires ' &
                      // 'UniformAccelerationOption', CONSOLE % ERROR )
          call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
          call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end if

      case ( 'CENTRAL_MASS' )

        GA % GravitationalConstant  =  CONSTANT % GRAVITATIONAL
        if ( present ( GravitationalConstantOption ) ) &
          GA % GravitationalConstant  =  GravitationalConstantOption

        if ( present ( CentralMassUnitOption ) ) &
          GA % CentralMassUnit  =  CentralMassUnitOption

        if ( present ( CentralMassOption ) ) then
          GA % CentralMass = CentralMassOption
          call Show ( GA % CentralMass, GA % CentralMassUnit, &
                      'CentralMass', GA % IGNORABILITY )
        else
          call Show ( 'GravitySolverType CENTRAL_MASS requires ' &
                      // 'CentralMassOption', CONSOLE % ERROR )
          call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
          call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end if

      case ( 'MULTIPOLE', 'MULTIPOLE_OLD_2', 'MULTIPOLE_OLD_1' )

        G => GA % Geometry_N ( )
        
        allocate ( GA % Storage_ASC )
        associate ( SA => GA % Storage_ASC )
        call SA % Initialize &
               ( GA, NameShort = 'PoissonStorage', &
                 iaSelectedOption = [ G % POTENTIAL ] )
        end associate !-- SA

        allocate ( GA % Gradient )
        associate ( Grad => GA % Gradient )
        call Grad % Initialize &
               ( 'GeometryGradient', [ G % nValues, 1 ] )
        end associate !-- Grad

      case default
        call Show ( 'GravitySolverType not recognized', CONSOLE % ERROR )
        call Show ( GA % GravitySolverType, 'GravitySolverType', &
                    CONSOLE % ERROR ) 
        call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select 

      nullify ( G )

    end select !-- GeometryType

  end subroutine Initialize
  
  
  subroutine AllocateDevice_G_ASC ( FA, AssociateVariablesOption )
  
    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      FA
    logical ( KDL ), intent ( in ), optional :: &
      AssociateVariablesOption
    
    call FA % AllocateDevice_ASC_Template ( AssociateVariablesOption )
    
    if ( allocated ( FA % Storage_ASC ) ) &
      call FA % Storage_ASC % AllocateDevice ( )
    
    if ( allocated ( FA % Gradient ) ) &
      call FA % Gradient % AllocateDevice ( )
    
  end subroutine AllocateDevice_G_ASC 


  subroutine UpdateHost_G_ASC ( FA )
  
    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      FA
    
    call FA % UpdateHost_ASC_Template ( )
    
    if ( allocated ( FA % Storage_ASC ) ) &
      call FA % Storage_ASC % UpdateHost ( )
    
    if ( allocated ( FA % Gradient ) ) &
      call FA % Gradient % UpdateHost ( )
    
  end subroutine UpdateHost_G_ASC


  function Geometry_G_CSL ( GA ) result ( G )

    class ( Geometry_ASC_Form ), intent ( in ) :: &
      GA
    class ( Geometry_G_Form ), pointer :: &
      G

    select type ( GC => GA % Chart )
    class is ( Geometry_CSL_Form )
      G => GC % Geometry_G ( )
    class default
      call Show ( 'Geometry Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Geometry_G_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GC

  end function Geometry_G_CSL


  function Geometry_N_CSL ( GA ) result ( G )

    class ( Geometry_ASC_Form ), intent ( in ) :: &
      GA
    class ( Geometry_N_Form ), pointer :: &
      G

    select type ( GC => GA % Chart )
    class is ( Geometry_CSL_Form )
      G => GC % Geometry_N ( )
    class default
      call Show ( 'Geometry Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Geometry_N_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- GC

  end function Geometry_N_CSL


  subroutine ComputeGravity &
               ( GA, FA, iBaryonMass, iBaryonDensity, iPressureOption )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Current_ASC_Template ), intent ( in ) :: &
      FA
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass, &
      iBaryonDensity
    integer ( KDI ), intent ( in ), optional :: &
      iPressureOption

    select case ( trim ( GA % GeometryType ) )
    case ( 'NEWTONIAN', 'NEWTONIAN_STRESS' )

    select case ( trim ( GA % GravitySolverType ) )
    case ( 'UNIFORM' )
      call ComputeGravityUniform ( GA )
    case ( 'CENTRAL_MASS' )
      call ComputeGravityCentralMass ( GA )
    case ( 'MULTIPOLE', 'MULTIPOLE_OLD_2', 'MULTIPOLE_OLD_1' )
      call ComputeGravityMultipole &
             ( GA, FA, iBaryonMass, iBaryonDensity, iPressureOption )
    case default
      call Show ( 'GravitySolverType not recognized', CONSOLE % ERROR )
      call Show ( GA % GravitySolverType, 'GravitySolverType', &
                  CONSOLE % ERROR ) 
      call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeGravity', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    end select !-- GeometryType

  end subroutine ComputeGravity


  impure elemental subroutine Finalize ( GA )

    type ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA

    if ( allocated ( GA % Poisson_ASC ) ) &
      deallocate ( GA % Poisson_ASC )
    if ( allocated ( GA % Gradient ) ) &
      deallocate ( GA % Gradient )
    if ( allocated ( GA % Storage_ASC ) ) &
      deallocate ( GA % Storage_ASC )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Template )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
    associate ( nValues => C % nProperCells + C % nGhostCells )

    allocate ( Geometry_CSL_Form :: FA % Chart )

    select type ( FC => FA % Chart )
    class is ( Geometry_CSL_Form )
      call FC % Initialize &
             ( C, FA % NameShort, FA % GeometryType, FA % UsePinnedMemory, &
               nValues, IgnorabilityOption = FA % IGNORABILITY )
    end select !-- FC

    end associate !-- nValues
    end select !-- C
    end select !-- A

  end subroutine SetField


  subroutine ComputeGravityUniform ( GA )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA

    class ( Geometry_N_Form ), pointer :: &
      G

    select type ( A => GA % Atlas )
    class is ( Atlas_SC_Form )

    G  =>  GA % Geometry_N ( )

    call ComputeGravityUniformKernel &
           ( G % Value ( :, G % POTENTIAL ), & 
             G % Value ( :, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
             G % Value ( :, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
             G % Value ( :, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
             G % Value ( :, G % CENTER_U ( 1 ) ), &
             G % Value ( :, G % CENTER_U ( 2 ) ), &
             G % Value ( :, G % CENTER_U ( 3 ) ), &
             GA % UniformAcceleration, A % nDimensions, &
             UseDeviceOption = G % AllocatedDevice )

    end select !-- A
    nullify ( G )

  end subroutine ComputeGravityUniform


  subroutine ComputeGravityCentralMass ( GA )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA

    class ( Geometry_N_Form ), pointer :: &
      G

    select type ( A => GA % Atlas )
    class is ( Atlas_SC_Form )

    G  =>  GA % Geometry_N ( )
    
    select case ( trim ( G % CoordinateSystem ) )
    case ( 'SPHERICAL' )

      call ComputeGravityCentralMassKernel &
             ( G % Value ( :, G % POTENTIAL ), & 
               G % Value ( :, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
               G % Value ( :, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
               G % Value ( :, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
               G % Value ( :, G % CENTER_U ( 1 ) ), &
               G % Value ( :, G % CENTER_U ( 2 ) ), &
               G % Value ( :, G % CENTER_U ( 3 ) ), &
               GA % GravitationalConstant, GA % CentralMass, &
               UseDeviceOption = G % AllocatedDevice )
    
    case default
      call Show ( 'CoordinateSystem not implemented', CONSOLE % ERROR )
      call Show ( G % CoordinateSystem, 'CoordinateSystem', &
                  CONSOLE % ERROR )
      call Show ( 'Geometry_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeGravityCentralMassKernel', 'subroutine', &
                   CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    
    end select !-- CoordinateSystem  

    end select !-- A
    nullify ( G )

  end subroutine ComputeGravityCentralMass


  subroutine ComputeGravityMultipole &
               ( GA, FA, iBaryonMass, iBaryonDensity, iPressureOption )

    class ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Current_ASC_Template ), intent ( in ) :: &
      FA
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass, &
      iBaryonDensity
    integer ( KDI ), intent ( in ), optional :: &
      iPressureOption

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      MaxDegree
    logical ( KDL ) :: &
      ComputeForce
    type ( StorageForm ), pointer :: &
      S, &
      S_P
    type ( Storage_ASC_Form ) :: &
      Storage_ASC_Pressure
    type ( GradientForm ) :: &
      GradientPressure
    class ( Geometry_N_Form ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      F

    select type ( A => GA % Atlas )
    class is ( Atlas_SC_Form )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )

    if ( .not. allocated ( GA % Poisson_ASC ) ) then

        MaxDegree = 10
        call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

        allocate ( GA % Poisson_ASC )
        associate ( PA => GA % Poisson_ASC )
        call PA % Initialize &
               ( A, SolverType = GA % GravitySolverType, &
                 MaxDegreeOption = MaxDegree, &
                 nEquationsOption = 1 )
        end associate !-- PA

    end if

    associate &
      ( PA => GA % Poisson_ASC, &
        SA => GA % Storage_ASC )

    S  =>  SA % Storage ( )
    G  =>  GA % Geometry_N ( )
    F  =>  FA % Current ( )

    call ComputeGravitySource &
           ( F % Value ( :, iBaryonMass ), &
             F % Value ( :, iBaryonDensity ), &
             GA % GravitationalConstant, &
             S % Value ( :, S % iaSelected ( 1 ) ), &
             UseDeviceOption = F % AllocatedDevice )
    
    call PA % Solve ( SA, SA )
    
    do iD = 1, C % nDimensions
      call GA % Gradient % Compute ( C, S, iDimension = iD )
      call Copy ( GA % Gradient % Output % Value ( :, 1 ), &
                  G % Value ( :, G % POTENTIAL_GRADIENT_D ( iD ) ), &
                  UseDeviceOption = G % AllocatedDevice )
    end do !-- iD

    if ( present ( iPressureOption ) ) then

      !-- Pressure gradient

      associate &
        ( SA_P => Storage_ASC_Pressure, &
          Grad_P => GradientPressure )

      call SA_P % Initialize &
             ( FA, NameShort = 'PressureStorage', &
               iaSelectedOption = [ iPressureOption ], &
               IgnorabilityOption = CONSOLE % INFO_3 )

      call Grad_P % Initialize &
             ( 'PressureGradient', [ F % nValues, 1 ] )

      S_P  =>  SA_P % Storage ( )

      do iD = 1, C % nDimensions
        call Grad_P % Compute ( C, S_P, iDimension = iD )
        call Copy ( Grad_P % Output % Value ( :, 1 ), &
                    G % Value ( :, G % PRESSURE_GRADIENT_D ( iD ) ), &
                    UseDeviceOption = G % AllocatedDevice )
      end do !-- iD

      end associate !-- SA_P, etc.

    end if !-- iPressureOption

    end associate !-- PA, etc.v
    end select !-- C
    end select !-- A
    nullify ( S, S_P, G, F )

  end subroutine ComputeGravityMultipole


end module Geometry_ASC__Form
