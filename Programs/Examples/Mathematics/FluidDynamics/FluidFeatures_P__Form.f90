module FluidFeatures_P__Form

  !-- FluidFeatures_Perfect__Form

  use Basics
  use Mathematics
  use FluidFeatures_Template
  use Fluid_P__Template

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_PERFECT  = 6, &
      N_VECTORS_PERFECT = 2

  type, public, extends ( FluidFeaturesTemplate ) :: FluidFeatures_P_Form
    integer ( KDI ) :: &
      N_FIELDS_PERFECT  = N_FIELDS_PERFECT, &
      N_VECTORS_PERFECT = N_VECTORS_PERFECT
    integer ( KDI ), dimension ( 3 ) :: &
      SHOCK_I   = 0, &
      CONTACT_I = 0
    real ( KDR ) :: &
      ShockThreshold
  contains
    procedure, public, pass :: &
      InitializeAllocate_P
    generic, public :: &
      Initialize => InitializeAllocate_P
    procedure, public, pass :: &
      Detect
    final :: &
      Finalize
    procedure, public, pass :: &
      SetOutput
  end type FluidFeatures_P_Form

    private :: &
      InitializeBasics, &
      Detect_CSL

contains


  subroutine InitializeAllocate_P &
               ( FF, Fluid, Grid, ShockThreshold, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    class ( VariableGroupForm ), intent ( in ) :: &
      Fluid
    class ( * ), intent ( in ) :: &
      Grid
    real ( KDR ), intent ( in ) :: &
      ShockThreshold
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), &
      optional :: &
        VectorIndicesOption

    type ( Integer_1D_Form ), dimension ( : ), allocatable :: &
      VectorIndices
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      VariableUnit
    character ( LDF ) :: &
      Name 
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      Vector

    call InitializeBasics &
           ( FF, Variable, Vector, Name, VariableUnit, VectorIndices, &
             VariableOption, VectorOption, NameOption, UnitOption, &
             VectorIndicesOption )

    call FF % InitializeTemplate &
           ( Fluid, Grid, nValues, VariableOption = Variable, &
             VectorOption = Vector, NameOption = Name, &
             ClearOption = ClearOption, UnitOption = VariableUnit, &
             VectorIndicesOption = VectorIndices )

    FF % ShockThreshold = ShockThreshold

  end subroutine InitializeAllocate_P


  subroutine Detect ( FF )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF

    integer ( KDI ) :: &
      iD
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      S_I_iD, &
      P, &
      V_iD

    select type ( F => FF % Fluid )
    class is ( Fluid_P_Template )

    select type ( Grid => FF % Grid )
    class is ( Chart_SL_Template )

      call Grid % SetVariablePointer &
             ( F % Value ( :, F % PRESSURE ), P )
      do iD = 1, Grid % nDimensions
        call Grid % SetVariablePointer &
               ( F % Value ( :, F % VELOCITY_U ( iD ) ), V_iD )
        call Grid % SetVariablePointer &
               ( FF % Value ( :, FF % SHOCK_I ( iD ) ), S_I_iD )
        call Detect_CSL &
               ( S_I_iD, Grid, P, V_iD, FF % ShockThreshold, iD, &
                 Grid % nGhostLayers ( iD ) )
      end do

    class default
      call Show ( 'Grid type not recognized', CONSOLE % ERROR )
      call Show ( 'Fluid_P__Template', 'module', CONSOLE % ERROR )
      call Show ( 'DetectFeaturesTemplate_P', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- Grid

    end select !-- F

    nullify ( S_I_iD, P, V_iD )

  end subroutine Detect


  impure elemental subroutine Finalize ( FF )

    type ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF

    call FF % FinalizeTemplate ( )
   
  end subroutine Finalize


  subroutine SetOutput ( FF, Output )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    type ( VariableGroupForm ), intent ( inout ) :: &
      Output

    call Output % Initialize ( FF )

  end subroutine SetOutput


  subroutine InitializeBasics &
               ( FF, Variable, Vector, Name, VariableUnit, VectorIndices, &
                 VariableOption, VectorOption, NameOption, &
                 VariableUnitOption, VectorIndicesOption )

    class ( FluidFeatures_P_Form ), intent ( inout ) :: &
      FF
    character ( LDL ), dimension ( : ), allocatable, intent ( out ) :: &
      Variable, &
      Vector
    character ( LDF ), intent ( out ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), allocatable, &
      intent ( out ) :: &
        VariableUnit
    !-- FIXME: intent(out) here caused ICE with Intel Compiler 15
    !          Temporarily set to intent(inout)
    !type ( Integer_1D_Form ), dimension ( : ), allocatable, &
    !  intent ( out ) :: &
    type ( Integer_1D_Form ), dimension ( : ), allocatable, &
      intent ( inout ) :: &
        VectorIndices
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      VariableUnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional :: &
      VectorIndicesOption

    integer ( KDI ) :: &
      iV  !-- iVector

    if ( FF % Type == '' ) &
      FF % Type = 'a FluidFeatures_P'

    Name = 'FluidFeatures'
    if ( present ( NameOption ) ) &
      Name = NameOption

    FF % IGNORABILITY = CONSOLE % INFO_4 
    call Show ( 'Initializing ' // trim ( FF % Type ), FF % IGNORABILITY )
    call Show ( Name, 'Name', FF % IGNORABILITY )

    !-- variable indices

    if ( FF % N_FIELDS == 0 ) &
      FF % N_FIELDS = FF % N_FIELDS_PERFECT

    FF % SHOCK_I   = [ 1, 2, 3 ]
    FF % CONTACT_I = [ 4, 5, 6 ]

    !-- variable names 

    if ( present ( VariableOption ) ) then
      allocate ( Variable ( size ( VariableOption ) ) )
      Variable = VariableOption
    else
      allocate ( Variable ( FF % N_FIELDS ) )
      Variable = ''
    end if

    Variable ( 1 : FF % N_FIELDS_PERFECT ) &
      = [ 'Shock_I_1  ', &
          'Shock_I_2  ', &
          'Shock_I_3  ', &
          'Contact_I_1', &
          'Contact_I_2', &
          'Contact_I_3' ]
          
    !-- units
    
    if ( present ( VariableUnitOption ) ) then
      allocate ( VariableUnit ( size ( VariableUnitOption ) ) )
      VariableUnit = VariableUnitOption
    else
      allocate ( VariableUnit ( FF % N_FIELDS ) )
      VariableUnit = UNIT % IDENTITY
    end if
    
    !-- vectors

    if ( FF % N_VECTORS == 0 ) &
      FF % N_VECTORS = FF % N_VECTORS_PERFECT

    if ( present ( VectorOption ) ) then
      allocate ( Vector ( size ( VectorOption ) ) )
      Vector = VectorOption
    else
      allocate ( Vector ( FF % N_VECTORS ) )
      Vector = ''
    end if

    Vector ( 1 : FF % N_VECTORS_PERFECT ) &
      = [ 'Shock_I  ', &
          'Contact_I' ]

    !-- vector indices

    if ( present ( VectorIndicesOption ) ) then
      allocate ( VectorIndices ( size ( VectorIndicesOption ) ) )
      do iV = FF % N_VECTORS_PERFECT + 1, size ( VectorIndices )
        call VectorIndices ( iV ) % Initialize ( VectorIndicesOption ( iV ) )
      end do
    else
      allocate ( VectorIndices ( FF % N_VECTORS ) )
    end if

    call VectorIndices ( 1 ) % Initialize ( FF % SHOCK_I )
    call VectorIndices ( 2 ) % Initialize ( FF % CONTACT_I )

  end subroutine InitializeBasics


  subroutine Detect_CSL ( S_I_iD, CSL, P, V_iD, ST, iD, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      S_I_iD
    class ( Chart_SL_Template ), intent ( in ) :: &
      CSL
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      P, &
      V_iD
    real ( KDR ), intent ( in ) :: &
      ST
    integer ( KDI ), intent ( in ) :: &
      iD, &
      oV

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &
      lV, uV
    real ( KDR ) :: &
      dP, &
      P_Min, &
      dLnP, &
      dV_iD, &
      SqrtTiny

    lV = 1
    where ( shape ( P ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( P ) > 1 )
      uV = shape ( P ) - oV
    end where
    uV ( iD ) = size ( P, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = -1
      
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    !$OMP parallel do private ( iV, jV, kV, iaVS )
    do kV = lV ( 3 ), uV ( 3 ) 
      do jV = lV ( 2 ), uV ( 2 )
        do iV = lV ( 1 ), uV ( 1 )

          iaVS = [ iV, jV, kV ] + iaS

          dP  =  abs ( P ( iV, jV, kV )  &
                       -  P ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) )
          P_Min  =  max ( min ( P ( iV, jV, kV ), &
                                P ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ), &
                          SqrtTiny )
          dLnP  =  dP / P_Min

          dV_iD  =  V_iD ( iV, jV, kV ) &
                    -  V_iD ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
 
          if ( dLnP > ST .and. dV_iD >= 0.0_KDR ) then
            S_I_iD ( iV, jV, kV )  =  1.0_KDR
          else
            S_I_iD ( iV, jV, kV )  =  0.0_KDR
          end if

        end do !-- iV
      end do !-- jV
    end do !-- kV
    !$OMP end parallel do
      
  end subroutine Detect_CSL


end module FluidFeatures_P__Form
