module Series_R_CC_C__Form

  !-- Series_Radiation_CentralCore_Collected__Form

  use Basics
  use Mathematics
  use Measures_R_CC_C__Form

  implicit none
  private

  type, public, extends ( Series_CS_1D_C_CS_Form ) :: Series_R_CC_C_Form
    integer ( KDI ) :: &
      iEnergy_F, &
      iEnergy_R, &
      iNumber_F, &
      iNumber_R
    type ( StorageForm ), allocatable :: &
      TotalChange, &
      Measures_CC
    class ( Measures_R_CC_C_Form ), pointer :: &
      Measures => null ( )
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    procedure, public, pass :: &
      Record
!     procedure, public, pass :: &
!       Restore
    final :: &
      Finalize
  end type Series_R_CC_C_Form


contains


  subroutine Initialize_R_CC_C &
               ( S, M, CS_1D, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
                 CommunicatorRank, nWrite, iCycle )

    class ( Series_R_CC_C_Form ), intent ( inout ) :: &
      S
    class ( Measures_R_CC_C_Form ), intent ( in ), target :: &
      M
    class ( CurrentSetForm ), dimension ( : ), intent ( in ) :: &
      CS_1D
    class ( CurrentSetForm ), intent ( in ) :: &
      CS
    type ( GridImageStreamForm ), intent ( in ) :: &
      GIS
    character ( * ), dimension ( : ), intent ( in ) :: &
      dT_Label
    type ( QuantityForm ), intent ( in ) :: &
      Unit_T
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      dT_Candidate
    real ( KDR ), intent ( in ), target :: &
      T
    integer ( KDI ), intent ( in ) :: &
      CommunicatorRank, &
      nWrite
    integer ( KDI ), intent ( in ), target :: &
      iCycle

    integer ( KDI ) :: &
      nVariables, &
      iI, &  !-- iIntegral
      iM  !-- iMeasure
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Unit, &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable, &
      SeriesName

    if ( S % Type == '' ) &
      S % Type = 'a Series_R_CC_C' 

    !-- Parent

    call S % Series_CS_1D_C_CS_Form % Initialize &
           ( CS_1D, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
             CommunicatorRank, nWrite, iCycle  )

    !-- TotalChange

    associate &
      ( TF  =>  S % TallyChange, &                    !-- Fluid
        TR  =>  S % TallyChange_1D ( 1 ) % Pointer )  !-- Radiation

    !-- Fluid indices
    do iI  =  1, TF % nIntegrals
      if ( trim ( TF % Variable ( iI ) )  ==  'TotalEnergy' ) then
        S % iEnergy_F  =  iI
      else if ( trim ( TF % Variable ( iI ) )  ==  'ElectronNumber' ) then
        S % iNumber_F  =  iI
      end if
    end do !-- iI

    !-- Radiation indices
    do iI  =  1, TR % nIntegrals
      if ( trim ( TR % Variable ( iI ) )  ==  'Energy' ) then
        S % iEnergy_R  =  iI
      else if ( trim ( TR % Variable ( iI ) )  ==  'Number' ) then
        S % iNumber_R  =  iI
      end if
    end do !-- iI

    nVariables = 2
    allocate ( Variable ( nVariables ) )
    allocate ( Unit ( nVariables ) )
    Variable ( 1 )  =  'Energy'
    Variable ( 2 )  =  'ElectronNumber'
        Unit ( 1 )  =  TF % Unit ( S % iEnergy_F )
        Unit ( 2 )  =  TF % Unit ( S % iNumber_F )

    allocate ( S % TotalChange )
    associate &
      ( STC  =>  S % TotalChange, &
        SB   =>  S % Basic )
    call STC % Initialize &
           ( [ SB % nValues, nVariables ], &
             VariableOption = Variable, UnitOption = Unit, &
             NameOption =  'TotalChange', &
             ClearOption = .true. )
    if ( allocated ( S % CurveImage ) ) then
      associate ( CI => S % CurveImage )
      call CI % AddStorage ( STC )
      end associate !-- CI
    end if
    end associate !-- STC, etc.

    end associate !-- TF, etc.

    !-- Measures

    S % Measures  =>  M

    allocate ( SeriesName ( M % nMeasures ) )
    allocate ( SeriesUnit ( M % nMeasures ) )
    do iM  =  1,  M % nMeasures
      SeriesName ( iM )  =  M % Name ( iM )
      SeriesUnit ( iM )  =  M % Unit ( iM )
    end do !-- iS

    allocate ( S % Measures_CC )
    associate &
      ( SM  =>  S % Measures_CC, &
        SB  =>  S % Basic )
    call SM % Initialize &
           ( [ SB % nValues, M % nMeasures ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption =  'Measures_CC', &
             ClearOption = .true. )
    if ( allocated ( S % CurveImage ) ) then
      associate ( CI => S % CurveImage )
      call CI % AddStorage ( SM )
      end associate !-- CI
    end if
    end associate !-- SM, etc.

  end subroutine Initialize_R_CC_C


  subroutine Record ( S )

    class ( Series_R_CC_C_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iN, &  !-- iNeutrinos
      iM     !-- iMeasure
    real ( KDR ) :: &
      Sign

    call S % Series_CS_1D_C_CS_Form % Record ( )

    !-- TotalChange

    associate &
      (  STCV  =>  S % TotalChange % Value, &
        iR     =>  S % iRecord )
    
    !-- Fluid
    associate ( TCFV  =>  S % TallyChange % Value )
    STCV ( iR, 1 )  =  TCFV ( S % iEnergy_F )
    STCV ( iR, 2 )  =  TCFV ( S % iNumber_F )
    end associate !-- TCFV

    !-- Radiation (neutrinos)
    do iN  =  1,  size ( S % TallyChange_1D )
      associate ( TCRV  =>  S % TallyChange_1D ( iN ) % Pointer % Value )
      select case ( iN )
        case ( 1 )
          Sign  =   1.0_KDR  !-- NEUTRINOS_E
        case ( 2 )
          Sign  =  -1.0_KDR  !-- NEUTRINOS_EB
        case default
          Sign  =   0.0_KDR  !-- NEUTRINOS_HL, etc.
      end select !-- iN
      STCV ( iR, 1 )  =  STCV ( iR, 1 )  +         TCRV ( S % iEnergy_R )
      STCV ( iR, 2 )  =  STCV ( iR, 2 )  +  Sign * TCRV ( S % iNumber_R )
      end associate !-- TCRV
    end do !-- iN

    end associate !-- STCV, etc.

    !-- Measures

    associate &
      (  SMV  =>  S % Measures_CC % Value, &
        iR    =>  S % iRecord, &
          MV  =>  S % Measures % Value, &
        nM    =>  S % Measures % nMeasures )
    do iM  =  1,  nM
      SMV ( iR, iM )  =  MV ( iM )
    end do !-- iM
    end associate !-- SMV, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( S )

    type ( Series_R_CC_C_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Measures_CC ) ) &
      deallocate ( S % Measures_CC )
    if ( allocated ( S % TotalChange ) ) &
      deallocate ( S % TotalChange )

    nullify ( S % Measures )

  end subroutine Finalize


end module Series_R_CC_C__Form
