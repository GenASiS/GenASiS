module InteractionsNeutrinos_BSLL_ASC_CSLD__Form

  !-- InteractionsNeutrinos_BundleSingleLevelLocal_AtlasSingleChart
  !   _ChartSingleLevelDistributed_Form

  use NULIBTABLE_INTERFACE, only: &
    NULIBTABLE_READER
  use Basics
  use Mathematics
  use Fluids
  use RadiationBasics
  use Interactions_OCO__Form
  use InteractionsNeutrinos_ASC__Form

  implicit none
  private

  type, public, extends ( Interactions_BSLL_ASC_CSLD_Template ) :: &
    InteractionsNeutrinos_BSLL_ASC_CSLD_Form
  contains
    procedure, public, pass :: &
      Set
    procedure, private, pass :: &
      AllocateField
  end type InteractionsNeutrinos_BSLL_ASC_CSLD_Form

contains


  subroutine Set ( IB, FA )

    class ( InteractionsNeutrinos_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    class ( Fluid_ASC_Form ), intent ( in ), target :: &
      FA

    integer ( KDI ) :: &
      iF  !-- iFiber
    logical ( KDL ) :: &
      Include_NES, &
      IncludePairs
    character ( LDF ) :: &
      Filename
    class ( GeometryFlatForm ), pointer :: &
      GF
    class ( Fluid_P_HN_Form ), pointer :: &
      Fluid
    class ( InteractionsTemplate ), pointer :: &
      I

    Filename = '../Parameters/NuLib_rho82_temp65_ye51_ng16_ns3_' &
               // 'Itemp65_Ieta61_version1.0_20190829.h5'

    select case ( trim ( IB % InteractionsType ) )
    case ( 'O_CONNOR_OTT' )
      call DelayFileAccess ( PROGRAM_HEADER % Communicator % Rank )
      call NULIBTABLE_READER &
           ( filename                   =  Filename, &
             include_Ielectron          =  Include_NES, &
             include_epannihil_kernels  =  IncludePairs )
    end select

    select type ( B => IB % Bundle )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    Fluid => FA % Fluid_P_HN ( )

    GF => B % GeometryFiber ( )
    associate &
      (    Energy => GF % Value ( :, GF % CENTER_U ( 1 ) ), &
        d3_Energy => GF % Value ( :, GF % VOLUME ) )

    do iF = 1, B % nFibers
      I => IB % Interactions ( iF )
      select type ( I )
      type is ( Interactions_OCO_Form )
        call I % Set &
               ( Fluid = Fluid, &
                 Energy = Energy, &
                 d3_Energy = d3_Energy, &
                 Include_NES = Include_NES, &
                 IncludePairs = IncludePairs, &
                 iBaseCell = B % iaBaseCell ( iF ) )
      end select !-- I
      nullify ( I )
    end do !-- iF

    end associate !-- Energy, etc.
    end select !-- B
    nullify ( GF )

  end subroutine Set


  subroutine AllocateField ( IB, iF )

    class ( InteractionsNeutrinos_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      IB
    integer ( KDI ), intent ( in ) :: &
      iF  !-- iFiber

    allocate ( InteractionsNeutrinos_ASC_Form :: &
                 IB % Fiber % Atlas ( iF ) % Element )

  end subroutine AllocateField


end module InteractionsNeutrinos_BSLL_ASC_CSLD__Form
