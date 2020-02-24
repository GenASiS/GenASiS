module TimeSeriesRadiationFluid_Form

  use Mathematics

  implicit none
  private

  type, public, extends ( TimeSeries_C_1D_C_Form ) :: &
    TimeSeriesRadiationFluidForm
  end type TimeSeriesRadiationFluidForm

end module TimeSeriesRadiationFluid_Form
