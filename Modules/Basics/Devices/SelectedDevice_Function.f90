module SelectedDevice_Function
  
  use OMP_LIB
  use iso_c_binding
  use Specifiers
  use Device_C
  use OffloadEnabled_Function

  implicit none
  private
  
  public :: &
    SelectedDevice
  
contains


  function SelectedDevice ( ErrorOption ) result ( SD )
    
    integer ( KDI ), intent ( out ), optional :: &
      ErrorOption
    integer ( KDI ) :: &
      SD
      
    integer ( KDI ) :: &
      iDevice_OMP, &
      iDevice_GPU, &
      Error
      
    iDevice_OMP = huge ( 1 )
    iDevice_GPU = huge ( 1 )
    SD          = huge ( 1 )
    
    if ( OffloadEnabled ( ) ) then
      iDevice_OMP = OMP_GET_DEFAULT_DEVICE ( )
      Error = GetDevice ( iDevice_GPU )
    
      if ( iDevice_OMP == iDevice_GPU ) then !-- Expect consistency
        SD = iDevice_OMP
      else
        Error = -1
      end if
    else
      Error = -2
    end if
    
    if ( present ( ErrorOption ) ) &
      ErrorOption = Error
      
  end function SelectedDevice


end module SelectedDevice_Function
