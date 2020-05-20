module nuc_eos_find

  use eosmodule, &
    only: nvars, nrho, ntemp, nye, nvars, t_max_hack
  use linterp, &
    only : intp3d_many, intp3d


contains

subroutine findthis(lr,lt,y,value,array,iVar,d1,d2,d3,logrho,logtemp,ye)

#ifdef ENABLE_OMP_OFFLOAD
  !$OMP declare target
#endif

  integer rip,rim
  integer tip,tim
  integer yip,yim

  real*8 lr,lt,y,value,d1,d2,d3
  real*8, dimension ( :, :, :, : ), intent(in) :: array
  integer, intent(in) :: iVar
  real*8, dimension ( : ), intent(in) :: &
    logrho, logtemp, ye
  
  real*8, dimension ( 1 ) :: &
    lr_in, lt_in, y_in, v_out
  
  lr_in ( 1 ) = lr
  lt_in ( 1 ) = lt
  y_in  ( 1 ) = y

! Ewald's interpolator           
!-- FIXME: need to pass logrho, logtemp, ye as args
  call intp3d ( lr_in,lt_in,y_in,v_out,1,array,nrho,ntemp,nye,iVar, &
                logrho,logtemp,ye,d1,d2,d3)
  value = v_out ( 1 )

end subroutine findthis


subroutine findall(lr, lt, y, ff, logrho, logtemp, ye, eos_table)

!  use eosmodule, only: &
!    nvars, nrho, ntemp, nye
!  use linterp, only: intp3d_many
!  implicit none

#ifdef ENABLE_OMP_OFFLOAD  
  !$OMP declare target
#endif
  
  real*8 ff(nvars)
  real*8 ffx(1,nvars)
  real*8 lr,lt,y
  real*8, dimension ( : ), intent(in) :: &
    logrho, logtemp, ye
  real*8, dimension( :, :, :, : ), intent ( in ) :: &
    eos_table
  integer i
  
  real*8, dimension ( 1 ) :: &
    lr_in, lt_in, y_in
  
  lr_in ( 1 ) = lr
  lt_in ( 1 ) = lt
  y_in  ( 1 ) = y
  
! Ewald's interpolator           
  call intp3d_many(lr_in, lt_in, y_in, ffx , 1, eos_table, &
       nrho,ntemp,nye,nvars,logrho,logtemp,ye)
  ff(:) = ffx(1,:)
  
end subroutine findall


subroutine findall_short(lr,lt,y,ff)

!  use eosmodule
!  implicit none

  real*8 ffx(8,1)
  real*8 ff(8)
  real*8 lr,lt,y
  integer i
  integer :: nvarsx = 8


! Ewald's interpolator           
!  call intp3d_many(lr,lt,y,ffx,1,alltables(:,:,:,1:8), &
!       nrho,ntemp,nye,nvarsx,logrho,logtemp,ye)
!  ff(:) = ffx(:,1)

end subroutine findall_short

subroutine findone(lr,lt,y,ff,index)

!  use eosmodule
!  implicit none

  real*8 ffx(1,1)
  real*8 ff(1)
  real*8 lr,lt,y
  integer index
  

! Ewald's interpolator           
!  call intp3d_many(lr,lt,y,ffx,1,alltables(:,:,:,index), &
!       nrho,ntemp,nye,1,logrho,logtemp,ye)
!  ff(:) = ffx(:,1)

end subroutine findone

end module nuc_eos_find
