      program test_fermi_dirac
      include 'implno.dek'

!..declare
      double precision eta,theta,fd(9),fdeta(9),fdtheta(9), &
                       fdeta2(9),fdtheta2(9),fdetadtheta(9)


!..evaluate some exact fermi-dirac functions and their derivatives
!..in physical terms eta is a chemical potential divided by kerg*temp
!..while theta is a rest mass energy divided by kerg*temp

      eta   = 0.4d0
      theta = 0.1d0

                      


      call dfermi(-0.5d0,eta,theta,fd(1),fdeta(1),fdtheta(1),fdeta2(1),fdtheta2(1),fdetadtheta(1))
      call dfermi(0.5d0,eta,theta,fd(2),fdeta(2),fdtheta(2),fdeta2(2),fdtheta2(2),fdetadtheta(2))
      call dfermi(1.0d0,eta,theta,fd(3),fdeta(3),fdtheta(3),fdeta2(3),fdtheta2(3),fdetadtheta(3))
      call dfermi(1.5d0,eta,theta,fd(4),fdeta(4),fdtheta(4),fdeta2(4),fdtheta2(4),fdetadtheta(4))
      call dfermi(2.0d0,eta,theta,fd(5),fdeta(5),fdtheta(5),fdeta2(5),fdtheta2(5),fdetadtheta(5))
      call dfermi(2.5d0,eta,theta,fd(6),fdeta(6),fdtheta(6),fdeta2(6),fdtheta2(6),fdetadtheta(6))
      call dfermi(3.0d0,eta,theta,fd(7),fdeta(7),fdtheta(7),fdeta2(7),fdtheta2(7),fdetadtheta(7))
      call dfermi(4.0d0,eta,theta,fd(8),fdeta(8),fdtheta(8),fdeta2(8),fdtheta2(8),fdetadtheta(8))
      call dfermi(5.0d0,eta,theta,fd(9),fdeta(9),fdtheta(9),fdeta2(9),fdtheta2(9),fdetadtheta(9))



!..say what we got
      write(6,110) 'fermi-dirac integrals at eta =',eta,' theta =',theta
 110  format(1x,a,1pe12.4,a,1pe12.4)


      write(6,111) 'fermi-dirac', &
                   'eta derivative','theta derivative'
 111  format(1x,t18,a,t32,a,t48,a)


      write(6,112) fd(1),fdeta(1),fdtheta(1),fd(2),fdeta(2),fdtheta(2), &
                   fd(3),fdeta(3),fdtheta(3),fd(4),fdeta(4),fdtheta(4), &
                   fd(5),fdeta(5),fdtheta(5),fd(6),fdeta(6),fdtheta(6), &
                   fd(7),fdeta(7),fdtheta(7),fd(8),fdeta(8),fdtheta(8), &
                   fd(9),fdeta(9),fdtheta(9)

 112  format(1x,'order -1/2 :',1p3e16.8,/, &
             1x,'order  1/2 :',1p3e16.8,/, &
             1x,'order    1 :',1p3e16.8,/, &
             1x,'order  3/2 :',1p3e16.8,/, &
             1x,'order    2 :',1p3e16.8,/, &
             1x,'order  5/2 :',1p3e16.8,/, &
             1x,'order    3 :',1p3e16.8,/, &
             1x,'order    4 :',1p3e16.8,/, &
             1x,'order    5 :',1p3e16.8)


      stop 'normal termination'
      end



!..routines for the brutally efficient quadrature integrations

      include 'fermi_dirac_quadrature.f90'
