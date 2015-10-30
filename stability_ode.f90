!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Stability diagrams for ODE solver
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 15, 2015
!-----------------------------------------------------------------------------!


program ode
implicit none
real*8 ::pi,theta,dtheta
integer::j,k,np,kmax
complex*16::z1,z2,ii

pi = 4.0d0*datan(1.0d0)
np = 360*10
dtheta = 2.0d0*pi/dfloat(np)
ii = dcmplx(0.0d0, 1.0d0)

kmax = 30 !for Newton Raphson Iteration for root finding

open(12, file="stability_curve_rk2.plt")
write(12,*)'variables ="Re","Im"'

open(13, file="stability_curve_rk4.plt")
write(13,*)'variables ="Re","Im"'

z1 = dcmplx(0.0d0, 0.0d0)
z2 = dcmplx(0.0d0, 0.0d0)
   
!for all angles
do j=0,np
theta = dfloat(j)*dtheta
	!compute roots of z
    do k = 1,kmax
    z1 = z1 - (1.0d0 + z1 + 0.5d0*z1*z1 - cdexp(ii*theta*2.0d0) )/(1.0d0 + z1)
    z2 = z2 - (1.0d0 + z2 + 0.5d0*z2*z2 + z2*z2*z2/6.0d0 + z2*z2*z2*z2/24.0d0 &
       - cdexp(ii*theta*4.0d0) )/(1.0d0 + z2 + 0.5d0*z2*z2 + z2*z2*z2/6.0d0)
    end do
    
    write(12,*) dreal(z1),dimag(z1)
    write(13,*) dreal(z2),dimag(z2)  
    
end do

close(12)
close(13)

end






