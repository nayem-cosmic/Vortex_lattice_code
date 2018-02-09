! 3d-VLM code for simple wing planforms with ground effect.
! This program is modified version of program written by Joe Katz, 1974.

! imax : maximum number of chordwise panels, jmax : maximum number of spanwise panels
! ib : number of chordwise panels, jb : number of spanwise panels, b : wing semi span
! c : chord length, ar : wing aspect ratio, ch : height above ground
! dxw : wake length, aLpha : angle of attack
! phi : dihedral angle, xLambda : sweep angle
! croot : chord length at wing root, ctip : chord length at wing tip
! qf : vortex ring corner points, qc : collocation points, ds : area vector
! qh : center points of panels, it is needed for pressure coefficient graphing
! a1 : auxiliary influence coefficient which is needed for calculation
! gamma : circulation, dL : partial lift force, dd : partial drag force
! dp : pressure difference, a : influence coefficient


module com
    integer, parameter :: imax=20, jmax=50, max=imax*jmax
    real, parameter :: pi=4.*atan(1.)
    integer :: ib,jb,ib1,ib2,jb1,isign
    real :: b,c,s,ar,ch,dxw,aLpha,phi,xLambda,croot,ctip,zm,p

    real :: qf(imax+1,jmax+1,3),qc(imax,jmax,3),qh(imax,jmax,2),ds(imax,jmax,4),a1(imax+1,jmax)
end module com

program main
    use com

    real :: gamma(imax,jmax),dL(imax,jmax),dd(imax,jmax),dp(imax,jmax),a(max,max), &
gamma1(max),dw(max),dLy(jmax),ddy(jmax),dLx(imax),ddx(imax)
    integer :: ip(max)

    open(10, file='outputdata/summary.txt')
    open(11, file='outputdata/mesh.txt')
    open(12, file='outputdata/gamma.txt')
    open(13, file='outputdata/dp_coeff.txt')
    open(14, file='outputdata/d_lift_span.txt')
    open(15, file='outputdata/d_drag_span.txt')
    open(16, file='outputdata/d_lift_chord.txt')
    open(17, file='outputdata/d_drag_chord.txt')
    open(18, file='outputdata/log.txt')

    ib=5
    jb=15
    ib1=ib+1
    ib2=ib+2
    jb1=jb+1
    croot=1.5
    ctip=.5
    b=1.415
    vt=5.0
    aLpha1=5.354
    aLpha=aLpha1*pi/180.0
    phi1=0.
    phi=phi1*pi/180.0
    xLambda1=53.54
    xLambda=xLambda1*pi/180
! If ch<=100, ground effect is counted
    ch=1000.
! zm : maximum camber for NACA profile, p : location of maximum camber 
    zm=0.0
    p=0.4
! constants
    dxw=100.0*b
    ro=1.0

! Write heading of files
    write(12,106)
    write(13,108) alpha1,xlambda1,phi1
    write(13,109)
    write(14,111) alpha1,xlambda1,phi1,ro,vt,ch,b/jb
    write(14,112)
    write(15,114) alpha1,xlambda1,phi1,ro,vt,ch,b/jb
    write(15,115)
    write(16,117) alpha1,xlambda1,phi1,ro,vt,ch,(croot+ctip)/2/ib
    write(16,118)
    write(17,120) alpha1,xlambda1,phi1,ro,vt,ch,(croot+ctip)/2/ib
    write(17,121)

! wing geometry
    call grid

! aerodynamic calculations
    do i=1,ib
        do j=1,jb
! gamma(i,j)=1.0 is required for influence matrix calculations
            gamma(i,j)=1.0
        end do
    end do

    write(*,*) "Induced coefficients calculation starts..."
! influence coefficients calculation
    k=0
    do i=1,ib
        do j=1,jb
            k=k+1
            isign=0
            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,ib
                do j1=1,jb
                    L=L+1
! a(k,L)-is the normal velocity component due to a unit vortex lattice
                    a(k,L)=a1(i1,j1)
                end do
            end do

! add influence of wing's other half
            isign=1
            call wing(qc(i,j,1),-qc(i,j,2),qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,ib
                do j1=1,jb
                    L=L+1
                    a(k,L)=a(k,L)+a1(i1,j1)
                end do
            end do

            if (ch.gt.100.0) goto 12
! add influence of mirror image (due to ground)
            isign=2
            call wing(qc(i,j,1),qc(i,j,2),-qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,ib
                do j1=1,jb
                    L=L+1
                    a(k,L)=a(k,L)+a1(i1,j1)
                end do
            end do

! add mirror image influence of wing's other half
            isign=3
            call wing(qc(i,j,1),-qc(i,j,2),-qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,ib
                do j1=1,jb
                    L=L+1
                    a(k,L)=a(k,L)+a1(i1,j1)
                end do
            end do

12          continue
            isign=0

! calculate wing geometrical downwash
            uinf=vt
            vinf=0.0
            winf=0.0
            dw(k)=-(uinf*ds(i,j,1)+vinf*ds(i,j,2)+winf*ds(i,j,3))
        end do
    end do
    write(*,*) "Induced coefficients calculated."
    write(18,*) "Induced coefficients calculated."

! Solution of the problem: dw(i) = a(i,j)*gamma(i)
    k1=ib*jb
    do k=1,k1
        gamma1(k)=dw(k)
    end do

    call decomp(k1,max,a,ip)
    call solver (k1,max,a,gamma1,ip)

! wing vortex lattice listing
    k=0
    do i=1,ib
        do j=1,jb
            k=k+1
            gamma(i,j)=gamma1(k)
        end do
    end do

! forces calculation
    fL=0.0
    fd=0.0
    fm=0.0
    que=0.5*ro*vt*vt

    do j=1,jb
        dLy(j)=0.
        ddy(j)=0.
        do i=1,ib
            if(i.eq.1) gammaij=gamma(i,j)
            if(i.gt.1) gammaij=gamma(i,j)-gamma(i-1,j)
            dym=qf(i,j+1,2)-qf(i,j,2)
            dL(i,j)=ro*vt*gammaij*dym

            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u1,v1,w1,0.0,i,j)
            call wing(qc(i,j,1),-qc(i,j,2),qc(i,j,3),gamma,u2,v2,w2,0.0,i,j)

            if(ch.gt.100.0) goto 194
            call wing(qc(i,j,1),qc(i,j,2),-qc(i,j,3),gamma,u3,v3,w3,0.0,i,j)
            call wing(qc(i,j,1),-qc(i,j,2),-qc(i,j,3),gamma,u4,v4,w4,0.0,i,j)

            goto 195

194         w3=0.0
            w4=0.0

195         wind=w1+w2-w3-w4

            dd(i,j)=-ro*gammaij*wind*dym
            dp(i,j)=dL(i,j)/ds(i,j,4)/que
            dLy(j)=dLy(j)+dL(i,j)
            ddy(j)=ddy(j)+dd(i,j)
            fL=fL+dL(i,j)
            fd=fd+dd(i,j)
            fm=fm+dL(i,j)*(qf(i,j,1)-0)

            write(12,107) i,j,gammaij
            write(13,110) qh(i,j,1),qh(i,j,2),dp(i,j)
        end do
        write(14,113) j,dLy(j)*jb/b
        write(15,116) j,ddy(j)*jb/b
    end do
    write(*,*) "Forces calculated."
    write(18,*) "Forces calculated."
    do i=1,ib
        dLx(i)=0.
        ddx(i)=0.
        do j=1,jb
            dLx(i)=dLx(i)+dL(i,j)
            ddx(i)=ddx(i)+dd(i,j)
        end do
        write(16,119) i,dLx(i)*ib/c
        write(17,122) i,ddx(i)*ib/c
    end do
    cL=fL/(que*s)
    cd=fd/(que*s)
    cm=fm/(que*s*c)
    write(*,*) "Coefficients calculated."
    write(18,*) "Coefficients calculated."

! Write statements
    write(10,101)
    write(10,102) 2.*b,croot,ctip,ib*jb,2*s,ar,aLpha1,xLambda1,phi1,vt,ro,ch
    write(10,103) cL, 2*fL, cm, cd
    write(*,103) cL, 2*fL, cm, cd
 
! Formats
101 format('Summary',/,7('-'),/)
102 format('Wing Span =',11x,f8.2,/,'Root Chord Length =',3x,f8.2,/,'Tip Chord Length =',4x,f8.2,&
/,'Total Number of',/, 'Panels (On Semi Span) =',i7,/,'Wing Area = ',10x,f8.2,/, &
'Aspect Ratio =',8x,f8.2,/,'Angle of Attack =',5x,f8.2,&
/,'Sweep Angle =',9x,f8.2,/,'Dihedral Angle =',6x,f8.2,/,'Free Stream Velocity =',f8.2, &
/,'Density of Medium =',3x,f8.2,/,'Height Above ground =',f9.2,/)
103 format('CL = ',f10.5,/,'FL = ',f10.2,/,'CM = ',f10.5,/,'CD = ',f10.5)
! 104,105 in grid
106 format(4x,'i',4x,'j',4x,'Gamma(i,j)')
107 format(2(i5),2x,f10.5)
108 format('Alpha:'f5.1,1x,'Lambda:',f5.1,1x,'Phi:',f5.1)
109 format(4x,'qh(i,j,1)',1x,'qh(i,j,2)',1x,'CdP(i,j) : Pressure difference coefficient of panel')
110 format(2(f10.3),2x,f10.5)
111 format('Alpha:'f5.1,1x,'Lambda:',f5.1,1x,'Phi:',f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'H.G.:',f6.1,/,'dy:',1x,f6.4)
112 format(4x,'j',4x,'dLy(j)*jb/b : Lift per span length')
113 format(i5,2x,f10.5)
114 format('Alpha:'f5.1,1x,'Lambda:',f5.1,1x,'Phi:',f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'H.G.:',f6.1,/,'dy:',1x,f6.4)
115 format(4x,'j',4x,'ddy(j)*jb/b : Drag per span length')
116 format(i5,2x,f10.5)
117 format('Alpha:'f5.1,1x,'Lambda:',f5.1,1x,'Phi:',f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'H.G.:',f6.1,/,'dx:',1x,f6.4)
118 format(4x,'i',4x,'dLx(i)*ib/c : Lift per chord length')
119 format(i5,2x,f10.5)
120 format('Alpha:'f5.1,1x,'Lambda:',f5.1,1x,'Phi:',f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'H.G.:',f6.1,/,'dx:',1x,f6.4)
121 format(4x,'i',4x,'ddx(i)*ib/c : Drag per chord length')
122 format(i5,2x,f10.5)

! Close opened files
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    write(*,*) "All opend files closed."
    write(18,*) "All opened files closed."
    write(*,*) "End of program."
    write(18,*) "End of program."

end program main

subroutine grid
    use com
    write(*,*) "Subroutine grid starts..."
    write(18,*) "Subroutine grid starts..."
    write(11,104)
104 format(5x,'xMesh',5x,'yMesh',5x,'zMesh')
105 format(3(f10.3))

! qf(i,j,1-3): wing fixed vortices location
    dy=b/jb
    do j=1,jb1
        yLe=dy*(j-1)
        xLe=0+yLe*tan(xLambda)
        xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
! ci : instantenous chord length
        ci=xte-xLe
        dx=ci/ib
        do i=1,ib1
            qf(i,j,1)=(xLe+dx*(i-0.75))*cos(aLpha)
            xc=dx*(i-0.75)
            if(xc.le.p*ci) then
                yc=zm*xc*(2*p-xc/ci)/p**2
            else
                yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
            end if
            qf(i,j,2)=yLe*cos(phi)-yc*tan(phi)
            qf(i,j,3)=yc*cos(alpha)*cos(phi)-qf(i,j,1)*tan(aLpha)+qf(i,j,2)*tan(phi)+ch
        end do
        qf(ib2,j,1)=xte+dxw
        qf(ib2,j,2)=qf(ib1,j,2)
        qf(ib2,j,3)=qf(ib1,j,3)
    end do

! Generating coordinate points for mesh plotting
    do j=1,jb1
        yLe=dy*(j-1)
        xLe=0+yLe*tan(xLambda)
        xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
        ci=xte-xLe
        dx=ci/ib
        if(mod(j,2).ne.0)then
            do i=1,ib1
                xme=(xLe+dx*(i-1))*cos(aLpha)
                xc=dx*(i-1)
                if(xc.le.p*ci) then
                    yc=zm*xc*(2*p-xc/ci)/p**2
                else
                    yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
                end if
                yme=yLe*cos(phi)-yc*tan(phi)
                zme=yc*cos(alpha)*cos(phi)-xme*tan(aLpha)+yme*tan(phi)
                write(11,105) xme, yme, zme
            end do
        else
            do i=ib1,1,-1
                xme=(xLe+dx*(i-1))*cos(aLpha)
                xc=dx*(i-1)
                if(xc.le.p*ci) then
                    yc=zm*xc*(2*p-xc/ci)/p**2
                else
                    yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
                end if
                yme=yLe*cos(phi)-yc*tan(phi)
                zme=yc*cos(alpha)*cos(phi)-xme*tan(aLpha)+yme*tan(phi)
                write(11,105) xme, yme, zme
            end do
        end if
    end do
    do i=1,ib1
        if(mod(jb,2).ne.0)then
            if(mod(i,2).ne.0)then
                do j=1,jb1
                    yLe=dy*(j-1)
                    xLe=0+yLe*tan(xLambda)
                    xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
                    ci=xte-xLe
                    dx=ci/ib
                    xme=(xLe+dx*(i-1))*cos(aLpha)
                    xc=dx*(i-1)
                    if(xc.le.p*ci) then
                        yc=zm*xc*(2*p-xc/ci)/p**2
                    else
                        yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
                    end if
                    yme=yLe*cos(phi)-yc*tan(phi)
                    zme=yc*cos(alpha)*cos(phi)-xme*tan(aLpha)+yme*tan(phi)
                    write(11,105) xme, yme, zme
                end do
            else
                do j=jb1,1,-1
                    yLe=dy*(j-1)
                    xLe=0+yLe*tan(xLambda)
                    xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
                    ci=xte-xLe
                    dx=ci/ib
                    xme=(xLe+dx*(i-1))*cos(aLpha)
                    xc=dx*(i-1)
                    if(xc.le.p*ci) then
                        yc=zm*xc*(2*p-xc/ci)/p**2
                    else
                        yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
                    end if
                    yme=yLe*cos(phi)-yc*tan(phi)
                    zme=yc*cos(alpha)*cos(phi)-xme*tan(aLpha)+yme*tan(phi)
                    write(11,105) xme, yme, zme
                end do
            end if
        else
            if(mod(i,2).eq.0)then
                do j=1,jb1
                    yLe=dy*(j-1)
                    xLe=0+yLe*tan(xLambda)
                    xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
                    ci=xte-xLe
                    dx=ci/ib
                    xme=(xLe+dx*(i-1))*cos(aLpha)
                    xc=dx*(i-1)
                    if(xc.le.p*ci) then
                        yc=zm*xc*(2*p-xc/ci)/p**2
                    else
                        yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
                    end if
                    yme=yLe*cos(phi)-yc*tan(phi)
                    zme=yc*cos(alpha)*cos(phi)-xme*tan(aLpha)+yme*tan(phi)
                    write(11,105) xme, yme, zme
                end do
            else
                do j=jb1,1,-1
                    yLe=dy*(j-1)
                    xLe=0+yLe*tan(xLambda)
                    xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
                    ci=xte-xLe
                    dx=ci/ib
                    xme=(xLe+dx*(i-1))*cos(aLpha)
                    xc=dx*(i-1)
                    if(xc.le.p*ci) then
                        yc=zm*xc*(2*p-xc/ci)/p**2
                    else
                        yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
                    end if
                    yme=yLe*cos(phi)-yc*tan(phi)
                    zme=yc*cos(alpha)*cos(phi)-xme*tan(aLpha)+yme*tan(phi)
                    write(11,105) xme, yme, zme
                end do
            end if
        end if 
    end do
    
    write(*,*) "Mesh generated."
    write(18,*) "Mesh generated."
! Generating center points of panels for contour graphing of pressure coefficients
    do j=1,jb
        yLe=dy*(j-0.5)
        xLe=0+yLe*tan(xLambda)
        xte=croot+(b*tan(xLambda)+ctip-croot)*yLe/b
        dx=(xte-xLe)/ib

        do i=1,ib1
            qh(i,j,1)=xLe+dx*(i-0.5)
            qh(i,j,2)=yLe
        end do
    end do

! Generating wing collocation points and panel area vectors
    do j=1,jb
        do i=1,ib
            qc(i,j,1)=(qf(i,j,1)+qf(i,j+1,1)+qf(i+1,j+1,1)+qf(i+1,j,1))/4
            qc(i,j,2)=(qf(i,j,2)+qf(i,j+1,2)+qf(i+1,j+1,2)+qf(i+1,j,2))/4
            qc(i,j,3)=(qf(i,j,3)+qf(i,j+1,3)+qf(i+1,j+1,3)+qf(i+1,j,3))/4

            call panel(qf(i,j,1),qf(i,j,2),qf(i,j,3),qf(i+1,j,1),qf(i+1,j,2), &
qf(i+1,j,3),qf(i,j+1,1),qf(i,j+1,2),qf(i,j+1,3),qf(i+1,j+1,1), &
qf(i+1,j+1,2),qf(i+1,j+1,3),ds(i,j,1),ds(i,j,2),ds(i,j,3),ds(i,j,4))
        end do
    end do
    write(*,*) "Collocation points generated."
    write(18,*) "Collocation points generated."

    s=0.5*(croot+ctip)*b
    c=s/b
    ar=2.0*b*b/s
    write(*,*) "End of subroutine grid."
    write(18,*) "End of subroutine grid."

    return
end subroutine grid

subroutine panel(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,v1,v2,v3,sp)

! calculation of panel area and vector
    a1=x2-x3
    a2=y2-y3
    a3=z2-z3

    b1=x4-x1
    b2=y4-y1
    b3=z4-z1

! normal vector
    x=a2*b3-a3*b2
    y=b1*a3-a1*b3
    z=a1*b2-a2*b1
    a=sqrt(x**2+y**2+z**2)

    v1=x/a
    v2=y/a
    v3=z/a

! Calculation of panel area. sp=0.5*|AxB|
    sp=a/2

! calculation of panel area
!    e1=x3-x1
!    e2=y3-y1
!    e3=z3-z1

!    f1=x2-x1
!    f2=y2-y1
!    f3=z2-z1

! normal areas (f*b+b*e)
!    s11=f2*b3-f3*b2
!    s12=b1*f3-f1*b3
!    s13=f1*b2-f2*b1
!    s21=b2*e3-b3*e2
!    s22=e1*b3-b1*e3
!    s23=b1*e2-b2*e1
!    sp=0.5*(sqrt(s11**2+s12**2+s13**2)+sqrt(s21**2+s22**2+s23**2))

    return
end subroutine panel

subroutine vortex(x,y,z,x1,y1,z1,x2,y2,z2,gamma,u,v,w)
    use com, only : pi
    rcut=1.0e-10

    r1r2x=(y-y1)*(z-z2)-(z-z1)*(y-y2)
    r1r2y=-((x-x1)*(z-z2)-(z-z1)*(x-x2))
    r1r2z=(x-x1)*(y-y2)-(y-y1)*(x-x2)

    square=r1r2x*r1r2x+r1r2y*r1r2y+r1r2z*r1r2z

    r1=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))
    r2=sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2))
    if(r1.Lt.rcut .or. r2.Lt.rcut .or. square.Lt.rcut) goto 1
    r0r1=(x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1)
    r0r2=(x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2)
    coef=gamma/(4.0*pi*square)*(r0r1/r1-r0r2/r2)
    u=r1r2x*coef
    v=r1r2y*coef
    w=r1r2z*coef
    goto 2
1   u=0.0
    v=0.0
    w=0.0
2   continue
    return
end subroutine vortex

subroutine wing(x,y,z,gamma,u,v,w,onoff,i1,j1)
    use com
    real :: gamma(imax,jmax)

    u=0.0
    v=0.0
    w=0.0
    do i=1,ib1
        do j=1,jb
            i3=i
            if(i.eq.ib1) i3=ib
            vortic=gamma(i3,j)
            if(onoff.Lt.0.1) goto 2
            call vortex(x,y,z,qf(i,j,1),qf(i,j,2),qf(i,j,3),qf(i,j+1,1),    &
            qf(i,j+1,2),qf(i,j+1,3),vortic,u1,v1,w1)
            call vortex(x,y,z,qf(i+1,j+1,1),qf(i+1,j+1,2),qf(i+1,j+1,3),qf(i+1,j,1),&
            qf(i+1,j,2),qf(i+1,j,3),vortic,u3,v3,w3)
2           call vortex(x,y,z,qf(i,j+1,1),qf(i,j+1,2),qf(i,j+1,3),qf(i+1,j+1,1), &
            qf(i+1,j+1,2),qf(i+1,j+1,3),vortic,u2,v2,w2)
            call vortex(x,y,z,qf(i+1,j,1),qf(i+1,j,2),qf(i+1,j,3),qf(i,j,1),&
            qf(i,j,2),qf(i,j,3),vortic,u4,v4,w4)

            u0=u2+u4+(u1+u3)*onoff
            v0=v2+v4+(v1+v3)*onoff
            w0=w2+w4+(w1+w3)*onoff

! magnitude of the influence co-efficient
            if(isign.eq.0)then
                a1(i,j)=u0*ds(i1,j1,1)+v0*ds(i1,j1,2)+w0*ds(i1,j1,3)
            else if(isign.eq.1)then
                a1(i,j)=u0*ds(i1,j1,1)-v0*ds(i1,j1,2)+w0*ds(i1,j1,3)
            else if(isign.eq.2)then
                a1(i,j)=u0*ds(i1,j1,1)+v0*ds(i1,j1,2)-w0*ds(i1,j1,3)
            else if(isign.eq.3)then
                a1(i,j)=u0*ds(i1,j1,1)-v0*ds(i1,j1,2)-w0*ds(i1,j1,3)
            end if

            if(i.eq.ib1) a1(ib,j)=a1(ib,j)+a1(ib1,j)

            u=u+u0
            v=v+v0
            w=w+w0
        end do
    end do

    return
end subroutine wing

subroutine decomp(n,ndim,a,ip)
    real :: a(ndim,ndim),t
    integer :: ip(ndim)
    write(*,*) "Subroutine decomp starts..."
    write(18,*) "Subroutine decomp starts..."
 
    ip(n)=1
    do k=1,n
        if(k.eq.n) goto 5
        kp1=k+1
        m=k
        do i=kp1,n
            if(abs(a(i,k)).gt.abs(a(m,k))) m=i
        end do
        ip(k)=m
        if(m.ne.k) ip(n)=-ip(n)
        t=a(m,k)
        a(m,k)=a(k,k)
        a(k,k)=t
        if(t.eq.0.e0) goto 5
        do i=kp1,n
            a(i,k)=-a(i,k)/t
        end do
        do j=kp1,n
            t=a(m,j)
            a(m,j)=a(k,j)
            a(k,j)=t
            if(t.eq.0.e0) goto 4
            do i=kp1,n
                a(i,j)=a(i,j)+a(i,k)*t
            end do
4       end do
5       if(a(k,k).eq.0.e0) ip(n)=0
    end do
    write(*,*) "Subroutine decomp ended."
    write(18,*) "Subroutine decomp ended."

    return
end subroutine  decomp

subroutine solver(n,ndim,a,b,ip)
    real :: a(ndim,ndim),b(ndim),t
    integer :: ip(ndim)
    write(*,*) "Subroutine solver starts..."
    write(18,*) "Subroutine solver starts..."
 
    if(n.eq.1) goto 9
    nm1=n-1
    do k=1,nm1
        kp1=k+1
        m=ip(k)
        t=b(m)
        b(m)=b(k)
        b(k)=t
        do i=kp1,n
            b(i)=b(i)+a(i,k)*t
        end do
    end do
    do kb=1,nm1
        km1=n-kb
        k=km1+1
        b(k)=b(k)/a(k,k)
        t=-b(k)
        do i=1,km1
            b(i)=b(i)+a(i,k)*t
        end do
    end do
9   b(1)=b(1)/a(1,1)
    write(*,*) "Subroutine solver ended."
    write(18,*) "Subroutine solver ended."

    return
end subroutine solver



!    subroutine cad
!    use para
!    integer:: i,j, k,n
!    integer:: i1,i2,i3,i4
!       integer :: k1, k2, k3, k4
!    integer:: n1(max),n2(max),n3(max),n4(max)
!    real:: x(max),y(max),z(max)
!
!    write(3,22)
! 22    format(3x,'    NODAL POINTS FROM THE CAD')
!
!    n=0
!    do k=1,JB
!    n=n+1
!      do i=1,IB
!        k1=i     + (n-1)*IB1
!        k2=i + 1 + (n-1)*IB1
!        k3=i + 1 + n* IB1
!        k4=i     + n*IB1
!
!        j = (n-1)*IB + i
!
!        n1(j)= k1
!        n2(j)= k2
!        n3(j)= k3
!        n4(j)= k4
!        write(3,32) j, k1,k2, k3, k4
! 32    format(i5,4i12)
!     end do
!    end do
!
!    write(3,33)
!33    format(/,3x,' No.'8x,'x'15x,'y'12x,'z')
!
!    n=0
!    do j=1,JB1
!    n=n+1
!      do i=1,IB1
!       k =(n-1)*IB1+i
!       x(k)=QF(I,J,1)
!       y(k)=QF(I,J,2)
!       z(k)=QF(I,J,3)
!       write(3,27) k, x(k), y(k), z(k)
! 27      format(i5,3f15.6)
!      end do
!    end do
!
!!    For body surface
!    open(2,file='body.scr',status='unknown')
!    write (2, '(A)') 'PFACE'
!    do i=1,max
!       call scrpoint(2, x(i), y(i), z(i))
!    end do
!    do j=1,IB*JB
!       i1=n1(j)
!       i2=n2(j)
!       i3=n3(j)
!       i4=n4(j)
!       call scrquadr(2, i1, i2, i3, i4)
!    end do
!    write(2,'(/)')
!    close(2)
!    return
!    end
!
!
!
!
!    subroutine scrpoint (un, x1, y1, z1)
!    implicit none
!    integer:: un
!    real:: x1, y1, z1
!    integer:: m, l
!    character line*80, lin*40
!
!6    format(2(f12.5, ','),f12.5)
!    write (lin, 6) x1, y1, z1
!
!    m = 0
!    do l = 1, 40
!     if (lin(l:l).ne.' ') then
!       m = m + 1
!       line(m:m) = lin(l:l)
!     end if
!    end do
!    write (un,'(a)') line(1:m)
!    return
!    end
!
!
!
!    subroutine scrquadr (un, p1, p2, p3, p4)
!    implicit none
!    integer:: un, p1, p2, p3, p4
!    integer:: m, l, pp(4), i
!    character lin*10, line*10
!
!6    format(i6)
!
!    pp(1) = p1
!    pp(2) = p2
!    pp(3) = p3
!    pp(4) = p4
!    write (un, '(1x)')
!    do i = 1, 4
!      write (lin, 6) pp(i)
!      m = 0
!      do l = 1, 6
!        if (lin(l:l).ne.' ') then
!        m = m + 1
!        line(m:m) = lin(l:l)
!        end if
!      end do
!      write (un,'(a)') line(1:m)
!    end do
!    return
!    end
!
!
!
!    SUBROUTINE HGRID
!    use para
!!    DATA GENERATION FOR NACA-0012 HYDROFOIL
!!    IB: NO. OF CHORDWISE PANELS
!!    JB: NO.OF SPANWISE PANELS
!!    QF(I,J,1 TO 3): WING FIXED VORTICES LOCATION
!
!    CROOT=1.0
!    CTIP=1.0
!    XTIP=0.0
!    ZTIP=0.0
!
!    WRITE(6,9)
!9    FORMAT(/,3X,'AIRFOIL COORDINATES',/,3X,19('='))
!
!    OPEN(5,FILE='NACA-0012.DAT')
!!    OPEN(5,FILE='REVERSE.DAT')
!
!    READ(5,*) IB1
!    WRITE(6,10) IB1
!
!    WRITE(6,11)
!    DO I=1,IB1
!      READ(5,*)   QF(I,1,1),QF(I,1,3)
!      WRITE(6,12) QF(I,1,1),QF(I,1,3)
!    END DO
!
!10    FORMAT(/,6X,'IB1=',I7)
!11    FORMAT(/,6X,'X',8X,'Z')
!
!12    FORMAT(3F10.4)
!100   FORMAT(/,3X,'I',4X,'J',6X,'QF(I,J,1)',5X,'QF(I,J,2)',5X,'QF(I,J,3)')
!110   FORMAT(/,3X,'I',4X,'J',7X,'QC(I,J,1)',5X,'QC(I,J,2)',5X,'QC(I,J,3)')
!120    FORMAT(2X,'J',12X,'Y',7X,'DXLE',6X,'DZLE',6X,'CHORD')
!
!200   FORMAT(I3,2X,I3,4F15.4)
!500   FORMAT(I3,5X,4F10.3)
!
!
!!    CALCULATE PANEL CORNER POINTS; QF(I,J,(X,Y,Z))
!    WRITE(7,120)
!
!    DO 3 J=1,JB1
!      Y=B/2.0/JB* FLOAT(J-1)
!      DXLE=XTIP*Y/B
!      DZLE=ZTIP*Y/B
!      CHORD=CROOT-(CROOT-CTIP)*Y/B
!      WRITE(7,500)J,Y,DXLE,DZLE,CHORD
!
!
!!     B- SEMI WING SPAN, DXLE-LOCAL SWEEP
!
!      WRITE(6,100)
!
!      DO I=1,IB1
!       QF(I,J,1)=QF(I,1,1)*CHORD+DXLE
!       QF(I,J,2)=Y
!       QF(I,J,3)=QF(I,1,3)*CHORD+DZLE
!
!       WRITE(6,200) I, J, QF(I,J,1), QF(I,J,2), QF(I,J,3)
!
!      END DO
!
!!     WAKE FAR FIELD POINTS (QF-IS IN BODY FRAME OF REFERENCE)
!
!      QF(IB2,J,1)=QF(IB1,J,1)+DXW
!      QF(IB2,J,2)=QF(IB1,J,2)
!      QF(IB2,J,3)=QF(IB1,J,3)
!
!      WRITE(6,200) I, J, QF(IB2,J,1), QF(IB2,J,2), QF(IB2,J,3)
!
!3    CONTINUE
!
!
!300    FORMAT(/,3X,'I',4X,'J',6X,'QC(I,J,1)',5X,'QC(I,J,2)',5X,'QC(I,J,3)')
!
!    WRITE(7,400)
!400    FORMAT(/2X,'I',4X,'J',8X,'DS(I,J,1)',6X,'DS(I,J,2)',6X,'DS(I,J,3)',5X,'DS(I,J,4)')
!
!!------> WING COLLOCATION POINTS
!
!      DO J=1,JB
!      WRITE(6,300)
!        DO I=1,IB1
!         QC(I,J,1)=(QF(I,J,1)+QF(I,J+1,1)+QF(I+1,J+1,1)+QF(I+1,J,1))/4.0
!         QC(I,J,2)=(QF(I,J,2)+QF(I,J+1,2)+QF(I+1,J+1,2)+QF(I+1,J,2))/4.0
!         QC(I,J,3)=(QF(I,J,3)+QF(I,J+1,3)+QF(I+1,J+1,3)+QF(I+1,J,3))/4.0
!
!         WRITE(6,200)I,J,QC(I,J,1),QC(I,J,2),QC(I,J,3)
!
!!------>    COMPUTATION OF CHORDWISE VECTORS DS(I,J,1-3), TANGENTIAL AND
!!------>    NORMAL VECTORS DS(I,J,4 TO 9),PANEL AREA DS(I,J,1-10)
!!------>    AND SOURCE STRENGTH (SIGMA)
!
!       CALL PANEL(QF(I,J,1),QF(I,J,2),QF(I,J,3),QF(I+1,J,1),QF(I+1,J,2),  &
!            QF(I+1,J,3),QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),QF(I+1,J+1,1),&
!            QF(I+1,J+1,2),QF(I+1,J+1,3),DS(I,J,1),DS(I,J,2),DS(I,J,3),    &
!            DS(I,J,4))
!          WRITE(7,200)I,J,DS(I,J,1),DS(I,J,2),DS(I,J,3),DS(I,J,4)
!
!        END DO
!        WRITE(6,*)
!        WRITE(7,*)
!      END DO
!
!!------> B-IS SEMI WING SPAN,C-ROOT CHORD,S-AREA
!
!      S=0.5*B*(CROOT+CTIP)
!      C=S/B
!      AR=B*B/S
!
!
!    WRITE(6,*) 'END OF DATA FROM SUBROUTINE HGRID'
!
!    RETURN
!    END
!
!
!
