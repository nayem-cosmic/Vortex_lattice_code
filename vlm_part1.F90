! 3d-VPM code for simple wing planforms with ground effect.
! This program is modified version of program written by Joe Katz, 1974.

! imax : maximum number of chordwise panels, jmax : maximum number of spanwise panels
! ib : number of chordwise panels, jb : number of spanwise panels, b : wing span
! c : chord length, ar : wing aspect ratio
! dxw : wake length, aLpha : angle of attack
! qf : vortex ring corner points, qc : collocation points, ds : area vector
! qh : center points of panels, it is needed for pressure coefficient graphing
! a1 : auxiliary influence coefficient which is needed for calculation
! gamma : circulation, dL : partial lift force, dd : partial drag force
! dp : pressure difference, a : influence coefficient
! zm : maximum camber, p : max camber position, t : thicknes of foil

module com
    integer, parameter :: imax=20, jmax=50, max=imax*jmax
    real, parameter :: pi=4.*atan(1.)
    integer :: ib,jb,d_ib,d_ib1,d_ib2,jb1  ! d_ib1 = ib*2+1 and similiar. This is for upper and lower surfaces panel number.
    real :: b,c,s,ar,vt,rho,dxw,aLpha,aLpha1,zm,p,t

    real :: qf(imax+1,jmax+1,3),qc(imax,jmax,3),qh(imax,jmax,2),ds(imax,jmax,4),a1(imax+1,jmax)
end module com

program main
    use com
    character :: ans
    integer :: inaca
    integer :: idig(4)

    write(*,*) "Program starts..."
    write(18,*) "Program starts..."
    open(10, file='outputdata/summary.txt')
    open(11, file='outputdata/mesh.txt')
    open(12, file='outputdata/gamma.txt')
    open(13, file='outputdata/dp_coeff.txt')
    open(14, file='outputdata/d_lift_span.txt')
    open(15, file='outputdata/d_drag_span.txt')
    open(16, file='outputdata/d_lift_chord.txt')
    open(17, file='outputdata/d_drag_chord.txt')
    open(18, file='outputdata/log.txt')
    open(19, file='inputdata/prev.txt')
    write(*,*) "Output data files opened."
    write(18,*) "Output data files opened."

    write(*,*) "What do you want to continue with?"
    write(*,*) "    (p) Previous data"
    write(*,*) "    (m) Manually assigned data in code"
    write(*,*) "    (c) Command line input(stdin)"
    read(*,*) ans
    if(ans.eq.'p') then
        read(19,*) zm,p,t,b,c,alpha1,vt,rho,jb,ib
    else if(ans.eq.'m') then
! Manual input
        ib=5
        jb=15
        b=1.415
        c=1.5
        alpha1=0.
        rho=1.
        vt=5.
        zm=0.04
        p=.4
        t=.12
    else
        write(*,*) "Enter 4-digit NACA foil number:"
1       read(*,*) inaca
! zm : maximum camber for NACA profile, p : location of maximum camber, t : thickness of foil
        if(inaca.ge.1000.and.inaca.le.9999.and.inaca.ne.0) then
            do k=1,4
                idig(k)=mod(inaca,10)
                inaca=inaca/10
            end do
            zm=idig(4)/100.
            p=idig(3)/10.
            t = mod(inaca,100)/100.
        else
            write(*,*) "Please enter a valid 4-digit number: "
            goto 1
        end if
        write(*,*) "Enter wing span:"
        read(*,*) b
        write(*,*) "Enter chord:"
        read(*,*) c
        write(*,*) "Enter angle of attack (in degree):"
        read(*,*) alpha1
        write(*,*) "Now enter free stream velocity:"
        read(*,*) vt
        write(*,*) "Enter density of the medium:"
        read(*,*) rho
        write(*,*) "Put spanwise panel number and chordwise panel number respectively:" 
! input ib and jb input for any surface: upper or lower. Total panel number at upper and total panel number at lower surface are equal.
2       read(*,*) jb,ib
        if(jb.gt.jmax.or.ib.gt.imax)then
            write(*,*) "Spanwise or chordwise array limit crossed! Please enter value in limit:(smax: 50, cmax: 20)"
            goto 2
        end if
    end if
    write(19,*) zm,p,t,b,c,alpha1,vt,rho,jb,ib

    d_ib = ib*2
    d_ib1=ib*2+1
    d_ib2=ib*2+2
    jb1=jb+1
    aLpha=aLpha1*pi/180.0
! constants
    dxw=100.0*b

    call vlm

! Close opened files
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    write(*,*) "End of program."
    write(18,*) "End of program."

    close(18)

end program

subroutine vlm
    use com
    real :: gamma(imax,jmax),dL(imax,jmax),dd(imax,jmax),dp(imax,jmax),a(max,max), &
gamma1(max),dw(max),dLy(jmax),ddy(jmax),dLx(imax),ddx(imax)
    integer :: ip(max)
    write(*,*) "Subroutine vlm starts..."
    write(18,*) "Subroutine vlm starts..."

! Write heading of files
    write(12,106)
    write(13,108) alpha1
    write(13,109)
    write(14,111) alpha1,rho,vt,b/jb
    write(14,112)
    write(15,114) alpha1,rho,vt,b/jb
    write(15,115)
    write(16,117) alpha1,rho,vt,c/ib
    write(16,118)
    write(17,120) alpha1,rho,vt,c/ib
    write(17,121)

! wing geometry
    call grid

! aerodynamic calculations
    do i=1,d_ib
        do j=1,jb
! gamma(i,j)=1.0 is required for influence matrix calculations
            gamma(i,j)=1.0
        end do
    end do

    write(*,*) "Induced coefficients calculation starts..."
! influence coefficients calculation
    k=0
    do i=1,d_ib
        do j=1,jb
            k=k+1
            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,d_ib
                do j1=1,jb
                    L=L+1
! a(k,L)-is the normal velocity component due to a unit vortex lattice
                    a(k,L)=a1(i1,j1)
                end do
            end do

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
    k1=d_ib*jb
    do k=1,k1
        gamma1(k)=dw(k)
    end do

    call decomp(k1,max,a,ip)
    call solver (k1,max,a,gamma1,ip)

! wing vortex lattice listing
    k=0
    do i=1,d_ib
        do j=1,jb
            k=k+1
            gamma(i,j)=gamma1(k)
        end do
    end do
 
! forces calculation
    write(*,*) "Force calculation starts..."
    write(18,*) "Force calculation starts..."
    fL=0.0
    fd=0.0
    fm=0.0
    que=0.5*rho*vt*vt

    do j=1,jb
        dLy(j)=0.
        ddy(j)=0.
        do i=1,d_ib
            if(i.eq.ib .or. i.eq.(ib+1)) gammaij=gamma(i,j)
            if(i.gt.(ib+1)) gammaij=gamma(i,j)-gamma(i-1,j)
            if(i.lt.ib) gammaij=gamma(i,j)-gamma(i+1,j)
            dym=qf(i,j+1,2)-qf(i,j,2)
            dL(i,j)=rho*vt*gammaij*dym

            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u1,v1,w1,0.0,i,j)

            wind=w1

            dd(i,j)=-rho*gammaij*wind*dym
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
    write(10,102) b,c,ib*jb,s,ar,aLpha1,vt,rho
    write(10,103) cL, fL, cm, cd
    write(*,103) cL, fL, cm, cd
 
! Formats
101 format('Summary',/,7('-'),/)
102 format('Wing Span =',11x,f8.2,/,'Chord Length =',8x,f8.2,&
/,'Total Number of Panels =',i6,/,'Wing Area = ',10x,f8.2,/, &
'Aspect Ratio =',8x,f8.2,/,'Angle of Attack =',5x,f8.2,&
/,'Free Stream Velocity =',f8.2, &
/,'Density of Medium =',3x,f8.2)
103 format('CL = ',f10.5,/,'FL = ',f10.2,/,'CM = ',f10.5,/,'CD = ',f10.5)
! 104,105 in grid
106 format(4x,'i',4x,'j',4x,'Gamma(i,j)')
107 format(2(i5),2x,f10.5)
108 format('Alpha:'f5.1)
109 format(4x,'qh(i,j,1)',1x,'qh(i,j,2)',1x,'CdP(i,j) : Pressure difference coefficient of panel')
110 format(2(f10.3),2x,f10.5)
111 format('Alpha:'f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'dy:',1x,f6.4)
112 format(4x,'j',4x,'dLy(j)*jb/b : Lift per span length')
113 format(i5,2x,f14.3)
114 format('Alpha:'f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'dy:',1x,f6.4)
115 format(4x,'j',4x,'ddy(j)*jb/b : Drag per span length')
116 format(i5,2x,f10.5)
117 format('Alpha:'f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'dx:',1x,f6.4)
118 format(4x,'i',4x,'dLx(i)*ib/c : Lift per chord length')
119 format(i5,2x,f14.3)
120 format('Alpha:'f5.1,/,'Density:',f6.1,1x,'V(inf):',f5.1,/,'dx:',1x,f6.4)
121 format(4x,'i',4x,'ddx(i)*ib/c : Drag per chord length')
122 format(i5,2x,f10.5)
    write(*,*) "Subroutine vlm ended."
    write(18,*) "Subroutine vlm ended."

    return
end subroutine vlm

subroutine grid
    use com
    write(*,*) "Subroutine grid starts..."
    write(18,*) "Subroutine grid starts..."
    write(11,104)
104 format(5x,'xMesh',5x,'yMesh',5x,'zMesh')
105 format(f10.4,2x,f10.4,2x,f10.4)

! qf(i,j,1-3): wing fixed vortices location
    dx=c/ib
    dy=b/jb
    do j=1,jb1
        yLe=dy*(j-1)
        do i=1,d_ib1
            if(i.le.ib) then ! following portion is for lower surface
                xc = (ib+1-i)*dx
                if(xc.lt.p*c) then
                    grad = 2*zm*(p-xc/c)/p**2
                    zc = zm*xc*(2*p-xc/c)/p**2
                else
                    grad = 2*zm*(p-xc/c)/(1-p)**2
                    zc = zm*(c-xc)*(1+xc/c-2*p)/(1-p)**2
                end if
                zt = 5*t*c*(.2969*sqrt(xc/c)-.1260*xc/c-.3516*(xc/c)**2+.2843*(xc/c)**3-.1015*(xc/c)**4)
                x_temp = xc+zt*sin(atan(grad)) ! setting temporary x to apply angle of attack
                z_temp = zc-zt*cos(atan(grad)) ! setting temporary z to apply angle of attack

                qf(i,j,1) = x_temp*cos(aLpha)-z_temp*sin(aLpha) ! applaying rotation matrix for angle of attack
                qf(i,j,2) = yLe ! angle of attack have no effect on angle of attack
                qf(i,j,3) = x_temp*sin(aLpha)+z_temp*cos(aLpha) ! applying rotation matrix
                
            else ! following portion is for upper surface vertex points
                xc = (i-ib-1)*dx
                if(xc.lt.p*c) then
                    grad = 2*zm*(p-xc/c)/p**2
                    zc = zm*xc*(2*p-xc/c)/p**2
                else
                    grad = 2*zm*(p-xc/c)/(1-p)**2
                    zc = zm*(c-xc)*(1+xc/c-2*p)/(1-p)**2
                end if
                zt = 5*t*c*(.2969*sqrt(xc/c)-.1260*xc/c-.3516*(xc/c)**2+.2843*(xc/c)**3-.1015*(xc/c)**4)
                x_temp = xc-zt*sin(atan(grad)) 
                z_temp = zc+zt*cos(atan(grad)) 

                qf(i,j,1) = x_temp*cos(aLpha)-z_temp*sin(aLpha) 
                qf(i,j,2) = yLe 
                qf(i,j,3) = x_temp*sin(aLpha)+z_temp*cos(aLpha) 
            end if

        end do
        qf(d_ib2,j,1)=xte+dxw
        qf(d_ib2,j,2)=qf(d_ib1,j,2)
        qf(d_ib2,j,3)=qf(d_ib1,j,3)
    end do
    write(*,*) "Ring vortex apex points generated."
    write(18,*) "Ring vortex apex points generated."

! Generating center points of panels for contour graphing of pressure coefficients
    do j=1,jb
        yLe=dy*(j-0.5)
        dx=c/ib

        do i=1,ib+1
            qh(i,j,1)=dx*(i-0.5)
            qh(i,j,2)=yLe
        end do
    end do

! Generating wing collocation points and panel area vectors
    do j=1,jb
        do i=1,d_ib
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

    s=c*b
    ar=b/c
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
    do i=1,d_ib1
        do j=1,jb
            i3=i
            if(i.eq.d_ib1) i3=ib
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
            a1(i,j)=u0*ds(i1,j1,1)+v0*ds(i1,j1,2)+w0*ds(i1,j1,3)

            if(i.eq.d_ib1) a1(ib,j)=a1(ib,j)+a1(d_ib1,j)-a1(1,j) ! Kutta condition

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
