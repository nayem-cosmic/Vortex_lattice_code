! 3d-VLM code for simple wing planforms with ground effect.
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


module com
    integer, parameter :: imax=20, jmax=50, max=imax*jmax
    real, parameter :: pi=4.*atan(1.)
    integer :: ib,jb,ib1,ib2,jb1,isign
    real :: b,c,s,ar,vt,rho,dxw,aLpha,zm,p,cL,cd,cm
    real :: qf(imax+1,jmax+1,3),qc(imax,jmax,3),qh(imax,jmax,2),ds(imax,jmax,4),a1(imax+1,jmax)
end module com

program main
    use com

    write(*,*) "Program starts..."
    write(21,*) "Program starts..."
    open(10, file='outputdata/alphacl.txt')
    open(13, file='outputdata/alphacm.txt')
    open(16, file='outputdata/alphacd.txt')
    open(21, file='outputdata/log.txt',access='append')
    open(22, file='inputdata/prev.txt')
    write(*,*) "Data files opened."
    write(21,*) "Data files opened."

    read(22,*) zm,p,b,c,alpha1,vt,rho,jb,ib
    aLpha=aLpha1*pi/180.0
    ib1=ib+1
    ib2=ib+2
    jb1=jb+1
! constants
    dxw=100.0*b
    write(10,100) b,c
100 format('Span:',f5.2,'Chord:',f5.2)
    write(10,*) "Alpha vs. CL"
    alpha=0
    do ic=0,40
        alpha=alpha+0.01
        call vlm
        write(10,*) alpha*180/pi,cL
    end do
    write(21,*) "Alpha vs. CL calculated."
    write(*,*) "Alpha vs. CL calculated."
    close(22)

    open(22, file='inputdata/prev.txt')
    read(22,*) zm,p,b,c,alpha1,vt,rho,jb,ib
    aLpha=aLpha1*pi/180.0
    write(13,103) b,c
103 format('Span:',f5.2,'Chord:',f5.2)
    write(13,*) "Alpha vs. CM"
    aLpha=0
    do ic=0,40
        aLpha=alpha+0.01
        call vlm
        write(13,*) alpha*180/pi,cm
    end do
    write(21,*) "Alpha vs. CM calculated."
    write(*,*) "Alpha vs. CM calculated."
    close(22)

    open(22, file='inputdata/prev.txt')
    read(22,*) zm,p,b,c,alpha1,vt,rho,jb,ib
    aLpha=aLpha1*pi/180.0
    write(16,106) b,c
106 format('Span:',f5.2,' Chord:',f5.2)
    write(16,*) "Alpha vs. CD"
    aLpha=0
    do ic=0,40
        aLpha=aLpha+0.01
        call vlm
        write(16,*) alpha*180/pi,cd
    end do
    write(21,*) "Alpha vs. CD calculated."
    write(*,*) "Alpha vs. CD calculated."
    close(22)

! Close opened files
    close(10)
    close(13)
    close(16)
    close(22)
    write(*,*) "End of program."
    write(21,*) "End of program."

    close(21)

end program

subroutine vlm
    use com
    real :: gamma(imax,jmax),dL(imax,jmax),dd(imax,jmax),dp(imax,jmax),a(max,max), &
gamma1(max),dw(max),dLy(jmax),ddy(jmax),dLx(imax),ddx(imax)
    integer :: ip(max)

! wing geometry
    call grid

! aerodynamic calculations
    do i=1,ib
        do j=1,jb
! gamma(i,j)=1.0 is required for influence matrix calculations
            gamma(i,j)=1.0
        end do
    end do

! influence coefficients calculation
    k=0
    do i=1,ib
        do j=1,jb
            k=k+1
            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,ib
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
    que=0.5*rho*vt*vt

    do j=1,jb
        dLy(j)=0.
        ddy(j)=0.
        do i=1,ib
            if(i.eq.1) gammaij=gamma(i,j)
            if(i.gt.1) gammaij=gamma(i,j)-gamma(i-1,j)
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
        end do
    end do
    cL=fL/(que*s)
    cd=fd/(que*s)
    cm=fm/(que*s*c)

    return
end subroutine vlm

subroutine grid
    use com

! qf(i,j,1-3): wing fixed vortices location
    dx=c/ib
    dy=b/jb
    do j=1,jb1
        yLe=dy*(j-1)
        do i=1,ib1
            qf(i,j,1)=(dx*(i-0.75))*cos(aLpha)
            xc=dx*(i-0.75)
            if(xc.le.p*c) then
                yc=zm*xc*(2*p-xc/c)/p**2
            else
                yc=zm*(c-xc)*(1+xc/c-2*p)/(1-p)**2
            end if
            qf(i,j,2)=yLe
            qf(i,j,3)=yc*cos(alpha)-(dx*(i-0.75))*sin(aLpha)
        end do
        qf(ib2,j,1)=xte+dxw
        qf(ib2,j,2)=qf(ib1,j,2)
        qf(ib2,j,3)=qf(ib1,j,3)
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

    s=c*b
    ar=b/c

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
            a1(i,j)=u0*ds(i1,j1,1)+v0*ds(i1,j1,2)+w0*ds(i1,j1,3)

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

    return
end subroutine  decomp

subroutine solver(n,ndim,a,b,ip)
    real :: a(ndim,ndim),b(ndim),t
    integer :: ip(ndim)
 
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

    return
end subroutine solver
