! 3d-vlm code for simple wing planforms with ground effect.
! This program is modified version of program written by Joe Katz, 1974.

! imax : maximum number of chordwise panels, jmax : maximum number of spanwise panels
! ib : number of chordwise panels, jb : number of spanwise panels, b : wing span
! c : chord length, ar : wing aspect ratio, ch : height above ground
! dxw : wake length, alpha : angle of attack
! qf : vortex ring corner points, qc : collocation points, ds : area vector
! xx(1) to xx(4): x-coordinates of the wing's four corner points
! a1 : auxiliary influence coefficient which is needed for calculation
! gamma : circulation, dl : partial lift force, dd : partial drag force
! dp : pressure difference, a : influence coefficient


module com
    integer, parameter :: imax=10, jmax=30, max=imax*jmax
    integer :: ib,jb,ib1,ib2,jb1,isign
    real :: b,c,s,ar,ch,dxw,alpha,pi=4.*atan(1.)
    real :: qf(imax+1,jmax+1,3),qc(imax,jmax,3),ds(imax,jmax,4),a1(imax+1,jmax),xx(4)
end module com

program main
    use com

    real :: gamma(imax,jmax),dl(imax,jmax),dd(imax,jmax),dp(imax,jmax),a(max,max),gamma1(max),dw(max),dly(jmax)
    integer :: ip(max)

    open(10, file='summary.txt')
    open(11, file='mesh.txt')


    ib=7
    jb=23
    ib1=ib+1
    ib2=ib+2
    jb1=jb+1
    xx(1)=0.0
    xx(2)=2.0
    xx(3)=2.0
    xx(4)=4.0
    b=13.0
    vt=5.0
    alpha1=5.0
    alpha=alpha1*pi/180.0
! If ch<=100, ground effect is counted
    ch= 1000 
! constants
    dxw=100.0*b
    ro=1.0

! wing geometry
    call grid
! call hgrid


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
            isign=0
            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u,v,w,1.0,i,j)
            l=0
            do i1=1,ib
                do j1=1,jb
                    l=l+1
! a(k,l)-is the normal velocity component due to a unit vortex lattice
                    a(k,l)=a1(i1,j1)
                end do
            end do

! add influence of wing's other half
            isign=1
            call wing(qc(i,j,1),-qc(i,j,2),qc(i,j,3),gamma,u,v,w,1.0,i,j)
            l=0
            do i1=1,ib
                do j1=1,jb
                    l=l+1
                    a(k,l)=a(k,l)+a1(i1,j1)
                end do
            end do

            if (ch.gt.100.0) goto 12
! add influence of mirror image (due to ground)
            isign=2
            call wing(qc(i,j,1),qc(i,j,2),-qc(i,j,3),gamma,u,v,w,1.0,i,j)
            l=0
            do i1=1,ib
                do j1=1,jb
                    l=l+1
                    a(k,l)=a(k,l)+a1(i1,j1)
                end do
            end do

! add mirror image influence of wing's other half
            isign=3
            call wing(qc(i,j,1),-qc(i,j,2),-qc(i,j,3),gamma,u,v,w,1.0,i,j)
            l=0
            do i1=1,ib
                do j1=1,jb
                    l=l+1
                    a(k,l)=a(k,l)+a1(i1,j1)
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
    fl=0.0
    fd=0.0
    fm=0.0
    que=0.5*ro*vt*vt

    do j=1,jb
        dly(j)=0.
        do i=1,ib
            if(i.eq.1) gammaij=gamma(i,j)
            if(i.gt.1) gammaij=gamma(i,j)-gamma(i-1,j)
            dym=qf(i,j+1,2)-qf(i,j,2)
            dl(i,j)=ro*vt*gammaij*dym

            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u1,v1,w1,0.0,i,j)
            call wing(qc(i,j,1),-qc(i,j,2),qc(i,j,3),gamma,u2,v2,w2,0.0,i,j)

            if(ch.gt.100.0) goto 194
            call wing(qc(i,j,1),qc(i,j,2),-qc(i,j,3),gamma,u3,v3,w3,0.0,i,j)
            call wing(qc(i,j,1),-qc(i,j,2),-qc(i,j,3),gamma,u4,v4,w4,0.0,i,j)

            goto 195

194         w3=0.0
            w4=0.0

195         wind=w1+w2-w3-w4

            dd(i,j)=-ro*dym*gammaij*wind
            dp(i,j)=dl(i,j)/ds(i,j,4)/que
            dly(j)=dly(j)+dl(i,j)
            fl=fl+dl(i,j)
            fd=fd+dd(i,j)
            fm=fm+dl(i,j)*(qf(i,j,1)-xx(1))
        end do
    end do
    cl=fl/(que*s)
    cd=fd/(que*s)
    cm=fm/(que*s*c)

! Write statements
    write(10,101)
    write(10,102) alpha1,b,c,s,ar,vt,ib,jb,ch
    write(10,103) cl, fl, cm, cd
 
! Formats
101 format(10x,'wing lift distribution calculation(with groun &
effect)',10x,54('-'))
102 format(10x,'alfa =',f10.2,8x,'b =', f10.2,8x,'c =',f13.2,&
10x,'s =',f10.2,8x,'ar = ',f10.2,8x,'v(inf) = ',f10.2,&
10x,'ib =',i10,8x,'jb =',i10,8x,'ch =',f16.2,/)
103 format(10x,'cl=',f10.4,2x,'fl=',f10.4,4x,'cm=',f10.4,3x,&
'cd=',f10.4)

! Close opened files
    close(10)
    close(11)

end program main

subroutine grid
    use com
    write(11,104)
104 format(5x,'xMesh',5x,'yMesh',5x,'zMesh')
105 format(3(f10.3))

! qf(i,j,1-3): wing fixed vortices location
    dy=b/jb
    do j=1,jb1
        yle=dy*(j-1)
        xle=xx(1)+(xx(2)-xx(1))*yle/b
        xte=xx(4)+(xx(3)-xx(4))*yle/b
        dx=(xte-xle)/ib

        do i=1,ib1
            qf(i,j,1)= (xle+dx*(i-0.75))*cos(alpha)
            qf(i,j,2)= yle
            qf(i,j,3)= -qf(i,j,1)*tan(alpha)+ch
        end do
        qf(ib2,j,1)=xte+dxw
        qf(ib2,j,2)=qf(ib1,j,2)
        qf(ib2,j,3)=qf(ib1,j,3)
    end do

! Generating coordinate points for mesh plotting
    do j=1,jb1
        yle=dy*(j-1)
        xle=xx(1)+(xx(2)-xx(1))*yle/b
        xte=xx(4)+(xx(3)-xx(4))*yle/b
        dx=(xte-xle)/ib
        if(mod(j,2).ne.0)then
            do i=1,ib1
                xme = (xle+dx*(i-1))*cos(alpha)
                yme = yle
                zme = -xme*tan(alpha)
                write(11,105) xme, yme, zme
            end do
        else
            do i=ib1,1,-1
                xme = (xle+dx*(i-1))*cos(alpha)
                yme = yle
                zme = -xme*tan(alpha)
                write(11,105) xme, yme, zme
            end do
        end if
    end do
    do i=1,ib1
        if(mod(jb,2).ne.0)then
            if(mod(i,2).ne.0)then
                do j=1,jb1
                    yle=dy*(j-1)
                    xle=xx(1)+(xx(2)-xx(1))*yle/b
                    xte=xx(4)+(xx(3)-xx(4))*yle/b
                    dx=(xte-xle)/ib
                    xme = (xle+dx*(i-1))*cos(alpha)
                    yme = yle
                    zme = -xme*tan(alpha)
                    write(11,105) xme, yme, zme
                end do
            else
                do j=jb1,1,-1
                    yle=dy*(j-1)
                    xle=xx(1)+(xx(2)-xx(1))*yle/b
                    xte=xx(4)+(xx(3)-xx(4))*yle/b
                    dx=(xte-xle)/ib
                    xme = (xle+dx*(i-1))*cos(alpha)
                    yme = yle
                    zme = -xme*tan(alpha)
                    write(11,105) xme, yme, zme
                end do
            end if
        else
            if(mod(i,2).eq.0)then
                do j=1,jb1
                    yle=dy*(j-1)
                    xle=xx(1)+(xx(2)-xx(1))*yle/b
                    xte=xx(4)+(xx(3)-xx(4))*yle/b
                    dx=(xte-xle)/ib
                    xme = (xle+dx*(i-1))*cos(alpha)
                    yme = yle
                    zme = -xme*tan(alpha)
                    write(11,105) xme, yme, zme
                end do
            else
                do j=jb1,1,-1
                    yle=dy*(j-1)
                    xle=xx(1)+(xx(2)-xx(1))*yle/b
                    xte=xx(4)+(xx(3)-xx(4))*yle/b
                    dx=(xte-xle)/ib
                    xme = (xle+dx*(i-1))*cos(alpha)
                    yme = yle
                    zme = -xme*tan(alpha)
                    write(11,105) xme, yme, zme
                end do
            end if
        end if 
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

    s=0.5*(xx(3)-xx(2)+xx(4)-xx(1))*b
    c=s/b
    ar=2.0*b*b/s

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
    if(r1.lt.rcut .or. r2.lt.rcut .or. square.lt.rcut) goto 1
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
            if(onoff.lt.0.1) goto 2
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
