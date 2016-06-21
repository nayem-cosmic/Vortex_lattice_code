! 3d-vlm code for simple wing planforms with ground effect(by joe katz,1974).

module com
    real :: qf(6,14,3),qc(4,13,3),ds(4,13,4),a1(5,13),x(4)
    real :: b,c,s,ar,ch,dxw,alpha,phi,zm,p
    integer :: ib,jb,isign
end module com

program vlm
    use com
    real :: gamma(4,13),dL(4,13),dd(4,13),dp(4,13),a(52,52),gamma1(52),gamma1j(5),dw(52),dLy(13)
    integer :: ip(52)
    pi=4.*atan(1.)

    ib=4
    jb=13
    x(1)=0.
    x(2)=0.
    x(3)=1.
    x(4)=1.
    b=3
    zm=.02
! zm is maximum camber for NACA profil.
    p=0.4
! p is the location of maximum camber.
    vt=5.0
    alpha=5.0
    alpha=alpha*pi/180.
    phi=0;
    phi=phi*pi/180.
    ch=1000.
! x(1) to x(4) are x-coordinates of the wing's four cornerpoints.
! b - wing span, vt - free stream speed, b - wing span,
! ch - height above ground

! constants
    dxw=100.0*b
    do  i=1,ib
        do  j=1,jb
! gamma(i,j)=1.0 is required for influence matrix calculations.
            gamma(i,j)=1.0
        end do
    end do
    ro=1000.
    ib1=ib+1
    ib2=ib+2
    jb1=jb+1

! wing geometry
    call grid
    write(6,101)
    write(6,102) alpha*180/pi,b,c,s,ar,vt,ib,jb,ch

! aerodynamic calculations
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
! a(k,L) - is the normal velocity component due to a unit vortex
! lattice.
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
            if(ch.gt.100.0) goto 12
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
! add mirror image influence of wing's other half.
            isign=3
            call wing(qc(i,j,1),-qc(i,j,2),-qc(i,j,3),gamma,u,v,w,1.0,i,j)
            L=0
            do i1=1,ib
                do j1=1,jb
                    L=L+1
                    a(k,L)=a(k,L)+a1(i1,j1)
                end do
            end do

! calculate wing geometrical downwash
12          uinf=vt
            vinf=0.0
            winf=0.0

! this is the general formulation for right hand side.
            dw(k)=-(uinf*ds(i,j,1)+vinf*ds(i,j,2)+winf*ds(i,j,3))
        end do
    end do
! solution of the problem: dw(i)=a(i,j)*gamma(i)
    k1=ib*jb
    do k=1,k1
         gamma1(k)=dw(k)
    end do

    call decomp(k1,52,a,ip)
    call solver(k1,52,a,gamma1,ip)
!           here   *   the same array size is required,
! as specified in the beginning of the code

! wing vortex lattice listing
    k=0
    do i=1,ib
        do j=1,jb
            k=k+1
            gamma(i,j)=gamma1(k)
        end do
    end do

! forces calculation
    fL=0.
    fd=0.
    fm=0.
    que=0.5*ro*vt*vt
    do j=1,jb
        dLy(j)=0.
        do i=1,ib
            if(i.eq.1) gammaij=gamma(i,j)
            if(i.gt.1) gammaij=gamma(i,j)-gamma(i-1,j)
            dym=qf(i,j+1,2)-qf(i,j,2)
            dL(i,j)=ro*vt*gammaij*dym

! induced drag calculation
            isign=4
            call wing(qc(i,j,1),qc(i,j,2),qc(i,j,3),gamma,u1,v1,w1,0.0,i,j)
            call wing(qc(i,j,1),-qc(i,j,2),qc(i,j,3),gamma,u2,v2,w2,0.0,i,j)
            if(ch.gt.100.0) goto 194
            call wing(qc(i,j,1),qc(i,j,2),-qc(i,j,3),gamma,u3,v3,w3,0.0,i,j)
            call wing(qc(i,j,1),-qc(i,j,2),-qc(i,j,3),gamma,u4,v4,w4,0.0,i,j)
            goto 195
194         w3=0.
            w4=0.
195         wind=w1+w2-w3-w4
            isign=0
! add influence of mirror image (ground).

            dd(i,j)=-ro*dym*gammaij*wind
            dp(i,j)=dL(i,j)/ds(i,j,4)/que
            dLy(j)=dLy(j)+dL(i,j)
            fL=fL+dL(i,j)
            fd=fd+dd(i,j)
            fm=fm+dL(i,j)*(qf(i,j,1)-x(1))
        end do
    end do

    cL=fL/(que*s)
    cd=fd/(que*s)
    cm=fm/(que*s*c)

! output
    write(6,104) cL,2*fL,cm,cd
    write(6,105)
    write(6,110)
    do j=1,jb
        do i=2,ib
            gamma1j(i)=gamma(i,j)-gamma(i-1,j)
        end do
        dLyj=dLy(j)/b*jb
        write(6,103) j,dLyj,dp(1,j),dp(2,j),dp(3,j),dp(4,j),gamma(1,j),gamma1j(2),gamma1j(3),gamma1j(4)
    end do
! end of program

! formats
101 format(20x,'Wing lift distribution calculation (with ground effect)',/,20x,56('-'))
102 format(10x,'Alpha:',f10.2,8x,'B: ',f10.2,8x,'C: ',f13.2,/,10x,'S: ', &
& f10.2,8x,'AR: ',f10.2,8x,'V(inf): ',f10.2,/,10x,'IB: ',i10,8x,'JB: ',i10,8x,'L.E. Height: ', f8.2,/)
103 format(i3,3x,f10.3,4(f12.3),3x,4(f12.3))
104 format(10X,'CL= ',f10.4,2x,'L= ',f10.4,4x,'CM= ',f10.7,3x,'CD= ',f10.4,/)
105 format(30X,'Distribution for One Half Side of the Wing',/,118('='))
110 format(3x,'J',8x,'DL',28x,'DCP',48x,'Gamma',/,118('=') &
& ,/,24x,'I=1',9x,'I=2',9x,'I=3',9x,'I=4',12x,'I=1',9x,'I=2',9x,'I=3',9x,'I=4',/,118('='))

    stop
end program vlm

subroutine grid
    use com,only : qf,qc,ds,x,b,c,s,ar,alpha,phi,ib,jb,ch,isign,dxw,zm,p
    pi=4.*atan(1.)
! x(1) - is root l.e., x(2) tip l.e., x(3) tip t.e., and x(4) is root t.e.
! ib: no. of chordwise boxes, jb: no. of spanwise boxes
    ib1=ib+1
    ib2=ib+2
    jb1=jb+1

! wing fixed vortices location ( qf(i,j,(x,y,z))...)
    dy=b/jb
    do j=1,jb1
        yLe=dy*(j-1)
        xLe=x(1)+(x(2)-x(1))*yLe/b
        xte=x(4)+(x(3)-x(4))*yLe/b
! xle and xte are l.e. and t.e. x-coordinates
        ci=xte-xLe
! ci instantenious chord length.
        dx=(xte-xLe)/ib
        do i=1,ib1
            qf(i,j,1)=(xLe+dx*(i-0.75))*cos(alpha)
            qf(i,j,2)=yLe*cos(phi)
            xc=dx*(i-1)
            if(xc.le.p*ci) then
                yc=zm*xc*(2*p-xc/ci)/p**2
            else
                yc=zm*(ci-xc)*(1+xc/ci-2*p)/(1-p)**2
            end if
            qf(i,j,3)=yc-qf(i,j,1)*tan(alpha)+qf(i,j,2)*tan(phi)+ch
        end do

! wake far field points
        qf(ib2,j,1)=xte+dxw
        qf(ib2,j,2)=qf(ib1,j,2)
        qf(ib2,j,3)=qf(ib1,j,3)
    end do

! wing collocation points
    do j=1,jb
        do i=1,ib
            qc(i,j,1)=(qf(i,j,1)+qf(i,j+1,1)+qf(i+1,j+1,1)+qf(i+1,j,1))/4
            qc(i,j,2)=(qf(i,j,2)+qf(i,j+1,2)+qf(i+1,j+1,2)+qf(i+1,j,2))/4
            qc(i,j,3)=(qf(i,j,3)+qf(i,j+1,3)+qf(i+1,j+1,3)+qf(i+1,j,3))/4

! computation of normal vectors!
            call panel(qf(i,j,1),qf(i,j,2),qf(i,j,3),qf(i+1,j,1),qf(i+1,j,2), &
& qf(i+1,j,3),qf(i,j+1,1),qf(i,j+1,2),qf(i,j+1,3),qf(i+1,j+1,1),&
& qf(i+1,j+1,2),qf(i+1,j+1,3),ds(i,j,1),ds(i,j,2),ds(i,j,3),ds(i,j,4))
        end do
    end do

! b -is semi span, c -av. chord, s - area
    s=0.5*(x(3)-x(2)+x(4)-x(1))*b
    c=s/b
    ar=2.*b*b/s
    return
end subroutine grid 

subroutine panel(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,c1,c2,c3,sp)
! calculation of panel area and normal vector.
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

    c1=x/a
    c2=y/a
    c3=z/a

! panel area 1/2*|a x b|
    sp = a/2

! calculation of panel area
!	e1=x3-x1
!	e2=y3-y1
!	e3=z3-z1
!	f1=x2-x1
!	f2=y2-y1
!	f3=z2-z1
! normal areas (f*b+b*e)
!	s11=f2*b3-f3*b2
!	s12=b1*f3-f1*b3
!	s13=f1*b2-f2*b1
!	s21=b2*e3-b3*e2
!	s22=e1*b3-b1*e3
!	s23=b1*e2-b2*e1
!	s=0.5*(sqrt(s11**2+s12**2+s13**2)+sqrt(s21**2+s22**2+s23**2))
    return
end subroutine panel

subroutine vortex(x,y,z,x1,y1,z1,x2,y2,z2,gamma,u,v,w)
! subroutine vortex calculates the induced velocity (u,v,w) at a poi
! (x,y,z) due to a vortex element vith strength gamma per unit length
! pointing to the direction (x2,y2,z2)-(x1,y1,z1).
    pi=4.*atan(1.)
    rcut=1.0e-10

! calculation of r1 x r2
    r1r2x=(y-y1)*(z-z2)-(z-z1)*(y-y2)
    r1r2y=-((x-x1)*(z-z2)-(z-z1)*(x-x2))
    r1r2z=(x-x1)*(y-y2)-(y-y1)*(x-x2)

! calculation of (r1 x r2 )**2
    square=r1r2x*r1r2x+r1r2y*r1r2y+r1r2z*r1r2z

! calculation of r0(r1/r(r1)-r2/r(r2))
    r1=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))
    r2=sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2))
    if((r1.lt.rcut).or.(r2.lt.rcut).or.(square.lt.rcut)) goto 1 
    r0r1=(x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1)
    r0r2=(x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2)
    coef=gamma/(4.0*pi*square)*(r0r1/r1-r0r2/r2)
    u=r1r2x*coef
    v=r1r2y*coef
    w=r1r2z*coef
    goto 2

! when point (x,y,z) lies on vortex element; its induced velocity is
1   u=0.
    v=0.
    w=0.
2   continue
    return
end subroutine vortex

subroutine wing(x,y,z,gamma,u,v,w,onoff,i1,j1)
    use com,only : ds,ib,jb,ch,isign,a1,qf
	real :: gamma(4,13)

! calculates induced velocity at a point (x,y,z), due to vorticity
! distribution gamma(i,j), of semi-configuration - in a wing fixed
! coordinate system.
    u=0
    v=0
    w=0
    ib1=ib+1
    do i=1,ib1
        do j=1,jb
! i3 is wake vortex counter
            i3=i
            if(i.eq.ib1) i3=ib
            vortic=gamma(i3,j)
            if(onoff.lt.0.1) goto 2
            call vortex(x,y,z,qf(i,j,1),qf(i,j,2),qf(i,j,3),qf(i,j+1,1),qf(i,j1+1,2),qf(i,j+1,3),vortic,u1,v1,w1)
            call vortex(x,y,z,qf(i+1,j+1,1),qf(i+1,j+1,2),qf(i+1,j+1,3),qf(i+1,j,1),qf(i+1,j,2),qf(i+1,j,3),vortic,u3,v3,w3)
2           call vortex(x,y,z,qf(i,j+1,1),qf(i,j+1,2),qf(i,j+1,3),qf(i+1,j+1,1),qf(i+1,j+1,2),qf(i+1,j+1,3),vortic,u2,v2,w2)
            call vortex(x,y,z,qf(i+1,j,1),qf(i+1,j,2),qf(i+1,j,3),qf(i,j,1),qf(i,j,2),qf(i,j,3),vortic,u4,v4,w4)

            u0=u2+u4+(u1+u3)*onoff
            v0=v2+v4+(v1+v3)*onoff
            w0=w2+w4+(w1+w3)*onoff

! influence coefficient at (x,y,z) for (i,j) panel
! when i .eq. ib1 a1(ib1,j)=a1(ib,j)+a1(ib1,j)
            if(isign.eq.0) then
                if(i.lt.ib1) then
                    a1(i,j)=u0*ds(i1,j1,1)+v0*ds(i1,j1,2)+w0*ds(i1,j1,3)
                else if(i.eq.ib1) then
                    a1(ib1,j)=u0*ds(ib,j,1)+v0*ds(ib,j,2)+w0*ds(ib,j,3)+u0*ds(ib1,j,1)+v0*ds(ib1,j,2)+w0*ds(ib1,j,3)
                end if

            else if(isign.eq.1) then
                if(i.lt.ib1) then
                    a1(i,j)=u0*ds(i1,j1,1)-v0*ds(i1,j1,2)+w0*ds(i1,j1,3)
                else if(i.eq.ib1) then
                    a1(ib1,j)=u0*ds(ib,j,1)-v0*ds(ib,j,2)+w0*ds(ib,j,3)+u0*ds(ib1,j,1)-v0*ds(ib1,j,2)+w0*ds(ib1,j,3)
                end if

            else if(isign.eq.2) then
                if(i.lt.ib1) then
                    a1(i,j)=u0*ds(i1,j1,1)+v0*ds(i1,j1,2)-w0*ds(i1,j1,3)
                else if(i.eq.ib1) then
                    a1(ib1,j)=u0*ds(ib,j,1)+v0*ds(ib,j,2)-w0*ds(ib,j,3)+u0*ds(ib1,j,1)+v0*ds(ib1,j,2)-w0*ds(ib1,j,3)
                end if

            else if(isign.eq.3) then
                if(i.lt.ib1) then
                    a1(i,j)=u0*ds(i1,j1,1)-v0*ds(i1,j1,2)-w0*ds(i1,j1,3)
                else if(i.eq.ib1) then
                    a1(ib1,j)=u0*ds(ib,j,1)-v0*ds(ib,j,2)-w0*ds(ib,j,3)+u0*ds(ib1,j,1)-v0*ds(ib1,j,2)-w0*ds(ib1,j,3)
                end if
            else if(isign.eq.4) then
! Do nothing. While caculating drag influence coefficient is not needed.
            end if

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

! matrix triangularization by gaussian elimination.
! n = order of matrix. ndim = declared dimension of array a.
! a = matrix to be triangularized.
! ip(k) , k .lt. n = index of k-th pivot row.

    ip(n) = 1
    do k=1,n
        if(k.eq.n) goto 5
        kp1=k+1
        m=k
        do  i=kp1, n
            if( abs(a(i,k)).gt.abs(a(m,k))) m=i
        end do
        ip(k) = m
        if(m.ne.k) ip(n) = -ip(n)
        t = a(m,k)
        a(m,k) = a(k,k)
        a(k,k) = t
        if(t.eq.0.e0) go to 5
        do i=kp1, n
            a(i,k) = -a(i,k)/t
        end do
        do j=kp1, n
            t = a(m,j)
            a(m,j) = a(k,j)
            a(k,j) = t
            if(t .eq. 0.e0) go to 4
            do i=kp1, n
                a(i,j) = a(i,j) + a(i,k)*t
            end do
4       end do
5       if(a(k,k) .eq. 0.e0) ip(n) = 0
    end do

    return
end subroutine decomp

subroutine solver(n,ndim,a,b,ip)

	real :: a(ndim,ndim), b(ndim), t
    integer :: ip(ndim)

! solution of linear system, a*x = b.
! n = order of matrix.
! ndim = declared dimension of the array a.
! b = right hand side vector.
! ip = pivot vector obtained from subroutine decomp.
! b = solution vector, x.

    if(n.eq.1) goto 9
    nm1=n-1
    do k=1,nm1
        kp1=k+1
        m = ip(k)
        t = b(m)
        b(m) = b(k)
        b(k) = t
        do i=kp1, n
            b(i) = b(i) + a(i,k)*t
        end do
    end do
    do kb=1,nm1
        km1=n-kb
        k=km1+1
        b(k) = b(k)/a(k,k)
        t = -b(k)
        do i=1,km1
            b(i) = b(i) + a(i,k)*t
        end do
    end do
9   b(1) = b(1)/a(1,1)


    return
end subroutine solver
