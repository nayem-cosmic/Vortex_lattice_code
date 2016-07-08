!      RECTANGULAR LIFTING SURFACE (VLM)
!      3D-VLM CODE FOR SIMPLE WING PLANFORMS WITH GROUND EFFECT
       module para
       integer,parameter::IN = 1
       integer,parameter:: mmax=22, nmax=22, max=mmax*nmax
       real,parameter:: PAY=3.141592654

          integer:: IB, JB, IB1, JB1, IB2
       real:: DS(mmax,nmax,4), XX(mmax),B,C,S,AR,SN1,CS1
       real:: ALPHA1, CROOT, CTIP, XTIP, ZTIP
       real:: QF(mmax,nmax,3),QC(mmax,nmax,3),DXW
       real:: CH,SIGN,A1(mmax,nmax)
       end module para


       PROGRAM MAIN
       use para
!       implicit none

       REAL:: GAMA(mmax,nmax),DL(mmax,nmax),DD(mmax,nmax),DP(mmax,nmax)
       REAL:: A(max,max),GAMA1(max),DW(max)
       REAL:: DLY(nmax)
       INTEGER:: IP(max)


       open(6,FILE='GRID.OUT')
       OPEN(7,FILE='CTLP.OUT')
       OPEN(8,FILE='INFLU.OUT')
       OPEN(9,FILE='VEL.OUT')
       OPEN(10,FILE='MATRX.OUT')
       OPEN(20, FILE='CP.DAT')


       open(3,FILE='CAD.DAT')

!      INPUT DATA FOR THE HYDROFOIL
!      IB: NO. OF CHORDWISE PANELS, JB: NO.OF SPANWISE PANELS
!      XX(1) TO XX(4): X-COORDINATES OF THE WING'S FOUR CORNER POINTS
!      B:  WING SPAN, VT: FREE STREAM SPEED, CH: HEIGHT ABOVE GROUND

       IB=4
       JB=13

       XX(1)=0.0
       XX(2)=1.91
       XX(3)=2.41
       XX(4)=1.5
       B=1.415
       VT=5.0
       ALPHA1=5.354
       CH= 1000        !100.0

!      CONSTANTS

       DXW=100.0*B

       RO=1.0
       ALPHA=ALPHA1*PAY/180.0
       SN1=SIN(ALPHA)
       CS1=COS(ALPHA)

       IB1=IB+1
       IB2=IB+2
       JB1=JB+1

!      WING GEOMETRY

       CALL GRID
!      CALL HGRID

       WRITE(6,101)
       WRITE(6,102) ALPHA1,B,C,S,AR,VT,IB,JB,CH


!      AERODYNAMIC CALCULATIONS
       DO I=1,IB
          DO J=1,JB
!         GAMA(I,J)=1.0 IS REQUIRED FOR INFLUENCE MATRIX CALCULATIONS
          GAMA(I,J)=1.0
          ENDDO
       ENDDO

!      INFLUENCE COEFFICIENTS CALCULATION

       WRITE(8,25)
       WRITE(9,25)

 25    FORMAT(2X'I',5X,'J',10X,'U',10X,'V',13X,'W')
 26    FORMAT(I3,3X,I3,3X,3(F10.3,3 X))

       K=0
       DO 15 I=1,IB
          DO 14 J=1,JB
             SIGN=0.0
             K=K+1

             CALL WING(QC(I,J,1),QC(I,J,2),QC(I,J,3),GAMA,U,V,W,1.0,I,J)
             WRITE(8,26)I,J,U,V,W

             L=0
             DO I1=1,IB
                DO J1=1,JB
                L=L+1
!               A(K,L)-IS THE NORMAL VELOCITY COMPONENT DUE TO A UNIT
!               VORTEX LATTICE
                A(K,L)=A1(I1,J1)
                ENDDO
             END DO

!            ADD INFLUENCE OF WING'S OTHER HALF
             CALL WING(QC(I,J,1),-QC(I,J,2),QC(I,J,3),GAMA,U,V,W,1.0,I,J)
             WRITE(9,26)I,J,U,V,W

             L=0
             DO I1=1,IB
                DO J1=1,JB
                   L=L+1
                   A(K,L)=A(K,L)+A1(I1,J1)
                 END DO
             END DO


             IF (IN.EQ.2)THEN

             IF (CH.GT.100.0) GOTO 12
!            ADD INFLUENCE OF MIRROR IMAGE (DUE TO GROUND)

             SIGN=10.0
             CALL WING(QC(I,J,1),QC(I,J,2),-QC(I,J,3),GAMA,U,V,W,1.0,I,J)

             L=0
             DO I1=1,IB
                DO J1=1,JB
                   L=L+1
                   A(K,L)=A(K,L)+A1(I1,J1)
                END DO
             END DO


!             ADD MIRROR IMAGE INFLUENCE OF WING'S OTHER HALF
             CALL WING(QC(I,J,1),-QC(I,J,2),-QC(I,J,3),GAMA,U,V,W,1.0,I,J)

             L=0
             DO I1=1,IB
                DO J1=1,JB
                   L=L+1
                   A(K,L)=A(K,L)+A1(I1,J1)
                END DO
             END DO

             SIGN=0.0
       12    CONTINUE

             ENDIF


!      CALCULATE WING GEOMETRICAL DOWNWASH

       UINF=VT
       VINF=0.0
       WINF=0.0

       DW(K)=-(UINF*DS(I,J,1)+VINF*DS(I,J,2)+WINF*DS(I,J,3))

14     CONTINUE
       WRITE(8,*)
       WRITE(9,*)

15     CONTINUE


!      SOLUTION OF THE PROBLEM: DW(I) = A(I,J)*GAMA(I)

       K1=IB*JB
       DO K=1,K1
          GAMA1(K)=DW(K)
       END DO

       CALL DECOMP(K1,max,A,IP)
       CALL SOLVER (K1,max,A,GAMA1,IP)


!      WING VORTEX LATTICE LISTING
       OPEN(11, FILE='GAMA.DAT')
       WRITE(11,60)
60     FORMAT(3X,'I',6X'J'6X,'GAMA(I,J)')

       K=0
       DO I=1,IB
         DO J=1,JB
            K=K+1
            GAMA(I,J)=GAMA1(K)
            WRITE(11,103)I,J,GAMA(I,J)
         END DO
       END DO

!      FORCES CALCULATION

       FL=0.0
       FD=0.0
       FM=0.0
       QUE=0.5*RO*VT*VT

       DO 20 J=1,JB
          DLY(J)=0.
          DO I=1,IB

             IF(I.EQ.1) GAMAIJ=GAMA(I,J)
             IF(I.GT.1) GAMAIJ=GAMA(I,J)-GAMA(I-1,J)

             DYM=QF(I,J+1,2)-QF(I,J,2)
             DL(I,J)=RO*VT*GAMAIJ*DYM

             CALL WING(QC(I,J,1),QC(I,J,2),QC(I,J,3),GAMA,U1,V1,W1,0.0,I,J)
             CALL WING(QC(I,J,1),-QC(I,J,2),QC(I,J,3),GAMA,U2,V2,W2,0.0,I,J)

             GO TO 194
             IF(CH.GT.100.0) GOTO 194
             CALL WING(QC(I,J,1),QC(I,J,2),-QC(I,J,3),GAMA,U3,V3,W3,0.0,I,J)
             CALL WING(QC(I,J,1),-QC(I,J,2),-QC(I,J,3),GAMA,U4,V4,W4,0.0,I,J)

             GOTO 195

194          W3=0.0
             W4=0.0

195          WIND=W1+W2-W3-W4

!             ADD INFLUENCE OF MIRROR (GROUND)

             ALFI=-WIND/VT
             DD(I,J)=RO*DYM*VT*GAMAIJ*ALFI

             DP(I,J)=DL(I,J)/DS(I,J,4)/QUE
             DLY(J)=DLY(J)+DL(I,J)
             FL=FL+DL(I,J)
             FD=FD+DD(I,J)
             FM=FM+DL(I,J)*(QF(I,J,1)-XX(1))
          END DO
20        CONTINUE

          CL=FL/(QUE*S)
          CD=FD/(QUE*S)
          CM=FM/(QUE*S*C)
            WRITE(6,104) CL,FL,CM,CD

          WRITE(20,110)
          DO J=1,JB
             DO I=1,IB
             WRITE(20,103) I, J, QC(I,J,1),QC(I,J,2),QC(I,J,3), DP(I,J), GAMA(I,J)
             END DO
             WRITE(20,*)
          END DO

!         END  PROGRAM
100       CONTINUE


!          CALL CAD
!         FORMATS

101       FORMAT(1X,//10X,'WING LIFT DISTRIBUTION CALCULATION(WITH GROUND    &
                 EFFECT)',/10X,54('-'))

102       FORMAT(1X,/,10X,'ALFA =',F10.2,8X,'B =', F10.2,8X,'C =',F13.2,     &
                    /,10X,'S =',F10.2,8X,'AR = ',F10.2,8X,'V(INF) = ',F10.2, &
                    /,10X,'IB =',I10,8X,'JB =',I10,8X,'CH =',F16.2,/)

103       FORMAT(1X,I3,3X,I3,5(F15.3))

104       FORMAT(/,10X,'CL=',F10.4,2X,'FL=',F10.4,4X,'CM=',F10.4,3X, &
                 'CD=',F10.4)

110       FORMAT(//,3X,82('='),/,3X,'I',5X,'J',7X,'QC(I,J,1)',7X, &
& 'QC(I,J,2)',7X,'QC(I,J,3)',6X,'DCP(I,J)',4X,'GAMA(I,J)',4X,/,3X,82('='))

112       FORMAT(1X,'QF(I=',I2,',J,X,Y,Z)= ',15(F6.1))

113       FORMAT(1X,110('='))

       END PROGRAM MAIN



       SUBROUTINE GRID
       use para
!      INPUT DATA FOR THE PLANAR WING SURFACE
!      QF(I,J,1-3): WING FIXED VORTICES LOCATION


100  FORMAT(I3,2X,I3,4F15.4)
110  FORMAT(2X,'I',4X,'J', 8X,'QF(I,J,1)', 6X,'QF(I,J,2)',5X,'QF(I,J,3)')
120     FORMAT(2X,'J',6X,'YLE',7X,'XLE',7X,'XTE',7X,'DX')
130  FORMAT( I3,4F10.3)

200  FORMAT(I3,2X,I3,4F15.4)


       WRITE(7,120)

       DY=B/JB

       DO 3 J=1,JB1
          YLE=DY*(J-1)
          XLE=XX(1)+(XX(2)-XX(1))*YLE/B
          XTE=XX(4)+(XX(3)-XX(4))*YLE/B
          DX=(XTE-XLE)/IB

          WRITE(7,130) J, YLE, XLE, XTE, DX

          IF (J.EQ.1) WRITE(6,110)

          DO I=1,IB1
             QF(I,J,1)= (XLE+DX*(I-0.75))*CS1
             QF(I,J,2)= YLE
             QF(I,J,3)= -QF(I,J,1)*SN1/CS1+CH

             WRITE(6,200)I,J,QF(I,J,1),QF(I,J,2),QF(I,J,3)

          END DO
            QF(IB2,J,1)=XTE+DXW
          QF(IB2,J,2)=QF(IB1,J,2)
          QF(IB2,J,3)=QF(IB1,J,3)

          WRITE(6,200) I, J, QF(IB2,J,1),QF(IB2,J,2),QF(IB2,J,3)
          WRITE(6,*)

3      CONTINUE

       WRITE(6,300)
300    FORMAT(/,2X,'I',4X,'J',8X,'QC(I,J,1)',6X,'QC(I,J,2)',5X,'QC(I,J,3)')

       WRITE(7,400)
400    FORMAT(/2X,'I',4X,'J',8X,'DS(I,J,1)',6X,'DS(I,J,2)',6X,'DS(I,J,3)',5X,'DS(I,J,4)')

       DO 4 J=1,JB
          DO I=1,IB
             QC(I,J,1)=(QF(I,J,1)+QF(I,J+1,1)+QF(I+1,J+1,1)+QF(I+1,J,1))/4
             QC(I,J,2)=(QF(I,J,2)+QF(I,J+1,2)+QF(I+1,J+1,2)+QF(I+1,J,2))/4
             QC(I,J,3)=(QF(I,J,3)+QF(I,J+1,3)+QF(I+1,J+1,3)+QF(I+1,J,3))/4

             WRITE(6,200)I,J,QC(I,J,1),QC(I,J,2),QC(I,J,3)


             CALL PANEL(QF(I,J,1),QF(I,J,2),QF(I,J,3),QF(I+1,J,1),QF(I+1,J,2), &
             QF(I+1,J,3),QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),QF(I+1,J+1,1),    &
             QF(I+1,J+1,2),QF(I+1,J+1,3),DS(I,J,1),DS(I,J,2),DS(I,J,3),DS(I,J,4))

             WRITE(7,100)I,J,DS(I,J,1),DS(I,J,2),DS(I,J,3),DS(I,J,4)

          END DO
          WRITE(6,*)
          WRITE(7,*)
4       CONTINUE

       S=0.5*(XX(3)-XX(2)+XX(4)-XX(1))*B
       C=S/B
       AR=2.0*B*B/S

       WRITE(6,*)' END OF DATA FROM SUBROUTINE GRID'

       RETURN
       END



       SUBROUTINE PANEL(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,V1,V2,V3,S)

!       CALCULATION OF PANEL AREA AND VECTOR
       A1=X2-X3
       A2=Y2-Y3
       A3=Z2-Z3

       B1=X4-X1
       B2=Y4-Y1
       B3=Z4-Z1

!      NORMAL VECTOR
       X=A2*B3-A3*B2
       Y=B1*A3-A1*B3
       Z=A1*B2-A2*B1
       A=SQRT(X**2+Y**2+Z**2)

       V1=X/A
       V2=Y/A
       V3=Z/A

!      CALCULATION OF PANEL AREA
       E1=X3-X1
       E2=Y3-Y1
       E3=Z3-Z1

       F1=X2-X1
       F2=Y2-Y1
       F3=Z2-Z1

!      NORMAL AREAS (F*B+B*E)
       S11=F2*B3-F3*B2
       S12=B1*F3-F1*B3
       S13=F1*B2-F2*B1
       S21=B2*E3-B3*E2
       S22=E1*B3-B1*E3
       S23=B1*E2-B2*E1
       S=0.5*(SQRT(S11**2+S12**2+S13**2)+SQRT(S21**2+S22**2+S23**2))
       RETURN
       END




       SUBROUTINE VORTEX(P,Y,Z,X1,Y1,Z1,X2,Y2,Z2,GAMA,U,V,W)

       PAY=3.141592654
       RCUT=1.0E-10

       R1R2X=(Y-Y1)*(Z-Z2)-(Z-Z1)*(Y-Y2)
       R1R2Y=-((P-X1)*(Z-Z2)-(Z-Z1)*(P-X2))
       R1R2Z=(P-X1)*(Y-Y2)-(Y-Y1)*(P-X2)

       SQUARE=R1R2X*R1R2X+R1R2Y*R1R2Y+R1R2Z*R1R2Z

       R1=SQRT((P-X1)*(P-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1))
       R2=SQRT((P-X2)*(P-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2))
       IF(R1.LT.RCUT) GOTO 1
       IF(R2.LT.RCUT) GOTO 1
       IF(SQUARE.LT.RCUT) GOTO 1
       R0R1=(X2-X1)*(P-X1)+(Y2-Y1)*(Y-Y1)+(Z2-Z1)*(Z-Z1)
       R0R2=(X2-X1)*(P-X2)+(Y2-Y1)*(Y-Y2)+(Z2-Z1)*(Z-Z2)
       COEF=GAMA/(4.0*PAY*SQUARE)*(R0R1/R1-R0R2/R2)
       U=R1R2X*COEF
       V=R1R2Y*COEF
       W=R1R2Z*COEF
       GOTO 2
1      U=0.0
       V=0.0
       W=0.0
2      CONTINUE
       RETURN
       END




       SUBROUTINE WING(P,Y,Z,GAMA,U,V,W,ONOFF,I1,J1)
       use para
       DIMENSION GAMA(mmax,nmax)

       WRITE(10,25)
 25    FORMAT(4X'I',4X,'J',6X,'U1',8X,'V1',8X,'W1',8X,'U2',8X,'V2',8X,'W2'    &
              8X,'U3',8X,'V3',8X,'W3'8X,'U4',8X,'V4',8X,'W4',8X,'A1')
 26    FORMAT(2I5,13F10.3)

       U=0.0
       V=0.0
       W=0.0
!       IB1=IB+1

       DO I=1,IB1
          DO J=1,JB
             I3=I
             IF(I.EQ.IB1) I3=IB
             VORTIC=GAMA(I3,J)
             IF(ONOFF.LT.0.1) GOTO 2
             CALL VORTEX(P,Y,Z,QF(I,J,1),QF(I,J,2),QF(I,J,3),QF(I,J+1,1),    &
                  QF(I,J+1,2),QF(I,J+1,3),VORTIC,U1,V1,W1)
             CALL VORTEX(P,Y,Z,QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),QF(I+1,J,1),&
                  QF(I+1,J,2),QF(I+1,J,3),VORTIC,U3,V3,W3)
2            CALL VORTEX(P,Y,Z,QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),QF(I+1,J+1,1), &
                     QF(I+1,J+1,2),QF(I+1,J+1,3),VORTIC,U2,V2,W2)
             CALL VORTEX(P,Y,Z,QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),QF(I,J,1),&
                  QF(I,J,2),QF(I,J,3),VORTIC,U4,V4,W4)

             U0=U2+U4+(U1+U3)*ONOFF
             V0=V2+V4+(V1+V3)*ONOFF
             W0=W2+W4+(W1+W3)*ONOFF

!            Magnitude of the influence co-efficient
              A1(I,J)=U0*DS(I1,J1,1)+V0*DS(I1,J1,2)+W0*DS(I1,J1,3)
              IF(SIGN.GE.1.0) A1(I,J)=U0*DS(I1,J1,1)+V0*DS(I1,J1,2)-W0*DS(I1,J1,3)

             IF(I.EQ.IB1) A1(IB,J)=A1(IB,J)+A1(IB1,J)

             U=U+U0
             V=V+V0
             W=W+W0

          WRITE(10,26)I,J,U1,V1,W1,U2,V2,W2,U3,V3,W3,U4,V4,W4,A1(I,J)

          END DO
       END DO
       RETURN
       END




       SUBROUTINE DECOMP(N,NDIM,A,IP)
       REAL:: A(NDIM,NDIM),T
       INTEGER:: IP(NDIM)

       IP(N)=1
       DO 6 K=1,N
          IF(K.EQ.N) GOTO 5
          KP1=K+1
          M=K
          DO I=KP1,N
             IF(ABS(A(I,K)).GT.ABS(A(M,K))) M=I
          END DO
          IP(K)=M
          IF(M.NE.K) IP(N)=-IP(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
          IF(T.EQ.0.E0) GOTO 5
          DO I=KP1,N
             A(I,K)=-A(I,K)/T
          END DO
          DO 4 J=KP1,N
             T=A(M,J)
             A(M,J)=A(K,J)
             A(K,J)=T
             IF(T.EQ.0.E0) GOTO 4
             DO I=KP1,N
                A(I,J)=A(I,J)+A(I,K)*T
             END DO
4         CONTINUE
5         IF(A(K,K).EQ.0.E0) IP(N)=0
6      CONTINUE
       RETURN
       END




       SUBROUTINE SOLVER(N,NDIM,A,B,IP)
       REAL:: A(NDIM,NDIM),B(NDIM),T
       INTEGER:: IP(NDIM)

       IF(N.EQ.1) GOTO 9
       NM1=N-1
       DO 7 K=1,NM1
          KP1=K+1
          M=IP(K)
          T=B(M)
          B(M)=B(K)
          B(K)=T
          DO I=KP1,N
             B(I)=B(I)+A(I,K)*T
          END DO
7      CONTINUE

       DO 8 KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
          DO I=1,KM1
             B(I)=B(I)+A(I,K)*T
          END DO
8      CONTINUE
9      B(1)=B(1)/A(1,1)
       RETURN
       END

!
!
!       subroutine cad
!       use para
!       integer:: i,j, k,n
!       integer:: i1,i2,i3,i4
!          integer :: k1, k2, k3, k4
!       integer:: n1(max),n2(max),n3(max),n4(max)
!       real:: x(max),y(max),z(max)
!
!       write(3,22)
! 22       format(3x,'    NODAL POINTS FROM THE CAD')
!
!       n=0
!       do k=1,JB
!       n=n+1
!         do i=1,IB
!           k1=i     + (n-1)*IB1
!           k2=i + 1 + (n-1)*IB1
!           k3=i + 1 + n* IB1
!           k4=i     + n*IB1
!
!           j = (n-1)*IB + i
!
!           n1(j)= k1
!           n2(j)= k2
!           n3(j)= k3
!           n4(j)= k4
!           write(3,32) j, k1,k2, k3, k4
! 32     format(i5,4i12)
!        end do
!       end do
!
!       write(3,33)
!33       format(/,3x,' No.'8x,'x'15x,'y'12x,'z')
!
!       n=0
!       do j=1,JB1
!       n=n+1
!         do i=1,IB1
!            k =(n-1)*IB1+i
!          x(k)=QF(I,J,1)
!          y(k)=QF(I,J,2)
!          z(k)=QF(I,J,3)
!          write(3,27) k, x(k), y(k), z(k)
! 27         format(i5,3f15.6)
!         end do
!       end do
!
!!      For body surface
!       open(2,file='body.scr',status='unknown')
!       write (2, '(A)') 'PFACE'
!       do i=1,max
!          call scrpoint(2, x(i), y(i), z(i))
!       end do
!       do j=1,IB*JB
!          i1=n1(j)
!          i2=n2(j)
!          i3=n3(j)
!          i4=n4(j)
!          call scrquadr(2, i1, i2, i3, i4)
!       end do
!       write(2,'(/)')
!       close(2)
!       return
!       end
!
!
!
!
!       subroutine scrpoint (un, x1, y1, z1)
!       implicit none
!       integer:: un
!       real:: x1, y1, z1
!       integer:: m, l
!       character line*80, lin*40
!
!6      format(2(f12.5, ','),f12.5)
!       write (lin, 6) x1, y1, z1
!
!       m = 0
!       do l = 1, 40
!        if (lin(l:l).ne.' ') then
!          m = m + 1
!          line(m:m) = lin(l:l)
!        end if
!       end do
!       write (un,'(a)') line(1:m)
!       return
!       end
!
!
!
!       subroutine scrquadr (un, p1, p2, p3, p4)
!       implicit none
!       integer:: un, p1, p2, p3, p4
!       integer:: m, l, pp(4), i
!       character lin*10, line*10
!
!6      format(i6)
!
!       pp(1) = p1
!       pp(2) = p2
!       pp(3) = p3
!       pp(4) = p4
!       write (un, '(1x)')
!       do i = 1, 4
!         write (lin, 6) pp(i)
!         m = 0
!         do l = 1, 6
!           if (lin(l:l).ne.' ') then
!             m = m + 1
!             line(m:m) = lin(l:l)
!           end if
!         end do
!         write (un,'(a)') line(1:m)
!       end do
!       return
!       end
!
!
!
!      SUBROUTINE HGRID
!      use para
!!     DATA GENERATION FOR NACA-0012 HYDROFOIL
!!     IB: NO. OF CHORDWISE PANELS
!!     JB: NO.OF SPANWISE PANELS
!!     QF(I,J,1 TO 3): WING FIXED VORTICES LOCATION
!
!      CROOT=1.0
!      CTIP=1.0
!      XTIP=0.0
!      ZTIP=0.0
!
!      WRITE(6,9)
!9     FORMAT(/,3X,'AIRFOIL COORDINATES',/,3X,19('='))
!
!      OPEN(5,FILE='NACA-0012.DAT')
!!     OPEN(5,FILE='REVERSE.DAT')
!
!      READ(5,*) IB1
!      WRITE(6,10) IB1
!
!      WRITE(6,11)
!      DO I=1,IB1
!         READ(5,*)   QF(I,1,1),QF(I,1,3)
!         WRITE(6,12) QF(I,1,1),QF(I,1,3)
!      END DO
!
!10    FORMAT(/,6X,'IB1=',I7)
!11    FORMAT(/,6X,'X',8X,'Z')
!
!12    FORMAT(3F10.4)
!100   FORMAT(/,3X,'I',4X,'J',6X,'QF(I,J,1)',5X,'QF(I,J,2)',5X,'QF(I,J,3)')
!110   FORMAT(/,3X,'I',4X,'J',7X,'QC(I,J,1)',5X,'QC(I,J,2)',5X,'QC(I,J,3)')
!120      FORMAT(2X,'J',12X,'Y',7X,'DXLE',6X,'DZLE',6X,'CHORD')
!
!200   FORMAT(I3,2X,I3,4F15.4)
!500   FORMAT(I3,5X,4F10.3)
!
!
!!     CALCULATE PANEL CORNER POINTS; QF(I,J,(X,Y,Z))
!      WRITE(7,120)
!
!      DO 3 J=1,JB1
!         Y=B/2.0/JB* FLOAT(J-1)
!         DXLE=XTIP*Y/B
!         DZLE=ZTIP*Y/B
!         CHORD=CROOT-(CROOT-CTIP)*Y/B
!         WRITE(7,500)J,Y,DXLE,DZLE,CHORD
!
!
!!        B- SEMI WING SPAN, DXLE-LOCAL SWEEP
!
!         WRITE(6,100)
!
!         DO I=1,IB1
!            QF(I,J,1)=QF(I,1,1)*CHORD+DXLE
!            QF(I,J,2)=Y
!            QF(I,J,3)=QF(I,1,3)*CHORD+DZLE
!
!            WRITE(6,200) I, J, QF(I,J,1), QF(I,J,2), QF(I,J,3)
!
!         END DO
!
!!        WAKE FAR FIELD POINTS (QF-IS IN BODY FRAME OF REFERENCE)
!
!         QF(IB2,J,1)=QF(IB1,J,1)+DXW
!         QF(IB2,J,2)=QF(IB1,J,2)
!         QF(IB2,J,3)=QF(IB1,J,3)
!
!         WRITE(6,200) I, J, QF(IB2,J,1), QF(IB2,J,2), QF(IB2,J,3)
!
!3     CONTINUE
!
!
!300    FORMAT(/,3X,'I',4X,'J',6X,'QC(I,J,1)',5X,'QC(I,J,2)',5X,'QC(I,J,3)')
!
!       WRITE(7,400)
!400    FORMAT(/2X,'I',4X,'J',8X,'DS(I,J,1)',6X,'DS(I,J,2)',6X,'DS(I,J,3)',5X,'DS(I,J,4)')
!
!!------> WING COLLOCATION POINTS
!
!         DO J=1,JB
!         WRITE(6,300)
!           DO I=1,IB1
!              QC(I,J,1)=(QF(I,J,1)+QF(I,J+1,1)+QF(I+1,J+1,1)+QF(I+1,J,1))/4.0
!              QC(I,J,2)=(QF(I,J,2)+QF(I,J+1,2)+QF(I+1,J+1,2)+QF(I+1,J,2))/4.0
!              QC(I,J,3)=(QF(I,J,3)+QF(I,J+1,3)+QF(I+1,J+1,3)+QF(I+1,J,3))/4.0
!
!              WRITE(6,200)I,J,QC(I,J,1),QC(I,J,2),QC(I,J,3)
!
!!------>      COMPUTATION OF CHORDWISE VECTORS DS(I,J,1-3), TANGENTIAL AND
!!------>      NORMAL VECTORS DS(I,J,4 TO 9),PANEL AREA DS(I,J,1-10)
!!------>      AND SOURCE STRENGTH (SIGMA)
!
!            CALL PANEL(QF(I,J,1),QF(I,J,2),QF(I,J,3),QF(I+1,J,1),QF(I+1,J,2),  &
!                 QF(I+1,J,3),QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),QF(I+1,J+1,1),&
!                 QF(I+1,J+1,2),QF(I+1,J+1,3),DS(I,J,1),DS(I,J,2),DS(I,J,3),    &
!                 DS(I,J,4))
!               WRITE(7,200)I,J,DS(I,J,1),DS(I,J,2),DS(I,J,3),DS(I,J,4)
!
!           END DO
!           WRITE(6,*)
!           WRITE(7,*)
!         END DO
!
!!------> B-IS SEMI WING SPAN,C-ROOT CHORD,S-AREA
!
!         S=0.5*B*(CROOT+CTIP)
!         C=S/B
!         AR=B*B/S
!
!
!     WRITE(6,*) 'END OF DATA FROM SUBROUTINE HGRID'
!
!     RETURN
!     END
!
!
!
