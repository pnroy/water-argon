
      SUBROUTINE H2OAr_3DPES(XPHI,TH,RR,V,iv)
!c******************************************************************************
!   H2O-Ar PES 
!** Subroutine to generate values of the 3D-MLR analyic vibrationally averaged
!   PESs for complexes formed between Ar and H2O, at given vibrationally quantum-
!   state-specific(v1,v2,v3) of H2O [paper,submitted 2015]  
!***********************************************************************
!** Input variables:
!  isAr-40 identifies the Ar isotope 
!  isO - 16 identifies the O isotope 
!  isH - 1 identifies the H isotope 
!  iv=1 for (v1,v2,v3)=(0,0,0); 
!  iv=2 for (v1,v2,v3)=(1,0,0); 
!  iv=3 for (v1,v2,v3)=(0,1,0); 
!  iv=4 for (v1,v2,v3)=(0,0,1). 
!  RR  - distance between Ar and H2O centre of mass in [Angst]
!  TH  - Jacobi angular coordinate 'theta' in degrees TH=0 corresponding Ar-X-O
!         X is the mess center of H2O, Ar at the O end 
!  PHI  - Jacobi angular coordinate 'phi' in degrees PHI=0 corresponding 
!  the angle between Ar 
!         and H-O-H plane
! PN changed angles to rad
!** Output:   V [cm-1]  is the calculated interaction energy '\Delta{V}'.
!-----------------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MXDATA,MXPARM,NPMAX,MMAX,MXDE,MXRE,MXPHI,isH,isO
      PARAMETER (MXDATA=5000, MXPARM=600, NPMAX=20, MMAX=16)
      PARAMETER (isp=4, MXDE=16,MXRE=16,MXPHI=38)
! isp--the number of CN,CS
      PARAMETER (isH=1,isO=16)
      REAL*8  YPHI,YTH,YR,XPHI,TH,RR,PHI(MXPARM)
      REAL*8  DM(MMAX),DMP(MMAX),DMPP(MMAX),ylm,RTP(MXDATA)
      REAL*8  YC,Re,De,phiINF,RTPp,RTPq,yp,yq,ype,yPOW,SUM,DSUM,XP,VLR,XPW
!-----------------------------------------------------------------------
      INTEGER NPHI(0:NPMAX),I,J,K,L,M,IP,IS,NP,NQ,NLDE,
     1 NLRE,NCN,NCS,MCM,NS,NL,NPOW,NPS,NTP0,MMN,M6,iv
      REAL*8 DEDE(0:ISP+1,MXDE),RERE(0:ISP+1,MXRE),PHIPHI(0:ISP+1,MXPHI)
     1  ,XCN(ISP),XRCN2CN(ISP),XRCN22CN(ISP),XCS(ISP),XRCS1CN(ISP),
     2  XRCS30CN(ISP),XRCS32CN(ISP),XRCMCN(ISP),
     2  XRCM2CN(ISP),XRCM40CN(ISP),XRCM42CN(ISP),Xua2ba(ISP),PI,
     #  XRCM22CN(ISP),XRCM44CN(ISP)
      REAL*8  PV(MXPARM),PD(MXPARM),PS(MXPARM),RMSR,Vasy,
     1  RREF,AREF,AREFp,AREFq,Rep,CN,RCN2CN,RCN22CN,RCS1CN,RCS30CN,
     2 RCMCN,RCS32CN,RCM22CN,RCM40CN,RCM42CN,
     2 RCM2CN,RCM44CN,VLRe,dVLRedRe,XRE,C6Sum,C7Sum,C8Sum,bDAMP,
     3 T0,V


!-----------------------------------------------------------------------

!c  input parameters which are the same for all isotopologues
      DATA NLDE/6/
      DATA NLRE/6/
      DATA NP/3/
      DATA NQ/3/
      DATA NS/4/
      DATA NL/3/
      DATA RREF/1.6d0/
      NPOW=NL
      DATA (NPHI(I),I=0,3)/6,4,4,2/
      DATA bDAMP/-3.08d0/
      DATA NCN/6/
      DATA NCS/7/
      DATA MCM/8/
!-----------------------------------------------------------------------
!  input initial parameters which are changed with different isotope case
      DATA (XCN(I),I=1,ISP)/0.288294408D+06,0.287327242D+06
     & ,0.301311296D+06,0.281305833D+06/
      DATA (XRCN2CN(I),I=1,ISP)/1.06557738,1.10187951,1.00693205
     & ,1.14323972/
      DATA (XRCN22CN(I),I=1,ISP)/0.649336216,0.686325455,
     & 0.617965547, 0.699238968/
      DATA (XRCS1CN(I),I=1,ISP)/-0.162854355D+01,-0.186338057D+01
     & ,-0.159303686D+01, -0.195694485D+01/
      DATA (XRCS30CN(I),I=1,ISP)/0.232054449,0.281908529,0.209086087
     & ,0.238174941/
      DATA (XRCS32CN(I),I=1,ISP)/-0.770046153,-0.863127347,-0.730141893
     & ,-0.856718817/
      DATA (XRCMCN(I),I=1,ISP)/0.199608450D+03,0.208107660D+03
     & ,0.189538198D+03,0.213898871D+03/
      DATA (XRCM2CN(I),I=1,ISP)/-0.616036922D+02,-0.642473017D+02,
     & -0.581790336D+02,-0.658002709D+02/
      DATA (XRCM22CN(I),I=1,ISP)/-0.376004518D+02,-0.384230884D+02,
     & -0.355114466D+02,-0.397432214D+02/
      DATA (XRCM40CN(I),I=1,ISP)/0.0,0.0,0.0,0.0/
      DATA (XRCM42CN(I),I=1,ISP)/0.0,0.0,0.0,0.0/
      DATA (XRCM44CN(I),I=1,ISP)/0.0,0.0,0.0,0.0/

!     isC=12, isF=16, isH=32, (v1,v2,v3)=(0,0,0)   opt-re-fitted
      DATA (DEDE(1,I),I=1,MXDE)/
     & 4.00110D+02, -4.05000D+00, 5.29000D+00,
     & 4.67700D+01, -1.54000D+00, -2.98000D+00,
     &-7.50000D-01, 1.47000D+00, 6.71000D+00,
     & 4.20000D-01, -4.60000D-01, 2.60000D-01,
     & 6.30000D-01, -8.50000D-01, 7.60000D-01,
     & 3.90000D-01/



      DATA (RERE(1,I),I=1,MXRE)/
     &  1.28607D+01, -2.07300D-01, 5.40000D-03,
     & -1.01400D-01, 5.99000D-02, -1.04500D-01,
     & -5.80000D-03, 9.80000D-03, 0.00000D+00,
     & -1.20000D-03, 1.30000D-03, -3.30000D-03,
     & -3.20000D-03, 4.50000D-03, -4.30000D-03,
     &  1.10000D-03/

    
      DATA (PHIPHI(1,I),I=1,MXPHI)/
     & -3.1340D+00, -3.2000D-01, 1.3200D-01,
     &  3.4900D-01, 9.0000D-02, -1.6600D-01,
     & -1.1000D-02, 2.1000D-02, 3.5000D-02,
     &  5.0000D-03, -3.0000D-03, 0.0000D+00,
     &  0.0000D+00, 4.0000D-03, -4.0000D-03,
     &  0.0000D+00, 3.9800D+00, 6.0000D-02,
     &  0.0000D+00, 0.0000D+00, 0.0000D+00,
     &  9.0000D-02, 0.0000D+00, -8.0000D-02,
     &  0.0000D+00, 2.2000D+00, 2.5000D-01,
     &  7.0000D-02, 0.0000D+00, 0.0000D+00,
     &  1.0000D-01, 0.0000D+00, -1.2800D-01,
     &  0.0000D+00, 3.5100D+00, 2.0000D-01,
     &  0.0000D+00, 1.3000D-01/

      
!     isO=16, isH=1, (v1,v2,v3)=(1,0,0)   opt-re-fitted
      DATA (DEDE(2,I),I=1,MXDE)/
     &  4.0691D+02, -1.3000D+00, 3.9900D+00,
     &  4.7610D+01, -1.3800D+00, -2.1000D+00,
     & -1.0800D+00, 1.7900D+00, 6.5700D+00,
     &  5.2000D-01, -4.6000D-01, 2.1000D-01,
     &  6.8000D-01, -1.0100D+00, 9.0000D-01,
     &  3.8000D-01/




      DATA (RERE(2,I),I=1,MXRE)/
     &  1.28819D+01, -2.41600D-01, 8.60000D-03,
     & -8.39000D-02, 6.40000D-02, -1.15700D-01,
     & -6.70000D-03, 1.05000D-02, 1.80000D-03,
     & -9.00000D-04, 1.10000D-03, -2.20000D-03,
     & -3.90000D-03, 5.50000D-03, -5.30000D-03,
     &  1.00000D-03/


      DATA (PHIPHI(2,I),I=1,MXPHI)/
     & -3.1600D+00, -3.3700D-01, 1.3800D-01,
     &  3.5200D-01, 9.6000D-02, -1.9300D-01,
     & -1.7000D-02, 3.0000D-02, 4.5000D-02,
     &  6.0000D-03, -3.0000D-03, 0.0000D+00,
     & -3.0000D-03, 4.0000D-03, -6.0000D-03,
     &  0.0000D+00, 4.0000D+00, 1.0000D-01,
     &  0.0000D+00, 0.0000D+00, 0.0000D+00,
     &  0.0000D+00, 0.0000D+00, -9.0000D-02,
     &  0.0000D+00, 2.2000D+00, 3.6000D-01,
     &  0.0000D+00, 2.0000D-01, 0.0000D+00,
     &  0.0000D+00, 2.0000D-02, -1.6100D-01,
     & -4.0000D-02, 3.5000D+00, 3.0000D-01,
     & -8.0000D-02, 4.0000D-01/



!     isO=16, isH=1, (v1,v2,v3)=(0,1,0)   opt-re-fitted
      DATA (DEDE(3,I),I=1,MXDE)/
     &  4.0256D+02, -4.4200D+00, 4.8800D+00,
     &  4.7880D+01, -2.2100D+00, -2.4400D+00,
     & -2.4000D-01, 9.3000D-01, 7.2100D+00,
     &  3.8000D-01, -4.5000D-01, 2.6000D-01,
     &  5.9000D-01, -7.9000D-01, 6.9000D-01,
     &  4.6000D-01/



      DATA (RERE(3,I),I=1,MXRE)/
     &  1.28615D+01, -2.07800D-01, 6.40000D-03,
     & -1.01300D-01, 5.97000D-02, -1.04500D-01,
     & -6.20000D-03, 1.03000D-02, -5.00000D-04,
     & -1.00000D-03, 9.00000D-04, -2.80000D-03,
     & -3.10000D-03, 4.50000D-03, -4.20000D-03,
     &  1.10000D-03/



      DATA (PHIPHI(3,I),I=1,MXPHI)/
     & -3.1140D+00, -3.2100D-01, 1.3400D-01,
     &  3.4400D-01, 8.5000D-02, -1.6100D-01,
     & -9.0000D-03, 2.1000D-02, 4.8000D-02,
     &  5.0000D-03, -4.0000D-03, 0.0000D+00,
     &  0.0000D+00, 4.0000D-03, -5.0000D-03,
     &  0.0000D+00, 3.9900D+00, 6.0000D-02,
     &  0.0000D+00, 0.0000D+00, 0.0000D+00,
     &  9.0000D-02, 0.0000D+00, -7.0000D-02,
     &  2.0000D-02, 2.2000D+00, 2.5000D-01,
     &  6.0000D-02, 1.0000D-01, 0.0000D+00,
     &  1.0000D-01, 0.0000D+00, -1.1700D-01,
     &  0.0000D+00, 3.4900D+00, 2.0000D-01,
     &  0.0000D+00, 2.5000D-01/



!     isO=16, isH=1, (v1,v2,v3)=(0,0,1)   opt-re-fitted
      DATA (DEDE(4,I),I=1,MXDE)/
     &  4.07980D+02, -2.64000D+00, 4.43000D+00,
     &  4.82200D+01, -9.90000D-01, -2.80000D+00,
     & -1.25000D+00, 2.06000D+00, 6.71000D+00,
     &  6.40000D-01, -6.20000D-01, 2.20000D-01,
     &  6.00000D-01, -8.70000D-01, 8.80000D-01,
     &  3.90000D-01/




      DATA (RERE(4,I),I=1,MXRE)/
     &  1.28777D+01, -2.42700D-01, 1.80000D-02,
     & -9.34000D-02, 5.63000D-02, -1.09700D-01,
     & -5.30000D-03, 9.80000D-03, 1.00000D-03,
     & -9.00000D-04, 1.20000D-03, -2.40000D-03,
     & -3.50000D-03, 4.90000D-03, -5.00000D-03,
     &  1.10000D-03/


      DATA (PHIPHI(4,I),I=1,MXPHI)/
     & -3.1790D+00, -3.3900D-01, 1.4700D-01,
     &  3.4300D-01, 7.1000D-02, -1.8700D-01,
     & -1.1000D-02, 2.8000D-02, 4.3000D-02,
     &  7.0000D-03, -5.0000D-03, 3.0000D-03,
     & -2.0000D-03, 5.0000D-03, -6.0000D-03,
     &  0.0000D+00, 3.9200D+00, 1.6000D-01,
     &  0.0000D+00, 0.0000D+00, -9.0000D-02,
     &  0.0000D+00, 0.0000D+00, -8.0000D-02,
     &  0.0000D+00, 2.1000D+00, 4.4000D-01,
     &  0.0000D+00, 2.0000D-01, -1.0000D-01,
     &  0.0000D+00, 0.0000D+00, -1.3900D-01,
     & -3.0000D-02, 3.5000D+00, 3.0000D-01,
     & -1.0000D-01, 3.8000D-01/



         call fct(40)


c-----------------------------------------------------------------------

      IF((isH.EQ.1).and.(isO.EQ.16).and.(iv.EQ.1)) isv= 1
      IF((isH.EQ.1).and.(isO.EQ.16).and.(iv.EQ.2))
     & isv= 2
      IF((isH.EQ.1).and.(isO.EQ.16).and.(iv.EQ.3))
     & isv= 3
      IF((isH.EQ.1).and.(isO.EQ.16).and.(iv.EQ.4))
     & isv= 4

      PI=DACOS(-1.0D0)
!... input initial parameters which are changed with different isotope case 
      CN=XCN(isv)
      RCN2CN=XRCN2CN(isv)
      RCN22CN=XRCN22CN(isv)
      RCS1CN=XRCS1CN(isv)
      RCS30CN=XRCS30CN(isv)
      RCS32CN=XRCS32CN(isv)
      RCMCN=XRCMCN(isv)
      RCM2CN=XRCM2CN(isv)
      RCM22CN=XRCM22CN(isv)
      RCM40CN=XRCM40CN(isv)
      RCM42CN=XRCM42CN(isv)
      RCM44CN=XRCM44CN(isv)
      IS=0
      IPDE=0
      DO L=0,NLDE
!c... first prepare well depth expansion parameters
          if (L.GT.6) then
              M6=6
          else
             M6=L
          endif
           DO M=0,M6,2
             IS=IS+1
             IPDE=IPDE+1
             PV(IS)=DEDE(isv,IPDE)
!c         write(13,123)DEDE,PV(IS)
123       format(2x,2f25.8)
           ENDDO
      ENDDO
      IPRE=0
      DO L=0,NLRE
!c... next prepare potential minimum position expansion parameters
          if (L.GT.6) then
              M6=6
          else
             M6=L
          endif
            DO M=0,M6,2
              IS=IS+1
              IPRE=IPRE+1
              PV(IS)=RERE(isv,IPRE)
            ENDDO
      ENDDO
        IPPHI=0
        NPS=0
        K=0
      DO J=0,NPOW
!c... then, prepare exponent coefficient \beta_i expansion parameters
         DO L=0,NPHI(J)
            if (L.GT.6) then
              M6=6
            else
             M6=L
            endif
              DO M=0,M6,2
                IS=IS+1
                IPPHI=IPPHI+1
                PV(IS)=PHIPHI(isv,IPPHI)
!c         write(14,124)IPPHI,PHIPHI,PV(IS)
124       format(2x,3f25.8)
              ENDDO
           ENDDO
       ENDDO
!c** Finally ... (if desired) read shift for definition of energy zero
      IS=IS+1
      PV(IS)=0.0D0
!cccc----------------------------------------------------------------ccccc
           PI=DACOS(-1.0D0)
           YR=RR
c           YPHI=XPHI*PI/180.D0
c           YTH=TH*PI/180.D0
           YPHI=XPHI
           YTH=TH

!c=======================================================================
!c  For the case of an  MLR_{p}  potential ...
!c-----------------------------------------------------------------------
!c caculate the derivative of the parameters of De not including
!!c the coefficient before, only the Legendre expansion.
        IP=0
        De=0.0d0
        DO L=0,NLDE        ! change with NDE  
          if (L.GT.6) then
              M6=6
          else
             M6=L
          endif
           DO M=0,M6,2
             IP=IP+1
             PD(IP)=ylm(L,M,YTH,YPHI)
             De=De+PV(IP)*PD(IP)
!c            WRITE(60,*)De,PV(IP),PD(IP)
            ENDDO
        ENDDO
!c caculate the derivative of the parameters of Re not including
!c the coefficient before, only the Associated Legendre expansion.

       Re=0.0d0
        DO L=0,NLRE        ! change with  NRE 
          if (L.GT.6) then
              M6=6
          else
             M6=L
          endif
           DO M=0,M6,2
             IP=IP+1
             PD(IP)=ylm(L,M,YTH,YPHI)
             Re=Re+PV(IP)*PD(IP)
            ENDDO
        ENDDO

       AREF= RREF*Re
       IF(RREF.LE.0.d0) AREF= Re
       AREFp= AREF**NP
       AREFq= AREF**NQ
       Rep= Re**NP
       C6Sum=0.0d0
       C6Sum=1.+ylm(2,0,YTH,YPHI)*RCN2CN+ylm(2,2,YTH,YPHI)*RCN22CN
       IF(bDAMP.GT.0.d0) CALL DAMPIG(Re,bDAMP,MMAX,DM,DMP,DMPP)
         T0= CN*C6Sum/Re**NCN
         IF(bDAMP.GT.0.d0) T0= T0*DM(NCN)
         VLRe= T0
         dVLRedRe= -NCN*T0/Re
         IF(bDAMP.GT.0.d0)
     &    dVLRedRe= dVLRedRe + T0*DMP(NCN)/DM(NCN)
       C7Sum=0.0d0
       C7Sum=ylm(1,0,YTH,YPHI)*RCS1CN+ylm(3,0,YTH,YPHI)*RCS30CN+
     &    ylm(3,2,YTH,YPHI)*RCS32CN
       C8Sum=0.0d0
       C8Sum=RCMCN+ylm(2,0,YTH,YPHI)*RCM2CN+ylm(2,2,YTH,YPHI)*RCM22CN+
     & ylm(4,0,YTH,YPHI)*RCM40CN+ylm(4,2,YTH,YPHI)*RCM42CN
     & +ylm(4,4,YTH,YPHI)*RCM44CN


       IF(MCM.GT.NCN) THEN
        MMN= MCM - NCN
       IF(NP.LE.MMN) MMN= 0
       IF(MMN.GT.0) THEN
!C  with and without damping function 
       IF(bDAMP.GT.0.d0) THEN
         VLRe= (CN/Re**NCN)*(C6Sum*DM(NCN)+C7Sum*DM(NCS)/Re
     &   +C8Sum*DM(MCM)/Re**MMN)
         dVLRedRe= dVLRedRe - NCS*CN*C7Sum*DM(NCS)/Re**(NCS+1)
     &   +CN*C7Sum*DMP(NCS)/Re**MCM - MCM*CN*C8Sum*DM(MCM)/Re**(MCM+1)
     &   +CN*C8Sum*DMP(MCM)/Re**MCM
         ELSE
         VLRe= (CN/Re**NCN)*(C6Sum+C7Sum/Re+C8Sum/Re**MMN)
         dVLRedRe= dVLRedRe -NCS*CN*C7Sum/Re**(NCS+1)
     &    -MCM*CN*C8Sum/Re**(MCM+1)
       ENDIF
       ENDIF
         phiINF= DLOG(2.d0*De/VLRe)
       ENDIF
       RTPp= YR**NP
       RTPq= YR**NQ
       yp= (RTPp - AREFp)/(RTPp + AREFp)
       yq= (RTPq - AREFq)/(RTPq + AREFq)
       ype= (RTPp - Rep)/(RTPp + Rep)
!c caculate the derivative of the parameters of PHI(N) not including
!c the coefficient before, the Legendre expansion and exponent expansion.

      NPOW= NS
       IF(YR.GE.Re) NPOW= NL
        yPOW= 1.d0 - yp
        SUM=0.0
        DSUM=0.0
        NPS=0
        DO J=0,NPOW
        PHI(J)=0.0d0
        DO L=0,NPHI(J)
          if (L.GT.6) then
              M6=6
          else
             M6=L
          endif
           DO M=0,M6,2
             IP=IP+1
             PHI(J)=PHI(J)+PV(IP)*ylm(L,M,YTH,YPHI)
             PD(IP)=ylm(L,M,YTH,YPHI)*yq**(J)
            ENDDO
        ENDDO
        SUM=SUM+PHI(J)*yq**(J)
        IF(RREF.LE.0.D0) DSUM=DSUM+yPOW*PHI(J)*(J)*yq**(J-1)
        ENDDO

        PD(IP+1)=1.0D0
        Vasy=PD(IP+1)*PV(IP+1)
        XP= SUM*yPOW+ phiINF*yp

        IF(bDAMP.GT.0.d0) CALL DAMPIG(YR,bDAMP,MMAX,DM,DMP,DMPP)
        IF(bDAMP.GT.0.d0)THEN
          VLR= CN*C6SUM*DM(NCN)/YR**NCN
          ELSE
          VLR= CN*C6Sum/YR**NCN
          ENDIF
        IF(MMN.GT.0) THEN
          IF(bDAMP.GT.0.d0)THEN
             VLR= (CN/YR**NCN)*(C6Sum*DM(NCN)+C7Sum*DM(NCS)/YR
     &       +C8Sum*DM(MCM)/YR**MMN)
             ELSE
             VLR= (CN/YR**NCN)*(C6Sum+C7Sum/YR+C8Sum/YR**MMN)
             ENDIF
       ENDIF
         XPW= DEXP(-XP*ype) * VLR/VLRe
         YC= De*(1.d0 - XPW)**2-De+Vasy 
99      format(1x,3F8.2,F15.5)
         V=YC
!         write(11,98) YTH, YPHI,YR, V
! commented by PN above
98      format(1x,3F8.2,F15.5)
         RETURN
         END

!C -------------------------------------------------------------------------
!c
!c compute the matrix of N!
!c
        subroutine fct(nmax)
        implicit real*8 (a-h,o-z)
        common/factorial/ fact(0:40)
!        common/factorial/ f(0:40)
!c
!        f(0) = 1.d0
        fact(0) = 1.d0
        do i=1,nmax
         fact(i) = fact(i-1)*i
!         f(i) = f(i-1)*i
        end do
        return
        end
!cccc----------------------------------------------------------------ccccc
        function ylm(L,M,th1,phi)
        implicit real*8 (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
        data pifact/12.566370614359d0/

        costh1=dcos(th1)
        call plmrb(p,costh1,L)

!c        Normal coeffcient builed in p(L,M)
        if(M.eq.0)then
        ylm=p(L,M)
        else
         if(M.gt.0)then
          ylm=sqrt(2.0)*p(L,M)*dcos(M*phi)
        else
          MM=-M
          ylm=sqrt(2.0)*(-1)**(MM)*p(L,MM)*dsin(-MM*phi)
        endif
        endif
        return
        end
!C ----------------------------------------------------------------------------
!c
!c Compute the set of associated Legendre polynomials P_lm
! for l=0,1,...,lmax, and m=0,1,...,l. First the standard
!c polynomials

!c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
!c
!c are computed, and then multiplied by
!c
!c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
!c
!c to get the P_lm polynomials....
!c
        subroutine plmrb(p,x,lmax)
        implicit real*8 (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
!c inverse of dsqrt(2Pi)
        data twopinv /0.3989422804014d0/
!c
!c starting value
!c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
!c
!c compute the diagonal elements
!c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
!c
!c compute P_lm along the columns with fixed m
!c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do
        end do
!c
!c Renormalize values...
!c
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do
!c
        return
        end





!c***********************************************************************
      SUBROUTINE DAMPIG(r,b,MMAX,DM,DMP,DMPP)
!c** Subroutine to generate values DM and the first and second radial
!c  derivatives DMP and DMPP of the Tang-Toennies/incomplete gamma 
!c  function damping function of all orders up to NMAX at radial distance
!c   'r' for damping parameter 'b'
!c***********************************************************************
      INTEGER MMAX,m,n
      REAL*8 b,r,br,XP,SSm,SSm1,SSm2,mfact,TK,DM(MMAX),DMP(MMAX),
     1       DMPP(MMAX),SM(-1:MMAX)
     
      br= b*r
      XP= DEXP(-br)
      SSm= 0.d0
      SSm1= 0.d0
      SSm2= 0.d0
      TK= 1.d0
      DO m=1,MMAX
         SSm2= SSm1
         SSm1= SSm
         SSm= SSm + TK
         DM(m)= 1.d0 - XP*SSm
         DMP(m)= b*XP*(SSm - SSm1)
         DMPP(m)= b**2 *XP* (2.d0*SSm1 - SSm - SSm2)
         TK= TK*br/DFLOAT(m)
      ENDDO
      IF(DM(MMAX).LT.1.0D-13) THEN
         mfact= 1.d0
         DO n=1,MMAX
            mfact= mfact*DFLOAT(n)
         ENDDO
         SSm= 0.d0
         TK= (br)**MMAX/mfact
         DO n=MMAX,MMAX+3
            SSm= SSm+TK
            TK= TK*br/DFLOAT(n+1)
         ENDDO
         SM(MMAX)= SSm
         TK= (br)**MMAX/mfact
         DO n=1,MMAX-1
            TK= TK*DFLOAT(MMAX+1-n)/br
            SM(MMAX-n)= TK+SM(MMAX-n+1)
         ENDDO
         SM(0)=1.d0 + SM(1)
         SM(-1)=SM(0)
         DO m=1,MMAX
            DM(m)= XP*SM(m)
            DMP(m)= -b*XP*(SM(m)-SM(m-1))
            DMPP(m)= b**2 *XP*(SM(m)-2.d0*SM(m-1)+SM(m-2))
         ENDDO
      ENDIF
      RETURN
      END
!c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12



