
! gfortran -g -fno-align-commons -o run_raytrace_old run_raytrace_old.for src/amsc.for src/asorp.for

! QNS is refrac
! HAI is final upper height
! HLA is lower height
! CFK is unit conversion to km
! CFM is unit conversion to m
! These variables are described in cards.txt!

! in raytracer:
! QHC is initial height above ground plane
! QHA is final height above ground plane

! in ata
! HRP is elevation of ground plane (converted from HPFI input)
! H1, H2 are antenna elevations

      program run_raytrace

        implicit real*8(A-H,O-Z)

        ! QNS surface refrac
        !
        !
        !
        !
        COMMON/RYTC/QNS,QHC,QHA,QHS,QQD

        DUM = 0.0
        HRP = 0.0 ! ground plane elev (km presumably)
        H1 = 4.8768 / 1000 ! km
        H2 = 0.3048 * 30000 / 1000 ! km
      HP2=H2-HRP
      HP1=H1-HRP
      ENS = 301 ! Maybe! They mat translate from ground to platform?
        ! Initialize stuff...
      write(6,*) ''
        QLIM = -1.56
      QNS=329.0
      QHC=H1
      QHA=HP2
      QHS=HRP
        DUM = RYTRA(DUM)
        RY = TRCRY(QLIM)
      DS0=QQD
      write(6,*) 'ry = ', ry
      write(6,*) 'qqd = ', qqd

      write(6,*) ''
        ZER = 0.0
      QNS=301.0
      QHC=0.0
      QHA=HP2
      QHS=HRP
        DUM = RYTRA(DUM)
        RY2 = TRCRY(ZER)
      DLSR=QQD
      write(6,*) 'ry2 = ', ry2
      write(6,*) 'qqd = ', qqd

      write(6,*) ''
      QNS=ENS
      QHC=ZER
      QHA=HP1
      QHS=HRP
      RY1=TRCRY(ZER)
      DLST=QQD
      write(6,*) 'ry1 = ', ry1
      write(6,*) 'qqd = ', qqd

      ! Made up....
      write(6,*) ''
      QNS=ENS
      QHC=HP1
      QHA=HP2
      QHS=HRP
      RAY=TRCRY(0.05d0)
      DLST=QQD
      write(6,*) 'ray = ', ray
      write(6,*) 'qqd = ', qqd

      ! Made up point downward....
      write(6,*) ''
      QNS=ENS
      QHC=HP1
      QHA=HP2
      QHS=HRP
      RAY=TRCRY(-0.05d0)
      DLST=QQD
      write(6,*) 'ray = ', ray
      write(6,*) 'qqd = ', qqd

      end program
