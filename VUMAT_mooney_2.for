       subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
		
		real*8 xJ, xP, U(3,3), U2dev(3,3), U4dev(3,3), C1, C2, D1 ,traceU2dev, traceU4dev , constant1, constant2
		real*8 sigma(3,3),xIdentity(3,3),lamda1_2,lamda2_2,lamda3_2,xB(3,3)
		
		C1 = props(1)
		C2 = props(2)
		D1 = props(3)
		
		sigma=0.
		
		xIdentity=0.0
		xIdentity(1,1)=1.0
		xIdentity(2,2)=1.0
		xIdentity(3,3)=1.0
		
		U2dev=0.0
		U4dev=0.0
		xB=0.0

		! loop through all blocks
		do km = 1, nblock

		U(1,1) = stretchNew(km,1)
		U(2,2) = stretchNew(km,2)
		U(3,3) = stretchNew(km,3)
		U(1,2) = stretchNew(km,4)
		U(2,1) = U(1,2)
		
		if (nshr .eq. 1) then
			U(2,3) = 0.0
			U(1,3) = 0.0
		else
			U(2,3) = stretchNew(km,5)
			U(1,3) = stretchNew(km,6)
		end if
		
		U(3,2)=U(2,3)
		U(3,1)=U(1,3)
		
		! calculate J
		xJ = (U(1,1)*(U(2,2)*U(3,3)-U(2,3)**2)) + (U(1,2)*(U(2,3)*U(1,3)-U(1,2)*U(3,3))) + (U(1,3)*(U(1,2)*U(2,3)-U(2,2)*U(1,3)))

		! Define the square of the deviatoric stretch tensor (B-bar)
		xB=matmul(U,U)
		U2dev=(xJ**(-2.0/3.0))*matmul(U,U)

		! Define the forth power of the deviatoric stretch tensor (B-bar * B-bar)
		U4dev=matmul(U2dev,U2dev)
				
		! Define some constants to simplify the stress expression
		
		lamda1_2=xB(1,1)
		lamda2_2=xB(2,2)
		lamda3_2=xB(3,3)
		
		I1 = (xJ**(-2.0/3.0))*(lamda1_2+lamda2_2+lamda3_2)  ! I-bar-1
		
		
		I2 =(xJ**(-4.0/3.0))*( lamda1_2*lamda2_2 + lamda2_2*lamda3_2 + lamda3_2*lamda1_2 )  ! I-bar-2
		
		
		constant1 = 2.0*(xJ-1.0)/D1
		constant2 = 2.0/xJ
		
		!new stress matrix
		
		
		sigma = constant1*xIdentity + constant2*((C1+I1*C2)*U2dev -C2*U4dev - (1.0/3.0)*(C1*I1 +2.0*C2*I2)*xIdentity)

		
		!The corotational stress
		stressNew(km,1) = sigma(1,1)
		stressNew(km,2) = sigma(2,2)
		stressNew(km,3) = sigma(3,3)
		stressNew(km,4) = sigma(1,2)
		
		if (nshr .eq. 3) then
			stressNew(km,5) = sigma(2,3)
			stressNew(km,6) = sigma(3,1)
		end if
		
		end do
		
      return
      end