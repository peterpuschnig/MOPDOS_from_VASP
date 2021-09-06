!!$******************************* ARUPS ************************************
!!$*
!!$*   input the WAVECAR file in binary format and the IBZKPT file from VASP
!!$*   output calculated ARUPS intensities for 3 different modes specified in
!!$*   the input file. 
!!$*
!!$*   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$*   for ifort.
!!$*
!!$*   with the help of the program Wavetrans version 2.0 - July 5, 2012 - by R. M. Feenstra and M. Widom
!!$*
!!$*
!!$*
!!$*   format of the input file which has to be named ARUPS.in
!!$*
!!$*   "string"                      ! path and name of output file 
!!$*   integerN                      ! number of k-points
!!$*   integerN integerN             ! first and last bandindex for summation over initial states
!!$*   realN                         ! Fermi-level (energy of HOMO) in [eV]
!!$*   realN                         ! Temperature in [eV] for Fermi-Dirac smearing, negative value means no FD smearing
!!$*   realN                         ! work function in [eV]
!!$*   realN realN                   ! energy range of incident photon [eV], only first number used for Energy vs k and kx-ky plot
!!$*   realN                         ! energy broadening [eV]
!!$*   realN realN                   ! broadening in k-space [1/A]
!!$*   integerN                      ! 0 = Energy versus k, 1 = kx-ky plot , 2 = CIS scan
!!$*   realN realN realN             ! binding energy range and step size in [eV], kx-ky plot and CIS scan use first binding energy only
!!$*   realN realN realN             ! range of k-vectors and step size in [1/A], Energy vs k and CIS scan use first vector only
!!$*   realN realN realN             ! unit vector for normal emission
!!$*   realN realN realN             ! unit vector for grazing emission
!!$*   realN realN                   ! damping parameter and z-position of surface, negative damping parameter means no damping
!!$*
!!$*
!!$*   note that the energy eigenvalues are complex, as provided in the
!!$*   WAVECAR file, but the imaginary part is zero (at least for all cases
!!$*   investigated thus far)
!!$*     
                                                                        
implicit real*8 (a-h, o-z)                                        
complex*8, allocatable :: coeff(:,:),coeffmol(:,:)
complex*16, allocatable :: cener(:),cenermol(:)
real*8, allocatable :: occ(:),igall(:,:),occmol(:),weights(:)
integer, allocatable :: hkl(:,:)
dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),wk(3),a1mol(3),a2mol(3),a3mol(3)


character (1000) :: fullfile,molfile
integer              :: kcount,nparwavecar,khi,numk,numE
integer              :: bandcount,bandlow,bandhi,bandcountmol,bandlowmol,bandhimol
integer              :: nkx,nek,ikx,iky,iek,iekstart,iekend
integer              :: iwindow,ephotonstart,ephotonend,nephk
integer              :: jg , hmil , kmil,  lmil , line
integer              :: jgz, hmilz, kmilz, lmilz      
integer              :: arupstype,isymmetry,nsymmetry
integer              :: ikxstart,ikxend,ikystart,ikyend
real*8               :: two_pi,ekintok,epsilon
real*8               :: efermi,temperature,fermifac,phiwork
real*8               :: ephotonmin,ephotonmax,ebroad,kxbroad,kybroad
real*8               :: vecz(3),vecx(3),vecy(3)  
real*8               :: kxlo,kxhi,kxstep
real*8               :: eklo,ekhi,ekstep,ekinmin,ekinmax
real*8  ,allocatable :: eigval(:,:,:)  
real*8               :: qvec(3),qweight,qplusG(3),qplusGx,qplusGy,qplusGz,qplusGlen
real*8               :: kfvec(3)
real*8               :: bvec(3,3)
real*8               :: deltae,deltak,ebind,efinal,kfinal
real*8               :: gz,gamma,z0,ccell
complex*16           :: fgz,cgccomplex(3),imag
complex*16 , allocatable:: matrixelement(:,:,:,:)
real*8               :: l1,l2,kxarr,kyarr,kpar,ekarr,pdos,gauss,ie   
logical              :: pureplanewave


     
!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)


data c/0.262465831d0/

!data c/0.26246582250210965422d0/ 
pi=4.*atan(1.)
      
!!$*   input


OPEN(UNIT=85,FILE='IBZKPT',STATUS='old',action='read', iostat=iost) 
if (iost.ne.0) write(80,*) 'IBZKPT open error - iostat =',iost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! read weights from IBZKPT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

read(85,*) 
read(85,*) numk
read(85,*)

allocate (weights(numk))

do line = 1,numk
   read(85,*) dummy, dummy, dummy, weights(line)
   weights(line) = weights(line)/2
end do

numk = sum(weights)
close(unit=85)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(UNIT=81,FILE='MOPDOS.in', STATUS='UNKNOWN')

      
two_pi    =  2.0d0*3.14159265d0

read(81,*) fullfile              ! name of output file
read(81,*) molfile               ! name of output file
read(81,*) khi                   ! number of k-points
read(81,*) bandlowmol,bandhimol     ! first and last bandindex of molecules
read(81,*) bandlow,bandhi     ! first and last bandindex of full system
read(81,*) efermi             ! Fermi-level (energy of HOMO) in [eV]
read(81,*) gaussfac        ! Temperature in [eV]
read(81,*) eklo,ekhi,ekstep   ! kinetic energy range (1 - 2) and step size (3) in [eV] 

OPEN(UNIT=80,FILE='MOPDOS.log',STATUS='UNKNOWN')
OPEN(UNIT=82,FILE='MOPDOS.out',STATUS='UNKNOWN')

WRITE(80,*) 'CALCULATION OF MOPDOS SPECTRA' 



nrecl=24
OPEN(UNIT=10,FILE=TRIM(ADJUSTL(fullfile)),access='direct',recl=nrecl, &
     iostat=iost,STATUS='old')
if (iost.ne.0) write(80,*) 'WAVECAR of full calculation open error - iostat =',iost





OPEN(UNIT=11,FILE=TRIM(ADJUSTL(molfile)),access='direct',recl=nrecl, &
     iostat=iost,STATUS='old')
if (iost.ne.0) write(80,*) '1 WAVECAR of mol calculation open error - iostat =',iost


read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)


read(unit=11,rec=1) xnreclmol,xnspinmol,xnprecmol
close(unit=11)


nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)


nreclmol=nint(xnreclmol)
nspinmol=nint(xnspinmol)
nprecmol=nint(xnprecmol)



if(nprec.eq.45210) then
   write(80,*) '*** error - full WAVECAR_double requires complex*16'
   stop
endif

if(nprec2.eq.45210) then
   write(80,*) '*** error - mol WAVECAR_double requires complex*16'
   stop
endif



write(80,*) 
write(80,*) ' full: record length  =',nrecl,' spins =',nspin, &
     ' prec flag ',nprec

write(80,*) 
write(80,*) ' mol: record length  =',nreclmol,' spins =',nspinmol, &
     ' prec flag ',nprecmol


OPEN(UNIT=10,FILE=TRIM(ADJUSTL(fullfile)),access='direct',recl=nrecl, &
     iostat=iost,STATUS='old')
if (iost.ne.0) write(80,*) 'WAVECAR of full calculation open error - iostat =',iost

OPEN(UNIT=11,FILE=TRIM(ADJUSTL(molfile)),access='direct',recl=nreclmol, &
     iostat=iost,STATUS='old')
if (iost.ne.0) write(80,*) '2 WAVECAR of mol calculation open error - iostat =',iost



read(unit=11,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)


nwkmol=nint(xnwk)
nbandmol=nint(xnband)
ecutmol = ecut
allocate(occmol(nbandmol))
allocate(cenermol(nbandmol))

a1mol = a1
a2mol = a2
a3mol = a3


write(80,*) 'no. k points =',nwkmol
write(80,*) 'no. bands =',nbandmol
write(80,*) 'max. energy =',sngl(ecut)
write(80,*) 'real space lattice vectors:'
write(80,*) 'a1 =',(sngl(a1(j)),j=1,3)
write(80,*) 'a2 =',(sngl(a2(j)),j=1,3)
write(80,*) 'a3 =',(sngl(a3(j)),j=1,3)
write(80,*) ' '


if(nbandmol .lt. bandhimol .or. nbandmol .lt. bandlowmol) then
  write(6,*) 'band range of mol for analysis exceeds number of bands in WAVECAR   stopping'
  stop
endif




read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)

nwk=nint(xnwk)
nband=nint(xnband)
allocate(occ(nband))
allocate(cener(nband))


if(nband .lt. bandhi) then
  write(6,*) 'band range of full system exceeds number of bands in WAVECAR'  
  write(6,*) 'reducing maximum bandnumber to NBANDS of full system'
!  bandlow = 1
  bandhi  = nband
endif




write(80,*) 'no. k points =',nwk
write(80,*) 'no. bands =',nband
write(80,*) 'max. energy =',sngl(ecut)
write(80,*) 'real space lattice vectors:'
write(80,*) 'a1 =',(sngl(a1(j)),j=1,3)
write(80,*) 'a2 =',(sngl(a2(j)),j=1,3)
write(80,*) 'a3 =',(sngl(a3(j)),j=1,3)
write(80,*) ' '


if(nwk .ne. nwkmol) then
  write(6,*) 'different number of k-points found in WAVECAR files   stopping'
  stop
endif

if(ecut .ne. ecutmol) then
  write(6,*) 'different energy cut off found in WAVECAR files   stopping'
  stop
endif



call vcross(vtmp,a2,a3)
Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
write(80,*) 'unit cell volume =',sngl(Vcell)
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
do j=1,3
   b1(j)=2.*pi*b1(j)/Vcell
   b2(j)=2.*pi*b2(j)/Vcell
   b3(j)=2.*pi*b3(j)/Vcell
enddo
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
write(80,*) 'reciprocal lattice vectors:'
write(80,*) 'b1 =',(sngl(b1(j)),j=1,3)
write(80,*) 'b2 =',(sngl(b2(j)),j=1,3)
write(80,*) 'b3 =',(sngl(b3(j)),j=1,3)


phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
npmax=min0(npmaxA,npmaxB,npmaxC)



allocate (igall(3,npmax))
allocate (hkl(3,npmax))
allocate (coeff(npmax,nband))
allocate (coeffmol(npmax,nbandmol))

allocate(matrixelement(nwk,nbandmol,nband,nspin))
allocate(eigval(nwk,nband,nspin))

irec=2
irec2=2
do isp=1,nspin
  do iwk=1,nwk                                   ! loop over q points in BZ
      irec=irec+1
      read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
           (cener(iband),occ(iband),iband=1,nband)
      nplane=nint(xnplane)
      write(80,*) 'k point #',iwk,'full  input no. of plane waves =', &
           nplane
      write(80,*) 'k value =',(sngl(wk(j)),j=1,3)
      do iband=1,nband
         irec=irec+1
         read(unit=10,rec=irec) (coeff(iplane,iband), &
              iplane=1,nplane)
      enddo

      irec2=irec2+1
      read(unit=11,rec=irec2) xnplane,(wk(i),i=1,3), &
           (cenermol(iband),occmol(iband),iband=1,nbandmol)
      nplanemol=nint(xnplane)
      write(80,*) 'k point #',iwk,'mol  input no. of plane waves =', &
           nplanemol
      write(80,*) 'k value =',(sngl(wk(j)),j=1,3)
      do iband=1,nbandmol
         irec2=irec2+1
         read(unit=11,rec=irec2) (coeffmol(iplane,iband), &
              iplane=1,nplane)
      enddo
         
      if( nplanemol .ne. nplane) then  
         write(6,*) 'different number of plane waves found in WAVECAR files   stopping'
         stop
      endif

      ncnt=0
      do ig3=0,2*nb3max
         ig3p=ig3
         if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
         do ig2=0,2*nb2max
            ig2p=ig2
            if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
            do ig1=0,2*nb1max
               ig1p=ig1
               if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
               do j=1,3
                  sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                       (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
               enddo
               gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
               etot=gtot**2/c
               if (etot.lt.ecut) then
                  ncnt=ncnt+1
               end if
            enddo
         enddo
      enddo
      if (ncnt.ne.nplane) then
         write(6,*) '*** error - computed no.',ncnt , ' != input no.', nplane
         stop
      endif
      if (ncnt.gt.npmax) then
         write(6,*) '*** error - plane wave count exceeds estimate'
         stop
      endif
   

      do bandcountmol = bandlowmol, bandhimol ! loop over bands mol  
!        pdos(iwk,bandcountmol) = 0.0d0
        do bandcount = bandlow , bandhi ! loop over bands full system
        

          eigval(iwk,bandcount,isp)   = real(cener(bandcount))             
          eigval(iwk,bandcount,isp)   = eigval(iwk,bandcount,isp)   - efermi
          matrixelement(iwk,bandcountmol,bandcount,isp) = (0.0d0,0.0d0)

!          eigval   = real(cener(bandcount))            
!          eigval   = eigval - efermi
!          matrixelement = (0.0d0,0.0d0)

          do  jg = 1, nplane   ! loop over G-vectors
            matrixelement(iwk,bandcountmol,bandcount,isp) = matrixelement(iwk,bandcountmol,bandcount,isp) &
                                                           +  conjg(coeffmol(jg,bandcountmol)) * coeff(jg,bandcount) 
          end do
          write(80,*)'spin= ',isp, 'kpoint= ', iwk, 'band= ',bandcountmol,bandcount, matrixelement(iwk,bandcountmol,bandcount,isp)
        end do        
      end do
      write(80,*) 'k-point number ', iwk, ' done' 
   end do
end do

numE = (ekhi -eklo)/ekstep +1



write(80,*) 'WAVECAR  read' 
write(80,*) 
write(80,*) 'calculating MOPDOS ' 


do isp=1,nspin! loop over spin
   write(82,*) '#   spin', isp
   do bandcountmol = bandlowmol, bandhimol ! loop over bands mol
      write(82,*) '#   bandnumber', bandcountmol
      do nie = 1,nume ! loop over energies
         ie = eklo + (nie-1)*ekstep
         pdos = 0.0d0
         do iwk = 1 , nwk ! kpoints
            do bandcount = bandlow, bandhi
               gauss = (1/(gaussfac*sqrt(two_pi)))*exp(-(ie-eigval(iwk,bandcount,isp))**2.0d0/(2.0d0*gaussfac**2.0d0))
               pdos = pdos + weights(iwk)* gauss * abs(matrixelement(iwk,bandcountmol,bandcount,isp))**2.0d0
!	       write(80,*)'k= ',iwk,'spin',isp, 'bandfull=  ', bandcount, 'bandmol=  ', bandcountmol, 'pdos=  ',pdos 
            end do            
         end do ! kpoints
	 pdos = pdos *((isp-2)*(-2)-1)
         write(82,*) ie, pdos
      end do ! end loop over energies
      write(82,*)
   end do ! loop over bands mol
end do ! loop over spin



write(80,*) 'Done with MOPDOS ' 

close(unit=10)
close(unit=11)
close(unit=80)
close(unit=81)
close(unit=82)

end program

subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross


