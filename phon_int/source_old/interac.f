!*     Phoninterac contains routines for calculation of redefined phonon
!*     interaction in p-h J-scheme
!*     last update 15.5.2012

      module phoninterac

      contains
      
      subroutine vinth(ipart,ifmx,isimax,lev)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_int.inc'
      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'


      double precision, dimension(:,:,:,:,:),allocatable ::
     *fp,fpn

      integer, dimension(:), allocatable :: jphon,ironp,iropp,
     *ironh,iroph,ndcamn,ndcamp,ius

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      type(level_typ),dimension(*) :: lev

      character*10 namer
      character*30 namefp,namefpn,namecp,namecn,namerpp,
     *namernp,namerph,namernh,namev
     
      jmin=0
      jmax=isimax

      ndlam=3200      
      ndamp=3200
      ifmxx=3200
      allocate (jphon(ifmxx))
      jphon=0

      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)
      ifmx=i


      allocate(ius(ifmx))
      ius=0


      open (3,file='0hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)



      open (3,file='1hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      write(*,*)' Number of 0hom+1hom phonons ',icount

      

      if (ipart.eq.1) then 
              namecp='1phonon/1f_cn.dat'
              namecn='1phonon/1f_cp.dat'
c              namefp='fmat_n.dat'
c              namefpn='fmat_pn.dat'
              namefp='f_n.dat'
              namefpn='f_pn.dat'
              namerpp='1phonon/1f_rnp.dat'
              namernp='1phonon/1f_rpp.dat'
              namerph='1phonon/1f_rnh.dat'
              namernh='1phonon/1f_rph.dat'
              namev='Vint_phon_h_n.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
              ihminp=ihpmn
              ipmaxp=ippmx
      endif

      if (ipart.eq.-1) then 
              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
c              namefp='fmat_p.dat'
c              namefpn='fmat_np.dat'
              namefp='f_p.dat'
              namefpn='f_np.dat'
              namerpp='1phonon/1f_rpp.dat'
              namernp='1phonon/1f_rnp.dat'
              namerph='1phonon/1f_rph.dat'
              namernh='1phonon/1f_rnh.dat'
              namev='Vint_phon_h_p.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
              ihminp=ihnmn
              ipmaxp=ipnmx

      endif

c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihmin,ipmax,ihmin,ihmax,fp)
c     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihminp,ipmaxp,ihmin,ipmax,fpn)
c     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
c     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      

      if (ipart.eq.1)
     *write(*,*)'Calculation of phonon neutron hole interaction'

      if (ipart.eq.-1)
     *write(*,*)'Calculation of phonon proton hole interaction'


       open(33,file=namernp,status='old'
     *,form='unformatted')

       open(34,file=namerpp,status='old'
     *,form='unformatted')

       open(43,file=namernh,status='old'
     *,form='unformatted')

       open(44,file=namerph,status='old'
     *,form='unformatted')


      open(5,file=namev,status='unknown'
     *,form='unformatted')

      do ig=1,ifmx

       if (ius(ig).eq.1) then 

       write(997,*)ig,ifmx


       jig=jphon(ig)

      write(5)ig

      call readro(33,ig,ronp,nronp)
      call readro(34,ig,ropp,nropp)
      call readro(43,ig,ronh,nronh)
      call readro(44,ig,roph,nroph)

    
      do ib=1,ifmx
       jib=jphon(ib)
      
      do isi=0,isimax 
      ifaz=(-1)**(jib+isi)
       
      allocate(ironp(nronp))
      ironp=0 
      allocate(iropp(nropp))
      iropp=0 
      allocate(ironh(nronh))
      ironh=0 
      allocate(iroph(nroph))
      iroph=0 

      call rosub(isi,ib,ropp,nropp,iropp,nropps)
      call rosub(isi,ib,ronp,nronp,ironp,nronps)     
      call rosub(isi,ib,roph,nroph,iroph,nrophs)
      call rosub(isi,ib,ronh,nronh,ironh,nronhs)     
    
     
      do ih1=ihmin,ihmax
       jh1=lev(ih1)%j 
      do ih2=ihmin,ihmax
       jh2=lev(ih2)%j
             
        vint=0.d0

        do ii=1,nropps
         i1=ropp(iropp(ii))%i1
         i2=ropp(iropp(ii))%i2
         vint=vint+dfloat(ifaz)*ropp(iropp(ii))%ro*fp(isi,i1,i2,ih1,ih2)
        enddo

        do ii=1,nronps
         i1=ronp(ironp(ii))%i1
         i2=ronp(ironp(ii))%i2
         vint=vint+dfloat(ifaz)
     **ronp(ironp(ii))%ro*fpn(isi,i1,i2,ih1,ih2)
        enddo

        do ii=1,nrophs
         i1=roph(iroph(ii))%i1
         i2=roph(iroph(ii))%i2
         vint=vint+0.5d0*dfloat(ifaz)*roph(iroph(ii))%ro
     **fp(isi,i1,i2,ih1,ih2)
        enddo

        do ii=1,nronhs
         i1=ronh(ironh(ii))%i1
         i2=ronh(ironh(ii))%i2
         vint=vint+dfloat(ifaz)
     **ronh(ironh(ii))%ro*fpn(isi,i1,i2,ih1,ih2)
        enddo


        factor=dfloat((-1)**(jib+(jh1+jh2)/2))*
     *((2*jib+1)*(2*jig+1)*(2*isi+1))**0.5d0

        do ii=1,ndcamn(ib)
          ipp=camn(ib,ii)%par
          ihp=camn(ib,ii)%hol
          campp=camn(ib,ii)%am
         do jj=1,ndcamp(ig)
           ip=camp(ig,jj)%par
           ih=camp(ig,jj)%hol
          if (ih.eq.ih2) then 
           jp=lev(ip)%j
           xsixj=sixj(2*jib,2*isi,2*jig,jh2,jp,jh1)          
           campn=camp(ig,jj)%am
           vint=vint+factor*campp*campn*xsixj*fpn(jib,ipp,ihp,ih1,ip)
          endif
         enddo
        enddo



             
c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.1) 
c     *write(997,102)isi,ig,ib,ih1,ih2,vint

c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.-1) 
c     *write(998,102)isi,ig,ib,ih1,ih2,vint


        if (dabs(vint).gt.xrotrunc) 
     *write(5)ib,isi,ih1,ih2,vint


        enddo ! loop ih2
      enddo ! loop ih1
      deallocate (ironp,ironh,iropp,iroph)
      enddo ! loop isi
      enddo ! loop ib
      endif
      write(5)0,0,0,0,0.d0      
      enddo ! loop ig      

      deallocate(camp,camn,jphon,fp,fpn,ronp,ronh,ropp,roph)

      write(5)10000000
      write(5)0,0,0,0,0.d0 

      close(33)
      close(34)
      close(43)
      close(44)
      close(5)
      return
      end subroutine vinth
!*************************************************************************
      subroutine vintp(ipart,ifmx,isimax,lev)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_int.inc'
      include 'input_phon_int.inc'
      include 'formats_phon_int.inc'


      double precision, dimension(:,:,:,:,:),allocatable ::
     *fp,fpn
      double precision, dimension(:,:,:,:,:,:), allocatable ::
     *sixjs

      integer, dimension(:), allocatable :: jphon,ironp,iropp,
     *ironh,iroph,ndcamn,ndcamp,ius

      type(amp_typ), dimension(:,:), allocatable :: camp,camn
      type(rho_typ), dimension(:), allocatable :: ronp,ropp,ronh,roph
      type(level_typ),dimension(*) :: lev

      character*10 namer
      character*30 namefp,namefpn,namecp,namecn,namerpp,
     *namernp,namerph,namernh,namev
     
      jmin=0
      jmax=isimax

      ndlam=3200      
      ndamp=3200
      ifmxx=3200
      allocate (jphon(ifmxx))
      jphon=0

      allocate(sixjs(0:jmax,0:jmax,0:jmax,0:jmax,
     *0:jmax,0:jmax))
      sixjs=0.0d0

      do i1=0,jmax
       do i2=0,jmax
        do i3=0,jmax
         do i4=1,jmax
          do i5=1,jmax
           do i6=1,jmax
       sixjs(i1,i2,i3,i4,i5,i6)=sixj(2*i1,2*i2,2*i3,i4,i5,i6)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      open (3,file='1phonon/1f_states.dat',
     *status='old',form='unformatted')

      do while (.not.eof(3))
       read(3)i,ipar,ijj,en
c       write(991,*)i,ipar,ijj,en
       jphon(i)=ijj
      enddo

      close(3)

      ifmx=i

      allocate(ius(ifmx))
      ius=0

      open (3,file='0hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)


      open (3,file='1hom_phon.dat',
     *status='old',form='formatted')

      do while (.not.eof(3))
       read(3,*)i,ee
       ius(i)=1
      enddo

      close(3)

      icount=0
      do i=1,ifmx
       if (ius(i).eq.1) icount=icount+1
      enddo

      write(*,*)' Number of 0hom+1hom phonons ',icount 

!      open (3,file='2hom_phon.dat',
!     *status='old',form='formatted')

!      do while (.not.eof(3))
!       read(3,*)i,ee
!       ius(i)=1
!      enddo

!      close(3)




      if (ipart.eq.1) then 
              namecp='1phonon/1f_cn.dat'
              namecn='1phonon/1f_cp.dat'
c              namefp='fmat_n.dat'
c              namefpn='fmat_pn.dat'
              namefp='f_n.dat'
              namefpn='f_pn.dat'
              namerpp='1phonon/1f_rnp.dat'
              namernp='1phonon/1f_rpp.dat'
              namerph='1phonon/1f_rnh.dat'
              namernh='1phonon/1f_rph.dat'
              namev='Vint_phon_p_n.dat'
              ipmin=ipnmn
              ipmax=ipnmx
              ihmin=ihnmn
              ihmax=ihnmx
              ihminp=ihpmn
              ipmaxp=ippmx
      endif

      if (ipart.eq.-1) then 
              namecp='1phonon/1f_cp.dat'
              namecn='1phonon/1f_cn.dat'
c              namefp='fmat_p.dat'
c              namefpn='fmat_np.dat'
              namefp='f_p.dat'
              namefpn='f_np.dat'
              namerpp='1phonon/1f_rpp.dat'
              namernp='1phonon/1f_rnp.dat'
              namerph='1phonon/1f_rph.dat'
              namernh='1phonon/1f_rnh.dat'
              namev='Vint_phon_p_p.dat'
              ipmin=ippmn
              ipmax=ippmx
              ihmin=ihpmn
              ihmax=ihpmx
              ihminp=ihnmn
              ipmaxp=ipnmx

      endif

c     loads F(p) or F(n) interaction
      call readfin(namefp,jmin,jmax,ihmin,ipmax,ipmin,ipmax,fp)
c     loads F(pn) 
      call readfin(namefpn,jmin,jmax,ihminp,ipmaxp,ihmin,ipmax,fpn)
c     loads Cph protons
      call readcam(namecp,ndlam,ndamp,camp,ndcamp)
c     loads Cph neutrons
      call readcam(namecn,ndlam,ndamp,camn,ndcamn)
      

      if (ipart.eq.1)
     *write(*,*)'Calculation of phonon neutron particle interaction'

      if (ipart.eq.-1)
     *write(*,*)'Calculation of phonon proton particle interaction'


       open(33,file=namernp,status='old'
     *,form='unformatted')

       open(34,file=namerpp,status='old'
     *,form='unformatted')

       open(43,file=namernh,status='old'
     *,form='unformatted')

       open(44,file=namerph,status='old'
     *,form='unformatted')


      open(5,file=namev,status='unknown'
     *,form='unformatted')

      do ig=1,ifmx
       if (ius(ig).eq.1) then 
       write(996,*)ig
       jig=jphon(ig)

      write(5)ig

      call readro(33,ig,ronp,nronp)
      call readro(34,ig,ropp,nropp)
      call readro(43,ig,ronh,nronh)
      call readro(44,ig,roph,nroph)

    
      do ib=1,ifmx
       jib=jphon(ib)
      
      do isi=0,isimax 
      ifaz=(-1)**(jib+isi)
       
      allocate(ironp(nronp))
      ironp=0 
      allocate(iropp(nropp))
      iropp=0 
      allocate(ironh(nronh))
      ironh=0 
      allocate(iroph(nroph))
      iroph=0 

      call rosub(isi,ib,ropp,nropp,iropp,nropps)
      call rosub(isi,ib,ronp,nronp,ironp,nronps)     
      call rosub(isi,ib,roph,nroph,iroph,nrophs)
      call rosub(isi,ib,ronh,nronh,ironh,nronhs)     
    
     
      do ip1=ipmin,ipmax
       jp1=lev(ip1)%j 
      do ip2=ipmin,ipmax
       jp2=lev(ip2)%j
             
        vint=0.d0
        if (iabs(jp1-jp2)/2.le.isi.and.(jp1+jp2)/2.ge.isi) then 

        do ii=1,nropps
         i1=ropp(iropp(ii))%i1
         i2=ropp(iropp(ii))%i2
         vint=vint+0.5d0*dfloat(ifaz)
     **ropp(iropp(ii))%ro*fp(isi,i1,i2,ip1,ip2)
        enddo

        do ii=1,nronps
         i1=ronp(ironp(ii))%i1
         i2=ronp(ironp(ii))%i2
         vint=vint+dfloat(ifaz)
     **ronp(ironp(ii))%ro*fpn(isi,i1,i2,ip1,ip2)
        enddo

        do ii=1,nrophs
         i1=roph(iroph(ii))%i1
         i2=roph(iroph(ii))%i2
         vint=vint+dfloat(ifaz)*roph(iroph(ii))%ro*fp(isi,i1,i2,ip1,ip2)
        enddo

        do ii=1,nronhs
         i1=ronh(ironh(ii))%i1
         i2=ronh(ironh(ii))%i2
         vint=vint+dfloat(ifaz)
     **ronh(ironh(ii))%ro*fpn(isi,i1,i2,ip1,ip2)
        enddo


        factor=dfloat((-1)**(isi-jig+(jp1+jp2)/2))*
     *((2*jib+1)*(2*jig+1)*(2*isi+1))**0.5d0

!        do ii=1,ndcamn(ib)
!          ipp=camn(ib,ii)%par
!          ihp=camn(ib,ii)%hol
!          campp=camn(ib,ii)%am
         do jj=1,ndcamp(ig)
           ip=camp(ig,jj)%par
!           ih=camp(ig,jj)%hol 
          if (ip.eq.ip1) then
           ih=camp(ig,jj)%hol
           jh=lev(ih)%j
           xsixj=sixjs(jib,isi,jig,jp1,jh,jp2)

           do ii=1,ndcamn(ib)
            ipp=camn(ib,ii)%par
            ihp=camn(ib,ii)%hol
            campp=camn(ib,ii)%am
!            jh=lev(ih)%j
!           xsixj=sixj(2*jib,2*isi,2*jig,jp1,jh,jp2)          
!            xsixj=sixjs(jib,isi,jig,jp1,jh,jp2)
            campn=camp(ig,jj)%am
            vint=vint-factor*campp*campn*xsixj*fpn(jib,ipp,ihp,ih,ip2)
 
           enddo
          endif
        enddo

        endif

c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.1) 
c     *write(897,102)isi,ig,ib,ip1,ip2,vint

c        if (dabs(vint).gt.xrotrunc.and.ipart.eq.-1) 
c     *write(898,102)isi,ig,ib,ip1,ip2,vint

        if (dabs(vint).gt.xrotrunc) 
     *write(5)ib,isi,ip1,ip2,vint


        enddo ! loop ih2
      enddo ! loop ih1
      deallocate (ironp,ironh,iropp,iroph)
      enddo ! loop isi
      enddo ! loop ib
      endif 
      write(5)0,0,0,0,0.d0 

      enddo ! loop ig      

      deallocate(camp,camn,jphon,fp,fpn,ronp,ronh,ropp,roph)

      write(5)10000000
      write(5)0,0,0,0,0.d0 

      close(33)
      close(34)
      close(43)
      close(44)
      close(5)
      return
      end subroutine vintp


************************************************************************
      subroutine readcam(fname,ndimi,ndimj,cam,ndcc)

      implicit double precision (a-h,o-z)

      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'

      type (amp_typ), dimension(:,:), allocatable :: cam
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      allocate(cam(ndimi,ndimj))
      allocate(ndcc(ndimi))
  
      open(2,file=fname,status='old',form='unformatted')

      ilam=0

      do while (.not.eof(2))
       ilam=ilam+1
       if (ilam.gt.ndimi) then 
               write(*,*)'Readcam: allocate bigger array in ndimi'
               stop
           endif
       read(2)ipar,ijj,ndc
         if (ndc.gt.ndimj) then 
            write(*,*)'Readcam: allocate bigger array in ndimj'
               stop
           endif
  
       read(2)(cam(ilam,i)%par,cam(ilam,i)%hol,cam(ilam,i)%am,i=1,ndc)
       ndcc(ilam)=ndc
      enddo

      close(2)

      return
      end subroutine readcam

***********************************************************************
      subroutine readro(ifile,ig,ron,ndgg)

      implicit double precision (a-h,o-z)

      include 'formats_phon_int.inc'
      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=30)fname

      rewind(ifile)

      ndro=10000000
      ndgg=0

      if (.not.allocated(ron)) allocate (ron(ndro))
  
       do igt=1,ig-1
         read(ifile)igggt,nggt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)itt
         read(ifile)xtt
       enddo


       read(ifile)igg,ndgg

       if (igg.ne.ig) then
               write(*,*)' Loaded Ig does not match !!! '
               stop
       endif
       
       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro'
                stop
       endif

       read(ifile)(ron(ii)%ilap,ii=1,ndgg)
       read(ifile)(ron(ii)%j,ii=1,ndgg)
       read(ifile)(ron(ii)%i1,ii=1,ndgg)
       read(ifile)(ron(ii)%i2,ii=1,ndgg)
       read(ifile)(ron(ii)%ro,ii=1,ndgg)


      return
      end subroutine readro
***********************************************************************
      subroutine rosub(j,ilam,rop,nrop,irop,nrops)

      implicit double precision (a-h,o-z)

      include 'types_phon_int.inc'

      type(rho_typ), dimension(:), allocatable :: rop
      integer, dimension(:), allocatable :: irop

      ii=0
      do i=1,nrop
      if (rop(i)%j.eq.j.and.rop(i)%ilap.eq.ilam) then
              ii=ii+1
              irop(ii)=i
             endif
 
      enddo

      nrops=ii

      end subroutine rosub
***********************************************************************

      subroutine readfin(fname,jmin,jmax,imin,imax,kmin,kmax,fpp)

      implicit double precision (a-h,o-z)

      include 'formats_phon_int.inc'

      double precision, dimension(:,:,:,:,:), allocatable ::
     *fpp

      character(len=30)fname

      allocate(fpp(jmin:jmax,
     *imin:imax,imin:imax,kmin:kmax,kmin:kmax))
      fpp=0.d0

c      open(2,file=fname,status='old',form='formatted')
      open(2,file=fname,status='old',form='unformatted')

      do while (.not.eof(2))
c       read(2,10)itt,ipt,ijt,i,j,k,l,vint
       read(2)itt,ipt,ijt,i,j,k,l,vint


c        if (ipt.ne.ipar) goto 11
      if (ijt.gt.jmax.or.ijt.lt.jmin) goto 11
      
      if (i.gt.imax.or.i.lt.imin
     *.or.j.gt.imax.or.j.lt.imin
     *.or.k.gt.kmax.or.k.lt.kmin
     *.or.l.gt.kmax.or.l.lt.kmin) goto 11


      fpp(ijt,i,j,k,l)=vint
  
 11   enddo
      close(2)


      return
      end subroutine readfin

     
      end module phoninterac 
      
      
