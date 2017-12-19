c        Reads in event and tries to find total energy and event topology

        include 'inputs/chargeReconstructionCommon.f'
        character*20 argncube, argnconfig, argmirror, argnevent
        tmp=rand(123456789)
c get command line arguments
        open(unit=6,file='/dev/null')        !write to null
        call getarg(1,argncube)
        if(argncube.eq.' ')then
                write(6,*)'call using: chargeReconstruction <ncube:1-25> <nconfig:1,2,3> <nmirror:0,1> <nevent>'
                stop
        endif
        call getarg(2,argnconfig)
        call getarg(3,argmirror)
        call getarg(4,argnevent)

        read(argncube,*)ncube
        ncube=min(ncube,25)
        read(argnconfig,*)nconfig
        read(argmirror,*)nmirror
        read(argnevent,*)nevent

c initialize pmts
        do i=1,ncube
        do j=1,ncube
        do k=1,6
        photon(i,j,k)=0
        enddo
        enddo
        enddo

c initialize energy depositions in cells
        do i=1,ncube
        do j=1,ncube
        do k=1,ncube
        ecell(i,j,k)=0
        enddo
        enddo
        enddo

c read in event data
        origtotal=0
        open(unit=4,file='outputs/trans-pmt.dat',err=999)
        do line=1,1000000000
                read(4,*,end=25)i,j,k,t
                photon(i,j,k)=photon(i,j,k)+t
                origtotal=origtotal+t
        enddo
25        close(unit=4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*******************************c
c THIS IS FOR 10 THOU GAPS ONLY c
c*******************************c

c map efficiency for collecting light in given tube based on options
        if(ncube.eq.15.and.nconfig.eq.1.and.nmirror.eq.0)then
c**************************************c
c THIS IS FOR 15x15x15 WITHOUT MIRRORS c
c**************************************c

c**** Warning these are for 2.25" cells

                efactor=2.306 !per keV (determined from center run)
                eff(1) = 1.000
                eff(2) = 1.092
                eff(3) = 1.201
                eff(4) = 1.325
                eff(5) = 1.460
                eff(6) = 1.606
                eff(7) = 1.756
                eff(8) = 1.927
                eff(9) = 2.132
                eff(10)= 2.364
                eff(11)= 2.665
                eff(12)= 2.931
                eff(13)= 3.294
                eff(14)= 3.705
                eff(15)= 4.310
                sb(1)  = 0.0034990792
                sb(2)  = 0.0043091066
                sb(3)  = 0.0057522124
                sb(4)  = 0.006761938
                sb(5)  = 0.0077572965
                sb(6)  = 0.0084495142
                sb(7)  = 0.0093345656
                sb(8)  = 0.0098884381
                sb(9)  = 0.0104318564
                sb(10) = 0.0106343284
                sb(11) = 0.0108695652
                sb(12) = 0.0108712413
                sb(13) = 0.0108318891
                sb(14) = 0.0107212476
                sb(15) = 0.0103174603
        elseif(ncube.eq.15.and.nconfig.eq.1.and.nmirror.eq.1)then
c***********************************c
c THIS IS FOR 15x15x15 WITH MIRRORS c
c***********************************c

c**** The values to the right of the ! are for 2.25" cells

                efactor= 1.937 !2.580 !per keV (determined from centroid of random events)
                eff(1) = 1.000 !1.000
                eff(2) = 1.073 !1.072
                eff(3) = 1.146 !1.142
                eff(4) = 1.218 !1.216
                eff(5) = 1.292 !1.297
                eff(6) = 1.364 !1.367
                eff(7) = 1.431 !1.423
                eff(8) = 1.496 !1.481
                eff(9) = 1.555 !1.549
                eff(10)= 1.611 !1.590
                eff(11)= 1.658 !1.656
                eff(12)= 1.697 !1.673
                eff(13)= 1.727 !1.727
                eff(14)= 1.750 !1.745
                eff(15)= 1.761 !1.752
                
                sb(1)  = 0.0051962002 !0.0052792116
                sb(2)  = 0.0064513714 !0.0065425264
                sb(3)  = 0.0075861941 !0.0080385852
                sb(4)  = 0.0086970257 !0.0094178082
                sb(5)  = 0.0097197868 !0.0105631659
                sb(6)  = 0.0106214973 !0.0118703882
                sb(7)  = 0.0114811642 !0.0131886477
                sb(8)  = 0.0122107174 !0.0140771637
                sb(9)  = 0.0128080200 !0.0150854235
                sb(10) = 0.0134752645 !0.016119403
                sb(11) = 0.0138968526 !0.0165954139
                sb(12) = 0.0143216038 !0.0169937206
                sb(13) = 0.0146469596 !0.0170178282
                sb(14) = 0.0147503340 !0.0168714169
                sb(15) = 0.0148441096 !0.0166461159
        elseif(ncube.eq.11.and.nconfig.eq.1.and.nmirror.eq.0)then
c**************************************c
c THIS IS FOR 11x11x11 WITHOUT MIRRORS c
c**************************************c

c**** Warning these are for 2.25" cells

                write(*,*) 'WARNING: Overall normalization not correctly set!'
                efactor= 2.23253  !per keV (determined from center run)
                eff(1) = 1.000
                eff(2) = 1.108
                eff(3) = 1.224
                eff(4) = 1.358
                eff(5) = 1.507
                eff(6) = 1.678
                eff(7) = 1.880
                eff(8) = 2.118
                eff(9) = 2.390
                eff(10)= 2.735
                eff(11)= 3.174
                sb(1)  = 0.0032146
                sb(2)  = 0.0039503
                sb(3)  = 0.0049281
                sb(4)  = 0.0057320
                sb(5)  = 0.0064991
                sb(6)  = 0.0070570
                sb(7)  = 0.0076542
                sb(8)  = 0.0079940
                sb(9)  = 0.0081792
                sb(10) = 0.0083524
                sb(11) = 0.0083220
        elseif (ncube.eq.11.and.nconfig.eq.1.and.nmirror.eq.1)then
c***********************************c
c THIS IS FOR 11x11x11 WITH MIRRORS c
c***********************************c

c**** Warning these are for 2.25" cells

                efactor= 2.6925  !per keV (determined from center run)
                eff(1) = 1.0000
                eff(2) = 1.0635
                eff(3) = 1.1253
                eff(4) = 1.1804
                eff(5) = 1.2301
                eff(6) = 1.2810
                eff(7) = 1.3256
                eff(8) = 1.3629
                eff(9) = 1.3917
                eff(10)= 1.4166
                eff(11)= 1.4278
                sb(1)  = 0.0050898
                sb(2)  = 0.0066348
                sb(3)  = 0.0081254
                sb(4)  = 0.0093376
                sb(5)  = 0.0102679
                sb(6)  = 0.0114255
                sb(7)  = 0.0120522
                sb(8)  = 0.0127374
                sb(9)  = 0.0132984
                sb(10) = 0.0134973
                sb(11) = 0.0140016
        elseif (ncube.eq.5.and.nconfig.eq.1.and.nmirror.eq.0)then
c***********************************c
c THIS IS FOR 5x5x5 WITHOUT MIRRORS c
c***********************************c

c**** Warning these are for 2.25" cells

                write(*,*) 'WARNING: Overall normalization not correctly set!'
                efactor= 2.61716  !per keV (determined from center run)
                eff(1) = 1.0000
                eff(2) = 1.1425
                eff(3) = 1.3176
                eff(4) = 1.5368
                eff(5) = 1.8286
                sb(1)  = 0.0019407
                sb(2)  = 0.0023125
                sb(3)  = 0.0028409
                sb(4)  = 0.0033491
                sb(5)  = 0.0036019
        elseif (ncube.eq.5.and.nconfig.eq.1.and.nmirror.eq.1)then
c********************************c
c THIS IS FOR 5x5x5 WITH MIRRORS c
c********************************c

c**** Warning these are for 2.25" cells

                efactor= 2.9949  !per keV (determined from center run)
                eff(1) = 1.0000
                eff(2) = 1.0370
                eff(3) = 1.0612
                eff(4) = 1.0839
                eff(5) = 1.0926
                sb(1)  = 0.0038928
                sb(2)  = 0.0049616
                sb(3)  = 0.0055874
                sb(4)  = 0.0060894
                sb(5)  = 0.0063010
        else
          write(*,*)'The reconstruction parameters not found! Exiting the program...'
          goto 999
        endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c******************************c
c THIS IS FOR 5 THOU GAPS ONLY c
c******************************c

c map efficiency for collecting light in given tube based on options
c        if(nconfig.eq.1.and.nmirror.eq.0)then
c                efactor=2.331 !per keV (determined from center run)
c                eff(1) =1.000
c                eff(2) =1.094
c                eff(3) =1.200
c                eff(4) =1.309
c                eff(5) =1.420
c                eff(6) =1.595
c                eff(7) =1.722
c                eff(8) =1.902
c                eff(9) =2.103
c                eff(10)=2.340
c                eff(11)=2.581
c                eff(12)=2.856
c                eff(13)=3.218
c                eff(14)=3.602
c                eff(15)=4.147
c                sb(1)  =0.001694473
c                sb(2)  =0.002494299
c                sb(3)  =0.003050688
c                sb(4)  =0.003411805
c                sb(5)  =0.004072566
c                sb(6)  =0.004573805
c                sb(7)  =0.004600539
c                sb(8)  =0.005205751
c                sb(9)  =0.005482456
c                sb(10) =0.005491153
c                sb(11) =0.005888291
c                sb(12) =0.005956813
c                sb(13) =0.006082215
c                sb(14) =0.005868545
c                sb(15) =0.005675676
c        elseif(nconfig.eq.1.and.nmirror.eq.1)then
c                efactor=2.646 !per keV (determined from centroid of random events)
c                eff(1) =1.000
c                eff(2) =1.063
c                eff(3) =1.134
c                eff(4) =1.204
c                eff(5) =1.280
c                eff(6) =1.344
c                eff(7) =1.397
c                eff(8) =1.463
c                eff(9) =1.526
c                eff(10)=1.578
c                eff(11)=1.617
c                eff(12)=1.645
c                eff(13)=1.681
c                eff(14)=1.713
c                eff(15)=1.717
c                sb(1)  =0.002773925
c                sb(2)  =0.003439803
c                sb(3)  =0.004323899
c                sb(4)  =0.005008347
c                sb(5)  =0.005769231
c                sb(6)  =0.006370416
c                sb(7)  =0.006782946
c                sb(8)  =0.007273342
c                sb(9)  =0.007762879
c                sb(10) =0.008205689
c                sb(11) =0.008594918
c                sb(12) =0.009125475
c                sb(13) =0.009324009
c                sb(14) =0.009303246
c                sb(15) =0.008928571
c        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c set event buffer to zero
        do i=1,1000
        nlogx(i)=0
        nlogy(i)=0
        nlogz(i)=0
        plog(i)=0
        enddo
        
c find which cells have light in all three principle axes

        ncell=0
        do i=1,ncube
        do j=1,ncube
        do k=1,ncube
                de(1)=photon(j,k,1)*eff(i)
                de(2)=photon(j,k,2)*eff(16-i)
                de(3)=photon(i,k,3)*eff(j)
                de(4)=photon(i,k,4)*eff(16-j)
                de(5)=photon(i,j,5)*eff(k)
                de(6)=photon(i,j,6)*eff(16-k)
                if((de(1)+de(2))*(de(3)+de(4))*(de(5)+de(6)).ne.0)then
                        ncell=ncell+1
                        nlogx(ncell)=i
                        nlogy(ncell)=j
                        nlogz(ncell)=k
                        plog(ncell)=min(min(de(1)+de(2),de(3)+de(4)),de(5)+de(6))/2.0
                        if(nmirror.eq.1) plog(ncell)=plog(ncell)*2
                endif
        enddo
        enddo
        enddo
        
c now sort using pointer array

        do i=1,1000
        nptr(i)=i
        enddo
        do i=1,999
        do j=1,999
        if(plog(nptr(j)).lt.plog(nptr(j+1)))then
                nsave=nptr(j)
                nptr(j)=nptr(j+1)
                nptr(j+1)=nsave
        endif
        enddo
        enddo

c sort the real arrays accordingly

        do ii=1,ncell
        i=nptr(ii)
        tnlogx(ii)=nlogx(i)
        tnlogy(ii)=nlogy(i)
        tnlogz(ii)=nlogz(i)
        tplog(ii)=plog(i)
        enddo
        do i=1,ncell
        nlogx(i)=tnlogx(i)
        nlogy(i)=tnlogy(i)
        nlogz(i)=tnlogz(i)
        plog(i)=tplog(i)
        enddo
        
c write output
        sum=0
        do i=1,ncell
                sum=sum+plog(i)
                write(6,"(3i5,2f10.1)")nlogx(i),nlogy(i),nlogz(i),plog(i),plog(i)*6/efactor
        enddo
        einit=sum*6/efactor
        chisqrinit = fchisqr(ncell)
        write(6,*)'Energy:',einit
        write(6,*)'chisq: ',chisqrinit
        
c minimize and repeat


        chiold=1e20
        write(6,"(a,$)")'Called gradls: '
        do ii=1,20
                do i=1,ncell
                        plog(i)=abs(plog(i)) ! avoid negative energies
                        dlog(i)= sqrt(plog(i)/1000)
                enddo
                write(6,"(i3,$)")ii
                call gradls(ncell)
                chinew=fchisqr(ncell)
                if(chinew.ge.chiold)goto 10
                chiold=chinew
        enddo

10        write(6,*)
        do i=1,ncell
                if(plog(i).lt.0.0)plog(i)=0.0 ! do not add any negative energies
        enddo

c repeat sort after chisqr minimization

        do i=1,1000
        nptr(i)=i
        enddo
        do i=1,999
        do j=1,999
        if(plog(nptr(j)).lt.plog(nptr(j+1)))then
                nsave=nptr(j)
                nptr(j)=nptr(j+1)
                nptr(j+1)=nsave
        endif
        enddo
        enddo

c sort the real arrays accordingly

        do ii=1,ncell
        i=nptr(ii)
        tnlogx(ii)=nlogx(i)
        tnlogy(ii)=nlogy(i)
        tnlogz(ii)=nlogz(i)
        tplog(ii)=plog(i)
        enddo
        do i=1,ncell
        nlogx(i)=tnlogx(i)
        nlogy(i)=tnlogy(i)
        nlogz(i)=tnlogz(i)
        plog(i)=tplog(i)
        enddo

c write output to screen
        sum=0
        do i=1,ncell
                sum=sum+plog(i)
                write(6,"(3i5,2f10.1)")nlogx(i),nlogy(i),nlogz(i),plog(i),plog(i)*6/efactor
        enddo
        efinal=sum*6/efactor
        chisqrfinal=fchisqr(ncell)
        write(6,*)'Energy:',efinal
        write(6,*)'chisq: ',chisqrfinal

c create a few nice metrics

c nearest wall
        nwallx=min(nlogx(1),16-nlogx(1))
        nwally=min(nlogy(1),16-nlogy(1))
        nwallz=min(nlogz(1),16-nlogz(1))
        nwall=min(nwallx,min(nwally,nwallz))

c find maximum cell
        nxm=nlogx(1)
        nym=nlogy(1)
        nzm=nlogz(1)
        pmax=plog(1)

c find total energyg
        sum=0
        do i=1,ncell
                sum=sum+plog(i)
        enddo

c add largest adjacent cell
        pnext=0
        do i=2,ncell
                if(abs(nlogx(i)-nxm).eq.1.and.nlogy(i).eq.nym.and.nlogz(i).eq.nzm)pnext=max(pnext,plog(i))
                if(abs(nlogy(i)-nym).eq.1.and.nlogx(i).eq.nxm.and.nlogz(i).eq.nzm)pnext=max(pnext,plog(i))
                if(abs(nlogz(i)-nzm).eq.1.and.nlogy(i).eq.nym.and.nlogx(i).eq.nxm)pnext=max(pnext,plog(i))
        enddo
        
c write to log file
        open(file='outputs/reconstruction.out',access='append',err=999,unit=7)
        write(7,'(A7,I4)')"#Event ", nevent        ! Event delimiter

c        write(7,*)nwall,pmax*6/efactor,pnext*6/efactor,sum*6/efactor        !Not intersted in these can get later in the analysis
        write(7,*)ncell        ! Number of cells hit
        do i=1,ncell           ! Loop over the hit cells and output
                sum=sum+plog(i)
                nx=nlogx(i)
                ny=nlogy(i)
                nz=nlogz(i)
                square=(nx-nxm)**2+(ny-nym)**2+(nz-nzm)**2
                radius=sqrt(square)
c                write(7,"(3i5,2f10.1)")nlogx(i),nlogy(i),nlogz(i),-1.0,plog(i)*6/efactor
                write(7,*)nlogx(i),nlogy(i),nlogz(i),-1.0,(plog(i)*6/efactor)/1000.0
        enddo
        close(unit=7)

999        stop
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        function fchisqr(ncell)
        include 'inputs/chargeReconstructionCommon.f'
        
c zero dummy array
        do i=1,ncube
        do j=1,ncube
        do k=1,6
        dummy(i,j,k)=0
        enddo
        enddo
        enddo

c load dummy array with predicted p.e. numbers
        do i=1,ncell
                nx=nlogx(i)
                ny=nlogy(i)
                nz=nlogz(i)
                dec=plog(i)

                tdec=dec/eff(nx)
                dummy(ny,nz,1)=dummy(ny,nz,1)+tdec
                        if(ny.gt.1)dummy(ny-1,nz,1)=dummy(ny-1,nz,1)+tdec*sb(nx)
                        if(ny.lt.15)dummy(ny+1,nz,1)=dummy(ny+1,nz,1)+tdec*sb(nx)
                        if(nz.gt.1)dummy(ny,nz-1,1)=dummy(ny,nz-1,1)+tdec*sb(nx)
                        if(nz.lt.15)dummy(ny,nz+1,1)=dummy(ny,nz+1,1)+tdec*sb(nx)
                tdec=dec/eff(ny)
                dummy(nx,nz,3)=dummy(nx,nz,3)+tdec
                        if(nx.gt.1)dummy(nx-1,nz,3)=dummy(nx-1,nz,3)+tdec*sb(ny)
                        if(nx.lt.15)dummy(nx+1,nz,3)=dummy(nx+1,nz,3)+tdec*sb(ny)
                        if(nz.gt.1)dummy(nx,nz-1,3)=dummy(nx,nz-1,3)+tdec*sb(ny)
                        if(nz.lt.15)dummy(nx,nz+1,3)=dummy(nx,nz+1,3)+tdec*sb(ny)
                tdec=dec/eff(nz)
                dummy(nx,ny,5)=dummy(nx,ny,5)+tdec
                        if(nx.gt.1)dummy(nx-1,ny,5)=dummy(nx-1,ny,5)+tdec*sb(nz)
                        if(nx.lt.15)dummy(nx+1,ny,5)=dummy(nx+1,ny,5)+tdec*sb(nz)
                        if(ny.gt.1)dummy(nx,ny-1,5)=dummy(nx,ny-1,5)+tdec*sb(nz)
                        if(ny.lt.15)dummy(nx,ny+1,5)=dummy(nx,ny+1,5)+tdec*sb(nz)
                if(nmirror.eq.0)then
                        tdec=dec/eff(16-nx)
                        dummy(ny,nz,2)=dummy(ny,nz,2)+tdec
                                if(ny.gt.1)dummy(ny-1,nz,2)=dummy(ny-1,nz,2)+tdec*sb(16-nx)
                                if(ny.lt.15)dummy(ny+1,nz,2)=dummy(ny+1,nz,2)+tdec*sb(16-nx)
                                if(nz.gt.1)dummy(ny,nz-1,2)=dummy(ny,nz-1,2)+tdec*sb(16-nx)
                                if(nz.lt.15)dummy(ny,nz+1,2)=dummy(ny,nz+1,2)+tdec*sb(16-nx)
                        tdec=dec/eff(16-ny)
                        dummy(nx,nz,4)=dummy(nx,nz,4)+tdec
                                if(nx.gt.1)dummy(nx-1,nz,4)=dummy(nx-1,nz,4)+tdec*sb(16-ny)
                                if(nx.lt.15)dummy(nx+1,nz,4)=dummy(nx+1,nz,4)+tdec*sb(16-ny)
                                if(nz.gt.1)dummy(nx,nz-1,4)=dummy(nx,nz-1,4)+tdec*sb(16-ny)


                                if(nz.lt.15)dummy(nx,nz+1,4)=dummy(nx,nz+1,4)+tdec*sb(16-ny)
                        tdec=dec/eff(16-nz)
                        dummy(nx,ny,6)=dummy(nx,ny,6)+tdec
                                if(nx.gt.1)dummy(nx-1,ny,6)=dummy(nx-1,ny,6)+tdec*sb(16-nz)
                                if(nx.lt.15)dummy(nx+1,ny,6)=dummy(nx+1,ny,6)+tdec*sb(16-nz)
                                if(ny.gt.1)dummy(nx,ny-1,6)=dummy(nx,ny-1,6)+tdec*sb(16-nz)
                                if(ny.lt.15)dummy(nx,ny+1,6)=dummy(nx,ny+1,6)+tdec*sb(16-nz)
                endif
        enddo

c find chisq
        fchisqr=0
        do k=1,6
        do i=1,ncube
        do j=1,ncube
        if(dummy(i,j,k).ne.0)then
                fchisqr=fchisqr+(photon(i,j,k)-dummy(i,j,k))**2
c                write(6,*)i,j,k,photon(i,j,k),dummy(i,j,k)
        endif
        enddo
        enddo
        enddo

c        write(6,*)"chisqr: ",fchisqr

        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine gradls(ncell)
        include 'inputs/chargeReconstructionCommon.f'
c        This routine is adapted from the one in Bevington.  

c        save just in case

        do i=1,ncell
        tlog(i)=plog(i)
        enddo

c        find gradient and normalize
        
        step=0.1        
5        chisq1=fchisqr(ncell)

        sum=0
        do i=1,ncell
                delta=step*dlog(i)
                plog(i)=plog(i)+delta
                grad(i)=chisq1-fchisqr(ncell)
                plog(i)=plog(i)-delta
                sum=sum+grad(i)**2
        enddo
        sum=sqrt(sum)
        if(sum.eq.0)then
                step=2*step
                goto 5
        endif
        do i=1,ncell
                grad(i)=dlog(i)*grad(i)/sum
        enddo

c        step until chisq reduces

        test=0
10        do i=1,ncell
        plog(i)=plog(i)+grad(i)
        enddo
        chisq2=fchisqr(ncell)
        if(chisq2.gt.chisq1)then
                test=test+1
                if(test.gt.10)then
                        write(6,*)"no decrease, revert"
                        do i=1,ncell
                        plog(i)=tlog(i)
                        enddo
                        chisqr=fchisqr(ncell)
                        return
                endif
                do i=1,ncell
                plog(i)=plog(i)-grad(i)
                grad(i)=grad(i)/2
                enddo
                goto 10
        endif

c        keep going until chisq stops decreasing

20        do i=1,ncell
        plog(i)=plog(i)+grad(i)
        enddo
        chisq3=fchisqr(ncell)
        if(chisq3.lt.chisq2)then
                chisq1=chisq2
                chisq2=chisq3
                goto 20
        endif

c        if the same, return

        if(chisq3.eq.chisq2)then
                chisqr=fchisqr(ncell)
                return
        endif

c        if it went up, then interpolate to find minimum

        delta=1.0/(1.0+(chisq1-chisq2)/(chisq3-chisq2))+0.5
        do i=1,ncell
        plog(i)=plog(i)-delta*grad(i)
        enddo
        chisqr=fchisqr(ncell)
        if(chisq2.gt.chisqr)then
                do i=1,ncell
                plog(i)=plog(i)+(delta-1.0)*grad(i)
                enddo
                chisqr=fchisqr(ncell)
        endif
        return
        end
