c call using: transport <ncube> <nconfig:1,2,3> <nmirror:0,1> <energy in MeV> <nxc> <nyc> <nzc> <ncombine:0-no guide,1,2,3> <"add" (optional)>
c Note that by addint 'add' to the end of the command line, it will add the
c results from the previous run.  Otherwise, it will just over-write them.
c mirror configurations:
c	0 - no mirrors
c	1 - mirror around PMTs (no guides), on planes 2,4,6 (with guides)

	character*20 argtr,arge,argt,argnx,argny,argnz
	character*100 argconfigfilename,argchargefilename,argtimefilename
	character*3 argadd
	real timezy(12)
        real cdfzy(12)
	include 'inputs/trans-common.f'
	dimension nssum(6)
	
	call srand(time())
c	call srand(1013)
	tmp=rand()

c get command line arguments
	open(unit=6,file='/dev/null')	!write to null
	call getarg(1,argtr)
	if(argtr.eq.' ')then
		write(6,*)'call using: trans <time response flag> <energy in MeV> <time in ps> <nxc> <nyc> <nzc> <config file> <"add" (optional)>'
		stop
	endif
	call getarg(2,arge)
	call getarg(3,argt)
	call getarg(4,argnx)
	call getarg(5,argny)
	call getarg(6,argnz)
	call getarg(7,argconfigfilename)
	call getarg(8,argadd)

        read(argtr,*)nscintimeresponse
	read(arge,*)energy
	read(argt,*)starttime


	read(argnx,*)rnxc
	read(argny,*)rnyc
	read(argnz,*)rnzc
	nxc=rnxc
	nyc=rnyc
	nzc=rnzc
	ntflag=0
	tcheck=rnxc+rnyc+rnzc-nxc-nyc-nzc
	if(tcheck.ne.0)ntflag=1

	open(unit=4,file=argconfigfilename)
	
	read(4,*)ncube
	read(4,*)nconfig
	read(4,*)nmirror
	read(4,*)ncombine
	read(4,*)npmtcouple
	read(4,*)argchargefilename
	read(4,*)argtimefilename

	close(unit=4)
	
	write(6,*)argchargefilename
	write(6,*)argtimefilename

        if(ncombine.eq.0)then
                ngflag=0
        else
                ngflag=0
                if(ncombine.gt.0)ngflag=1
        endif

c user adjustable variables
c detector has one corner at 0,0,0 and fills the positive octant
c cubes are of unit dimension to facilitate coding; this means
c that attenuation lengths afre cell path lengths times cell
c dimensions.
c All dimensions are in cm.

c read in parameters from input file
	open(unit=4,file='inputs/trans.dat')
	
	read(4,*)pmtreflect
	read(4,*)ntbin
	read(4,*)ntprint
	read(4,*)celldim
	read(4,*)pmtr
	read(4,*)tfilm
	read(4,*)glength
	read(4,*)gztrans
	if(glength.le.gztran)then
		write(6,*)'Guide geometry error:'
		write(6,*)'glength: ',glength
		write(6,*)'gztrans:   ',gztrans
		stop
	endif
	read(4,*)agrcoef
	read(4,*)tgfilm

	read(4,*)pcyield
	read(4,*)qe
	read(4,*)relativeyield
	read(4,*)rscint
	read(4,*)rshield
	read(4,*)ascint
	read(4,*)agfilm
	read(4,*)agbuffer

	read(4,*)afilm1
	read(4,*)nfilms1
	read(4,*)rgap1
	read(4,*)gapwidth1
	read(4,*)afilm2
	read(4,*)nfilms2
	read(4,*)rgap2
	read(4,*)gapwidth2
	read(4,*)afilm3
	read(4,*)nfilms3
	read(4,*)rgap3
	read(4,*)gapwidth3

	close(unit=4)

c to simplify comparisons, have some predefined configurations
c nconfig=1
	if(nconfig.eq.1)then
	afilm=afilm1
	nfilms=nfilms1
	rgap=rgap1
	gapwidth=gapwidth1

c nconfig=2
	elseif(nconfig.eq.2)then
	afilm=afilm2
	nfilms=nfilms2
	rgap=rgap2
	gapwidth=gapwidth2

c nconfig=3
	elseif(nconfig.eq.3)then
	afilm=afilm3
	nfilms=nfilms3
	rgap=rgap3
	gapwidth=gapwidth3

	endif

c initialize common parameters
	pi=3.141592654
	epsilon=1.0e-5
	finity=1.0e5
	clight=0.030 !cm/ps

c initialize counters
	nalost=0		!how many lost due to scint absorption
	nflost=0		!how many lost in film
	ntrapped=0		!how many were completely trapped
	nerror=0		!how many had errors within the detector region
	nglost=0                !how many were lost in channel to pmt transition
	
	ngstart=0		!how many started into guide
	ngflost=0		!how many were lost in guide film
	ngblost=0		!how many were lost in guide fill
	ngtlost=0		!how many were lost because of too many bounces in guide
	ngrlost=0		!how many were lost due to imperfect reflection of VM2000
	ngerror=0		!how many were lost due to geometry errors in the guide region
	ngreflect=0		!how many were reflected back into detector
	nsuccess=0		!how many were detected in PMTs
	nwithincrit=0		!how many were being 'channeled' into guide
	nsc=0			!how many 'channeled' reached PMTs

c initialize array
	do i=1,ncube
	do j=1,ncube
	do k=1,6
	nphoton(i,j,k)=0
	enddo
	enddo
	enddo
	do i=1,ncube
	do j=1,ncube
	do k=1,6
	do l=1,200
	ntp(i,j,k,l)=0
	enddo
	enddo
	enddo
	enddo

c	Calculate number of photons to transport
	num=pcyield*relativeyield*energy
c	num=num+dnormal()/sqrt(1.0*num) fix here
c if using guides, calculate their plane definitions
	if(ngflag.eq.1)call getplanes

c need next line when determining if photon gets into pmt from square guide
c (remember that cell dimensions are in units of celldim)
	fsq=(pmtr/(celldim))**2 !radius squared for pmts

c determine critical angle for total internal reflection (within detector)
	tcrit=0
	ccrit=1
	if(rgap.le.rscint)then
		tcrit=asin(rgap/rscint)
		ccrit=cos(tcrit)
	endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c MAIN LOOP
	do ntrack=1,num
	if(rand().gt.qe)goto 20 !did it satisfy quantum efficiency?
	nint=0 !used to find 'trapped' photons
	
	emissiontime=0.0
	if(nscintimeresponse.eq.0)then
	        emissiontime=-2200.0*log(rand())
	elseif(nscintimeresponse.eq.1)then
	        timezy(1)=0.0           !Odd name to ensure that I don't re-use a variable here
	        timezy(2)=5000.0
	        timezy(3)=10000.0
	        timezy(4)=20000.0
	        timezy(5)=30000.0
	        timezy(6)=40000.0

                cdfzy(1)=0.0000
                cdfzy(2)=0.6020
                cdfzy(3)=0.9272
                cdfzy(4)=0.9934
                cdfzy(5)=0.9994
                cdfzy(6)=1.0000

                rnzy=rand()
                nlowzy=0
                nhighzy=1
                do izy=1,5
                        if(cdfzy(izy).le.rnzy.and.rnzy.lt.cdfzy(izy+1))then
                                nlowzy=izy
                                nhighzy=izy+1
                        endif
                enddo
                if(nlowzy.eq.0)goto 20

                p1xzy=timezy(nlowzy)
                p1yzy=cdfzy(nlowzy)
                p2xzy=timezy(nhighzy)
                p2yzy=cdfzy(nhighzy)

                slopezy=(p2yzy-p1yzy)/(p2xzy-p1xzy)
                bzy=p1yzy-slopezy*p1xzy
                emissiontime=(rnzy-bzy)/slopezy
        elseif(nscintimeresponse.eq.2)then
                timezy(1)=0.0
                timezy(2)=5000.0
                timezy(3)=10000.0
                timezy(4)=20000.0
                timezy(5)=30000.0
                timezy(6)=40000.0
                timezy(7)=50000.0
                timezy(8)=60000.0
                timezy(9)=70000.0
                timezy(10)=80000.0
                timezy(11)=90000.0
                timezy(12)=100000.0

                cdfzy(1)=0.0000
                cdfzy(2)=0.4866
                cdfzy(3)=0.7494
                cdfzy(4)=0.8370
                cdfzy(5)=0.8905
                cdfzy(6)=0.9197
                cdfzy(7)=0.9392
                cdfzy(8)=0.9562
                cdfzy(9)=0.9708
                cdfzy(10)=0.9830
                cdfzy(11)=0.9927
                cdfzy(12)=1.0000

                rnzy=rand()
                nlowzy=0
                nhighzy=1
                do izy=1,11
                        if(cdfzy(izy).le.rnzy.and.rnzy.lt.cdfzy(izy+1))then
                                nlowzy=izy
                                nhighzy=izy+1
                        endif
                enddo
                if(nlowzy.eq.0)goto 20

                p1xzy=timezy(nlowzy)
                p1yzy=cdfzy(nlowzy)
                p2xzy=timezy(nhighzy)
                p2yzy=cdfzy(nhighzy)

                slopezy=(p2yzy-p1yzy)/(p2xzy-p1xzy)
                bzy=p1yzy-slopezy*p1xzy
                emissiontime=(rnzy-bzy)/slopezy
	elseif(nscintimeresponse.eq.3)then
                timezy(1)=0.0
                timezy(2)=5000.0
                timezy(3)=10000.0
                timezy(4)=20000.0
                timezy(5)=30000.0
                timezy(6)=40000.0
                timezy(7)=50000.0
                timezy(8)=60000.0
                timezy(9)=70000.0
                timezy(10)=80000.0
                timezy(11)=90000.0
                timezy(12)=100000.0

                cdfzy(1)=0.0000
                cdfzy(2)=0.4032
                cdfzy(3)=0.6210
                cdfzy(4)=0.7056
                cdfzy(5)=0.7621
                cdfzy(6)=0.8065
                cdfzy(7)=0.8468
                cdfzy(8)=0.8831
                cdfzy(9)=0.9173
                cdfzy(10)=0.9476
                cdfzy(11)=0.9758
                cdfzy(12)=1.0000

                rnzy=rand()
                nlowzy=0
                nhighzy=1
                do izy=1,11
                        if(cdfzy(izy).le.rnzy.and.rnzy.lt.cdfzy(izy+1))then
                                nlowzy=izy
                                nhighzy=izy+1
                        endif
                enddo
                if(nlowzy.eq.0)goto 20

                p1xzy=timezy(nlowzy)
                p1yzy=cdfzy(nlowzy)
                p2xzy=timezy(nhighzy)
                p2yzy=cdfzy(nhighzy)

                slopezy=(p2yzy-p1yzy)/(p2xzy-p1xzy)
                bzy=p1yzy-slopezy*p1xzy
                emissiontime=(rnzy-bzy)/slopezy
        else
                emissiontime=0.0
	endif
	
	dtpath=starttime+emissiontime !used to track time of arrival for photons

c find random direction and then direction cosines
	theta=acos(2*rand()-1)
	phi=rand()*2*pi
	a=sin(theta)*cos(phi)
	b=sin(theta)*sin(phi)
	c=cos(theta)
	abc=sqrt(a*a+b*b+c*c)
	a=a/abc
	b=b/abc
	c=c/abc

c pick random location just within cell
	x=nxc-1+epsilon+(1-2*epsilon)*rand()
	y=nyc-1+epsilon+(1-2*epsilon)*rand()
	z=nzc-1+epsilon+(1-2*epsilon)*rand()
	if(ntflag.eq.1)then
		x=rnxc-1
		y=rnyc-1
		z=rnzc-1
	endif

c return here to transport to next intersection	

10	call transport(x,y,z,a,b,c,np,theta,d)
c Used to find which plane will get struck; the angle between
c that plane's normal and the incoming photon; and the distance
c the photon will have to travel to reach the struck plane.
c Don't actually transport photon yet, since we want to know
c if it will reflect or not, so we can transport just a little
c short or a little beyond the exact plane location.

cc	rphoton=sqrt((x-ncube/2.0)**2+(y-ncube/2.0)**2+(z-ncube/2.0)**2)
cc	write(6,'(i6,f10.2)')ntrack,rphoton
	if(d.lt.0)then
c		write(6,'(a,i6,6f11.7,i3)')"D neg ",ntrack,x,y,z,a,b,c,np
		nerror=nerror+1
		goto 20
	endif
	xold=x
	yold=y
	zold=z
	x=x+a*d
	y=y+b*d
	z=z+c*d
	dtpath=dtpath+d*celldim*rscint/clight

c reflect if necessary
	nrflg=0
	if(rgap.ne.rscint)then
		if(theta.gt.tcrit)then !theta came from the transport call
			nrflg=1
		elseif(rand().gt.coef(theta,rscint,rgap))then
			nrflg=1
		endif
	endif
	if(nrflg.eq.1)then
		if(np.eq.1.or.np.eq.2)then
			if(z-int(z).gt.gapwidth/celldim)a=-a
		elseif(np.eq.3.or.np.eq.4)then
			if(x-int(x).gt.gapwidth/celldim)b=-b
		elseif(np.eq.5.or.np.eq.6)then
			if(y-int(y).gt.gapwidth/celldim)c=-c
		endif
	endif

c Move epsilon along new direction

	x=x+a*epsilon
	y=y+b*epsilon
	z=z+c*epsilon
	
c check if still within detector at all
c	test=ncube/2.0
c	if(abs(x-test).gt.test+10*epsilon)write(6,'(i6,9f9.5)')ntrack,xold,yold,zold,x,y,z,a,b,c
c	if(abs(y-test).gt.test+10*epsilon)write(6,'(i6,9f9.5)')ntrack,xold,yold,zold,x,y,z,a,b,c
c	if(abs(z-test).gt.test+10*epsilon)write(6,'(i6,9f9.5)')ntrack,xold,yold,zold,x,y,z,a,b,c

c increment the counter for the number of intersections
	nint=nint+1

c check for various ways the photon might have been lost

c were there too many intersections with the cell walls?
	if(nint.gt.400)then
		ntrapped=ntrapped+1
c		write(6,'(9f10.5)')x,y,z,a,b,c,a/b,b/c,c/a
		goto 20
	endif

c was the photon lost by attenuation in scintillator?
	acoef=exp(-d*celldim/ascint)
	if(rand().gt.acoef)then
		nalost=nalost+1
		goto 20
	endif

c was the photon lost by attenuation in the film?
	if(nfilms.eq.1)then
		if(nrflg.eq.0)then !did not pass through film if reflected
			fcoef=exp(-tfilm/cos(theta)/afilm)
			if(rand().gt.fcoef)then
				nflost=nflost+1
				goto 20
			endif
		endif
	elseif(nfilms.eq.2)then
		fcoef=exp(-2*(tfilm/cos(theta))/afilm) !always goes through two film layers
		if(rand().gt.fcoef)then
			nflost=nflost+1
			goto 20
		endif
	endif		
	
c check if passed a detector side (remember, we are epsilon off from any plane)



	nhitside=0
	if(nrflg.eq.0)then !if it reflected from this wall, then clearly it didn't go into the guide or PMT
		if(x.lt.0.and.a.lt.0)then !passing through the x=0 plane, heading towards negative x
			nplane=1
			x1=y
			y1=z
			nrefl=1
			nhitside=1
		elseif(x.gt.ncube.and.a.gt.0)then
			nplane=2
			x1=y
			y1=z
			nrefl=1
			nhitside=1
		elseif(y.lt.0.and.b.lt.0)then
			nplane=3
			x1=x
			y1=z
			nrefl=2
			nhitside=1
		elseif(y.gt.ncube.and.b.gt.0)then
			nplane=4
			x1=x
			y1=z
			nrefl=2
			nhitside=1
		elseif(z.lt.0.and.c.lt.0)then
			nplane=5
			x1=x
			y1=y
			nrefl=3
			nhitside=1
		elseif(z.gt.ncube.and.c.gt.0)then
			nplane=6
			x1=x
			y1=y
			nrefl=3
			nhitside=1
		endif
	endif
c if nhitside is 1, log this photon (or transport it through guide - which may in fact reflect it)
	if(nhitside.eq.1)then
c find out if it was being 'channeled' (ie: within the critical angle of planes perpindicular to struck plane)
		nwc=0
		if(nplane.le.2)then
			if(abs(b).lt.ccrit.and.abs(c).lt.ccrit)nwc=1
		elseif(nplane.le.4)then
			if(abs(a).lt.ccrit.and.abs(c).lt.ccrit)nwc=1
		elseif(nplane.le.6)then
			if(abs(a).lt.ccrit.and.abs(b).lt.ccrit)nwc=1
		endif

		nx1=int(x1)+1
		ny1=int(y1)+1
		nx1=max(1,min(nx1,ncube))
		ny1=max(1,min(ny1,ncube))
c check if actually made it to a pmt
		if(ngflag.eq.0)then
			if(nwc.eq.1)nwithincrit=nwithincrit+1 !number entering which are 'channeled'
			rsq=((x1-int(x1)-0.5)**2+(y1-int(y1)-0.5)**2)
			if(rsq.gt.fsq)then !not within PMT radius
				if(nmirror.eq.0)then
					nglost=nglost+1
				else
					if(rand().lt.agrcoef)then
						if(nplane.eq.1.or.nplane.eq.2)then
							a=-a
						elseif(nplane.eq.3.or.nplane.eq.4)then
							b=-b
						elseif(nplane.eq.5.or.nplane.eq.6)then
							c=-c
						endif
						x=x+a*2*epsilon
						y=y+b*2*epsilon
						z=z+c*2*epsilon
						goto 10
					else
						nglost=nglost+1
					endif
				endif
			else
				if(nmirror.eq.0.or.(nplane/2.ne.nplane/2.0))then !no mirror, or only for planes 1,3,5
					nphoton(nx1,ny1,nplane)=nphoton(nx1,ny1,nplane)+1
					ndt=max(1.0,min(dtpath/ntbin,1.0*ntprint))
					ntp(nx1,ny1,nplane,ndt)=ntp(nx1,ny1,nplane,ndt)+1 !log time of arrival for this photon
					if(rand().lt.pmtreflect)then !counted, but also reflected
						if(nplane.eq.1.or.nplane.eq.2)then
							a=-a
						elseif(nplane.eq.3.or.nplane.eq.4)then
							b=-b
						elseif(nplane.eq.5.or.nplane.eq.6)then
							c=-c
						endif
						x=x+a*2*epsilon
						y=y+b*2*epsilon
						z=z+c*2*epsilon
						goto 10
					endif
				elseif(nmirror.eq.1.and.(nplane/2.eq.nplane/2.0))then !mirror, only for planes 2,4,6
					if(rand().lt.agrcoef)then
						if(nplane.eq.1.or.nplane.eq.2)then
							a=-a
						elseif(nplane.eq.3.or.nplane.eq.4)then
							b=-b
						elseif(nplane.eq.5.or.nplane.eq.6)then
							c=-c
						endif
						x=x+a*2*epsilon
						y=y+b*2*epsilon
						z=z+c*2*epsilon
						goto 10
					else
						nglost=nglost+1
					endif
				endif
			endif
			goto 20
		else
			if(nmirror.eq.1.and.nplane/2.eq.nplane/2.0)then !only for planes 2,4,6
				if(rand().lt.agrcoef)then !did it actually reflect, or get absorbed
					if(nplane.eq.2)then
						a=-a
					elseif(nplane.eq.4)then
						b=-b
					elseif(nplane.eq.6)then
						c=-c
					endif
					x=x+a*2*epsilon
					y=y+b*2*epsilon
					z=z+c*2*epsilon
					goto 10
				else
					nglost=nglost+1 !was absorbed
					goto 20
				endif
			endif
			call planerefract(a,b,c,nrefl,rscint,rshield,nflag) !check if it reflected from the detector-guide transition
			if(nflag.eq.1)then !reflected, so move a bit back into detector
				x=x+a*2*epsilon
				y=y+b*2*epsilon
				z=z+c*2*epsilon
				goto 10
			endif
			ngstart=ngstart+1
			if(nwc.eq.1)nwithincrit=nwithincrit+1 !number entering which are 'channeled'
			call d2g(x,y,z,a,b,c,xg,yg,zg,ag,bg,cg,nplane) !switch to guide frame
c			if(ntrack.eq.27400.or.ntrack.eq.32728.or.ntrack.eq.46658)write(6,'(a,i6,12f9.5,i5)')'Bad launch: ',ntrack,x,y,z,a,b,c,xg,yg,zg,ag,bg,cg,nplane
			if(cg.le.0)write(6,'(a,i6,12f9.5,i5)')'Bad launch: ',ntrack,x,y,z,a,b,c,xg,yg,zg,ag,bg,cg,nplane
			call guide(xg,yg,zg,ag,bg,cg,nreturn,nx1,ny1,nplane) !transport through light guide
c				if(cg.gt.0.and.nreturn.eq.1)write(6,*)'cg returned positive: ',xy,yg,zg,ag,by,cg
			if(nreturn.eq.0)then !lost in guide
				nglost=nglost+1
				goto 20
			elseif(nreturn.eq.1)then !transported into and back from guide
				call g2d(xg,yg,zg,ag,bg,cg,x,y,z,a,b,c,nplane) !switch back to detector frame
				if(npmtreturn.eq.0)ngstart=ngstart-1 !this guy got a second chance
				if(nwc.eq.1.and.npmtreturn.eq.0)nwithincrit=nwithincrit-1
				if(npmtreturn.eq.1)then !made it to pmt, but was ALSO reflected
					nphoton(nx1,ny1,nplane)=nphoton(nx1,ny1,nplane)+1
					if(nwc.eq.1)nsc=nsc+1 !was 'channeled' and made it
				endif
				goto 10
			elseif(nreturn.eq.2)then !made it to PMT
				nphoton(nx1,ny1,nplane)=nphoton(nx1,ny1,nplane)+1
				if(nwc.eq.1)nsc=nsc+1 !was 'channeled' and made it
			endif
			goto 20
		endif
	endif
	goto 10
20	enddo
c done with main loop

c prepare PMT hit output

c first add previous results if requested
	if(argadd.eq.'add')then
c		open(unit=4,file='outputs/trans-pmt.dat',err=999)
		open(unit=4,file=argchargefilename,err=999)
		do line=1,1000000000
			read(4,*,end=25)i,j,k,t
			nt=t
			nphoton(i,j,k)=nphoton(i,j,k)+nt
		enddo
25		close(unit=4)
	endif
c save results to file
	nsum=0                  !how many actually reached a detector face
c	open(unit=7,file='outputs/trans-pmt.dat')
	open(unit=7,file=argchargefilename)
	do i=1,ncube
	do j=1,ncube
	do k=1,6
	nt=nphoton(i,j,k)
	nsum=nsum+nt
	if(nt.ne.0)write(7,*)i,j,k,nt
	enddo
	enddo
	enddo
	close(unit=7)

c find number in each vertex-aligned pmt and then sum them
	nvertex(1)=nphoton(nyc,nzc,1)
	nvertex(2)=nphoton(nyc,nzc,2)
	nvertex(3)=nphoton(nxc,nzc,3)
	nvertex(4)=nphoton(nxc,nzc,4)
	nvertex(5)=nphoton(nxc,nyc,5)
	nvertex(6)=nphoton(nxc,nyc,6)
	nchannel=0
	do i=1,6
	nchannel=nchannel+nvertex(i)
	enddo

c calculate percentages lost
	pcts=100.0*nalost/num
	pctf=100.0*nflost/num
	pctt=100.0*ntrapped/num
	pctg=100.0*nglost/num

c output major results

	do nsk=1,6
		nsidesum=0
		do nsi=1,ncube
		do nsj=1,ncube
			nsidesum=nsidesum+nphoton(nsi,nsj,nsk)
		enddo
		enddo
		nssum(nsk)=nsidesum
	enddo
	
	write(6,*)'*************************************************************'
	write(6,*)
	do k=1,6
	write(6,*)
	write(6,"(a,i3,a,i8)")"Side ",k,"   Sum: ",nssum(k)
	write(6,*)
	do j=ncube,1,-1
	write(6,'(3x,i3,4x,65i6)')j,(nphoton(i,j,k),i=1,ncube)
	enddo
	write(6,*)
	write(6,'(10x,65i6)')(i,i=1,ncube)
	enddo

	write(6,*)
	write(6,'(a,3i3,f7.3,4i3,1x,a)')'Command line parameters: ',ncube,nconfig,nmirror,energy,nxc,nyc,nzc,ncombine,argadd

	write(6,*)'---GEOMETRY---'
	write(6,'(a,i2,a,i2,a,i2)')' dimension           ',ncube,' x ',ncube,' x ',ncube
	write(6,'(a,f7.3)')        ' cell dimension      ',celldim
	write(6,'(a,f7.3)')        ' photocathode radius ',pmtr
	write(6,'(a,f7.3)')        ' pnt reflection      ',pmtreflect
	write(6,'(a,f7.3)')        ' t film              ',tfilm
	write(6,*)                  'film layers         ',nfilms
	write(6,'(a,f7.1)')        ' guide length        ',glength
	write(6,'(a,f7.1)')        ' guide transition    ',gztrans
	write(6,'(a,f7.3)')        ' guide VM200 refl    ',agrcoef
	write(6,'(a,f7.3)')        ' guide t film        ',tgfilm
	write(6,'(a,f7.3)')        ' gapwidth            ',gapwidth
	write(6,*)'---FLUID---'
	write(6,*)                  'PC light yield      ',int(pcyield)
	write(6,'(a,f7.3)')        ' q.e.                ',qe
	write(6,'(a,f7.3)')        ' relative yield      ',relativeyield
	write(6,'(a,f7.3)')        ' index (scint)       ',rscint
	write(6,'(a,f7.3)')        ' index (gap)         ',rgap
	write(6,'(a,f7.3)')        ' index (shield)      ',rshield
	write(6,'(a,f7.3)')        ' abs (scint)         ',ascint


	write(6,'(a,f7.3)')        ' abs (film)          ',afilm
	write(6,'(a,f7.3)')        ' abs (guide film)    ',agfilm
	write(6,'(a,f7.1)')        ' abs (shield)        ',agbuffer
	write(6,*)'---EVENT---'
	write(6,*)                  'x start             ',nxc
	write(6,*)                  'y start             ',nyc
	write(6,*)                  'z start             ',nzc
	write(6,'(a,f7.3)')        ' Energy (MeV)        ',energy
	write(6,*)                  'number to track     ',num
	write(6,*)                  'n channeled @ PMTs  ',nwithincrit
	write(6,*)                  'n error             ',nerror
	write(6,*)                  '% abs scint         ',int(pcts)
	write(6,*)                  '% abs film          ',int(pctf)
	write(6,*)                  '% loss det to pmt   ',int(pctg)
	write(6,*)                  '% too many bounces  ',int(pctt)
	write(6,'(a,3i8)')'                                                       unchanneled, channeled, total: ',nsum-nchannel,nchannel,nsum
c output guide information
	if(ngflag.eq.1)then
		write(6,*)'---GUIDE RESULTS---'
		write(6,*)'into guides (-refl): ',ngstart
		write(6,*)'reflected:           ',ngreflect
		write(6,*)'PMT hits:            ',nsuccess
		write(6,*)'film loss:           ',ngflost
		write(6,*)'buffer loss:         ',ngblost
		write(6,*)'trap loss:           ',ngtlost
		write(6,*)'refl loss:           ',ngrlost
		write(6,*)'error loss:          ',ngerror
		write(6,*)'channeled:           ',nwithincrit
		write(6,*)'% ch success:        ',int(100.0*nsc/nwithincrit)
		pctuc=100.0*(nsuccess-nsc)/(ngstart-nwithincrit)
		if(pctuc.ge.0.and.pctuc.le.100)write(6,*)'% unch success:      ',int(pctuc)
	endif

c output pmt timing informationi
c first add previous results if requested
	if(argadd.eq.'add')then
c		open(unit=4,file='outputs/trans-time.dat',err=999)
		open(unit=4,file=argtimefilename,err=999)
		do line=1,1000000000
			read(4,*,end=30)i,j,k,l,nt
			ntp(i,j,k,l)=ntp(i,j,k,l)+nt
		enddo
30		close(unit=4)
	endif
c save results to file
	ntimecheck=0
	open(unit=7,file=argtimefilename)
	do i=1,ncube
	do j=1,ncube
	do k=1,6
	do l=1,ntprint
	nt=ntp(i,j,k,l)
	ntimecheck=ntimecheck+nt
	if(nt.ne.0)write(7,*)i,j,k,l,nt
	enddo
	enddo
	enddo
	enddo
	close(unit=7)
c print timing results
	write(6,*)
	write(6,*)"Timing in source PMT channel (in ",ntbin," ps bins)"
	write(6,*)
	write(6,'("    ps    ",6i10)')(i,i=1,6)
	write(6,*)
	do j=1,ntprint
	write(6,'(5x,i5,6i10)')j*ntbin,ntp(nyc,nzc,1,j),ntp(nyc,nzc,2,j),ntp(nxc,nzc,3,j),ntp(nxc,nzc,4,j),ntp(nxc,nyc,5,j),ntp(nxc,nyc,6,j)
	enddo
	write(6,*)
	write(6,*)"Time integral check: ",ntimecheck
	write(6,*)'*************************************************************'

	stop

999	write(6,*)'previous results do not exist for adding'

	stop
	end

c*******************************************************************************

	function coef(ti,rhigh,rlow)

c Calculates the transmitted fraction for unpolarized light passing through a film
c (thus two sides) at an angle of 'ti'

	coef=0.0
	if(rhigh.ne.rlow)then
		tt=asin(rhigh*sin(ti)/rlow)
		coef=1-(sin(tt-ti)/sin(tt+ti))**2/2+(tan(tt-ti)/tan(tt+ti))**2/2
		coef=coef*coef !two sides of foil
	endif
	return
	end

c*******************************************************************************

	subroutine transport(x,y,z,a,b,c,nplane,theta,shortest)
	include 'inputs/trans-common.f'

c assume unit cells (scale absorption, etc, to appropriate dimensions)
c ray starts at x,y,z; direction cosines a,b,c
c determines which face is hit first, and returns 'shortest' as the distance to that face
c face 1 x=i
c face 2 x=i+1
c face 3 y=j
c face 4 y=j+1
c face 5 z=k
c face 6 z=k+1

c normalize direction cosines
	abc=sqrt(a*a+b*b+c*c)
	a=a/abc
	b=b/abc
	c=c/abc

c find face with the shortest distance along ray
	nplane=0
	shortest=finity ! longer than any path

c check x-planes
	if(a.eq.0)then
	elseif(a.gt.0)then
		xp=int(x+1)
		dist=(xp-x)/a
		if(dist.lt.shortest)then
			nplane=2
			ct=a
			shortest=dist
		endif
	else
		xp=int(x)
		if(x.eq.xp)xp=xp-1
		dist=(xp-x)/a
		if(dist.lt.shortest)then
			nplane=1
			ct=-a
			shortest=dist
		endif
	endif

c check y-planes
	if(b.eq.0)then
	elseif(b.gt.0)then
		yp=int(y+1)
		dist=(yp-y)/b
		if(dist.lt.shortest)then
			nplane=4
			ct=b
			shortest=dist
		endif
	else
		yp=int(y)
		if(y.eq.yp)yp=yp-1
		dist=(yp-y)/b
		if(dist.lt.shortest)then
			nplane=3
			ct=-b
			shortest=dist
		endif
	endif

c check z-planes
	if(c.eq.0)then
	elseif(c.gt.0)then
		zp=int(z+1)
		dist=(zp-z)/c
		if(dist.lt.shortest)then
			nplane=6
			ct=c
			shortest=dist
		endif
	else
		zp=int(z)
		if(z.eq.zp)zp=zp-1
		dist=(zp-z)/c
		if(dist.lt.shortest)then
			nplane=5
			ct=-c
			shortest=dist
		endif
	endif

c returns angle between ray and struck plane normal
	theta=acos(ct)
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine guide(x,y,z,a,b,c,nreturn,nx1,ny1,nplane)
c	incident ray at x,y,z 
c	direction cosines are a,b,c
c	nreturn = 0 -lost in guide
c	nreturn = 1 -reflected back from guide
c	nreturn = 2 -successful transport to pmt
	include 'inputs/trans-common.f'
c	write(6,'(a,6f12.5)')"Entering Guide with: ",x,y,z,a,b,c
	npmtreturn=0
	nreturn=0
	nbounces=0
	iminold=0

10	continue

c	write(6,'(i6,10x,f10.2)')ntrack,z

c find first plane hit and project new location
	d=finity
	imin=0
	do i=1,ngplanes
	top=plane(4,i)*(plane(1,i)-x)+plane(5,i)*(plane(2,i)-y)+plane(6,i)*(plane(3,i)-z)
	bottom=plane(4,i)*a+plane(5,i)*b+plane(6,i)*c
	if(bottom.ne.0)then
		u=top/bottom
		if(u.gt.0.and.u.lt.d)then
			d=u
			imin=i
		endif
	endif
	enddo
	if(iminold.eq.imin.or.imin.eq.0.or.d.eq.0)then
		ngerror=ngerror+1
		return
	endif
c	if(iminold.eq.imin)then
c		write(6,*)
c		write(6,*)"Hit same plane in guide twice in a row: ",imin
c		call reflect(x,y,z,xn,yn,zn,plane(4,imin),plane(5,imin),plane(6,imin),at,bt,ct)
c		write(6,*)x,y,z,xn,yn,zn,plane(4,imin),plane(5,imin),plane(6,imin),at,bt,ct
c		write(6,*)
c	endif
	iminold=imin
c	if(imin.eq.0.or.d.eq.0)then
c		write(6,*)"Did not find a valid intersection with guide plane."
c		do i=1,18
c		write(6,'(i3,6f12.5)')i,(plane(j,i),j=1,6)
c		enddo
c		write(6,*)
c		write(6,'(6f12.5)')x,y,z,a,b,c
c		write(6,*)"Top,Bottom,u,d: ",top,bottom,u,d
c		write(6,*)
c	endif
	dtpath=dtpath+d*celldim*rshield/clight
	xold=x
	yold=y
	zold=z
	x=x+d*a
	y=y+d*b
	z=z+d*c
	zhalf=(glength/celldim)/2.0
	if(
     $		(abs(x).gt.0.5+epsilon).or.
     $		(abs(y).gt.ncombine/2.0+epsilon).or.
     $		(abs(z-zhalf).gt.zhalf+epsilon)
     $	)then
		write(6,'(a,i6,a,i4,a,6f10.6,3f6.2)')"Photon ",ntrack," on bounce ",nbounces," went out of guide.",xold,yold,zold,x,y,z,a,b,c
		ngerror=ngerror+1
		return
	endif
	gpath=gpath+d*celldim

c check if too many bounces
	nbounces=nbounces+1
	if(nbounces.gt.100)then
		ngtlost=ngtlost+1
		return
	endif

c was the photon lost by attenuation in buffer?
	acoef=exp(-d*celldim/agbuffer)
	if(rand().gt.acoef)then
c	write(6,*)d, celldim, agbuffer, acoef,gpath
		ngblost=ngblost+1
		return
	endif

c was the photon lost by attenuation in the guide film?
	cosphoton=a*plane(4,imin)+b*plane(5,imin)+c*plane(6,imin) !for not-normal path through film
	acoef=exp(-2*(tgfilm/cosphoton)*celldim/agfilm)
	if(rand().gt.acoef)then
		ngflost=ngflost+1
		return
	endif

c check if reached pmt plane (ngpmtp)
	if(imin.eq.ngpmtp)then
	
		rfinal=sqrt(x*x+y*y)    ! Check if the photon hits inside the pmt radius
		if(rfinal.gt.pmtr/celldim)then
			write(6,*)"Reached PMT plane, but outside of PMT radius."
			ngerror=ngerror+1
			return
		endif

		reflectioncoeff=0.0     ! Get the reflection coefficient
                if(pmtreflect.ne.0.or.npmtcouple.eq.1)then
                        reflectioncoeff=pmtreflect
                else
                        thetazwy=acos(c)
                        thetacritical=asin(1.0/rshield)
                        if(thetazwy.ge.thetacritical)then
                                reflectioncoeff=1.0
                        else
                                reflectioncoeff=1.0-coef(thetazwy,rshield,1.0)
                        endif
                endif
		
		if(rand().gt.reflectioncoeff)then ! Check the reflection
		        nsuccess=nsuccess+1
		        ndt=max(1.0,min(dtpath/ntbin,1.0*ntprint)) !which time bin was this in
		        ntp(nx1,ny1,nplane,ndt)=ntp(nx1,ny1,nplane,ndt)+1 !log time of arrival for this photon

			nreturn=2
			return
		else
			if(pmtreflect.ne.0)npmtreturn=1
			else npmtreturn=0
			
			c=-c
			x=x+a*epsilon
			y=y+b*epsilon
			z=z+c*epsilon
			goto 10
		endif
	endif

c check if it came back to the entry plane (ngdp)
	if(imin.eq.ngdp)then
c	check for geometry error
		if(
     $			(abs(x).gt.0.5+epsilon).or.
     $			(abs(y).gt.ncombine/2.0+epsilon)
     $		)then
			write(6,'(a,i6,a,6f12.5)')"Photon ",ntrack," returning outside of guide",x,y,z,a,b,c
			ngerror=ngerror+1
			return
		endif
c	refract back into scintillator, or else reflect back into guide
		nflag=0
		if(rshield.ne.rscint)call planerefract(a,b,c,3,rshield,rscint,nflag)
c	back into guide
		if(nflag.eq.1)then
			x=x+a*epsilon
			y=y+b*epsilon
			z=z+c*epsilon
			goto 10
		endif
c 	refracted, and now on into detector, but it immediately hits the guide plane on the side of the detector where it might reflect
		nrflg=0
		if(rgap.ne.rscint)then
			if(acos(abs(c)).gt.tcrit)then
				nrflg=1
				c=-c
			elseif(rand().gt.coef(theta,rscint,rgap))then
				nrflg=1
				c=-c
			endif
		endif
c 	if back into guide
		if(nrflg.eq.1)then
			call planerefract(a,b,c,3,rscint,rshield,nflag)
			if(nflag.eq.1)c=-c !simply ignore yet another reflection
			x=x+a*epsilon
			y=y+b*epsilon
			z=z+c*epsilon
			goto 10
		endif
c	on into detector
		x=x+a*epsilon
		y=y+b*epsilon
		z=z+c*epsilon
		nreturn=1
		if(npmtreflect.eq.0)ngreflect=ngreflect+1
		return
	endif

c check if it actually reflected from the guide sides
	ntirflg=0
	ntaperrgn=0
	if(z.gt.gztrans/celldim)then
		ntaperrgn=1
		tcrit=asin(1.0/rshield) !assuming air gap light-guides
	else
		tcrit=asin(rgap/rshield) !assuming same gap as inside detector for light-guides
	endif
	theta=acos(-(a*plane(4,imin)+b*plane(5,imin)+c*plane(6,imin)))
	if(theta.gt.tcrit)ntirflg=1
	if((ntirflg.eq.0).and.(rand().gt.agrcoef.or.ntaperrgn.eq.0))then !no total reflection, and not reflected by VM2000
		ngrlost=ngrlost+1
		return
	endif

c made it so far, so determine reflection direction and then move epsilon along new direction

	call reflect(xold,yold,zold,x,y,z,plane(4,imin),plane(5,imin),plane(6,imin),a,b,c)
c	write(6,*)xold,yold,zold,x,y,z,plane(4,imin),plane(5,imin),plane(6,imin),a,b,c
	x=x+a*epsilon
	y=y+b*epsilon
	z=z+c*epsilon

c go to next reflection
	goto 10

	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine reflect(x1,y1,z1,x2,y2,z2,u,v,w,ur,vr,wr)
c	Finds direction of reflected ray going from p1 to p2, where the reflecting
c	surface has a unit normal given by u (u,v,w).  Returns ur,vr,wr, leaving from p2.
c	method:
c	projection of (p1-p2) onto u gives symmetry point of reflection s
c	project line from p1 to s, and continue for twice the length to get outgoing point o
c	normalize the vector o - p2 to find ur
c	in summary: ur = [p1+2{p2+[(p1-p2)dot(u)]u-p1}]-p2
c	or reducing: ur = p2-p1+2[(p1-p2)dot(u)]u  and then normalize
	dot=(x1-x2)*u+(y1-y2)*v+(z1-z2)*w
	ur=x2-x1+2*dot*u
	vr=y2-y1+2*dot*v
	wr=z2-z1+2*dot*w
	uvw=sqrt(ur*ur+vr*vr+wr*wr)
	ur=ur/uvw
	vr=vr/uvw
	wr=wr/uvw
	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine getplanes
c	this routine defines a point and unit normal for 16 planes
c	which transition from a rectangular opening to a circular one.
c	Finally, it also defines the PMT plane, and the entry plane
	include 'inputs/trans-common.f'

c determine basic geometry
	rx=0.5
	ry=ncombine/2.0
	rr=pmtr/celldim
	rh=glength/celldim
	rt=gztrans/celldim
	dt=pi/6.0

c define number of planes and special ones
	ngplanes=22
	ngdp=18
	ngpmtp=17

c first, do the triangles which have a side along the circular opening
	do i=0,11
	if(i.le.2)then
		p3x=rx
		p3y=ry
		p3z=rt
	elseif(i.le.5)then
		p3x=-rx
		p3y=ry
		p3z=rt
	elseif(i.le.8)then
		p3x=-rx
		p3y=-ry
		p3z=rt
	elseif(i.le.11)then
		p3x=rx
		p3y=-ry
		p3z=rt
	endif
	p1x=rr*cos(i*dt)
	p1y=rr*sin(i*dt)
	p1z=rh
	p2x=rr*cos((i+1)*dt)
	p2y=rr*sin((i+1)*dt)
	p2z=rh
	
	pax=p2x-p1x
	pay=p2y-p1y
	paz=p2z-p1z
	pbx=p3x-p1x
	pby=p3y-p1y
	pbz=p3z-p1z
	call cross(pax,pay,paz,pbx,pby,pbz,cx,cy,cz)
	xyz=sqrt(cx*cx+cy*cy+cz*cz)
	cx=cx/xyz
	cy=cy/xyz
	cz=cz/xyz
c	write(6,'(6f15.8)')p1x,p1y,p1z,cx,cy,cz
	plane(1,i+1)=p1x
	plane(2,i+1)=p1y
	plane(3,i+1)=p1z
	plane(4,i+1)=cx
	plane(5,i+1)=cy
	plane(6,i+1)=cz
	enddo
c now do the triangles which align with the sides of the rectangular opening
	s1a=rx-rr
	s1b=rh-rt
	s=sqrt(s1a*s1a+s1b*s1b)
	s1a=s1a/s
	s1b=s1b/s

	s2a=ry-rr
	s2b=rh-rt
	s=sqrt(s2a*s2a+s2b*s2b)
	s2a=s2a/s
	s2b=s2b/s

	plane(1,13)=rr
	plane(2,13)=0
	plane(3,13)=rh
	plane(4,13)=-s1b
	plane(5,13)=0
	plane(6,13)=-s1a

	plane(1,14)=-rr
	plane(2,14)=0
	plane(3,14)=rh
	plane(4,14)=s1b
	plane(5,14)=0
	plane(6,14)=-s1a

	plane(1,15)=0
	plane(2,15)=rr
	plane(3,15)=rh
	plane(4,15)=0
	plane(5,15)=-s2b
	plane(6,15)=-s2a

	plane(1,16)=0
	plane(2,16)=-rr
	plane(3,16)=rh
	plane(4,16)=0
	plane(5,16)=s2b
	plane(6,16)=-s2a

c PMT plane
	plane(1,17)=0
	plane(2,17)=0
	plane(3,17)=rh
	plane(4,17)=0
	plane(5,17)=0
	plane(6,17)=-1

c Entry plane
	plane(1,18)=0
	plane(2,18)=0
	plane(3,18)=0
	plane(4,18)=0
	plane(5,18)=0
	plane(6,18)=1

c Rectangular extension planes
	plane(1,19)=rx
	plane(2,19)=0
	plane(3,19)=0
	plane(4,19)=-1
	plane(5,19)=0
	plane(6,19)=0
	
	plane(1,20)=-rx
	plane(2,20)=0
	plane(3,20)=0
	plane(4,20)=1
	plane(5,20)=0
	plane(6,20)=0
	
	plane(1,21)=0
	plane(2,21)=ry
	plane(3,21)=0
	plane(4,21)=0
	plane(5,21)=-1
	plane(6,21)=0
	
	plane(1,22)=0
	plane(2,22)=-ry
	plane(3,22)=0
	plane(4,22)=0
	plane(5,22)=1
	plane(6,22)=0
	
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine cross(ax,ay,az,bx,by,bz,cx,cy,cz)
c	calculate the cross product C=AxB
	cx=ay*bz-az*by
	cy=az*bx-ax*bz
	cz=ax*by-ay*bx
	return
	end
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine d2g(xd,yd,zd,ad,bd,cd,xg,yg,zg,ag,bg,cg,nplane)
	include 'inputs/trans-common.f'

	zg=epsilon
	if(nplane.eq.1)then
		nyoffset=yd/ncombine
		nzoffset=zd
		xg=zd-nzoffset
		yg=yd-nyoffset*ncombine
		ag=cd
		bg=bd
		cg=-ad
	elseif(nplane.eq.2)then
		nyoffset=yd
		nzoffset=zd/ncombine
		xg=yd-nyoffset
		yg=zd-nzoffset*ncombine
		ag=bd
		bg=cd
		cg=ad
	elseif(nplane.eq.3)then
		nxoffset=xd
		nzoffset=zd/ncombine
		xg=xd-nxoffset
		yg=zd-nzoffset*ncombine
		ag=ad
		bg=cd
		cg=-bd
	elseif(nplane.eq.4)then
		nxoffset=xd/ncombine
		nzoffset=zd
		xg=zd-nzoffset
		yg=xd-nxoffset*ncombine
		ag=cd
		bg=ad
		cg=bd
	elseif(nplane.eq.5)then
		nxoffset=xd/ncombine
		nyoffset=yd
		xg=yd-nyoffset
		yg=xd-nxoffset*ncombine
		ag=bd
		bg=ad
		cg=-cd
	elseif(nplane.eq.6)then
		nxoffset=xd
		nyoffset=yd/ncombine
		xg=xd-nxoffset
		yg=yd-nyoffset*ncombine
		ag=ad
		bg=bd
		cg=cd
	endif
	xg=xg-0.5
	yg=yg-ncombine/2.0
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine g2d(xg,yg,zg,ag,bg,cg,xd,yd,zd,ad,bd,cd,nplane)
c	place epsilon inside of detector, since we already handled the
c	outer detector planes in the guide section.

	include 'inputs/trans-common.f'
	xg=xg+0.5
	yg=yg+ncombine/2.0
	if(nplane.eq.1)then
		xd=epsilon
		yd=yg+nyoffset*ncombine
		zd=xg+nzoffset
		ad=-cg
		bd=bg
		cd=ag
	elseif(nplane.eq.2)then
		xd=ncube-epsilon
		yd=xg+nyoffset
		zd=yg+nzoffset*ncombine
		ad=cg
		bd=ag
		cd=bg
	elseif(nplane.eq.3)then
		xd=xg+nxoffset
		yd=epsilon
		zd=yg+nzoffset*ncombine
		ad=ag
		bd=-cg
		cd=bg		
	elseif(nplane.eq.4)then
		xd=yg+nxoffset*ncombine
		yd=ncube-epsilon
		zd=xg+nzoffset
		ad=bg
		bd=cg
		cd=ag		
	elseif(nplane.eq.5)then
		xd=yg+nxoffset*ncombine
		yd=xg+nyoffset
		zd=epsilon
		ad=bg
		bd=ag
		cd=-cg		
	elseif(nplane.eq.6)then
		xd=xg+nxoffset
		yd=yg+nyoffset*ncombine
		zd=ncube-epsilon
		ad=ag
		bd=bg
		cd=cg
	endif
	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine planerefract(u,v,w,nplane,ri,ro,nflag)
c	photon coming with direction cosines u,v,w intersects plane
c	perpindicular to nplane (1,2,3 for x,y,z). Either it gets reflected
c	setting nflag=1, and returns the new u,v,w, or else it gets
c	transmitted and the new directions are returned in u,v,w, with
c	nflag=0.  ri and ro are the incident and outgoing refractive indices.
c	It is assumed that the incident ray will indeed strike the plane,
c	and this is not checked.


	if(ri.eq.ro)then
		nflag=0
		return
	endif
	nflag=0
	su=u/abs(u)
	sv=v/abs(v)
	sw=w/abs(w)
c first rotate so reflection is always from z plane
	if(nplane.eq.1)then
		xi=v
		yi=w
		zi=u
	elseif(nplane.eq.2)then
		xi=w
		yi=u
		zi=v
	elseif(nplane.eq.3)then
		xi=u
		yi=v
		zi=w
	endif
	ai=acos(abs(zi))

c check first for total internal reflection
	if(ri.gt.ro)then
		crit=asin(ro/ri)
		if(ai.gt.crit)then
			nflag=1
			if(nplane.eq.1)then
				u=-u
			elseif(nplane.eq.2)then
				v=-v
			elseif(nplane.eq.3)then
				w=-w
			endif
		return
		endif
	endif

c next check for Fresnel reflection
	ao=asin(ri*sin(ai)/ro)
	if(rand().lt.(sin(ao-ai)/sin(ao+ai))**2/2+(tan(ao-ai)/tan(ao+ai))**2/2)then
		nflag=1
		if(nplane.eq.1)then
			u=-u
		elseif(nplane.eq.2)then
			v=-v
		elseif(nplane.eq.3)then
			w=-w
		endif
		return
	endif

c refract
	phi=atan(abs(yi/xi))
	xo=sin(ao)*cos(phi)
	yo=sin(ao)*sin(phi)
	zo=cos(ao)

c switch back to original frame
	if(nplane.eq.1)then
		u=zo
		v=xo
		w=yo
	elseif(nplane.eq.2)then
		u=yo
		v=zo
		w=xo
	elseif(nplane.eq.3)then
		u=xo
		v=yo
		w=zo
	endif
	uvw=sqrt(u*u+v*v+w*w)
	u=su*abs(u)/uvw
	v=sv*abs(v)/uvw
	w=sw*abs(w)/uvw
	nflag=0

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	function dnormal()
c using the Box-Muller method
	pi=3.141592654
	dnormal=sqrt(-2.0*alog(rand()))*cos(2*pi*rand())
	return
	end
