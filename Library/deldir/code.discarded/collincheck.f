C Output from Public domain Ratfor, version 1.03
      subroutine collincheck(nadj,madj,npts,x,y,ntot,eps)
      implicit double precision(a-h,o-z)
      dimension x(-3:ntot), y(-3:ntot)
      dimension nadj(-3:ntot,0:madj)
      logical collin, changed
      nerror = -1
      changed = .false.
23000 continue
      do23003 j = 1,npts 
      nj = nadj(j,0)
      do23005 k = 1,nj 
      k1 = nadj(j,k)
      call succ(k2,j,k1,nadj,madj,ntot,nerror)
      if(nerror .gt. 0)then
      call intpr("Error number =",-1,nerror,1)
      call rexit("Error in succ, called from collincheck.")
      endif
      call crossutil(j,k1,k2,x,y,ntot,eps,collin)
      if(collin)then
      changed = .true.
      sd1 = (x(k1) - x(j))**2 + (y(k1) - y(j))**2
      sd2 = (x(k2) - x(j))**2 + (y(k2) - y(j))**2
      if(sd1 .lt. sd2)then
      kr = k2
      else
      kr = k1
      endif
      call delet(j,kr,nadj,madj,ntot,nerror)
      if(nerror .gt. 0)then
      call intpr("Error number =",-1,nerror,1)
      call rexit("Error in collincheck.")
      endif
      goto 23006
      endif
23005 continue
23006 continue
      if(changed)then
      goto 23004
      endif
23003 continue
23004 continue
23001 if(.not.(.not.changed))goto 23000
23002 continue
      return
      end
