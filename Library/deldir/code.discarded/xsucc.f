C Output from Public domain Ratfor, version 1.03
      subroutine xsucc(ksc,i,j,nadj,madj,ntot,nerror)
      implicit double precision(a-h,o-z)
      dimension nadj(-3:ntot,0:madj)
      dimension junk(20)
      nerror = -1
      n = nadj(i,0)
      if(n.eq.0)then
      nerror = 9
      return
      endif
      do23002 k = 1,n 
      junk(k) = nadj(i,k)
23002 continue
23003 continue
      call intpr("i =",-1,i,1)
      call intpr("adj. list of i:",-1,junk,n)
      do23004 k = 1,n 
      if(j.eq.nadj(i,k))then
      kp = k+1
      if(kp.gt.n)then
      kp = 1
      endif
      ksc = nadj(i,kp)
      call intpr("k =",-1,k,1)
      call intpr("kp =",-1,kp,1)
      call intpr("ksc =",-1,ksc,1)
      call intpr("junk(k) =",-1,junk(k),1)
      call intpr("junk(kp) =",-1,junk(kp),1)
      return
      endif
23004 continue
23005 continue
      nerror = 10
      return
      end
