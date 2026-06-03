subroutine collincheck(nadj,madj,npts,x,y,ntot,eps)

# Collinearity check --- experimental.  Runs through the adjacency
# list to see if any of the putative triangles in the triangulation
# that has so far been created are "degenerate", i.e. are actually
# just three points lying on a straight line.

implicit double precision(a-h,o-z)
dimension x(-3:ntot), y(-3:ntot)
dimension nadj(-3:ntot,0:madj)
logical collin, changed

nerror  = -1
changed = .false.
repeat {
    do j = 1,npts {
       nj = nadj(j,0)
       do k = 1,nj {
           k1 = nadj(j,k)
           call succ(k2,j,k1,nadj,madj,ntot,nerror)
           if(nerror > 0) {
               call intpr("Error number =",-1,nerror,1)
               call rexit("Error in succ, called from collincheck.")
           }

# Check whether triangle j, k1, k2 is really a triangle.
           call crossutil(j,k1,k2,x,y,ntot,eps,collin)

# If collinear, remove the triangle from the mix.
           if(collin) {
               changed = .true.
# First determine which of k1 and k2 is closer to j.  It
# *should* be k1, but y'never know in these chaotic
# circumstances.
           sd1 = (x(k1) - x(j))**2 + (y(k1) - y(j))**2
           sd2 = (x(k2) - x(j))**2 + (y(k2) - y(j))**2
           if(sd1 < sd2) {
               kr = k2
           } else {
               kr = k1
           }
# Delete kr ("r" for "remove") from the adjacency list of j and j
# from the adjacency list of kr.
               call delet(j,kr,nadj,madj,ntot,nerror)
               if(nerror > 0) {
                   call intpr("Error number =",-1,nerror,1)
                   call rexit("Error in collincheck.")
               }
               break
           }
       }
       if(changed) break
    }
}
until(!changed)
return
end
