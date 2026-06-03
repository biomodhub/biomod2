# Fast version of trigraf assuming that 
#           ie(m) < je(m)
#           ie[is in ascending order
#           je is in ascending order within ie,
# that is, je(ie==i) is in ascending order for each fixed i.
# Code adapted from C code in trigraf.c, from the spatstat package
# by Adrian Baddeley.
#

subroutine trigraf(nv, ne, ie, je, nt, it, jt, kt)
#
# nv --- number of points being triangulated.
# ne --- number of triangle edges
# ie and je --- vectors of indices of ends of each edge
# nt --- number of triangles assumed to be at most ne
# it, jt, kt --- vectors of indices of vertices of triangles
#
integer firstedge, lastedge;
dimension ie(1), je(1), it(1), jt(1), kt(1)

# Initialise output.
nt = 1
lastedge = 0
while(lastedge < ne) {
# Consider next vertex i.
# The edges (i,j) with i < j appear contiguously in the edge list.
 firstedge = lastedge + 1
 i = ie(firstedge)
 do m = firstedge+1,ne {
  if ( ie(m) != i ) break
 }
  lastedge = m-1
# Consider each pair j, k of neighbours of i, where i < j < k. 
# Scan entire edge list to determine whether j, k are joined by an edge.
# If so, save triangle (i,j,k) 
  if(lastedge > firstedge) {
   do mj = firstedge,lastedge-1 {
    j = je(mj)
    do mk = firstedge+1,lastedge {
     k = je(mk)
# Run through edges to determine whether j, k are neighbours.
     do m = 1,ne {
      if(ie(m) >= j) break
     }
      while(m <= ne & ie(m) == j) {
       if(je(m) == k) {
# Add (i, j, k) to list of triangles.
        it(nt) = i;
        jt(nt) = j;
        kt(nt) = k;
        nt = nt+1
       }
      m = m+1
      }
     }
    }
   }
  }

return
end
