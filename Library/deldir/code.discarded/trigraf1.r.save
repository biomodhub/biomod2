#
# trigraf.r
# Code adapted from C code in trigraf.c, from the spatstat package
# by Adrian Baddeley.
#

subroutine trigraf(nv, ne, ie, je, nt, it, jt, kt, scratch)
#
# nv --- number of points being triangulated.
# ne --- number of triangle edges
# ie and je --- vectors of indices of ends of each edge
# nt --- number of triangles assumed to be at most ne
# it, jt, kt --- vectors of indices of vertices of triangles
# scratch --- integer vector of lenght at least ne.
#
integer scratch(1)
dimension ie(1), je(1), it(1), jt(1), kt(1)
do i = 1,nv {
# Find triangles involving vertex 'i' in which 'i' is the
# lowest-numbered vertex.
# First, find vertices j > i connected to i.
    nj = 1
    do m = 1, ne {
      if(ie[m] == i) {
	j = je[m]
	if(j > i) {
	  jj[nj] = j
          nj = nj+1
	}
      } else if(je[m] == i) {
	j = ie[m];
	if(j > i) {
	  jj[nj] = j
          nj = nj+1
	}
      }
    }

# Determine which pairs of vertices j, k are joined by an edge;
# and save triangles (i,j,k).

    if(nj > 1) {
# Sort jj in ascending order
      do mj = 1,nj {
	j = jj[mj]
        do mk = mj+1,nj {
	  k = jj[mk]
	  if(k < j) {
# Swap.
	    jj[mk] = j
	    jj[mj] = k
	    j = k
	  }
	}
      }
      do mj = 1,nj {
	j = jj[mj]
        do mk = mj+1,nj {
	  k = jj[mk];
	  if(j != k) {
# Run through edges to determine whether j, k are neighbours.
	    for(m = 0; m < Ne; m++) {
            do m = 1,ne {
	      if((ie[m] == j & je[m] == k)
		 | (ie[m] == k & je[m] == j)) {
# Add (i, j, k) to list of triangles.
		it[nt] = i
		jt[nt] = j
		kt[nt] = k
                nt = nt+1
	      }
	    }
	  }
	}
      }
    }
  }
}
