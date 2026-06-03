/*

  trigraf.c

  $Revision: 1.2 $     $Date: 2009/02/20 19:24:30 $

  trigraf()  Form list of all triangles in a planar graph, given list of edges

  trigrafS() Faster version when input data are sorted.

*/

void trigraf(nv, ne, ie, je, nt, it, jt, kt, scratch)
     /* inputs */
     int *nv;         /* number of graph vertices */
     int *ne;         /* number of edges */
     int *ie, *je;    /* vectors of indices of ends of each edge */ 
     /* scratch area */
     int *scratch;    /* integer vector, at least 'ne' in length */
     /* output */
     int *nt;              /* number of triangles (assumed <= ne) */
     int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
{
  int Nv, Ne, Nt;
  int Nj, m, i, j, k, mj, mk;
  int *jj;
  
  Nv = *nv;
  Ne = *ne;

  /* initialise storage */
  jj = scratch;
  Nt = 0;

  for(i=0; i < Nv; i++) {
    /* Find triangles involving vertex 'i'
       in which 'i' is the lowest-numbered vertex */

    /* First, find vertices j > i connected to i */
    Nj = 0;
    for(m = 0; m < Ne; m++) {
      if(ie[m] == i) {
	j = je[m];
	if(j > i) {
	  jj[Nj] = j;
	  Nj++;
	}
      } else if(je[m] == i) {
	j = ie[m];
	if(j > i) {
	  jj[Nj] = j;
	  Nj++;
	}
      }
    }

    /* 
       Determine which pairs of vertices j, k are joined by an edge;
       save triangles (i,j,k) 
    */

    if(Nj > 1) {
      /* Sort jj in ascending order */
      for(mj = 0; mj < Nj-1; mj++) {
	j = jj[mj];
	for(mk = mj+1; mk < Nj; mk++) {
	  k = jj[mk];
	  if(k < j) {
	    /* swap */
	    jj[mk] = j;
	    jj[mj] = k;
	    j = k;
	  }
	}
      }
      for(mj = 0; mj < Nj-1; mj++) {
	j = jj[mj];
	for(mk = mj+1; mk < Nj; mk++) {
	  k = jj[mk];
	  if(j != k) {
	    /* Run through edges to determine whether j, k are neighbours */
	    for(m = 0; m < Ne; m++) {
	      if((ie[m] == j && je[m] == k)
		 || (ie[m] == k && je[m] == j)) {
  	        /* add (i, j, k) to list of triangles */
		it[Nt] = i;
		jt[Nt] = j;
		kt[Nt] = k;
		Nt++;
	      }
	    }
	  }
	}
      }
    }
  }

  *nt = Nt;

}


/* faster version of trigraf() 
   assuming that 
            ie[m] < je[m]
            ie[] is in ascending order
            je[] is in ascending order within ie[],
	          that is, je[ie[]=i] is in ascending order for each fixed i
*/

void trigrafS(nv, ne, ie, je, nt, it, jt, kt)
     /* inputs */
     int *nv;         /* number of graph vertices */
     int *ne;         /* number of edges */
     int *ie, *je;    /* vectors of indices of ends of each edge */ 
     /* output */
     int *nt;              /* number of triangles (assumed <= ne) */
     int *it, *jt, *kt;    /* vectors of indices of vertices of triangles */ 
{
  int Nv, Ne, Nt;
  int m, i, j, k, mj, mk;
  int firstedge, lastedge;
  
  Nv = *nv;
  Ne = *ne;

  /* initialise output */
  Nt = 0;

  lastedge = -1;
  while(lastedge + 1 < Ne) {
    /* 
       Consider next vertex i.
       The edges (i,j) with i < j appear contiguously in the edge list.
    */
    firstedge = lastedge + 1;
    i = ie[firstedge]; 
    for(m= firstedge+1; m < Ne && ie[m] == i; m++)
      ;
    lastedge = m-1;
    /* 
       Consider each pair j, k of neighbours of i, where i < j < k. 
       Scan entire edge list to determine whether j, k are joined by an edge.
       If so, save triangle (i,j,k) 
    */
    if(lastedge > firstedge) {
      for(mj = firstedge; mj < lastedge; mj++) {
	j = je[mj];
	for(mk = firstedge+1; mk <= lastedge; mk++) {
	  k = je[mk];
	  /* Run through edges to determine whether j, k are neighbours */
	  for(m = 0; m < Ne && ie[m] < j; m++) 
	    ;
	  while(m < Ne && ie[m] == j) {
	    if(je[m] == k) {
	      /* add (i, j, k) to list of triangles */
	      it[Nt] = i;
	      jt[Nt] = j;
	      kt[Nt] = k;
	      Nt++;
	    }
	    m++;
	  }
	}
      }
    }
  }

  *nt = Nt;

}


  
  
  
