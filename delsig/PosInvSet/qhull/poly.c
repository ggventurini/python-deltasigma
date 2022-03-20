/* poly.c -- implements polygons and simplices

   see README, poly.h and qhull.h
   
   copyright (c) 1993-1995, The Geometry Center

   infrequent code is in poly2.c 
        (all but top 50 and their callers 12/3/95)
*/

#include "qhull_a.h"

/*======== functions in alphabetical order ==========*/

/*-------------------------------------------------
-appendfacet- appends facet to end of qh facet_list,
  updates qh facet_list, facet_tail, newfacet_list, facet_next
  increments qh numfacets
  assumes qh facet_list/facet_tail is defined (createsimplex)
*/
void qh_appendfacet(facetT *facet) {
  facetT *tail= qh facet_tail;

  if (tail == qh newfacet_list)
    qh newfacet_list= facet;
  if (tail == qh facet_next)
    qh facet_next= facet;
  facet->previous= tail->previous;
  facet->next= tail;
  if (tail->previous)
    tail->previous->next= facet;
  else
    qh facet_list= facet;
  tail->previous= facet;
  qh num_facets++;
  trace4(("qh_appendfacet: append f%d to facet_list\n", facet->id));
} /* appendfacet */


/*-------------------------------------------------
-appendvertex- appends vertex to end of qh vertex_list,
  updates qh vertex_list, vertex_tail, newvertex_list
  increments qh num_vertices
  sets vertex->newlist
  assumes qh vertex_list/vertex_tail is defined (createsimplex)
*/
void qh_appendvertex (vertexT *vertex) {
  vertexT *tail= qh vertex_tail;

  if (tail == qh newvertex_list)
    qh newvertex_list= vertex;
  vertex->newlist= True;
  vertex->previous= tail->previous;
  vertex->next= tail;
  if (tail->previous)
    tail->previous->next= vertex;
  else
    qh vertex_list= vertex;
  tail->previous= vertex;
  qh num_vertices++;
  trace4(("qh_appendvertex: append v%d to vertex_list\n", vertex->id));
} /* appendvertex */


/*-------------------------------------------------
-attachnewfacets- attach horizon facets to new facets in qh newfacet_list
  only needed for qh ONLYgood
  newfacets have neighbor and ridge links to horizon but not vice versa
  will set NEWfacets
returns:
  horizon facets linked to new facets 
     ridges changed from visible facets to new facets
     simplicial ridges deleted
  qh visible_list, no ridges valid
     ->f.replace is a newfacet (if any)
*/
void qh_attachnewfacets (void ) {
  facetT *newfacet= NULL, *neighbor, **neighborp, *horizon, *visible;
  ridgeT *ridge, **ridgep;

  qh NEWfacets= True;
  trace3(("qh_attachnewfacets: delete interior ridges\n"));
  qh visit_id++;
  FORALLvisible_facets {
    visible->visitid= qh visit_id;
    if (visible->ridges) {
      FOREACHridge_(visible->ridges) {
	neighbor= otherfacet_(ridge, visible);
	if (neighbor->visitid == qh visit_id
	    || (!neighbor->visible && neighbor->simplicial)) {
	  if (!neighbor->visible)  /* delete ridge for simplicial horizon */
	    qh_setdel (neighbor->ridges, ridge);
	  qh_setfree (&(ridge->vertices)); /* delete on 2nd visit */
	  qh_memfree (ridge, sizeof(ridgeT));
	}
      }
      SETfirst_(visible->ridges)= NULL;
    }
    SETfirst_(visible->neighbors)= NULL;
  }
  trace1(("qh_attachnewfacets: attach horizon facets to new facets\n"));
  FORALLnew_facets {
    horizon= SETfirst_(newfacet->neighbors);
    if (horizon->simplicial) {
      visible= NULL;
      FOREACHneighbor_(horizon) {   /* may have more than one horizon ridge */
	if (neighbor->visible) {
	  if (visible) {
	    if (qh_setequal_skip (newfacet->vertices, 0, horizon->vertices,
				  SETindex_(horizon->neighbors, neighbor))) {
	      visible= neighbor;
	      break;
	    }
	  }else
	    visible= neighbor;
	}
      }
      if (visible) {
	visible->f.replace= newfacet;
	qh_setreplace (horizon->neighbors, visible, newfacet);
      }else {
	mexPrintf( "qhull internal error (qh_attachnewfacets): couldn't find visible facet for horizon f%d of newfacet f%d\n",
		 horizon->id, newfacet->id);
	qh_errexit2 (qh_ERRqhull, horizon, newfacet);
      }
    }else { /* non-simplicial, with a ridge for newfacet */
      FOREACHneighbor_(horizon) {    /* may hold for many new facets */
	if (neighbor->visible) {
	  neighbor->f.replace= newfacet;
	  qh_setdelnth (horizon->neighbors,
			SETindex_(horizon->neighbors, neighbor));
	  neighborp--; /* repeat */
	}
      }
      qh_setappend (&horizon->neighbors, newfacet);
      ridge= SETfirst_(newfacet->ridges);
      if (ridge->top == horizon)
	ridge->bottom= newfacet;
      else
	ridge->top= newfacet;
      }
  } /* newfacets */
  if (qh PRINTstatistics) {
    FORALLvisible_facets {
      if (!visible->f.replace) 
	zinc_(Zinsidevisible);
    }
  }
} /* attachnewfacets */

/*-------------------------------------------------
-checkflipped- checks facet orientation to interior point
  tests against 0 if !allerror since tested against DISTround before
returns:
  False if flipped orientation (sets facet->flipped)
  distance if non-NULL
*/
boolT qh_checkflipped (facetT *facet, realT *distp, boolT allerror) {
  realT dist;

  if (facet->flipped && !distp)
    return False;
  zzinc_(Zdistcheck);
  qh_distplane(qh interior_point, facet, &dist);
  if (distp)
    *distp= dist;
  if ((allerror && dist > -qh DISTround)|| (!allerror && dist >= 0.0)) {
    facet->flipped= True;
    zzinc_(Zflippedfacets);
    trace0(("qh_checkflipped: facet f%d is flipped, distance= %6.12g during p%d\n",
              facet->id, dist, qh furthest_id));
    return False;
  }
  return True;
} /* checkflipped */

/*-------------------------------------------------
-delfacet- removes facet from facet_list and frees up its memory
   assumes vertices and ridges already freed
*/
void qh_delfacet(facetT *facet) {
  void **freelistp;

  trace5(("qh_delfacet: delete f%d\n", facet->id));
  if (facet == qh tracefacet)
    qh tracefacet= NULL;
  qh_removefacet(facet);
  qh_memfree_(facet->normal, qh normal_size, freelistp);
  if (qh CENTERtype == qh_ASvoronoi) {   /* uses macro calls */
    qh_memfree_(facet->center, qh center_size, freelistp);
  }else /* AScentrum */ {
    qh_memfree_(facet->center, qh normal_size, freelistp);
  }
  qh_setfree(&(facet->neighbors));
  if (facet->ridges)
    qh_setfree(&(facet->ridges));
  qh_setfree(&(facet->vertices));
  if (facet->outsideset)
    qh_setfree(&(facet->outsideset));
  if (facet->coplanarset)
    qh_setfree(&(facet->coplanarset));
  qh_memfree_(facet, sizeof(facetT), freelistp);
} /* delfacet */


/*-------------------------------------------------
-deletevisible- delete visible facets and vertices
    ridges already deleted
    horizon facets do not reference facets on qh visible_list
    new facets in qh newfacet_list
returns:
    deletes each facet and removes from facetlist
    uses qh visit_id;
    qh visible_list empty (== qh newfacet_list)
*/
void qh_deletevisible (/*qh visible_list*/) {
  facetT *visible, *nextfacet;
  vertexT *vertex, **vertexp;
  int numvisible= 0, numdel= qh_setsize(qh del_vertices);

  trace1(("qh_deletevisible: delete %d visible facets and %d vertices\n",
         qh num_visible, numdel));
  for (visible= qh visible_list; visible && visible->visible; 
                visible= nextfacet) { /* deleting current */
    nextfacet= visible->next;        
    numvisible++;
    qh_delfacet(visible);
  }
  if (numvisible != qh num_visible) {
    mexPrintf( "qhull internal error (qh_deletevisible): qh num_visible %d is not number of visible facets %d\n",
             qh num_visible, numvisible);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  qh num_visible= 0;
  zadd_(Zvisfacettot, numvisible);
  zmax_(Zvisfacetmax, numvisible);
  zadd_(Zdelvertextot, numdel);
  zmax_(Zdelvertexmax, numdel);
  FOREACHvertex_(qh del_vertices) 
    qh_delvertex (vertex);
  qh_settruncate (qh del_vertices, 0);
} /* deletevisible */

/*-------------------------------------------------
-facetintersect- return vertices for intersection of two simplicial facets
  may include 1 prepended entry (if more, need to settemppush)
returns:
  returns set of hull_dim-1 + optional extra
  returns skipped index for each test and checks for exactly one
notes:
  does not need settemp since set in quick memory
  see also qh_vertexintersect and qh_vertexintersect_new
  use qh_setnew_delnthsorted to get nth ridge (no skip information)
*/
setT *qh_facetintersect (facetT *facetA, facetT *facetB,
			 int *skipA,int *skipB, int prepend) {
  setT *intersect;
  int dim= qh hull_dim, i, j;
  facetT **neighborsA, **neighborsB;

  neighborsA= SETaddr_(facetA->neighbors, facetT);
  neighborsB= SETaddr_(facetB->neighbors, facetT);
  i= j= 0;
  if (facetB == *neighborsA++)
    *skipA= 0;
  else if (facetB == *neighborsA++)
    *skipA= 1;
  else if (facetB == *neighborsA++)
    *skipA= 2;
  else {
    for (i= 3; i < dim; i++) {
      if (facetB == *neighborsA++) {
        *skipA= i;
        break;
      }
    }
  }
  if (facetA == *neighborsB++)
    *skipB= 0;
  else if (facetA == *neighborsB++)
    *skipB= 1;
  else if (facetA == *neighborsB++)
    *skipB= 2;
  else {
    for (j= 3; j < dim; j++) {
      if (facetA == *neighborsB++) {
        *skipB= j;
        break;
      }
    }
  }
  if (i >= dim || j >= dim) {
    mexPrintf( "qhull internal error (qh_facetintersect): f%d or f%d not in others neighbors\n",
            facetA->id, facetB->id);
    qh_errexit2 (qh_ERRqhull, facetA, facetB);
  }
  intersect= qh_setnew_delnthsorted (facetA->vertices, qh hull_dim, *skipA, prepend);
  trace4(("qh_facetintersect: f%d skip %d matches f%d skip %d\n",
	  facetA->id, *skipA, facetB->id, *skipB));
  return(intersect);
} /* facetintersect */

/*----------------------------------------
-gethash- return hashvalue for a set with firstindex and skipelem
  assumes at least firstindex+1 elements
  sum of elements does badly in high d
  assumes skipelem is NULL, in set, or part of hash
*/
unsigned qh_gethash (int hashsize, setT *set, int size, int firstindex, void *skipelem) {
  void **elemp= SETelemaddr_(set, firstindex, void);
  ptr_intT hash, elem;
  int i;

  switch (size-firstindex) {
  case 1:
    hash= (ptr_intT)(*elemp) - (ptr_intT) skipelem;
    break;
  case 2:
    hash= (ptr_intT)(*elemp) + (ptr_intT)elemp[1] - (ptr_intT) skipelem;
    break;
  case 3:
    hash= (ptr_intT)(*elemp) + (ptr_intT)elemp[1] + (ptr_intT)elemp[2]
      - (ptr_intT) skipelem;
    break;
  case 4:
    hash= (ptr_intT)(*elemp) + (ptr_intT)elemp[1] + (ptr_intT)elemp[2]
      + (ptr_intT)elemp[3] - (ptr_intT) skipelem;
    break;
  case 5:
    hash= (ptr_intT)(*elemp) + (ptr_intT)elemp[1] + (ptr_intT)elemp[2]
      + (ptr_intT)elemp[3] + (ptr_intT)elemp[4] - (ptr_intT) skipelem;
    break;
  case 6:
    hash= (ptr_intT)(*elemp) + (ptr_intT)elemp[1] + (ptr_intT)elemp[2]
      + (ptr_intT)elemp[3] + (ptr_intT)elemp[4]+ (ptr_intT)elemp[5]
      - (ptr_intT) skipelem;
    break;
  default:
    hash= 0;
    i= 3;
    do {     /* this is about 10% in 10-d */
      if ((elem= (ptr_intT)*elemp++) != (ptr_intT)skipelem) {
        hash ^= (elem << i) + (elem >> (32-i));
	i += 3;
	if (i >= 32)
	  i -= 32;
      }
    }while(*elemp);
    break;
  }
  hash %= (ptr_intT) hashsize;
  /* hash= 0;   for debugging purposes */
  return hash;
} /* gethash */

/*-------------------------------------------------
-makenewfacet- creates a toporient? facet from vertices and apex
   modifies vertices 
returns:
    adds newfacet to qh facet_list 
       facet->neighbor= horizon, but not vice versa
    facet->vertices= vertices= apex+vertices
    newvertex_list updated
*/
facetT *qh_makenewfacet(setT *vertices, boolT toporient,facetT *horizon) {
  facetT *newfacet;
  vertexT *vertex, **vertexp;

  FOREACHvertex_(vertices) {
    if (!vertex->newlist) {
      qh_removevertex (vertex);
      qh_appendvertex (vertex);
    }
  }
  newfacet= qh_newfacet();
  newfacet->vertices= vertices;
  newfacet->toporient= toporient;
  qh_setappend(&(newfacet->neighbors), horizon);
  qh_appendfacet(newfacet);
  return(newfacet);
} /* makenewfacet */


/*--------------------------------------------
-makenewplanes- make new hyperplanes for facets
  ->f.samecycle defined for ->mergehorizon facets
returns:
  all facets have hyperplanes or are marked for merging
    doesn't create plane if horizon is coplanar (will merge)
*/
void qh_makenewplanes (void /* newfacet_list */) {
  facetT *newfacet;

  FORALLnew_facets {
    if (!newfacet->mergehorizon)
      qh_setfacetplane (newfacet);  
  }
} /* makenewplanes */

/*---------------------------------------------
-makenew_nonsimplicial- make new facets for ridges of visible facets
  qh visit_id if visible has already been seen
  attaches new facets if !qh ONLY good
  assumes all 'seen' flags false
returns:
  newfacet or NULL, bumps numnew as needed
  marks ridge neighbors for simplicial visible
  if (qh ONLYgood)
    ridges on newfacet, horizon, and visible
  else
    ridge and neighbors between newfacet and horizon
    visible facet's ridges are deleted    
*/
#ifndef qh_NOmerge
facetT *qh_makenew_nonsimplicial (facetT *visible, vertexT *apex, int *numnew) {
  void **freelistp;
  ridgeT *ridge, **ridgep;
  facetT *neighbor, *newfacet= NULL, *samecycle;
  setT *vertices;
  boolT toporient;

  FOREACHridge_(visible->ridges) {
    neighbor= otherfacet_(ridge, visible);
    if (neighbor->visible) {
      if (!qh ONLYgood) {
        if (neighbor->visitid == qh visit_id) {
          qh_setfree (&(ridge->vertices));  /* delete on 2nd visit */
	  qh_memfree_(ridge, sizeof(ridgeT), freelistp);
	}
      }
    }else {  /* neighbor is an horizon facet */
      toporient= (ridge->top == visible);
      vertices= qh_setnew (qh hull_dim); /* makes sure this is quick */
      qh_setappend (&vertices, apex);
      qh_setappend_set (&vertices, ridge->vertices);
      newfacet= qh_makenewfacet(vertices, toporient, neighbor);
      (*numnew)++;
      if (neighbor->coplanar) {
	newfacet->mergehorizon= True;
        if (!neighbor->seen) {
          newfacet->f.samecycle= newfacet;
          neighbor->f.newcycle= newfacet;
        }else {
          samecycle= neighbor->f.newcycle;
          newfacet->f.samecycle= samecycle->f.samecycle;
          samecycle->f.samecycle= newfacet;
	}
      }
      if (qh ONLYgood) {
        if (!neighbor->simplicial)
 	  qh_setappend(&(newfacet->ridges), ridge);
      }else {  /* qh_attachnewfacets */
        if (neighbor->seen) {
	  if (neighbor->simplicial) {
	    mexPrintf( "qhull internal error (qh_makenew_nonsimplicial): simplicial f%d sharing two ridges with f%d\n", 
	           neighbor->id, visible->id);
	    qh_errexit2 (qh_ERRqhull, neighbor, visible);
	  }
	  qh_setappend (&(neighbor->neighbors), newfacet);
	}else
          qh_setreplace (neighbor->neighbors, visible, newfacet);
        if (neighbor->simplicial) {
          qh_setdel (neighbor->ridges, ridge);
          qh_setfree (&(ridge->vertices)); 
	  qh_memfree (ridge, sizeof(ridgeT));
	}else {
 	  qh_setappend(&(newfacet->ridges), ridge);
 	  if (toporient)
 	    ridge->top= newfacet;
 	  else
 	    ridge->bottom= newfacet;
 	}
      trace4(("qh_makenew_nonsimplicial: created facet f%d from v%d and r%d of horizon f%d\n",
	    newfacet->id, apex->id, ridge->id, neighbor->id));
      }
    }
    neighbor->seen= True;        
  } /* for each ridge */
  if (!qh ONLYgood)
    SETfirst_(visible->ridges)= NULL;
  return newfacet;
} /* makenew_nonsimplicial */
#else /* qh_NOmerge */
facetT *qh_makenew_nonsimplicial (facetT *visible, vertexT *apex, int *numnew) {
  return NULL;
}
#endif /* qh_NOmerge */

/*---------------------------------------------
-makenew_simplicial- make new facets for simplicial facet
  nop if neighbor->seen or neighbor->visible (see makenew_nonsimplicial)
  attaches new facets if !qh ONLY good
returns:
  newfacet or NULL, bumps numnew as needed
  if (!qh ONLYgood)
    neighbors between newfacet and horizon
*/
facetT *qh_makenew_simplicial (facetT *visible, vertexT *apex, int *numnew) {
  facetT *neighbor, **neighborp, *newfacet= NULL;
  setT *vertices;
  boolT flip, toporient;
  int horizonskip, visibleskip;

  FOREACHneighbor_(visible) {
    if (!neighbor->seen && !neighbor->visible) {
      vertices= qh_facetintersect(neighbor,visible, &horizonskip, &visibleskip, 1);
      SETfirst_(vertices)= apex;
      flip= ((horizonskip & 0x1) ^ (visibleskip & 0x1));
      if (neighbor->toporient)         
	toporient= horizonskip & 0x1;
      else
	toporient= (horizonskip & 0x1) ^ 0x1;
      newfacet= qh_makenewfacet(vertices, toporient, neighbor);
      (*numnew)++;
      if (neighbor->coplanar && (qh PREmerge || qh MERGEexact)) {
#ifndef qh_NOmerge
	newfacet->f.samecycle= newfacet;
	newfacet->mergehorizon= True;
#endif
      }
      if (!qh ONLYgood)
        SETelem_(neighbor->neighbors, horizonskip)= newfacet;
      trace4(("qh_makenew_simplicial: create facet f%d top %d from v%d and horizon f%d skip %d top %d and visible f%d skip %d, flip? %d\n",
	    newfacet->id, toporient, apex->id, neighbor->id, horizonskip,
	      neighbor->toporient, visible->id, visibleskip, flip));
    }
  }
  return newfacet;
} /* makenew_simplicial */

/*-------------------------------------------------
-matchneighbor- match subridge of newfacet with neighbor or add to hash_table
  ridge is newfacet->vertices w/o newskip vertex
returns:
  duplicate ridges are unmatched and marked by qh_DUPLICATEridge
notes:
  do not allocate memory (need to free hash_table cleanly)
  similar to matchduplicates
  uses linear hash chains
*/
void qh_matchneighbor (facetT *newfacet, int newskip, int hashsize, int *hashcount) {
  boolT newfound= False;   /* True, if new facet is already in hash chain */
  boolT same, ismatch;
  unsigned hash, scan;
  facetT *facet, *matchfacet;
  int skip, matchskip;

  hash= qh_gethash (hashsize, newfacet->vertices, qh hull_dim, 1, 
                     SETelem_(newfacet->vertices, newskip));
  trace4(("qh_matchneighbor: newfacet f%d skip %d hash %d hashcount %d\n",
	  newfacet->id, newskip, hash, *hashcount));
  zinc_(Zhashlookup);
  for (scan= hash; (facet= SETelem_(qh hash_table, scan)); 
       scan= (++scan >= hashsize ? 0 : scan)) {
    if (facet == newfacet) {
      newfound= True;
      continue;
    }
    zinc_(Zhashtests);
    if (qh_matchvertices (1, newfacet->vertices, newskip, facet->vertices, &skip, &same)) {
      if (SETelem_(newfacet->vertices, newskip) == 
          SETelem_(facet->vertices, skip)) {
        mexPrintf( "qhull precision error: Vertex sets are the same for f%d and f%d.  Can not force output.\n",
          facet->id, newfacet->id);
        qh_errexit2 (qh_ERRprec, facet, newfacet);
      }
      ismatch= (same == (newfacet->toporient ^ facet->toporient));
      matchfacet= SETelem_(facet->neighbors, skip);
      if (ismatch && !matchfacet) {
        SETelem_(facet->neighbors, skip)= newfacet;
        SETelem_(newfacet->neighbors, newskip)= facet;
        (*hashcount)--;
        trace4(("qh_matchneighbor: f%d skip %d matched with new f%d skip %d\n",
           facet->id, skip, newfacet->id, newskip));
        return;
      }
      if (!qh PREmerge && !qh MERGEexact) {
	mexPrintf( "qhull precision error: facets f%d, f%d and f%d meet at a ridge with more than 2 neighbors.  Can not continue.\n",
		 facet->id, newfacet->id, getid_(matchfacet));
	qh_errexit2 (qh_ERRprec, facet, newfacet);
      }
      SETelem_(newfacet->neighbors, newskip)= qh_DUPLICATEridge;
      newfacet->dupridge= True;
      if (!newfacet->normal)
	qh_setfacetplane (newfacet);
      qh_addhash (newfacet, qh hash_table, hashsize, hash);
      (*hashcount)++;
      if (!facet->normal)
	qh_setfacetplane (facet);
      if (matchfacet != qh_DUPLICATEridge) {
	SETelem_(facet->neighbors, skip)= qh_DUPLICATEridge;
	facet->dupridge= True;
	if (!facet->normal)
	  qh_setfacetplane (facet);
	if (matchfacet) {
	  matchskip= qh_setindex (matchfacet->neighbors, facet);
	  SETelem_(matchfacet->neighbors, matchskip)= qh_DUPLICATEridge;
	  matchfacet->dupridge= True;
	  if (!matchfacet->normal)
	    qh_setfacetplane (matchfacet);
	  qh_addhash (matchfacet, qh hash_table, hashsize, hash);
	  *hashcount += 2;
	}
      }
      trace4(("qh_matchneighbor: new f%d skip %d duplicates ridge for f%d skip %d matching f%d ismatch %d at hash %d\n",
	   newfacet->id, newskip, facet->id, skip, 
	   (matchfacet == qh_DUPLICATEridge ? -2 : getid_(matchfacet)), 
	   ismatch, hash));
      return; /* end of duplicate ridge */
    }
  }
  if (!newfound) 
    SETelem_(qh hash_table, scan)= newfacet;  /* same as qh_addhash */
  (*hashcount)++;
  trace4(("qh_matchneighbor: no match for f%d skip %d at hash %d\n",
           newfacet->id, newskip, hash));
} /* matchneighbor */


/*-------------------------------------------------
-matchnewfacets- match newfacets in newfacet_list to their newfacet neighbors
  newfacets already have neighbor[0] (horizon facet)
  assumes qh hash_table is NULL
  vertex->neighbors has not been updated yet
returns:
  qh newfacet_list with full neighbor sets
    get vertices with nth neighbor by deleting nth vertex
  if PREmerge/MERGEexact or FORCEoutput 
    all facets check for flipped (also prevents point partitioning)
  if duplicate ridges and PREmerge/MERGEexact
    facet->dupridge set
    missing neighbor links identifies extra ridges to be merging
notes:
  do not allocate memory after hash_table (need to free it cleanly)
*/
void qh_matchnewfacets (void) {
  int numnew=0, hashcount=0, newskip;
  facetT *newfacet, *neighbor;
  int dim= qh hull_dim, hashsize, neighbor_i, neighbor_n;
  setT *neighbors;
#ifndef qh_NOtrace
  int facet_i, facet_n, numfree= 0;
  facetT *facet;
#endif
  
  trace1(("qh_matchnewfacets: match neighbors for new facets.\n"));
  FORALLnew_facets {
    numnew++;
    {  /* inline qh_setzero (newfacet->neighbors, 1, qh hull_dim); */
      neighbors= newfacet->neighbors;
      neighbors->e[neighbors->maxsize].i= dim+1; /*may be overwritten*/
      memset ((char *)SETelemaddr_(neighbors, 1, void), 0, dim * SETelemsize);
    }    
  }
  qh_newhashtable (numnew*(qh hull_dim-1)); /* twice what is normally needed,
                                     but every ridge could be DUPLICATEridge */
  hashsize= qh_setsize (qh hash_table);
  FORALLnew_facets {
    for (newskip=1; newskip<qh hull_dim; newskip++) /* furthest/horizon already matched */
      qh_matchneighbor (newfacet, newskip, hashsize, &hashcount);
#if 0   /* use the following to trap hashcount errors */
    {
      int count= 0, k;
      facetT *facet, *neighbor;

      count= 0;
      FORALLfacet_(qh newfacet_list) {  /* newfacet already in use */
	for (k=1; k<qh hull_dim; k++) {
	  neighbor= SETelem_(facet->neighbors, k);
	  if (!neighbor || neighbor == qh_DUPLICATEridge)
	    count++;
	}
	if (facet == newfacet)
	  break;
      }
      if (count != hashcount) {
	mexPrintf( "qh_matchnewfacets: after adding facet %d, hashcount %d != count %d\n",
		 newfacet->id, hashcount, count);
	qh_errexit (qh_ERRqhull, newfacet, NULL);
      }
    }
#endif  /* end of trap code */
  }
  if (hashcount) {
    FORALLnew_facets {
      if (newfacet->dupridge) {
        FOREACHneighbor_i_(newfacet) {
          if (neighbor == qh_DUPLICATEridge) {
            qh_matchduplicates (newfacet, neighbor_i, hashsize, &hashcount);
         	    /* this may report MERGEfacet */
	  }
        }
      }
    }
  }
  if (hashcount) {
    mexPrintf( "qhull internal error (qh_matchnewfacets): %d neighbors did not match up\n",
        hashcount);
    qh_printhashtable (qh ferr);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
#ifndef qh_NOtrace
  if (qh IStracing >= 2) {
    FOREACHfacet_i_(qh hash_table) {
      if (!facet)
        numfree++;
    }
    mexPrintf( "qh_matchnewfacets: %d new facets, %d unused hash entries .  hashsize %d\n",
	     numnew, numfree, qh_setsize (qh hash_table));
  }
#endif /* !qh_NOtrace */
  qh_setfree (&qh hash_table);
  if (qh PREmerge || qh MERGEexact) {
    if (qh IStracing >= 4)
      qh_printfacetlist (qh newfacet_list, NULL, qh_ALL);
    FORALLnew_facets {
      if (newfacet->normal)
	qh_checkflipped (newfacet, NULL, qh_ALL);
    }
  }else if (qh FORCEoutput)
    qh_checkflipped_all (qh newfacet_list);  /* prints warnings for flipped */
} /* matchnewfacets */

    
/*----------------------------------------
-matchvertices- tests whether vertices match with a single skip
  starts match at firstindex since all new facets have a common vertex
  assumes skipA is in A and both sets are the same size
returns:
  skip index
  sets same iff vertices have the same orientation
*/
boolT qh_matchvertices (int firstindex, setT *verticesA, int skipA, 
       setT *verticesB, int *skipB, boolT *same) {
  vertexT **elemAp, **elemBp, **skipBp=NULL, **skipAp;

  elemAp= SETelemaddr_(verticesA, firstindex, vertexT);
  elemBp= SETelemaddr_(verticesB, firstindex, vertexT);
  skipAp= SETelemaddr_(verticesA, skipA, vertexT);
  do if (elemAp != skipAp) {
    while (*elemAp != *elemBp++) {
      if (skipBp)
        return False;
      skipBp= elemBp;  /* one extra like FOREACH */
    }
  }while(*(++elemAp));
  if (!skipBp)
    skipBp= ++elemBp;
  *skipB= SETindex_(verticesB, skipB);
  *same= !(((ptr_intT)skipA & 0x1) ^ ((ptr_intT)*skipB & 0x1));
  trace4(("qh_matchvertices: matched by skip %d (v%d) and skip %d (v%d) same? %d\n",
	  skipA, (*skipAp)->id, *skipB, (*(skipBp-1))->id, *same));
  return (True);
} /* matchvertices */

/*----------------------------------------
-newfacet- creates and allocates space for a facet
returns:
    all fields initialized or cleared (NULL)
    preallocates neighbors
*/
facetT *qh_newfacet(void) {
  facetT *facet;
  void **freelistp;
  
  qh_memalloc_(sizeof(facetT), freelistp, facet, facetT);
  memset ((char *)facet, 0, sizeof(facetT));
  if (qh facet_id == qh tracefacet_id)
    qh tracefacet= facet;
  facet->id= qh facet_id++;
  facet->neighbors= qh_setnew(qh hull_dim);
#if !qh_COMPUTEfurthest
  facet->furthestdist= 0.0;
#endif
#if qh_MAXoutside
  if (qh FORCEoutput && qh APPROXhull)
    facet->maxoutside= qh MINoutside;
  else
    facet->maxoutside= qh DISTround;
#endif
  facet->simplicial= True;
  facet->good= True;
  facet->newfacet= True;
  trace4(("qh_newfacet: created facet f%d\n", facet->id));
  return (facet);
} /* newfacet */


/*----------------------------------------
-newridge- creates and allocates space for a ridge
*/
ridgeT *qh_newridge(void) {
  ridgeT *ridge;
  void **freelistp;

  qh_memalloc_(sizeof(ridgeT), freelistp, ridge, ridgeT);
  memset ((char *)ridge, 0, sizeof(ridgeT));
  zinc_(Ztotridges);
  if (qh ridge_id == 0xFFFFFF) {
    mexPrintf("\
qhull warning: more than %d ridges.  Id field overflows and two ridges\n\
may have the same identifier.  Otherwise output ok.\n", 0xFFFFFF);
  }
  ridge->id= qh ridge_id++;     
  trace4(("qh_newridge: created ridge r%d\n", ridge->id));
  return (ridge);
} /* newridge */


/*-------------------------------------------------
-pointid- return id for a point, -3 if null, -2 if interior, or -1 if not known
notes:
  alternative code:
  unsigned long id;
  id= ((unsigned long)point - (unsigned long)qh first_point)/qh normal_size;
*/
int qh_pointid (pointT *point) {
  long offset, id;

  if (!point)
    return -3;
  offset= point - qh first_point;
  id= offset / qh hull_dim;
  if (id >= qh num_points || id < 0) {
    if (point == qh interior_point)
      id= -2;
    else if ((id= qh_setindex (qh other_points, point)) != -1)
      id += qh num_points;
  }
  return (int) id;
} /* pointid */
  
/*-------------------------------------------------
-removefacet- unlinks facet from qh facet_list,
updates qh facet_list .newfacet_list .facet_next visible_list

decrements qh num_facets
*/
void qh_removefacet(facetT *facet) {
  facetT *next= facet->next, *previous= facet->previous;
  
  if (facet == qh newfacet_list)
    qh newfacet_list= next;
  if (facet == qh facet_next)
    qh facet_next= next;
  if (facet == qh visible_list)
    qh visible_list= next; 
  if (previous) {
    previous->next= next;
    next->previous= previous;
  }else {  /* 1st facet in qh facet_list */
    qh facet_list= next;
    qh facet_list->previous= NULL;
  }
  qh num_facets--;
  trace4(("qh_removefacet: remove f%d from facet_list\n", facet->id));
} /* removefacet */


/*-------------------------------------------------
-removevertex- unlinks vertex from qh vertex_list,
updates qh vertex_list .newvertex_list 

decrements qh num_vertices
*/
void qh_removevertex(vertexT *vertex) {
  vertexT *next= vertex->next, *previous= vertex->previous;
  
  if (vertex == qh newvertex_list)
    qh newvertex_list= next;
  if (previous) {
    previous->next= next;
    next->previous= previous;
  }else {  /* 1st vertex in qh vertex_list */
    qh vertex_list= vertex->next;
    qh vertex_list->previous= NULL;
  }
  qh num_vertices--;
  trace4(("qh_removevertex: remove v%d from vertex_list\n", vertex->id));
} /* removevertex */


/*-----------------------------------------------
-updatevertices - update vertex neighbors and delete interior vertices
  if qh VERTEXneighbors, update neighbors for each vertex
  interior vertices added to qh del_vertices for later partitioning
*/
void qh_updatevertices (void) {
  facetT *newfacet= NULL, *neighbor, **neighborp, *visible;
  vertexT *vertex, **vertexp;

  trace3(("qh_updatevertices: delete interior vertices and update vertex->neighbors\n"));
  if (qh VERTEXneighbors) {
    FORALLvertex_(qh newvertex_list) {
      FOREACHneighbor_(vertex) {
	if (neighbor->visible) 
	  SETref_(neighbor)= NULL;
      }
      qh_setcompact (vertex->neighbors);
    }
    FORALLnew_facets {
      FOREACHvertex_(newfacet->vertices)
        qh_setappend (&vertex->neighbors, newfacet);
    }
    FORALLvisible_facets {
      FOREACHvertex_(visible->vertices) {
        if (!vertex->newlist && !vertex->deleted) {
  	  FOREACHneighbor_(vertex) { /* this can happen under merging */
	    if (!neighbor->visible)
	      break;
	  }
	  if (neighbor)
	    qh_setdel (vertex->neighbors, visible);
	  else {
	    vertex->deleted= True;
	    qh_setappend (&qh del_vertices, vertex);
	    trace2(("qh_updatevertices: delete vertex p%d (v%d) in f%d\n",
		  qh_pointid(vertex->point), vertex->id, visible->id));
  	  }
        }
      }
    }
  }else {  /* !VERTEXneighbors */
    FORALLvisible_facets {
      FOREACHvertex_(visible->vertices) {
        if (!vertex->newlist && !vertex->deleted) {
          vertex->deleted= True;
	  qh_setappend (&qh del_vertices, vertex);
	  trace2(("qh_updatevertices: delete vertex p%d (v%d) in f%d\n",
		  qh_pointid(vertex->point), vertex->id, visible->id));
  	}
      }
    }
  }
} /* updatevertices */



