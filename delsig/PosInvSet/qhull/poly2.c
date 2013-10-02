/* poly2.c -- implements polygons and simplices

   see README, poly.h and qhull.h
   
   copyright (c) 1993-1995, The Geometry Center

   frequently used code is in poly.c
*/

#include "qhull_a.h"

/*======== functions in alphabetical order ==========*/



/*-----------------------------------------------
-addhash- add hash element to linear hash table if not already there
*/
void qh_addhash (void* newelem, setT *hashtable, int hashsize, unsigned hash) {
  unsigned scan;
  void *elem;

  for (scan= hash; (elem= SETelem_(hashtable, scan)); 
       scan= (++scan >= hashsize ? 0 : scan)) {
    if (elem == newelem)
      break;
  }
  /* loop terminates because qh_HASHfactor >= 1.1 by qh_initbuffers */
  if (!elem)
    SETelem_(hashtable, scan)= newelem;
} /* addhash */

/*-----------------------------------------------
-check_bestdist- check that points are within max_outside of the nearest facet
  if ONLYgood, ignores !good facets
  see: check_maxout
*/
void qh_check_bestdist (void) {
  boolT waserror= False, isoutside, unassigned;
  facetT *facet, *bestfacet, *errfacet1= NULL, *errfacet2= NULL;
  facetT *facetlist; 
  realT dist, maxoutside;
  pointT *point;
  int numpart, facet_i, facet_n, notgood= 0, notverified= 0;
  setT *facets;

  maxoutside= fmax_(qh max_outside, qh DISTround);
  maxoutside += 2 * qh DISTround;
  /* 1 DISTround to actual point and another DISTround to computed point */
  if (qh RANDOMdist) /* repeated computations can differ by 2*distround */
    maxoutside += qh DISTround;
  trace1((qh ferr, "qh_check_bestdist: check that all points are within %2.2g of best facet\n", maxoutside));
  facets= qh_pointfacet (/*qh facet_list*/);
  if (!qh_QUICKhelp && qh PRINTprecision)
    fprintf (qh ferr, "\n\
qhull output completed.  Verifying that %d points are\n\
below %2.2g of the nearest %sfacet.\n",
	     qh_setsize(facets), maxoutside, (qh ONLYgood ?  "good " : ""));
  FOREACHfacet_i_(facets) {  /* for each point with facet assignment */
    if (facet)
      unassigned= False;
    else {
      unassigned= True;
      facet= qh facet_list;
    }
    point= qh_point(facet_i);
    if (point == qh GOODpointp)
      continue;
    bestfacet= qh_findbest (point, facet, qh_ALL, False,
			    &dist, &isoutside, &numpart);
    /* occurs after statistics reported */
    if (dist > maxoutside) {
      if (qh ONLYgood && !bestfacet->good 
	  && !((bestfacet= qh_findgooddist (point, bestfacet, &dist, &facetlist))
	       && dist > maxoutside))
	notgood++;
      else {
	waserror= True;
	fprintf(qh ferr, "qhull precision error: point p%d is outside facet f%d, distance= %6.8g maxoutside= %6.8g\n", 
		facet_i, bestfacet->id, dist, maxoutside);
	errfacet2= errfacet1;
	errfacet1= bestfacet;		    
      }
    }else if (unassigned && dist < -qh MAXcoplanar)
      notverified++;
  }
  qh_settempfree (&facets);
  if (notverified && !qh DELAUNAY && !qh_QUICKhelp && qh PRINTprecision) 
    fprintf(qh ferr, "%d points were well inside the hull.  If the hull contains\n\
a lens-shaped component, these points were not verified.  Use\n\
options 'Qi Tv' to verify all points.\n", notverified); 
  if (waserror)
    qh_errexit2 (qh_ERRprec, errfacet1, errfacet2);
} /* check_bestdist */

/*-----------------------------------------------
-check_maxout- updates max_outside by checking all points against bestfacet
  updates facet->maxoutside via findbest
  if printing min_vertex, it is updated to the current vertices
  if ONLYgood, ignores !good facets
  see check_bestdist
notes:
  may not need to check near-inside points if KEEPcoplanar 
     (since coplanar is now MAXcoplanar instead of -DISTround)
*/
#ifndef qh_NOmerge
void qh_check_maxout (void) {
  facetT *facet, *bestfacet, *neighbor, **neighborp, *facetlist;
  realT dist, maxoutside, minvertex;
  pointT *point, **pointp;
  int numpart, facet_i, facet_n, notgood= 0;
  setT *facets, *vertices;
  vertexT *vertex;

  maxoutside= minvertex= 0;
  trace1((qh ferr, "qh_check_maxout: determine actual maxoutside and minoutside\n"));
  if (qh VERTEXneighbors 
  && (qh PRINTsummary || qh PRINTout[0] == qh_PRINTnone || qh PRINTstatistics
      || qh PRINTout[0] == qh_PRINTsummary || qh TRACElevel)) {
    vertices= qh_pointvertex (/*qh facet_list*/);
    FORALLvertices {
      FOREACHneighbor_(vertex) {
	zinc_(Zdistvertex);  /* distance also computed by main loop below */
	qh_distplane (vertex->point, neighbor, &dist);
	minimize_(minvertex, dist);
	if (-dist > qh TRACEdist || dist > qh TRACEdist 
	    || neighbor == qh tracefacet || vertex == qh tracevertex)
	  fprintf (qh ferr, "qh_check_maxout: p%d (v%d) is %.2g from f%d\n",
		   qh_pointid (vertex->point), vertex->id, dist, neighbor->id);
      }
    }
    if (qh MERGING) {
      wmin_(Wminvertex, qh min_vertex);
    }
    qh min_vertex= minvertex;
    qh_settempfree (&vertices);
  }
  facets= qh_pointfacet (/*qh facet_list*/);
  FOREACHfacet_i_(facets) {     /* for each point with facet assignment */
    if (facet) { 
      zinc_(Ztotcheck);
      point= qh_point(facet_i);
      if (point == qh GOODpointp)
	continue;
      bestfacet= qh_findbest (point, facet, qh_ALL,
			    False, &dist, NULL, &numpart);
      zadd_(Zcheckpart, numpart);
      if (bestfacet && dist > maxoutside) {
        if (qh ONLYgood && !bestfacet->good 
        && !((bestfacet= qh_findgooddist (point, bestfacet, &dist, &facetlist))
             && dist > maxoutside))
          notgood++;
        else
	  maxoutside= dist;
      }
      if (dist > qh TRACEdist || (bestfacet && bestfacet == qh tracefacet))
	fprintf (qh ferr, "qh_check_maxout: p%d is %.2g above f%d\n",
		   qh_pointid (point), dist, bestfacet->id);
    }
  }
  qh_settempfree (&facets);
  wval_(Wmaxout)= maxoutside - qh max_outside;
  wmax_(Wmaxoutside, qh max_outside);
  qh max_outside= maxoutside;
  if (qh KEEPnearinside && !qh KEEPcoplanar && !qh KEEPinside) {
    FORALLfacets {
      if (facet->coplanarset) 
        qh_setfree( &facet->coplanarset);
    }
  }else if (qh KEEPnearinside) {
    numpart= 0;
    FORALLfacets {   /* could be combined with qh_findbest */
      if (facet->coplanarset) {
        FOREACHpoint_(facet->coplanarset) {
  	  numpart++;
	  qh_distplane (point, facet, &dist); 
	  if (dist < -qh MAXcoplanar) {
	    if (!qh KEEPinside)
              SETref_(point)= NULL;
          }else if (!qh KEEPcoplanar)
            SETref_(point)= NULL;
        }
	qh_setcompact (facet->coplanarset);
      }
    }
    zadd_(Zcheckpart, numpart);
  }
  trace1((qh ferr, "qh_check_maxout: maxoutside %2.2g, minvertex %2.2g, outside of not good %d\n",
       maxoutside, minvertex, notgood));
} /* check_maxout */
#else /* qh_NOmerge */
void qh_check_maxout (void) {
}
#endif

/*----------------------------------------
-check_output- performs the checks at the end of qhull algorithm
  does not check points (may take a long time)
*/
void qh_check_output (void) {
  int i;

  if (qh STOPcone)
    return;
  if (qh VERIFYoutput | qh IStracing | qh CHECKfrequently) {
    qh_checkpolygon (qh facet_list);
    qh_checkflipped_all (qh facet_list);
    qh_checkconvex (qh facet_list, qh_ALGORITHMfault);
  }else if (!qh MERGING && qh_newstats (qhstat precision, &i)) {
    qh_checkflipped_all (qh facet_list);
    qh_checkconvex (qh facet_list, qh_ALGORITHMfault);
  }
} /* check_output */



/*-------------------------------------------------------------
-check_point- check that point is not outside facet
  if maxerror, doesn't report an error
*/
void qh_check_point (pointT *point, facetT *facet, realT *maxoutside, facetT **errfacet1, facetT **errfacet2) {
  realT dist;

  /* occurs after statistics reported */
  qh_distplane(point, facet, &dist);
  if (dist > *maxoutside) {
    *errfacet2= *errfacet1;
    *errfacet1= facet;
    fprintf(qh ferr, "qhull precision error: point p%d is outside facet f%d, distance= %6.8g maxoutside= %6.8g\n", 
	      qh_pointid(point), facet->id, dist, *maxoutside);
  }
} /* qh_check_point */


/*-------------------------------------------------
-check_points- checks that all points are inside all facets
     uses findbest if lots of points
     ignores flipped facets
notes:
  maxoutside includes 2 DISTrounds, one for the computed
  distances in qh_check_points
  qh_printafacet and qh_printsummary needs only one DISTround
*/
void qh_check_points (void) {
  facetT *facet, *errfacet1= NULL, *errfacet2= NULL;
  realT total, maxoutside;
  pointT *point, **pointp, *pointtemp;
  boolT testouter;

  maxoutside= fmax_(qh max_outside, qh DISTround);  /* agrees with qh_printafacet, qh_printsummary */
  maxoutside += 2* qh DISTround;
  /* 1 DISTround to actual point and another DISTround to computed point */
  if (qh RANDOMdist) /* repeated computations can differ by 2*distround */
    maxoutside += qh DISTround;
  trace1((qh ferr, "qh_check_points: check all points below %2.2g of all facet planes\n",
	  maxoutside));
  if (qh num_good)   /* miss counts other_points and !good facets */
     total= (float) qh num_good * qh num_points;
  else
     total= (float) qh num_facets * qh num_points;
  if (total >= qh_VERIFYdirect 
  && (!qh MERGING || qh SKIPcheckmax || qh ZEROall_ok)) {
    if (!qh_QUICKhelp && qh SKIPcheckmax && qh MERGING)
      fprintf (qh ferr, "\n\
qhull input warning: merging without checking outer planes ('Q5').\n\
Verify may report that a point is outside of a facet.\n");
    qh_check_bestdist();
  }else {
    if (qh MERGING && qh_MAXoutside && !qh SKIPcheckmax)
      testouter= True;  /* agrees with qh_printafacet() */
    else
      testouter= False;
    if (!qh_QUICKhelp) {
      if (qh MERGEexact || qh SKIPcheckmax || qh NOnearinside)
	fprintf (qh ferr, "\n\
qhull input warning: exact merge ('Qx'), no outer plane check ('Q5'), or\n\
no processing of near-inside points ('Q8').  Verify may report that a point\n\
is outside of a facet.\n");
    }
    if (qh PRINTprecision) {
      if (testouter)
	fprintf (qh ferr, "\n\
Output completed.  Verifying that all points are below outer planes of\n\
all %sfacets.  Will make %2.0f distance computations.\n", 
	      (qh ONLYgood ?  "good " : ""), total);
      else
	fprintf (qh ferr, "\n\
Output completed.  Verifying that all points are below %2.2g of\n\
all %sfacets.  Will make %2.0f distance computations.\n", 
	      maxoutside, (qh ONLYgood ?  "good " : ""), total);
    }
    FORALLfacets {
      if (!facet->good && qh ONLYgood)
        continue;
      if (facet->flipped)
        continue;
      if (testouter) {
#if qh_MAXoutside
	maxoutside= facet->maxoutside + 2* qh DISTround;
	/* 1 DISTround to actual point and another to computed point */
#endif
      }
      FORALLpoints {
	if (point != qh GOODpointp)
	  qh_check_point (point, facet, &maxoutside, &errfacet1, &errfacet2);
      }
      FOREACHpoint_(qh other_points) {
	if (point != qh GOODpointp)
	  qh_check_point (point, facet, &maxoutside, &errfacet1, &errfacet2);
      }
    }
    if (errfacet1)
      qh_errexit2(qh_ERRprec, errfacet1, errfacet2);
  }
} /* check_points */


/*-------------------------------------------------
-checkconvex- check that each ridge in facetlist is convex
returns:
    counts Zconcaveridges and Zcoplanarridges
    errors if concaveridge or if merging an coplanar ridge
note:
    if not merging, tests vertices for neighboring simplicial facets
    else if ZEROcentrum, tests vertices for neighboring simplicial facets
    else tests centrums of neighboring facets
*/
void qh_checkconvex(facetT *facetlist, int fault) {
  facetT *facet, *neighbor, **neighborp, *errfacet1=NULL, *errfacet2=NULL;
  vertexT *vertex;
  realT dist;
  pointT *centrum;
  boolT waserror= False, tempcentrum= False, allsimplicial;
  int neighbor_i;

  trace1((qh ferr, "qh_checkconvex: check all ridges are convex\n"));
  zzval_(Zconcaveridges)= 0;
  zzval_(Zcoplanarridges)= 0;
  FORALLfacet_(facetlist) {
    if (facet->flipped) {
      fprintf (qh ferr, "qhull precision error: f%d is flipped (interior point is outside)\n",
	       facet->id);
      errfacet1= facet;
      waserror= True;
      continue;
    }
    if (qh MERGING && (!qh ZEROcentrum || !facet->simplicial))
      allsimplicial= False;
    else {
      allsimplicial= True;
      neighbor_i= 0;
      FOREACHneighbor_(facet) {
        vertex= SETelem_(facet->vertices, neighbor_i++);
	if (!neighbor->simplicial) {
	  allsimplicial= False;
	  continue;
	}
        qh_distplane (vertex->point, neighbor, &dist);
        if (dist > -qh DISTround) {
	  if (fault == qh_DATAfault) {
	    fprintf (qh ferr, "qhull precision error: initial simplex is not convex. Distance=%.2g\n", dist);
	    qh_errexit(qh_ERRsingular, NULL, NULL);
	  }
          if (dist > qh DISTround) {
            zzinc_(Zconcaveridges);
            fprintf (qh ferr, "qhull precision error: f%d is concave to f%d, since p%d (v%d) is %6.4g above\n",
              facet->id, neighbor->id, qh_pointid(vertex->point), vertex->id, dist);
            errfacet1= facet;
            errfacet2= neighbor;
            waserror= True;
          }else if (qh ZEROcentrum) {
            if (dist > 0) {     /* qh_checkzero checks that dist < - qh DISTround */
              zzinc_(Zcoplanarridges); 
              fprintf (qh ferr, "qhull precision error: f%d is clearly not convex to f%d, since p%d (v%d) is %6.4g above\n",
                facet->id, neighbor->id, qh_pointid(vertex->point), vertex->id, dist);
              errfacet1= facet;
              errfacet2= neighbor;
              waserror= True;
	    }
	  }else {              
            zzinc_(Zcoplanarridges);
            trace0((qh ferr, "qhull precision error: f%d may be coplanar to f%d, since p%d (v%d) is within %6.4g during p%d\n",
              facet->id, neighbor->id, qh_pointid(vertex->point), vertex->id, dist, qh furthest_id));
          }
        }
      }
    }
    if (!allsimplicial) {
      if (qh CENTERtype == qh_AScentrum) {
        if (!facet->center)
          facet->center= qh_getcentrum (facet);
        centrum= facet->center;
      }else {
        centrum= qh_getcentrum(facet);
        tempcentrum= True;
      }
      FOREACHneighbor_(facet) {
	if (qh ZEROcentrum && facet->simplicial && neighbor->simplicial)
	  continue;
        zzinc_(Zdistconvex);
        qh_distplane (centrum, neighbor, &dist);
        if (dist > qh DISTround) {
          zzinc_(Zconcaveridges);
          fprintf (qh ferr, "qhull precision error: f%d is concave to f%d.  Centrum of f%d is %6.4g above f%d\n",
            facet->id, neighbor->id, facet->id, dist, neighbor->id);
          errfacet1= facet;
          errfacet2= neighbor;
          waserror= True;
        }else if (dist >= 0.0) {   /* if arithmetic always rounds the same,
				     can test against centrum radius instead */
          zzinc_(Zcoplanarridges);
          fprintf (qh ferr, "qhull precision error: f%d is coplanar or concave to f%d.  Centrum of f%d is %6.4g above f%d\n",
            facet->id, neighbor->id, facet->id, dist, neighbor->id);
	  errfacet1= facet;
	  errfacet2= neighbor;
	  waserror= True;
        }
      }
      if (tempcentrum)
        qh_memfree(centrum, qh normal_size);
    }
  }
  if (waserror && !qh FORCEoutput)
    qh_errexit2 (qh_ERRprec, errfacet1, errfacet2);
} /* checkconvex */


/*-------------------------------------------------
-checkfacet- checks for consistency errors in facet (called from merge.c)
    vertex ids are inverse sorted
    unless newmerge, at least hull_dim neighbors and vertices (exactly if simplicial)
    if non-simplicial, at least as many ridges as neighbors
    neighbors are not duplicated
    ridges are not duplicated
    in 3-d, ridges=verticies
    (hull_dim-1) ridge vertices
    neighbors are reciprocated
    ridge neighbors are facet neighbors and a ridge for every neighbor
    simplicial neighbors match facetintersect
    vertex intersection matches vertices of common ridges 
    vertex neighbors and facet vertices agree
    all ridges have distinct vertex sets
  sets waserror if any error occurs
  uses neighbor->seen
*/
void qh_checkfacet(facetT *facet, boolT newmerge, boolT *waserrorp) {
  facetT *neighbor, **neighborp, *errother=NULL;
  ridgeT *ridge, **ridgep, *errridge= NULL, *ridge2;
  vertexT *vertex, **vertexp;
  unsigned previousid= INT_MAX;
  int numneighbors, numvertices, numridges=0, numRvertices=0;
  boolT waserror= False;
  int skipA, skipB, ridge_i, ridge_n, i;
  setT *intersection;

  if (facet->visible) {
    fprintf (qh ferr, "qhull internal error (qh_checkfacet): facet f%d is on the visible_list\n",
      facet->id);
    qh_errexit (qh_ERRqhull, facet, NULL);
  }
  if (!facet->normal) {
    fprintf (qh ferr, "qhull internal error (qh_checkfacet): facet f%d does not have  a normal\n",
      facet->id);
    waserror= True;
  }
  qh_setcheck (facet->vertices, "vertices for f", facet->id);
  qh_setcheck (facet->ridges, "ridges for f", facet->id);
  qh_setcheck (facet->outsideset, "outsideset for f", facet->id);
  qh_setcheck (facet->coplanarset, "coplanarset for f", facet->id);
  qh_setcheck (facet->neighbors, "neighbors for f", facet->id);
  FOREACHvertex_(facet->vertices) {
    if (vertex->deleted) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): deleted vertex v%d in f%d\n", vertex->id, facet->id);
      qh_errprint ("ERRONEOUS", NULL, NULL, NULL, vertex);
      waserror= True;
    }
    if (vertex->id >= previousid) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): vertices of f%d are not in descending id order at v%d\n", facet->id, vertex->id);
      waserror= True;
      break;
    }
    previousid= vertex->id;
  }
  numneighbors= qh_setsize(facet->neighbors);
  numvertices= qh_setsize(facet->vertices);
  numridges= qh_setsize(facet->ridges);
  if (facet->simplicial) {
    if (numvertices+numneighbors != 2*qh hull_dim 
    && !facet->degenerate && !facet->redundant) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): for simplicial facet f%d, #vertices %d + #neighbors %d != 2*qh hull_dim\n", 
                facet->id, numvertices, numneighbors);
      qh_setprint (qh ferr, "", facet->neighbors);
      waserror= True;
    }
  }else { /* non-simplicial */
    if (!newmerge 
    &&(numvertices < qh hull_dim || numneighbors < qh hull_dim)
    && !facet->degenerate && !facet->redundant) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): for facet f%d, #vertices %d or #neighbors %d < qh hull_dim\n",
         facet->id, numvertices, numneighbors);
       waserror= True;
    }
    if (numridges < numneighbors
    ||(qh hull_dim == 3 && numvertices != numridges && !qh NEWfacets)
    ||(qh hull_dim == 2 && numridges + numvertices + numneighbors != 6)) {
      if (!facet->degenerate && !facet->redundant) {
	fprintf(qh ferr, "qhull internal error (qh_checkfacet): for facet f%d, #ridges %d < #neighbors %d or (3-d) != #vertices %d or (2-d) not all 2\n",
	    facet->id, numridges, numneighbors, numvertices);
	waserror= True;
      }
    }
  }
  FOREACHneighbor_(facet) {
    if (neighbor == qh_MERGEridge || neighbor == qh_DUPLICATEridge) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): facet f%d still has a MERGE or DUP neighbor\n", facet->id);
      qh_errexit (qh_ERRqhull, facet, NULL);
    }
    neighbor->seen= True;
  }
  FOREACHneighbor_(facet) {
    if (!qh_setin(neighbor->neighbors, facet)) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): facet f%d has neighbor f%d, but f%d does not have neighbor f%d\n",
	      facet->id, neighbor->id, neighbor->id, facet->id);
      errother= neighbor;
      waserror= True;
    }
    if (!neighbor->seen) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): facet f%d has a duplicate neighbor f%d\n",
	      facet->id, neighbor->id);
      errother= neighbor;
      waserror= True;
    }    
    neighbor->seen= False;
  }
  FOREACHridge_(facet->ridges) {
    qh_setcheck (ridge->vertices, "vertices for r", ridge->id);
    ridge->seen= False;
  }
  FOREACHridge_(facet->ridges) {
    if (ridge->seen) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): facet f%d has a duplicate ridge r%d\n",
	      facet->id, ridge->id);
      errridge= ridge;
      waserror= True;
    }    
    ridge->seen= True;
    numRvertices= qh_setsize(ridge->vertices);
    if (numRvertices != qh hull_dim - 1) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): ridge between f%d and f%d has %d vertices\n", 
                ridge->top->id, ridge->bottom->id, numRvertices);
      errridge= ridge;
      waserror= True;
    }
    neighbor= otherfacet_(ridge, facet);
    neighbor->seen= True;
    if (!qh_setin(facet->neighbors, neighbor)) {
      fprintf(qh ferr, "qhull internal error (qh_checkfacet): for facet f%d, neighbor f%d of ridge r%d not in facet\n",
           facet->id, neighbor->id, ridge->id);
      errridge= ridge;
      waserror= True;
    }
  }
  if (!facet->simplicial) {
    FOREACHneighbor_(facet) {
      if (!neighbor->seen) {
        fprintf(qh ferr, "qhull internal error (qh_checkfacet): facet f%d does not have a ridge for neighbor f%d\n",
	      facet->id, neighbor->id);
	errother= neighbor;
        waserror= True;
      }
      intersection= qh_vertexintersect_new(facet->vertices, neighbor->vertices);
      qh_settemppush (intersection);
      FOREACHvertex_(facet->vertices) {
	vertex->seen= False;
	vertex->seen2= False;
      }
      FOREACHvertex_(intersection)
	vertex->seen= True;
      FOREACHridge_(facet->ridges) {
	if (neighbor != otherfacet_(ridge, facet))
	    continue;
	FOREACHvertex_(ridge->vertices) {
	  if (!vertex->seen) {
	    fprintf (qh ferr, "qhull internal error (qh_checkfacet): vertex v%d in r%d not in f%d intersect f%d\n",
  	          vertex->id, ridge->id, facet->id, neighbor->id);
	    qh_errexit (qh_ERRqhull, facet, ridge);
	  }
	  vertex->seen2= True;
	}
      }
      if (!newmerge) {
	FOREACHvertex_(intersection) {
	  if (!vertex->seen2) {
	    if (qh IStracing >=3 || !qh MERGING) {
	      fprintf (qh ferr, "qhull precision error (qh_checkfacet): vertex v%d in f%d intersect f%d but\n\
 not in a ridge.  This is ok under merging.  Last point was p%d\n",
		     vertex->id, facet->id, neighbor->id, qh furthest_id);
	      if (!qh FORCEoutput && !qh MERGING) {
		qh_errprint ("ERRONEOUS", facet, neighbor, NULL, vertex);
		if (!qh MERGING)
		  qh_errexit (qh_ERRqhull, NULL, NULL);
	      }
	    }
	  }
	}
      }      
      qh_settempfree (&intersection);
    }
  }else { /* simplicial */
    FOREACHneighbor_(facet) {
      if (neighbor->simplicial) {    
	skipA= SETindex_(facet->neighbors, neighbor);
	skipB= qh_setindex (neighbor->neighbors, facet);
	if (!qh_setequal_skip (facet->vertices, skipA, neighbor->vertices, skipB)) {
	  fprintf (qh ferr, "qhull internal error (qh_checkfacet): facet f%d skip %d and neighbor f%d skip %d do not match \n",
		   facet->id, skipA, neighbor->id, skipB);
	  errother= neighbor;
	  waserror= True;
	}
      }
    }
  }
  if (qh hull_dim < 5 && (qh IStracing > 2 || qh CHECKfrequently)) {
    FOREACHridge_i_(facet->ridges) {           /* expensive */
      for (i= ridge_i+1; i < ridge_n; i++) {
	ridge2= SETelem_(facet->ridges, i);
	if (qh_setequal (ridge->vertices, ridge2->vertices)) {
	  fprintf (qh ferr, "qh_checkfacet: ridges r%d and r%d have the same vertices\n",
		  ridge->id, ridge2->id);
	  errridge= ridge;
	  waserror= True;
	}
      }
    }
  }
  if (waserror) {
    qh_errprint("ERRONEOUS", facet, errother, errridge, NULL);
    *waserrorp= True;
  }
} /* checkfacet */


/*-------------------------------------------------
-checkflipped_all- checks orientation of facets in list against interior point
*/
void qh_checkflipped_all (facetT *facetlist) {
  facetT *facet;
  boolT waserror= False;
  realT dist;

  if (facetlist == qh facet_list)
    zzval_(Zflippedfacets)= 0;
  FORALLfacet_(facetlist) {
    if (facet->normal && !qh_checkflipped (facet, &dist, !qh_ALL)) {
      fprintf(qh ferr, "qhull precision error: facet f%d is flipped, distance= %6.12g\n",
	      facet->id, dist);
      if (!qh FORCEoutput) {
	qh_errprint("ERRONEOUS", facet, NULL, NULL, NULL);
	waserror= True;
      }
    }
  }
  if (waserror) {
    fprintf (qh ferr, "\n\
A flipped facet occurs when its distance to the interior point is\n\
greater than %2.2g, the maximum round-off error.\n", -qh DISTround);
    qh_errexit(qh_ERRprec, NULL, NULL);
  }
} /* checkflipped_all */

/*-------------------------------------------------
-checkpolygon- checks the correctness of the structure
  check num_facets and num_vertices if qh facet_list
  call with either qh facet_list or qh newfacet_list
*/
void qh_checkpolygon(facetT *facetlist) {
  facetT *facet;
  vertexT *vertex, **vertexp, *vertexlist;
  int numfacets= 0, numvertices= 0, numridges= 0;
  int totvneighbors= 0, totvertices= 0;
  boolT waserror= False, nextseen= False, visibleseen= False;
  
  trace1((qh ferr, "qh_checkpolygon: check all facets from f%d\n", facetlist->id));
  if (facetlist != qh facet_list || qh ONLYgood)
    nextseen= True;
  FORALLfacet_(facetlist) {
    if (facet == qh visible_list)
      visibleseen= True;
    if (!facet->visible) {
      if (!nextseen) {
	if (facet == qh facet_next)
	  nextseen= True;
	else if (qh_setsize (facet->outsideset)) {
	  fprintf (qh ferr, "qhull internal error (qh_checkpolygon): f%d has outside set before qh facet_next\n",
		   facet->id);
	  qh_errexit (qh_ERRqhull, facet, NULL);
	}
      }
      numfacets++;
      qh_checkfacet(facet, False, &waserror);
    }
  }
  if (qh visible_list && !visibleseen && facetlist == qh facet_list) {
    fprintf (qh ferr, "qhull internal error (qh_checkpolygon): visible list f%d no longer on facet list\n", qh visible_list->id);
    qh_printlists();
    qh_errexit (qh_ERRqhull, qh visible_list, NULL);
  }
  if (facetlist == qh facet_list)
    vertexlist= qh vertex_list;
  else if (facetlist == qh newfacet_list)
    vertexlist= qh newvertex_list;
  else
    vertexlist= NULL;
  FORALLvertex_(vertexlist) {
    vertex->seen= False;
    vertex->visitid= 0;
  }  
  FORALLfacet_(facetlist) {
    if (facet->visible)
      continue;
    if (facet->simplicial)
      numridges += qh hull_dim;
    else
      numridges += qh_setsize (facet->ridges);
    FOREACHvertex_(facet->vertices) {
      vertex->visitid++;
      if (!vertex->seen) {
	vertex->seen= True;
	numvertices++;
	if (qh_pointid (vertex->point) == -1) {
	  fprintf (qh ferr, "qhull internal error (qh_checkpolygon): unknown point %p for vertex v%d first_point %p\n",
		   vertex->point, vertex->id, qh first_point);
	  waserror= True;
	}
      }
    }
  }
  qh vertex_visit += numfacets;
  if (facetlist == qh facet_list) {
    if (numfacets != qh num_facets - qh num_visible) {
      fprintf(qh ferr, "qhull internal error (qh_checkpolygon): actual number of facets is %d, cumulative facet count is %d\n",
	      numfacets, qh num_facets- qh num_visible);
      waserror= True;
    }
    qh vertex_visit++;
    if (qh VERTEXneighbors) {
      FORALLvertices {
	qh_setcheck (vertex->neighbors, "neighbors for v", vertex->id);
	if (vertex->deleted)
	  continue;
	totvneighbors += qh_setsize (vertex->neighbors);
      }
      FORALLfacet_(facetlist)
	totvertices += qh_setsize (facet->vertices);
      if (totvneighbors != totvertices) {
	fprintf(qh ferr, "qhull internal error (qh_checkpolygon): vertex neighbors inconsistent.  Totvneighbors %d, totvertices %d\n",
		totvneighbors, totvertices);
	waserror= True;
      }
    }
    if (numvertices != qh num_vertices - qh_setsize(qh del_vertices)) {
      fprintf(qh ferr, "qhull internal error (qh_checkpolygon): actual number of vertices is %d, cumulative vertex count is %d\n",
	      numvertices, qh num_vertices - qh_setsize(qh del_vertices));
      waserror= True;
    }
    if (qh hull_dim == 2 && numvertices != numfacets) {
      fprintf (qh ferr, "qhull internal error (qh_checkpolygon): #vertices %d != #facets %d\n",
        numvertices, numfacets);
      waserror= True;
    }
    if (qh hull_dim == 3 && numvertices + numfacets - numridges/2 != 2) {
      fprintf (qh ferr, "qhull internal error (qh_checkpolygon): #vertices %d + #facets %d - #edges %d != 2\n",
        numvertices, numfacets, numridges/2);
      waserror= True;
    }
  }
  if (waserror) 
    qh_errexit(qh_ERRqhull, NULL, NULL);
} /* checkpolygon */


/*-------------------------------------------------
-checkvertex- check vertex for consistency
notes:
  neighbors checked efficiently in checkpolygon
*/
void qh_checkvertex (vertexT *vertex) {
  boolT waserror= False;
  facetT *neighbor, **neighborp, *errfacet=NULL;
  int size;

  if (qh_pointid (vertex->point) == -1) {
    fprintf (qh ferr, "qhull internal error (qh_checkvertex): unknown point id %p\n", vertex->point);
    waserror= True;
  }
  if (vertex->id >= qh vertex_id) {
    fprintf (qh ferr, "qhull internal error (qh_checkvertex): unknown vertex id %d\n", vertex->id);
    waserror= True;
  }
  if (!waserror && !vertex->deleted) {
    if ((size= qh_setsize (vertex->neighbors))) {
      FOREACHneighbor_(vertex) {
        if (!qh_setin (neighbor->vertices, vertex)) {
          fprintf (qh ferr, "qhull internal error (qh_checkvertex): neighbor f%d does not contain v%d\n", neighbor->id, vertex->id);
	  errfacet= neighbor;
	  waserror= True;
	}
      }
    }
  }
  if (waserror) {
    qh_errprint ("ERRONEOUS", NULL, NULL, NULL, vertex);
    qh_errexit (qh_ERRqhull, errfacet, NULL);
  }
} /* checkvertex */
  
/*-------------------------------------------------
-clearcenters- clear old data from facet->center
  sets new centertype
  nop if CENTERtype is the same
*/
void qh_clearcenters (qh_CENTER type) {
  facetT *facet;
  
  if (qh CENTERtype != type) {
    FORALLfacets {
      if (qh CENTERtype == qh_ASvoronoi){
        if (facet->center) {
          qh_memfree (facet->center, qh center_size);
          facet->center= NULL;
        }
      }else /* qh CENTERtype == qh_AScentrum */ {
        if (facet->center) {
          qh_memfree (facet->center, qh normal_size);
	  facet->center= NULL;
        }
      }
    }
    qh CENTERtype= type;
  }
  trace2((qh ferr, "clearcenters: switched to center type %d\n", type));
} /* clearcenters */

/*----------------------------------------
-createsimplex- creates a simplex from a set of vertices
returns:
    initializes qh facet_list to the simplex
*/
void qh_createsimplex(setT *vertices) {
  facetT *facet= NULL, *newfacet;
  boolT toporient= True;
  int vertex_i, vertex_n, nth;
  setT *newfacets= qh_settemp (qh hull_dim+1);
  vertexT *vertex;
  
  qh facet_list= qh newfacet_list= qh facet_tail= qh_newfacet();
  qh vertex_list= qh newvertex_list= qh vertex_tail= qh_newvertex(NULL);
  FOREACHvertex_i_(vertices) {
    newfacet= qh_newfacet();
    newfacet->vertices= qh_setnew_delnthsorted (vertices, vertex_n,
						vertex_i, 0);
    newfacet->toporient= toporient;
    qh_appendfacet(newfacet);
    newfacet->newfacet= True;
    qh_appendvertex (vertex);
    qh_setappend (&newfacets, newfacet);
    toporient ^= True;
  }
  FORALLnew_facets {
    nth= 0;
    FORALLfacet_(qh newfacet_list) {
      if (facet != newfacet) 
        SETelem_(newfacet->neighbors, nth++)= facet;
    }
    qh_settruncate (newfacet->neighbors, qh hull_dim);
  }
  qh_settempfree (&newfacets);
  trace1((qh ferr, "qh_createsimplex: created simplex\n"));
} /* createsimplex */

/*-------------------------------------------------
-delridge- deletes ridge from data structures it belongs to and frees up the
memory occupied by it
notes:
  in merge.c, caller sets vertex->delridge for each vertex
  also freed in qh_freeqhull
*/
void qh_delridge(ridgeT *ridge) {
  void **freelistp;
  
  qh_setdel(ridge->top->ridges, ridge);
  qh_setdel(ridge->bottom->ridges, ridge);
  qh_setfree(&(ridge->vertices));
  qh_memfree_(ridge, sizeof(ridgeT), freelistp);
} /* delridge */


/*-------------------------------------------------
-delvertex- deletes a vertex and frees its memory
  assumes vertex->adjacencies have been updated if needed
  unlinks for vertex_list
*/
void qh_delvertex (vertexT *vertex) {

  if (vertex == qh tracevertex)
    qh tracevertex= NULL;
  qh_removevertex (vertex);
  qh_setfree (&vertex->neighbors);
  qh_memfree(vertex, sizeof(vertexT));
} /* delvertex */


/*----------------------------------------
-facet3vertex- return temporary set of 3-d vertices
  in qh_ORIENTclock order
*/
setT *qh_facet3vertex (facetT *facet) {
  ridgeT *ridge, *firstridge;
  vertexT *vertex;
  int cntvertices, cntprojected=0;
  setT *vertices;

  cntvertices= qh_setsize(facet->vertices);
  vertices= qh_settemp (cntvertices);
  if (facet->simplicial) {
    if (cntvertices != 3) {
      fprintf (qh ferr, "qhull internal error (qh_facet3vertex): only %d vertices for simplicial facet f%d\n", 
                  cntvertices, facet->id);
      qh_errexit(qh_ERRqhull, facet, NULL);
    }
    qh_setappend (&vertices, SETfirst_(facet->vertices));
    if (facet->toporient ^ qh_ORIENTclock)
      qh_setappend (&vertices, SETsecond_(facet->vertices));
    else
      qh_setaddnth (&vertices, 0, SETsecond_(facet->vertices));
    qh_setappend (&vertices, SETelem_(facet->vertices, 2));
  }else {
    ridge= firstridge= SETfirst_(facet->ridges);   /* no infinite */
    while ((ridge= qh_nextridge3d (ridge, facet, &vertex))) {
      qh_setappend (&vertices, vertex);
      if (++cntprojected > cntvertices || ridge == firstridge)
        break;
    }
    if (!ridge || cntprojected != cntvertices) {
      fprintf (qh ferr, "qhull internal error (qh_facet3vertex): ridges for facet %d don't match up.  got at least %d\n", 
                  facet->id, cntprojected);
      qh_errexit(qh_ERRqhull, facet, ridge);
    }
  }
  return vertices;
} /* facet3vertex */


/*-------------------------------------------------
-findfacet- find facet that is furthest below a point 
  if 'facet' defined, starts search at 'facet'
  if NULL 'facet', starts at first non-flipped non-upperDelaunay facet
  searches neighbors of coplanar and most flipped facets
  does not search upper envelope of Delaunay triangulations
returns:
  best facet
  distance to facet
  isoutside if point is outside of the hull
  numpart counts the number of distance tests
notes:
  this works for Delaunay triangulations.
  this works if the initial facet is not above the point.
  warning: a lens-shaped hull may mistakenly return a facet
      that is above the point
  uses qh visit_id, qh searchset
*/
facetT *qh_findfacet (pointT *point, facetT *facet, 
           realT *dist, boolT *isoutside, int *numpart) {
  facetT *bestfacet;
  
  if (!facet)
    facet= qh_nonupper( qh facet_list);
  bestfacet= qh_findbest (point, facet, qh_ALL, False, dist, isoutside, numpart);
  return bestfacet;
} /* findfacet */ 
 
/*-------------------------------------------------
-findgood- identify good facets for qh ONLYgood
  GOODvertex>0 - facet includes point as vertex
    if !match, returns goodhorizon
    inactive if qh MERGING
  GOODpoint - facet is visible or coplanar (>0) or not visible (<0) 
  GOODthreshold - facet->normal matches threshold
    if !goodhorizon and !match, selects facet with closest angle
       and sets GOODclosest
returns:
  number of new, good facets found
  determins facet->good
  may update GOODclosest
notes:
  findgood_all further reduces the good region
*/
int qh_findgood (facetT *facetlist, int goodhorizon) {
  facetT *facet, *bestfacet;
  realT angle, bestangle, dist;
  int  numgood=0;

  if (qh GOODclosest) {
    bestfacet= qh GOODclosest;
    qh_inthresholds (bestfacet->normal, &bestangle);
  }else {
    bestfacet= NULL;
    bestangle= -REALmax;
  }
  FORALLfacet_(facetlist) {
    if (facet->good)
      numgood++;
  }
  if (qh GOODvertex>0 && !qh MERGING) {
    FORALLfacet_(facetlist) {
      if (!qh_isvertex (qh GOODvertexp, facet->vertices)) {
        facet->good= False;
        numgood--;
      }
    }
  }
  if (qh GOODpoint && numgood) {
    FORALLfacet_(facetlist) {
      if (facet->good && facet->normal) {
        zinc_(Zdistgood);
        qh_distplane (qh GOODpointp, facet, &dist);
        if ((qh GOODpoint > 0) ^ (dist > 0.0)) {
          facet->good= False;
          numgood--;
        }
      }
    }
  }
  if (qh GOODthreshold && (numgood || goodhorizon)) {
    FORALLfacet_(facetlist) {
      if (facet->good && facet->normal) {
        if (!qh_inthresholds (facet->normal, &angle)) {
          facet->good= False;
          numgood--;
          angle= fabs_(angle);
          if (angle > bestangle) {
            bestangle= angle;
            bestfacet= facet;
          }
        }
      }
    }
    if (!numgood && bestfacet && bestfacet != qh GOODclosest) {
      if (qh GOODclosest)
	qh GOODclosest->good= False;
      qh GOODclosest= bestfacet;
      bestfacet->good= True;
      numgood++;
      trace2((qh ferr, "qh_findgood: f%d is closest (%2.2g) to thresholds\n", 
           bestfacet->id, bestangle));
      return numgood;
    }else if (numgood && qh GOODclosest)
      qh GOODclosest->good= False;
  }
  zadd_(Zgoodfacet, numgood);
  trace2((qh ferr, "qh_findgood: found %d good facets\n", numgood));
  if (!numgood && qh GOODvertex>0 && !qh MERGING) 
    return goodhorizon;
  return numgood;
} /* findgood */

/*-------------------------------------------------
-findgood_all- apply other constraints for good facets (used by qh PRINTgood)
  GOODvertex - facet includes (>0) or doesn't include (<0) point as vertex
    if last good facet, prints warning and continues
  SPLITthresholds- facet->normal matches threshold, or if none, the closest one
  calls findgood if !ONLYgood
  nop if good not used
returns:
  clears facet->good if not good
  sets qh num_good
notes:
  this is like findgood but more restrictive
*/
void qh_findgood_all (facetT *facetlist) {
  facetT *facet, *bestfacet=NULL;
  realT angle, bestangle= REALmax;
  int  numgood=0, startgood;

  if (!qh GOODvertex && !qh GOODthreshold && !qh GOODpoint 
  && !qh SPLITthresholds)
    return;
  if (!qh ONLYgood)
    qh_findgood (qh facet_list, 0);
  FORALLfacet_(facetlist) {
    if (facet->good)
      numgood++;
  }
  if (qh GOODvertex <0 || (qh GOODvertex > 0 && qh MERGING)) {
    FORALLfacet_(facetlist) {
      if (facet->good && ((qh GOODvertex > 0) ^ !!qh_isvertex (qh GOODvertexp, facet->vertices))) {
        if (!--numgood) {
          fprintf (qh ferr, "qhull warning: good vertex p%d does not match last good facet f%d.  Ignored.\n",
             qh_pointid(qh GOODvertexp), facet->id);
          return;
        }
        facet->good= False;
      }
    }
  }
  startgood= numgood;
  if (qh SPLITthresholds) {
    FORALLfacet_(facetlist) {
      if (facet->good) {
        if (!qh_inthresholds (facet->normal, &angle)) {
          facet->good= False;
          numgood--;
          angle= fabs_(angle);
          if (angle < bestangle) {
            bestangle= angle;
            bestfacet= facet;
          }
        }
      }
    }
    if (!numgood) {
      bestfacet->good= True;
      numgood++;
      trace0((qh ferr, "qh_findgood_all: f%d is closest (%2.2g) to thresholds\n", 
           bestfacet->id, bestangle));
      return;
    }
  }
  qh num_good= numgood;
  trace0((qh ferr, "qh_findgood_all: %d good facets remain out of %d facets\n",
        numgood, startgood));
} /* findgood_all */

/*-------------------------------------------------
-infiniteloop- report infinite loop error due to facet
*/
void qh_infiniteloop (facetT *facet) {

  fprintf (qh ferr, "qhull internal error (qh_infiniteloop): potential infinite loop detected\n");
  qh_errexit (qh_ERRqhull, facet, NULL);
} /* qh_infiniteloop */

/*--------------------------------------------------
-initialhull- constructs the initial hull as a qh hull_dim simplex of vertices
*/
void qh_initialhull(setT *vertices) {
  facetT *facet, *firstfacet;
  realT dist;
#ifndef qh_NOtrace
  int k;
#endif

  qh_createsimplex(vertices);  /* qh facet_list */
  qh interior_point= qh_getcenter(vertices);
  firstfacet= qh facet_list;
  qh_setfacetplane(firstfacet);
  zinc_(Znumvisibility); /* needs to be in printsummary */
  qh_distplane(qh interior_point, firstfacet, &dist);
  if (dist > 0) {  
    FORALLfacets
      facet->toporient ^= True;
  }
  FORALLfacets
    qh_setfacetplane(facet);
  FORALLfacets {
    if (!qh_checkflipped (facet, NULL, qh_ALL)) {/* due to axis-parallel facet */
      trace1((qh ferr, "qh_initialhull: initial orientation incorrect.  Correct all facets\n"));
      facet->flipped= False;
      FORALLfacets {
	facet->toporient ^= True;
	qh_orientoutside (facet);
      }
      break;
    }
  }
  FORALLfacets {
    if (!qh_checkflipped (facet, NULL, !qh_ALL)) {  /* can happen with 'R0.1' */
      fprintf (qh ferr, "qhull precision error: facet %d is coplanar with the interior point\n",
                   facet->id);
      qh_errexit (qh_ERRsingular, facet, NULL);
    }
  }
  zzval_(Zprocessed)= qh hull_dim+1;
  qh_checkpolygon (qh facet_list);
  qh_checkconvex(qh facet_list,   qh_DATAfault);
#ifndef qh_NOtrace
  if (qh IStracing >= 1) {
    fprintf(qh ferr, "qh_initialhull: simplex constructed, interior point:");
    for (k=0; k<qh hull_dim; k++) 
      fprintf (qh ferr, " %6.4g", qh interior_point[k]);
    fprintf (qh ferr, "\n");
  }
#endif
} /* initialhull */

/*-------------------------------------------------
-initialvertices- determines a non-singular set of initial vertices
  maxpoints are not unique
returns:
  temporary set of dim+1 vertices in descending order by vertex id
notes:
  unless qh ALLpoints, uses maxpoints as long as determinate is non-zero
  picks random points if qh RANDOMoutside && !ALLpoints
  if dim >= qh_INITIALmax, uses min/max x and max points with non-zero determinants
*/
setT *qh_initialvertices(int dim, setT *maxpoints, pointT *points, int numpoints) {
  pointT *point, **pointp;
  setT *vertices, *simplex, *tested;
  realT randr, det;
  int index, point_i, point_n, k;
  boolT nearzero= False;
  
  vertices= qh_settemp (dim + 1);
  simplex= qh_settemp (dim+1);
  if (qh ALLpoints) 
    qh_maxsimplex (dim, NULL, points, numpoints, &simplex);
  else if (qh RANDOMoutside) {
    while (qh_setsize (simplex) != dim+1) {
      randr= qh_RANDOMint;
      randr= randr/(qh_RANDOMmax+1);
      index= floor(qh num_points * randr);
      point= qh_point (index);
      qh_setunique (&simplex, point);
    }
  }else if (qh hull_dim >= qh_INITIALmax) {
    tested= qh_settemp (dim+1);
    qh_setappend (&simplex, SETfirst_(maxpoints));   /* max and min X coord */
    qh_setappend (&simplex, SETsecond_(maxpoints));
    qh_maxsimplex (fmin_(qh_INITIALsearch, dim), maxpoints, points, numpoints, &simplex);
    k= qh_setsize (simplex);
    FOREACHpoint_i_(maxpoints) { 
      if (point_i & 0x1) {     /* first pick up max. coord. points */
      	if (!qh_setin (simplex, point) && !qh_setin (tested, point)){
	  det= qh_detsimplex(point, simplex, k, &nearzero);
          if (nearzero)
            qh_setappend (&tested, point);
          else {
            qh_setappend (&simplex, point);
            if (++k == dim)  /* use search for last point */
	      break;
	  }
	}
      }
    }
    while (k != dim && (point= (pointT*)qh_setdellast (maxpoints))) {
      if (!qh_setin (simplex, point) && !qh_setin (tested, point)){
        det= qh_detsimplex (point, simplex, k, &nearzero);
        if (nearzero)
          qh_setappend (&tested, point);
        else {
          qh_setappend (&simplex, point);
          k++;
	}
      }
    }
    index= 0;
    while (k != dim && (point= qh_point (index++))) {
      if (!qh_setin (simplex, point) && !qh_setin (tested, point)){
        det= qh_detsimplex (point, simplex, k, &nearzero);
        if (!nearzero){
          qh_setappend (&simplex, point);
          k++;
	}
      }
    }
    qh_settempfree (&tested);
    qh_maxsimplex (dim, maxpoints, points, numpoints, &simplex);
  }else
    qh_maxsimplex (dim, maxpoints, points, numpoints, &simplex);
  FOREACHpoint_(simplex) 
    qh_setaddnth (&vertices, 0, qh_newvertex(point)); /* descending order */
  qh_settempfree (&simplex);
  return vertices;
} /* initialvertices */


/*-------------------------------------------------
-isvertex- returns vertex if point is in vertex set, else returns NULL
*/
vertexT *qh_isvertex (pointT *point, setT *vertices) {
  vertexT *vertex, **vertexp;

  FOREACHvertex_(vertices) {
    if (vertex->point == point)
      return vertex;
  }
  return NULL;
} /* isvertex */

/*-------------------------------------------------
-makenewfacets- make new facets from point and qh visible_list
returns:
  qh newfacet_list= list of new facets with hyperplanes and ->newfacet
  qh newvertex_list= list of vertices in new facets with ->newlist set
  if (qh ONLYgood)
    newfacets reference horizon facets, but not vice versa
    ridges reference non-simplicial horizon ridges, but not vice versa
    does not change existing facets
  otherwise
    NEWfacets set
    newfacets attached to horizon facets and ridges
    visible->r.replace is corresponding new facet
*/
vertexT *qh_makenewfacets (pointT *point /*visible_list*/) {
  facetT *visible, *newfacet= NULL, *newfacet2= NULL, *neighbor, **neighborp;
  vertexT *apex;
  int numnew=0;

  qh newfacet_list= qh facet_tail;
  qh newvertex_list= qh vertex_tail;
  apex= qh_newvertex(point);
  qh_appendvertex (apex);  
  qh visit_id++;
  if (!qh ONLYgood)
    qh NEWfacets= True;
  FORALLvisible_facets {
    FOREACHneighbor_(visible) 
      neighbor->seen= False;
    if (visible->ridges) {
      visible->visitid= qh visit_id;
      newfacet2= qh_makenew_nonsimplicial (visible, apex, &numnew);
    }
    if (visible->simplicial)
      newfacet= qh_makenew_simplicial (visible, apex, &numnew);
    if (!qh ONLYgood) {
      if (newfacet2)  /* newfacet is null if all ridges defined */
        newfacet= newfacet2;
      if (newfacet)
      	visible->f.replace= newfacet;
      else
        zinc_(Zinsidevisible);
      SETfirst_(visible->neighbors)= NULL;
    }
  }
  trace1((qh ferr, "qh_makenewfacets: created %d new facets from point p%d to horizon\n",
	  numnew, qh_pointid(point)));
  if (qh IStracing >= 4)
    qh_printfacetlist (qh newfacet_list, NULL, qh_ALL);
  return apex;
} /* makenewfacets */

/*-------------------------------------------------
-matchduplicates- match duplicate ridges in hashtable
  marked with ->dupridge and qh_DUPLICATEridge
  picks match with worst merge (min distance apart)
returns:
  updates hashcount
  similar to qh_matchneighbor
*/
#ifndef qh_NOmerge
void qh_matchduplicates (facetT *atfacet, int atskip, int hashsize, int *hashcount) {
  boolT same, ismatch;
  unsigned hash, scan;
  facetT *facet, *newfacet, *maxmatch= NULL, *maxmatch2= NULL, *nextfacet;
  int skip, newskip, nextskip= 0, maxskip= 0, maxskip2= 0, makematch;
  realT maxdist= -REALmax, mindist, dist2, low, high;

  hash= qh_gethash (hashsize, atfacet->vertices, qh hull_dim, 1, 
                     SETelem_(atfacet->vertices, atskip));
  trace2((qh ferr, "qh_matchduplicates: find duplicate matches for f%d skip %d hash %d hashcount %d\n",
	  atfacet->id, atskip, hash, *hashcount));
  for (makematch= 0; makematch < 2; makematch++) {
    qh visit_id++;
    for (newfacet= atfacet, newskip= atskip; newfacet; newfacet= nextfacet, newskip= nextskip) {
      zinc_(Zhashlookup);
      nextfacet= NULL;
      newfacet->visitid= qh visit_id;
      for (scan= hash; (facet= SETelem_(qh hash_table, scan)); 
	   scan= (++scan >= hashsize ? 0 : scan)) {
	if (!facet->dupridge || facet->visitid == qh visit_id)
	  continue;
	zinc_(Zhashtests);
	if (qh_matchvertices (1, newfacet->vertices, newskip, facet->vertices, &skip, &same)) {
	  ismatch= (same == (newfacet->toporient ^ facet->toporient));
	  if (SETelem_(facet->neighbors, skip) != qh_DUPLICATEridge) {
	    if (!makematch) {
	      fprintf (qh ferr, "qhull internal error (qh_matchduplicates): missing dupridge at f%d skip %d for new f%d skip %d hash %d\n",
		     facet->id, skip, newfacet->id, newskip, hash);
	      qh_errexit2 (qh_ERRqhull, facet, newfacet);
	    }
	  }else if (ismatch && makematch) {
	    if (SETelem_(newfacet->neighbors, newskip) == qh_DUPLICATEridge) {
	      SETelem_(facet->neighbors, skip)= newfacet;
	      SETelem_(newfacet->neighbors, newskip)= qh_MERGEridge;
	      *hashcount -= 2; /* removed two unmatched facets */
	      trace4((qh ferr, "qh_matchduplicates: duplicate f%d skip %d matched with new f%d skip %d merge\n",
		    facet->id, skip, newfacet->id, newskip));
	    }
	  }else if (ismatch) {
	    mindist= qh_getdistance (facet, newfacet, &low, &high);
	    dist2= qh_getdistance (newfacet, facet, &low, &high);
	    minimize_(mindist, dist2);
	    if (mindist > maxdist) {
	      maxdist= mindist;
	      maxmatch= facet;
	      maxskip= skip;
	      maxmatch2= newfacet;
	      maxskip2= newskip;
	    }
	    trace3((qh ferr, "qh_matchduplicates: duplicate f%d skip %d new f%d skip %d at dist %2.2g, max is now f%d f%d\n",
		    facet->id, skip, newfacet->id, newskip, mindist, 
		    maxmatch->id, maxmatch2->id));
	  }else { /* !ismatch */
	    nextfacet= facet;
	    nextskip= skip;
	  }
	}
	if (makematch && !facet 
        && SETelem_(facet->neighbors, skip) == qh_DUPLICATEridge) {
	  fprintf (qh ferr, "qhull internal error (qh_matchduplicates): no MERGEridge match for duplicate f%d skip %d at hash %d\n",
		     newfacet->id, newskip, hash);
	  qh_errexit (qh_ERRqhull, newfacet, NULL);
	}
      }
    } /* end of for each new facet at hash */
    if (!makematch) {
      if (!maxmatch) {
	fprintf (qh ferr, "qhull internal error (qh_matchduplicates): no maximum match at duplicate f%d skip %d at hash %d\n",
		     atfacet->id, atskip, hash);
	qh_errexit (qh_ERRqhull, atfacet, NULL);
      }
      SETelem_(maxmatch->neighbors, maxskip)= maxmatch2;
      SETelem_(maxmatch2->neighbors, maxskip2)= maxmatch;
      *hashcount -= 2; /* removed two unmatched facets */
      zzinc_(Zmultiridge);
      trace1((qh ferr, "qh_matchduplicates: duplicate f%d skip %d matched with new f%d skip %d keep\n",
	      maxmatch->id, maxskip, maxmatch2->id, maxskip2));
      if (qh IStracing >= 4)
	qh_errprint ("DUPLICATED/MATCH", maxmatch, maxmatch2, NULL, NULL);
    }
  }
} /* matchduplicates */

/*----------------------------------------
-nearvertex- return nearest vertex to point
returns:
  vertex and distance
*/
vertexT *qh_nearvertex (facetT *facet, pointT *point, realT *bestdistp) {
  realT bestdist= REALmax, dist;
  vertexT *bestvertex= NULL, *vertex, **vertexp;

  FOREACHvertex_(facet->vertices) {
    dist= qh_pointdist (vertex->point, point, -qh hull_dim);
    if (dist < bestdist) {
      bestdist= dist;
      bestvertex= vertex;
    }
  }
  *bestdistp= sqrt (bestdist);
  return bestvertex;
} /* nearvertex */

/*-------------------------------------------------
-newhashtable- returns size of qh hash_table of at least newsize slots
  assumes qh hash_table is NULL
  qh_HASHfactor determines the number of extra slots
*/
int qh_newhashtable(int newsize) {
  int size;

  size= ((newsize+1)*qh_HASHfactor) | 0x1;  /* odd number */
  while (True) { 
    if ((size%3) && (size%5))
      break;
    size += 2;
    /* loop terminates because there is an infinite number of primes */
  }
  qh hash_table= qh_setnew (size);
  qh_setzero (qh hash_table, 0, size);
  return size;
} /* newhashtable */

/*----------------------------------------
-newvertex- creates and allocates space for a vertex
*/
vertexT *qh_newvertex(pointT *point) {
  vertexT *vertex;

  zinc_(Ztotvertices);
  vertex= (vertexT *)qh_memalloc(sizeof(vertexT));
  memset ((char *) vertex, 0, sizeof (vertexT));
  if (qh vertex_id == 0xFFFFFF) {
    fprintf(qh ferr, "qhull input error: more than %d vertices.  Id field overflows and two vertices\n\
may have the same identifier.  Vertices not sorted correctly.\n", 0xFFFFFF);
    qh_errexit(qh_ERRinput, NULL, NULL);
  }
  if (qh vertex_id == qh tracevertex_id)
    qh tracevertex= vertex;
  vertex->id= qh vertex_id++;
  vertex->point= point;
  trace4((qh ferr, "qh_newvertex: vertex p%d (v%d) created\n", qh_pointid(vertex->point), 
	  vertex->id));
  return (vertex);
} /* newvertex */

/*----------------------------------------
-nextridge3d- return next ridge and vertex for a 3d facet
  in qh_ORIENTclock order
  n^2 implementation to trace all ridges
  be sure to stop on any 2nd visit
*/
ridgeT *qh_nextridge3d (ridgeT *atridge, facetT *facet, vertexT **vertexp) {
  vertexT *atvertex, *vertex, *othervertex;
  ridgeT *ridge, **ridgep;

  if ((atridge->top == facet) ^ qh_ORIENTclock)
    atvertex= SETsecond_(atridge->vertices);
  else
    atvertex= SETfirst_(atridge->vertices);
  FOREACHridge_(facet->ridges) {
    if (ridge == atridge)
      continue;
    if ((ridge->top == facet) ^ qh_ORIENTclock) {
      othervertex= SETsecond_(ridge->vertices);
      vertex= SETfirst_(ridge->vertices);
    }else {
      vertex= SETsecond_(ridge->vertices);
      othervertex= SETfirst_(ridge->vertices);
    }
    if (vertex == atvertex) {
      if (vertexp)
        *vertexp= othervertex;
      return ridge;
    }
  }
  return NULL;
} /* nextridge3d */
#else /* qh_NOmerge */
void qh_matchduplicates (facetT *atfacet, int atskip, int hashsize, int *hashcount) {
}
ridgeT *qh_nextridge3d (ridgeT *atridge, facetT *facet, vertexT **vertexp) {

  return NULL;
}
#endif /* qh_NOmerge */
  
/*----------------------------------------
-nonupper- return first facet without upperdelaunay or flipped
*/
facetT *qh_nonupper (facetT *facetlist) {
  facetT *facet;

  FORALLfacet_(facetlist) {
    if (!facet->upperdelaunay && !facet->flipped)
      return facet;
  }
  fprintf(qh ferr, "qhull errror (qh_nonupper): all facets from f%d are upper-Delaunay or flipped\n",
    getid_(facetlist)); 
  qh_errexit( qh_ERRqhull, facetlist, NULL);
  return NULL; /* avoid warning */
} /* nonupper */

/*-------------------------------------------------
-point- return point for a point id, or NULL if unknown
alternative code:
    return ((pointT *)((unsigned long)qh first_point
           + (unsigned long)((id)*qh normal_size)));
*/
pointT *qh_point (int id) {

  if (id < 0)
    return NULL;
  if (id < qh num_points)
    return qh first_point + id * qh hull_dim;
  id -= qh num_points;
  if (id < qh_setsize (qh other_points))
    return SETelem_(qh other_points, id);
  return NULL;
} /* point */
  
/*-------------------------------------------------
-point_add- access function for pointfacet and pointvertex
*/
void qh_point_add (setT *set, pointT *point, void *elem) {
  int id, size;

  SETreturnsize_(set, size);
  if ((id= qh_pointid(point)) < 0)
    fprintf (qh ferr, "qhull internal warning (point_add): unknown point %p id %d\n", 
      point, id);
  else if (id >= size) {
    fprintf (qh ferr, "qhull internal errror (point_add): point p%d is out of bounds (%d)\n",
	     id, size);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }else
    SETelem_(set, id)= elem;
} /* point_add */


/*-------------------------------------------------
-pointfacet- return temporary set of facets indexed by point id
  for vertices, coplanarset, and outsideset
  access with FOREACHfacet_i_(facets) and SETelem_(facets, i)
  NULL if no facet for point (inside)
     this will include qh GOODpointp
*/
setT *qh_pointfacet (void /*qh facet_list*/) {
  int numpoints= qh num_points + qh_setsize (qh other_points);
  setT *facets;
  facetT *facet;
  vertexT *vertex, **vertexp;
  pointT *point, **pointp;
  
  facets= qh_settemp (numpoints);
  qh_setzero (facets, 0, numpoints);
  qh vertex_visit++;
  FORALLfacets {
    FOREACHvertex_(facet->vertices) {
      if (vertex->visitid != qh vertex_visit) {
        vertex->visitid= qh vertex_visit;
        qh_point_add (facets, vertex->point, facet);
      }
    }
    FOREACHpoint_(facet->coplanarset) 
      qh_point_add (facets, point, facet);
    FOREACHpoint_(facet->outsideset) 
      qh_point_add (facets, point, facet);
  }
  return facets;
} /* pointfacet */

/*-------------------------------------------------
-pointvertex- return temporary set of vertices indexed by point id
  access with FOREACHvertex_i_(vertices) and SETelem_(vertices, i)
  NULL if no vertex for point
     this will include qh GOODpointp
*/
setT *qh_pointvertex (void /*qh facet_list*/) {
  int numpoints= qh num_points + qh_setsize (qh other_points);
  setT *vertices;
  vertexT *vertex;
  
  vertices= qh_settemp (numpoints);
  qh_setzero (vertices, 0, numpoints);
  FORALLvertices 
    qh_point_add (vertices, vertex->point, vertex);
  return vertices;
} /* pointvertex */


/*-------------------------------------------------
-prependfacet- prepend facet to the start of a facetlist
  increments qh numfacets
  updates facetlist, qh facet_list, facet_next
notes:
  be careful of prepending since it can lose a pointer.
  e.g., can lose _next by deleting and then prepending before _next
*/
void qh_prependfacet(facetT *facet, facetT **facetlist) {
  facetT *prevfacet, *list= *facetlist;
  

  trace4((qh ferr, "qh_prependfacet: prepend f%d before f%d\n",
	  facet->id, list->id));
  prevfacet= list->previous;
  facet->previous= prevfacet;
  if (prevfacet)
    prevfacet->next= facet;
  list->previous= facet;
  facet->next= *facetlist;
  if (qh facet_list == list)  /* this may change *facetlist */
    qh facet_list= facet;
  if (qh facet_next == list)
    qh facet_next= facet;
  *facetlist= facet;
  qh num_facets++;
} /* prependfacet */


/*-----------------------------------------
-printhashtable- print hash table
  not in I/O to avoid bringing io.c in
*/
void qh_printhashtable(FILE *fp) {
  facetT *facet, *neighbor;
  int id, facet_i, facet_n, neighbor_i= 0, neighbor_n= 0;
  vertexT *vertex, **vertexp;

  FOREACHfacet_i_(qh hash_table) {
    if (facet) {
      FOREACHneighbor_i_(facet) {
        if (!neighbor || neighbor == qh_MERGEridge || neighbor == qh_DUPLICATEridge) 
          break;
      }
      if (neighbor_i == neighbor_n)
        continue;
      fprintf (fp, "hash %d f%d ", facet_i, facet->id);
      FOREACHvertex_(facet->vertices)
        fprintf (fp, "v%d ", vertex->id);
      fprintf (fp, "\n neighbors:");
      FOREACHneighbor_i_(facet) {
	if (neighbor == qh_MERGEridge)
	  id= -3;
	else if (neighbor == qh_DUPLICATEridge)
	  id= -2;
	else
	  id= getid_(neighbor);
        fprintf (fp, " %d", id);
      }
      fprintf (fp, "\n");
    }
  }
} /* printhashtable */
     

/*-------------------------------------------------
-printlists- print out facet and vertex list for debugging (without 'f/v' tags)
*/
void qh_printlists (void) {
  facetT *facet;
  vertexT *vertex;
  
  fprintf (qh ferr, "qh_printlists: facets:");
  FORALLfacets 
    fprintf (qh ferr, " %d", facet->id);
  fprintf (qh ferr, "\n  new facets %d visible facets %d next facet for addpoint %d\n  vertices (new %d):",
     getid_(qh newfacet_list), getid_(qh visible_list), getid_(qh facet_next),
     getid_(qh newvertex_list));
  FORALLvertices
    fprintf (qh ferr, " %d", vertex->id);
  fprintf (qh ferr, "\n");
} /* printlists */
  
/*-------------------------------------------------
-resetlists- reset newvertex_list, newfacet_list, visible_list
  clears num_visible
  if stats, maintains statistics
  visible_list is empty if qh_deletevisible was called
*/
void qh_resetlists (boolT stats /*qh newvertex_list newfacet_list visible_list*/) {
  vertexT *vertex;
  facetT *newfacet, *visible;
  int totnew=0, totver=0;
  
  if (stats) {
    FORALLvertex_(qh newvertex_list)
      totver++;
    FORALLnew_facets 
      totnew++;
    zadd_(Zvisvertextot, totver);
    zmax_(Zvisvertexmax, totver);
    zadd_(Znewfacettot, totnew);
    zmax_(Znewfacetmax, totnew);
  }
  FORALLvertex_(qh newvertex_list)
    vertex->newlist= False;
  qh newvertex_list= NULL;
  FORALLnew_facets
    newfacet->newfacet= False;
  qh newfacet_list= NULL;
  FORALLvisible_facets {
    visible->f.replace= NULL;
    visible->visible= False;
  }
  qh visible_list= NULL;
  qh num_visible= 0;
  qh NEWfacets= False;
} /* resetlists */

/*-------------------------------------------------
-setvoronoi_all- compute Voronoi centers for all facets
returns:
  facet->center is the Voronoi center
  includes upperDelaunay facets if !qh ATinfinity (i.e., 'Qu')
notes:
  unused/untested code: please email barber@tiac.net if this works ok for you
  Use FORALLvertices to locate the vertex for a point.  
  Use FOREACHneighbor_(vertex) to visit the Voronoi centers for a Voronoi cell.
*/
void qh_setvoronoi_all (void) {
  facetT *facet;

  qh_clearcenters (qh_ASvoronoi);
  qh_vertexneighbors();
  
  FORALLfacets {
    if (!facet->normal || !facet->upperdelaunay || !qh ATinfinity) {
      if (!facet->center)
        facet->center= qh_facetcenter (facet->vertices);
    }
  }
} /* setvoronoi_all */

/*-------------------------------------------------
-vertexintersect- intersects two vertex sets (inverse id ordered)
  temporary set vertexsetA is replaced by the intersection
     must be at top of stack
  could overwrite vertexsetA if currently too slow
*/
void qh_vertexintersect(setT **vertexsetA,setT *vertexsetB) {
  setT *intersection;

  intersection= qh_vertexintersect_new (*vertexsetA, vertexsetB);
  qh_settempfree (vertexsetA);
  *vertexsetA= intersection;
  qh_settemppush (intersection);
} /* vertexintersect */

/*-------------------------------------------------
-vertexintersect_new- intersects two vertex sets (inverse id ordered)
returns:
  a new set
*/
setT *qh_vertexintersect_new (setT *vertexsetA,setT *vertexsetB) {
  setT *intersection= qh_setnew (qh hull_dim - 1);
  vertexT **vertexA= SETaddr_(vertexsetA, vertexT); 
  vertexT **vertexB= SETaddr_(vertexsetB, vertexT); 

  while (*vertexA && *vertexB) {
    if (*vertexA  == *vertexB) {
      qh_setappend(&intersection, *vertexA);
      vertexA++; vertexB++;
    }else {
      if ((*vertexA)->id > (*vertexB)->id)
        vertexA++;
      else
        vertexB++;
    }
  }
  return intersection;
} /* vertexintersect_new */

/*-------------------------------------------
-vertexneighhbors- for each vertex in hull, determine facet neighbors
  nop if VERTEXneighbors
  assumes all vertex->neighbors are NULL
returns:
  sets qh VERTEXneighbors, qh_addpoint() will maintain them
*/
void qh_vertexneighbors (void /*qh facet_list*/) {
  facetT *facet;
  vertexT *vertex, **vertexp;

  if (qh VERTEXneighbors)
    return;
  trace1((qh ferr, "qh_vertexneighbors: determing neighboring facets for each vertex\n"));
  qh vertex_visit++;
  FORALLfacets {
    if (facet->visible)
      continue;
    FOREACHvertex_(facet->vertices) {
      if (vertex->visitid != qh vertex_visit) {
        vertex->visitid= qh vertex_visit;
        vertex->neighbors= qh_setnew (qh hull_dim);
      }
      qh_setappend (&vertex->neighbors, facet);
    }
  }
  qh VERTEXneighbors= True;
} /* vertexneighbors */

/*-------------------------------------------------
-vertexsubset- returns True if vertexsetA is a subset of vertexsetB, False
  otherwise; relies on vertexsets being sorted;
an empty set is a subset of any other set
*/
boolT qh_vertexsubset(setT *vertexsetA, setT *vertexsetB) {
  vertexT **vertexA= (vertexT **) SETaddr_(vertexsetA, vertexT);
  vertexT **vertexB= (vertexT **) SETaddr_(vertexsetB, vertexT);

  while (True) {
    if (!*vertexA)
      return True;
    if (!*vertexB)
      return False;
    if ((*vertexA)->id > (*vertexB)->id)
      return False;
    if (*vertexA  == *vertexB)
      vertexA++;
    vertexB++; 
  }
} /* vertexsubset */
