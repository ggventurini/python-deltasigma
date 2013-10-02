/* geom2.c -- infrequently used geometric routines of qhull
   
   see README and geom.h

   copyright (c) 1993-1995 The Geometry Center        

   frequently used code goes into geom.c
*/
   
#include "qhull_a.h"
   


/*-------------------------------------------------
-crossproduct- of 2 dim vectors, C= A x B
    from Glasner, Graphics Gems I, p. 639
    NOTE: only defined for dim==3
*/
void qh_crossproduct (int dim, realT vecA[3], realT vecB[3], realT vecC[3]){

  if (dim == 3) {
    vecC[0]=   det2_(vecA[1], vecA[2],
		     vecB[1], vecB[2]);
    vecC[1]= - det2_(vecA[0], vecA[2],
		     vecB[0], vecB[2]);
    vecC[2]=   det2_(vecA[0], vecA[1],
		     vecB[0], vecB[1]);
  }
} /* vcross */

/*-------------------------------------------------
-determinant- compute signed determinant of a square matrix
  rows= row vectors
  uses qh NEARzero to test for degenerate matrices
    this does look right, probably no easy way of doing it
returns:
  determinant
  overwrites rows and the matrix
  nearzero set if degenerate, else cleared
   if dim == 2 or 3
     nearzero iff determinant < qh NEARzero[dim-1]  (not quite correct)
   if dim >= 4
     nearzero iff diagonal[k] < qh NEARzero[k]
*/
realT qh_determinant (realT **rows, int dim, boolT *nearzero) {
  realT det=0;
  int i;
  boolT sign= False;

  *nearzero= False;
  if (dim < 2) {
    fprintf (qh ferr, "qhull internal error (qh_determinate): only implemented for dimension >= 2\n");
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }else if (dim == 2) {
    det= det2_(rows[0][0], rows[0][1],
		 rows[1][0], rows[1][1]);
    if (fabs_(det) < qh NEARzero[1])  /* not really correct, what should this be? */
      *nearzero= True;
  }else if (dim == 3) {
    det= det3_(rows[0][0], rows[0][1], rows[0][2],
		 rows[1][0], rows[1][1], rows[1][2],
		 rows[2][0], rows[2][1], rows[2][2]);
    if (fabs_(det) < qh NEARzero[2])  /* not really correct, what should this be? */
      *nearzero= True;
  }else {	
    qh_gausselim(rows, dim, dim, &sign, nearzero);  /* if nearzero, diagonal still ok*/
    det= 1.0;
    for (i= dim; i--; )
      det *= (rows[i])[i];
    if (sign)
      det= -det;
  }
  return det;
} /* determinant */

/*-------------------------------------------------
-detsimplex- compute determinant of a simplex with point apex and base points
   uses qh gm_matrix/qh gm_row (assumes they're big enough)
   uses dim coordinates of point and vertex->point
returns:
   signed determinant and nearzero from qh_determinant
*/
realT qh_detsimplex(pointT *apex, setT *points, int dim, boolT *nearzero) {
  pointT *coorda, *coordp, *gmcoord, *point, **pointp;
  coordT **rows;
  int k,  i=0;
  realT det;

  zinc_(Zdetsimplex);
  gmcoord= qh gm_matrix;
  rows= qh gm_row;
  FOREACHpoint_(points) {
    if (i == dim)
      break;
    rows[i++]= gmcoord;
    coordp= point;
    coorda= apex;
    for (k= dim; k--; )
      *(gmcoord++)= *coordp++ - *coorda++;
  }
  if (i < dim) {
    fprintf (qh ferr, "qhull internal error (qh_detsimplex): #points %d < dimension %d\n", 
               i, dim);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  det= qh_determinant (rows, dim, nearzero);
  trace2((qh ferr, "qh_detsimplex: det=%2.2g for point p%d, dim %d, nearzero? %d\n",
	  det, qh_pointid(apex), dim, *nearzero)); 
  return det;
} /* detsimplex */

/*--------------------------------------------------
-divzero -- divide by a number that's nearly zero
  mindenom1= minimum denominator for dividing into 1.0
returns:
  zerodiv and 0.0 if it would overflow
*/
realT qh_divzero (realT numer, realT denom, realT mindenom1, boolT *zerodiv) {
  realT temp, numerx, denomx;
  

  if (numer < mindenom1 && numer > -mindenom1) {
    numerx= fabs_(numer);
    denomx= fabs_(denom);
    if (numerx < denomx) {
      *zerodiv= False;
      return numer/denom;
    }else {
      *zerodiv= True;
      return 0.0;
    }
  }
  temp= denom/numer;
  if (temp > mindenom1 || temp < -mindenom1) {
    *zerodiv= False;
    return numer/denom;
  }else {
    *zerodiv= True;
    return 0.0;
  }
} /* divzero */
  

/*-------------------------------------------------
-facetarea- return area for a facet
  if non-simplicial, uses centrum to triangulate facet.  Sums the projected areas.
  assumes facet->normal exists
  if (qh DELAUNAY),
     computes projected area instead for last coordinate
*/
realT qh_facetarea (facetT *facet) {
  vertexT *apex;
  pointT *centrum;
  realT area= 0.0;
  ridgeT *ridge, **ridgep;

  if (facet->simplicial) {
    apex= SETfirst_(facet->vertices);
    area= qh_facetarea_simplex (qh hull_dim, apex->point, facet->vertices, 
                    apex, facet->toporient, facet->normal, &facet->offset);
  }else {
    if (qh CENTERtype == qh_AScentrum)
      centrum= facet->center;
    else
      centrum= qh_getcentrum (facet);
    FOREACHridge_(facet->ridges) 
      area += qh_facetarea_simplex (qh hull_dim, centrum, ridge->vertices, 
                 NULL, (ridge->top == facet),  facet->normal, &facet->offset);
    if (qh CENTERtype != qh_AScentrum)
      qh_memfree (centrum, qh normal_size);
  }
  trace4((qh ferr, "qh_facetarea: f%d area %2.2g\n", facet->id, area)); 
  return area;
} /* facetarea */

/*-------------------------------------------------
-facetarea_simplex- return area for a simplex defined by an apex,
        a base of vertices, an orientation, and a unit normal
  if simplicial facet, notvertex is defined and it is skipped in vertices
  if (qh DELAUNAY),
     computes projected area instead for last coordinate
notes:
  computes area of simplex projected to plane [normal,offset]
  returns 0 if vertex too far below plane (qh WIDEfacet)
  uses qh gm_matrix/gm_row and qh hull_dim
  helper function for qh_facetarea
*/
realT qh_facetarea_simplex (int dim, coordT *apex, setT *vertices, 
        vertexT *notvertex,  boolT toporient, coordT *normal, realT *offset) {
  pointT *coorda, *coordp, *gmcoord;
  coordT **rows, *normalp;
  int k,  i=0;
  realT area, dist;
  vertexT *vertex, **vertexp;
  boolT nearzero;

  gmcoord= qh gm_matrix;
  rows= qh gm_row;
  FOREACHvertex_(vertices) {
    if (vertex == notvertex)
      continue;
    rows[i++]= gmcoord;
    coorda= apex;
    coordp= vertex->point;
    normalp= normal;
    if (notvertex) {
      for (k= dim; k--; )
	*(gmcoord++)= *coordp++ - *coorda++;
    }else {
      dist= *offset;
      for (k= dim; k--; )
	dist += *coordp++ * *normalp++;
      if (dist < -qh WIDEfacet) {
	zinc_(Znoarea);
	return 0.0;
      }
      coordp= vertex->point;
      normalp= normal;
      for (k= dim; k--; )
	*(gmcoord++)= (*coordp++ - dist * *normalp++) - *coorda++;
    }
  }
  if (i != dim-1) {
    fprintf (qh ferr, "qhull internal error (qh_facetarea_simplex): #points %d != dim %d -1\n", 
               i, dim);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  rows[i]= gmcoord;
  if (qh DELAUNAY) {
    for (i= 0; i < dim-1; i++)
      rows[i][dim-1]= 0.0;
    for (k= dim; k--; )
      *(gmcoord++)= 0.0;
    rows[dim-1][dim-1]= 1.0;
  }else {
    normalp= normal;
    for (k= dim; k--; )
      *(gmcoord++)= *normalp++;
  }
  zinc_(Zdetsimplex);
  area= qh_determinant (rows, dim, &nearzero);
  if (toporient)
    area= -area;
  area *= qh AREAfactor;
  trace4((qh ferr, "qh_facetarea_simplex: area=%2.2g for point p%d, toporient %d, nearzero? %d\n",
	  area, qh_pointid(apex), toporient, nearzero)); 
  return area;
} /* facetarea_simplex */

/*-------------------------------------------------
-facetcenter- return Voronoi center for a facet's vertices
*/
pointT *qh_facetcenter (setT *vertices) {
  setT *points= qh_settemp (qh_setsize (vertices));
  vertexT *vertex, **vertexp;
  pointT *center;
  
  FOREACHvertex_(vertices) 
    qh_setappend (&points, vertex->point);
  center= qh_voronoi_center (qh hull_dim-1, points);
  qh_settempfree (&points);
  return center;
} /* facetcenter */

/*-------------------------------------------------
-findbestsharp- find best facet on newfacet_list
  skips visited facets (qh visit_id) on qh newfacet_list
  skips upperdelaunay facets unless point is outside
returns:
  true if could be an acute angle (facets in different quadrants)
  returns bestfacet and distance to facet
  increments numpart by number of distance tests
  for qh_findbest
*/
boolT qh_findbestsharp (pointT *point, facetT **bestfacet, 
           realT *bestdist, int *numpart) {
  facetT *facet;
  realT dist;
  boolT issharp = False;
  int *quadrant, k;
  
  quadrant= (int*)qh_memalloc (qh hull_dim * sizeof(int));
  FORALLfacet_(qh newfacet_list) {
    if (facet == qh newfacet_list) {
      for (k= qh hull_dim; k--; )
      	quadrant[ k]= (facet->normal[ k] > 0);
    }else if (!issharp) {
      for (k= qh hull_dim; k--; ) {
        if (quadrant[ k] != (facet->normal[ k] > 0)) {
          issharp= True;
          break;
        }
      }
    }
    if (facet->visitid != qh visit_id) {
      qh_distplane (point, facet, &dist);
      (*numpart)++;
      if (dist > *bestdist) {
      	if (!facet->upperdelaunay || dist > qh MINoutside) {
 	  *bestdist= dist;
	  *bestfacet= facet;
	}
      }
    }
  }
  qh_memfree( quadrant, qh hull_dim * sizeof(int));
  return issharp;
} /* findbestsharp */
  
/*-------------------------------------------------
-findgooddist- find best good facet visible for point from facetA
  assumes facetA is visible from point
  uses qh visit_id
returns:
  furthest distance to good facet, if any
  list of good, visible facets (and some other visible facets)
     at end of qh facet_list
*/
facetT *qh_findgooddist (pointT *point, facetT *facetA, realT *distp, 
               facetT **facetlist) {
  realT bestdist= -REALmax, dist;
  facetT *neighbor, **neighborp, *bestfacet=NULL, *facet;
  boolT goodseen= False;  

  if (facetA->good) {
    zinc_(Zcheckpart);  /* calls from check_bestdist occur after print stats */
    qh_distplane (point, facetA, &bestdist);
    bestfacet= facetA;
    goodseen= True;
  }
  qh_removefacet (facetA);
  qh_appendfacet (facetA);
  *facetlist= facetA;
  facetA->visitid= ++qh visit_id;
  FORALLfacet_(*facetlist) {
    FOREACHneighbor_(facet) {
      if (neighbor->visitid == qh visit_id)
        continue;
      neighbor->visitid= qh visit_id;
      if (goodseen && !neighbor->good)
        continue;
      zinc_(Zcheckpart); 
      qh_distplane (point, neighbor, &dist);
      if (dist > 0) {
        qh_removefacet (neighbor);
        qh_appendfacet (neighbor);
        if (neighbor->good) {
          goodseen= True;
          if (dist > bestdist) {
            bestdist= dist;
            bestfacet= neighbor;
          }
        }
      }
    }
  }
  if (bestfacet) {
    *distp= bestdist;
    trace2((qh ferr, "qh_findgooddist: p%d is %2.2g above good facet f%d\n",
      qh_pointid(point), bestdist, bestfacet->id));
    return bestfacet;
  }
  trace4((qh ferr, "qh_findgooddist: no good facet for p%d above f%d\n", 
      qh_pointid(point), facetA->id));
  return NULL;
}  /* findgooddist */
    
/*------------------------------------------
-getarea- get area of all facets in facetlist, collect statistics

notes:
  could compute outer volume by expanding facet area by rays from interior
#if qh_MAXoutside  // a perpendicular projection underestimates badly
      qh totoutvol += (-dist + facet->maxoutside + qh DISTround) 
                            * area/ qh hull_dim;
#endif
*/
void qh_getarea (facetT *facetlist) {
  realT area;
  realT dist;
  facetT *facet;

  if (qh REPORTfreq)
    fprintf (qh ferr, "computing area of each facet and volume of the convex hull\n");
  else 
    trace1((qh ferr, "qh_getarea: computing volume and area for each facet\n"));
  qh totarea= qh totvol= 0.0;
  FORALLfacet_(facetlist) {
    if (!facet->normal)
      continue;
    facet->f.area= area= qh_facetarea (facet);
    facet->isarea= True;
    if (!qh DELAUNAY) {
      qh_distplane (qh interior_point, facet, &dist);
      qh totvol += -dist * area/ qh hull_dim;
    }
    qh totarea += area;
    if (qh PRINTstatistics) {
      wadd_(Wareatot, area);
      wmax_(Wareamax, area);
      wmin_(Wareamin, area);
    }
  }
} /* getarea */

/*-------------------------------------------------
-gram_schmidt- implements Gram-Schmidt orthogonalization by rows
   overwrites rows[dim][dim]
returns:
   false if gets a zero norm
notes:
   see Golub & van Loan Algorithm 6.2-2
   overflow due to small divisors not handled
*/
boolT qh_gram_schmidt(int dim, realT **row) {
  realT *rowi, *rowj, norm;
  int i, j, k;
  
  for(i=0; i < dim; i++) {
    rowi= row[i];
    for (norm= 0.0, k= dim; k--; rowi++)
      norm += *rowi * *rowi;
    norm= sqrt(norm);
    wmin_(Wmindenom, norm);
    if (norm == 0.0)  /* either 0 or overflow due to sqrt */
      return False;
    for(k= dim; k--; )
      *(--rowi) /= norm;  
    for(j= i+1; j < dim; j++) {
      rowj= row[j];
      for(norm= 0.0, k=dim; k--; )
	norm += *rowi++ * *rowj++;
      for(k=dim; k--; )
	*(--rowj) -= *(--rowi) * norm;
    }
  }
  return True;
} /* gram_schmidt */



	
/*--------------------------------------------------
-inthresholds- return True if normal within qh lower_/upper_threshold
returns:
  angle cos to a threshold border (may be NULL, invalid if qh SPLITthresholds)
*/
boolT qh_inthresholds (coordT *normal, realT *angle) {
  boolT within= True;
  int k;

  if (angle)
    *angle= 0.0;
  for(k= 0; k < qh hull_dim; k++) {
    if (qh lower_threshold[k] > -REALmax/2) {
      if (normal[k] < qh lower_threshold[k])
        within= False;
      if (angle)
        *angle += normal[k] * qh lower_threshold[k];
    }
    if (qh upper_threshold[k] < REALmax/2) {
      if (normal[k] > qh upper_threshold[k])
        within= False;
      if (angle)
        *angle += normal[k] * qh upper_threshold[k];
    }
  }
  return within;
} /* inthresholds */
    

/*--------------------------------------------------
-maxabsval -- return pointer to maximum absolute value of a dim vector
   returns NULL if dim==0
*/
realT *qh_maxabsval (realT *normal, int dim) {
  realT maxval= -REALmax;
  realT *maxp= NULL, *colp, absval;
  int k;

  for (k= dim, colp= normal; k--; colp++) {
    absval= fabs_(*colp);
    if (absval > maxval) {
      maxval= absval;
      maxp= colp;
    }
  }
  return maxp;
} /* maxabsval */


/*-------------------------------------------------
-maxmin- collects the maximum and minimum points of input into a set
  determines maximum roundoff errors
returns:
  returns a temporary set, without qh GOODpoint
  points are not unique
*/
setT *qh_maxmin(pointT *points, int numpoints, int dimension) {
  int k;
  realT maxsum= 0.0, maxcoord, temp, maxdistsum;
  realT maxneg= REALmax, maxpos= -REALmax;
  pointT *minimum, *maximum, *point, *pointtemp;
  setT *set;

  if (REALmin < REALepsilon && REALmin < REALmax && REALmin > -REALmax
  && REALmax > 0.0 && -REALmax < 0.0)
    ; /* all ok */
  else {
    fprintf (qh ferr, "qhull error: floating point constants in user.h are wrong\n\
REALepsilon %g REALmin %g REALmax %g -REALmax %g\n",
	     REALepsilon, REALmin, REALmax, -REALmax);
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  set= qh_settemp(2*dimension);
  for(k= 0; k < dimension; k++) {
    if (points == qh GOODpointp)
      minimum= maximum= points + qh hull_dim;
    else
      minimum= maximum= points;
    FORALLpoint_(points, numpoints) {
      if (point == qh GOODpointp)
	continue;
      if (maximum[k] < point[k])
	maximum= point;
      else if (minimum[k] > point[k])
	minimum= point;
    }
    maxcoord= fmax_(maximum[k], -minimum[k]);
    if (qh GOODpointp) {
      temp= fmax_(qh GOODpointp[k], -qh GOODpointp[k]);
      maximize_(maxcoord, temp);
    }
    maximize_(qh maxmaxcoord, maxcoord);
    maxsum += maxcoord;
    maximize_(maxpos, maximum[k]);
    minimize_(maxneg, minimum[k]);
    qh_setappend (&set, maximum);
    qh_setappend (&set, minimum);
    /* calculation of qh NEARzero is based on error formula 4.4-13 of
       Golub & van Loan, authors say n^3 can be ignored and 10 be used in
       place of rho */
    qh NEARzero[k]= 80 * maxsum * REALepsilon;
  }
  /* calculate roundoff error according to
     Lemma 3.2-1 of Golub and van Loan "Matrix Computation"
     use sqrt(dim) since one vector is normalized */
  maxdistsum= sqrt (qh hull_dim) * qh maxmaxcoord;
  if (!qh SETroundoff) {
    qh DISTround= REALepsilon * (qh hull_dim * maxdistsum * 1.01
			   	       + qh maxmaxcoord);  /* for offset */
    if (qh RANDOMdist)
      qh DISTround += qh RANDOMfactor * qh maxmaxcoord;
    qh_option ("Error-roundoff", NULL, &qh DISTround);
  }
  qh MINdenom= qh MINdenom_1 * qh maxmaxcoord;
  qh MINdenom_1_2= sqrt (qh MINdenom_1 * qh hull_dim) ;  /* if will be normalized */
  qh MINdenom_2= qh MINdenom_1_2 * qh maxmaxcoord;
                                              /* for inner product */
  qh ANGLEround= 1.01 * qh hull_dim * REALepsilon;
  if (qh RANDOMdist)
    qh ANGLEround += qh RANDOMfactor;
  if (qh premerge_cos < REALmax/2) {
    qh premerge_cos -= qh ANGLEround;
    if (qh RANDOMdist) 
      qh_option ("Angle-premerge-with-random", NULL, &qh premerge_cos);
  }
  if (qh postmerge_cos < REALmax/2) {
    qh postmerge_cos -= qh ANGLEround;
    if (qh RANDOMdist)
      qh_option ("Angle-postmerge-with-random", NULL, &qh postmerge_cos);
  }
  qh premerge_centrum += 2 * qh DISTround;    /*2 for centrum and distplane()*/
  qh postmerge_centrum += 2 * qh DISTround;
  if (qh RANDOMdist && (qh MERGEexact || qh PREmerge))
    qh_option ("Centrum-premerge-with-random", NULL, &qh premerge_centrum);
  if (qh RANDOMdist && qh POSTmerge)
    qh_option ("Centrum-postmerge-with-random", NULL, &qh postmerge_centrum);
  qh_option ("_max-coord", NULL, &maxpos);
  qh_option ("_min-coord", NULL, &maxneg);
  { /* compute ONEmerge, max vertex offset for merging simplicial facets */
    realT maxangle= 1.0, maxrho;
    
    minimize_(maxangle, qh premerge_cos);
    minimize_(maxangle, qh postmerge_cos);
    /* max diameter * sin theta + DISTround for vertex to its hyperplane */
    qh ONEmerge= sqrt (qh hull_dim) * (maxpos - maxneg) *
      sqrt (1.0 - maxangle * maxangle) + qh DISTround;  
    maxrho= qh hull_dim * qh premerge_centrum + qh DISTround;
    maximize_(qh ONEmerge, maxrho);
    maxrho= qh hull_dim * qh postmerge_centrum + qh DISTround;
    maximize_(qh ONEmerge, maxrho);
    if (qh MERGING)
      qh_option ("_one-merge", NULL, &qh ONEmerge);
  }
  qh NEARinside= qh ONEmerge * qh_RATIOnearinside; /* only used if qh KEEPnearinside */
  if (qh KEEPnearinside)
    qh_option ("_near-inside", NULL, &qh NEARinside);
  qh MINnorm= sqrt( REALepsilon) * qh maxmaxcoord; /* FIXUP: what is correct?*/
  if (qh hull_dim <= 4) /* used in qh_sethyperplane_det */
    qh_option ("_min-norm", NULL, &qh MINnorm);
  if (qh MINvisible > REALmax/2) {
    if (!qh MERGING)
      qh MINvisible= qh DISTround;
    else if (qh hull_dim <= 3)
      qh MINvisible= qh premerge_centrum;
    else
      qh MINvisible= qh_COPLANARratio * qh premerge_centrum;
    if (qh APPROXhull && qh MINvisible > qh MINoutside)
      qh MINvisible= qh MINoutside;
    qh_option ("Visible-distance", NULL, &qh MINvisible);
  }
  if (qh MAXcoplanar > REALmax/2) {
    qh MAXcoplanar= qh MINvisible;
    qh_option ("U-coplanar-distance", NULL, &qh MAXcoplanar);
  }
  if (!qh APPROXhull) {             /* user may specify qh MINoutside */
    qh MINoutside= 2 * qh MINvisible;
    if (qh premerge_cos < REALmax/2) 
      maximize_(qh MINoutside, (1- qh premerge_cos) * qh maxmaxcoord);
    qh_option ("Width-outside", NULL, &qh MINoutside);
  }
  qh WIDEfacet= qh MINoutside;
  maximize_(qh WIDEfacet, qh_WIDEcoplanar * qh MAXcoplanar); 
  maximize_(qh WIDEfacet, qh_WIDEcoplanar * qh MINvisible); 
  qh_option ("_wide-facet", NULL, &qh WIDEfacet);
  if (qh MINvisible > qh MINoutside + 3 * REALepsilon 
  && !qh BESToutside && !qh FORCEoutput)
    fprintf (qh ferr, "qhull input warning: minimum visibility V%.2g is greater than \nminimum outside W%.2g.  Flipped facets are likely.\n",
	     qh MINvisible, qh MINoutside);
  qh max_vertex= qh DISTround;
  qh min_vertex= -qh DISTround;
  if (qh IStracing >=1)
    qh_printpoints (qh ferr, "qh_maxmin: found the max and min points (by dim):", set);
  /* numeric constants reported in printsummary */
  return(set);
} /* maxmin */


/*-------------------------------------------------
-maxsimplex- determines maximum simplex for a set of points 
  assumes at least pointsneeded points in points
  skips qh GOODpointp (assumes that it isn't in maxpoints)
  starts from points already in simplex
returns:
  temporary set of dim+1 points
notes:
  maximizes determinate for x,y,z,w, etc.
  uses maxpoints as long as determinate is clearly non-zero
*/
void qh_maxsimplex (int dim, setT *maxpoints, pointT *points, int numpoints, setT **simplex) {
  pointT *point, **pointp, *pointtemp, *maxpoint, *minx=NULL, *maxx=NULL;
  boolT nearzero, maxnearzero= False;
  int k, sizinit;
  realT maxdet= -REALmax, det, mincoord= REALmax, maxcoord= -REALmax;

  sizinit= qh_setsize (*simplex);
  if (sizinit < 2) {
    if (qh_setsize (maxpoints) >= 2) {
      FOREACHpoint_(maxpoints) {
	
        if (maxcoord < point[0]) {
          maxcoord= point[0];
          maxx= point;
        }
	if (mincoord > point[0]) {
          mincoord= point[0];
          minx= point;
        }
      }
    }else {
      FORALLpoint_(points, numpoints) {
	if (point == qh GOODpointp)
	  continue;
        if (maxcoord < point[0]) {
	  maxcoord= point[0];
          maxx= point;
        }
	if (mincoord > point[0]) {
          mincoord= point[0];
          minx= point;
	}
      }
    }
    qh_setunique (simplex, minx);
    if (qh_setsize (*simplex) < 2)
      qh_setunique (simplex, maxx);
    sizinit= qh_setsize (*simplex);
    if (sizinit < 2) {
      if (zzval_(Zsetplane) > qh hull_dim+1) {
	fprintf (qh ferr, "qhull precision error (qh_maxsimplex for voronoi_center):\n%d points with the same x coordinate.\n",
		 qh_setsize(maxpoints)+numpoints);
	qh_errexit (qh_ERRprec, NULL, NULL);
      }else {
	fprintf (qh ferr, "qhull input error: input is less than %d-dimensional since it has the same x coordinate\n", qh hull_dim);
	qh_errexit (qh_ERRinput, NULL, NULL);
      }
    }
  }
  for(k= sizinit; k < dim+1; k++) {
    maxpoint= NULL;
    maxdet= -REALmax;
    FOREACHpoint_(maxpoints) {
      if (!qh_setin (*simplex, point)) {
        det= qh_detsimplex(point, *simplex, k, &nearzero);
        if ((det= fabs_(det)) > maxdet) {
	  maxdet= det;
          maxpoint= point;
	  maxnearzero= nearzero;
        }
      }
    }
    if (!maxpoint || maxnearzero) {
      zinc_(Zsearchpoints);
      if (!maxpoint) {
        trace0((qh ferr, "qh_maxsimplex: searching all points for %d-th initial vertex.\n", k));
      }else {
        trace0((qh ferr, "qh_maxsimplex: searching all points for %d-th initial vertex, better than p%d det %2.2g\n",
		k+1, qh_pointid(maxpoint), maxdet));
      }
      FORALLpoint_(points, numpoints) {
	if (point == qh GOODpointp)
	  continue;
        if (!qh_setin (*simplex, point)) {
          det= qh_detsimplex(point, *simplex, k, &nearzero);
          if ((det= fabs_(det)) > maxdet) {
	    maxdet= det;
            maxpoint= point;
	    maxnearzero= nearzero;
	  }
        }
      }
    } /* !maxpoint */
    if (!maxpoint) {
      fprintf (qh ferr, "qhull internal error (qh_maxsimplex): not enough points available\n");
      qh_errexit (qh_ERRqhull, NULL, NULL);
    }
    qh_setappend(simplex, maxpoint);
    trace1((qh ferr, "qh_maxsimplex: selected point p%d for %d`th initial vertex, det=%2.2g\n",
	    qh_pointid(maxpoint), k, maxdet));
  } /* k */ 
} /* maxsimplex */

/*--------------------------------------------------
-minabsval -- return min absolute value of a dim vector
*/
realT qh_minabsval (realT *normal, int dim) {
  realT minval= 0;
  realT maxval= 0;
  realT *colp;
  int k;

  for (k= dim, colp= normal; k--; colp++) {
    maximize_(maxval, *colp);
    minimize_(minval, *colp);
  }
  return fmax_(maxval, -minval);
} /* minabsval */


/*--------------------------------------------------
-mindiff -- return index of min abs. difference of two vectors
*/
int qh_mindiff (realT *vecA, realT *vecB, int dim) {
  realT mindiff= REALmax, diff;
  realT *vecAp= vecA, *vecBp= vecB;
  int k, mink;

  for (k= 0; k<dim; k++) {
    diff= *vecAp++ - *vecBp++;
    diff= fabs_(diff);
    if (diff < mindiff) {
      mindiff= diff;
      mink= k;
    }
  }
  return k;
} /* mindiff */



/*-------------------------------------------
-orientoutside- make facet outside oriented via qh interior_point
  returns True if reversed orientation.
*/
boolT qh_orientoutside (facetT *facet) {
  int k;
  realT dist;

  qh_distplane (qh interior_point, facet, &dist);
  if (dist > 0) {
    for (k= qh hull_dim; k--; )
      facet->normal[k]= -facet->normal[k];
    facet->offset= -facet->offset;
    return True;
  }
  return False;
} /* orientoutside */

/*-------------------------------------------
-pointdist- distance between points
  returns distance squared if 'dim' is negative
*/
coordT qh_pointdist(pointT *point1, pointT *point2, int dim) {
  coordT dist, diff;
  int k;
  
  dist= 0.0;
  for (k= (dim > 0 ? dim : -dim); k--; ) {
    diff= *point1++ - *point2++;
    dist += diff * diff;
  }
  if (dim > 0)
    return(sqrt(dist));
  return dist;
} /* pointdist */


/*-------------------------------------------------
-printmatrix- print matrix given by row vectors
  print a vector by (fp, "", &vect, 1, len)
*/
void qh_printmatrix (FILE *fp, char *string, realT **rows, int numrow, int numcol) {
  realT *rowp;
  realT r; /*bug fix*/
  int i,k;

  fprintf (fp, "%s\n", string);
  for (i= 0; i<numrow; i++) {
    rowp= rows[i];
    for (k= 0; k<numcol; k++)
      fprintf (fp, "%6.3g ", r=*rowp++);
    fprintf (fp, "\n");
  }
} /* printmatrix */

  
/*-------------------------------------------------
-printpoints- print pointids for a set of points starting at index 
  prints string and 'p' if defined
*/
void qh_printpoints (FILE *fp, char *string, setT *points) {
  pointT *point, **pointp;

  if (string) {
    fprintf (fp, "%s", string);
    FOREACHpoint_(points) 
      fprintf (fp, " p%d", qh_pointid(point));
    fprintf (fp, "\n");
  }else {
    FOREACHpoint_(points) 
      fprintf (fp, " %d", qh_pointid(point));
    fprintf (fp, "\n");
  }
} /* printpoints */

  
/*-------------------------------------------------
-projectinput- project input points using qh DELAUNAY and qh lower_bound/upper_bound
  input points in qh first_point, num_points, input_dim
     if POINTSmalloc, will free old point array
  if low[k]=high[k]= 0, removes dimension k 
     checks that hull_dim agrees with input_dim, PROJECTinput, and DELAUNAY
  if DELAUNAY 
    projects points to paraboloid
returns:
  new point array in first_point of qh hull_dim coordinates
  sets POINTSmalloc
  lowbound/highbound is also projected
*/
void qh_projectinput (void) {
  int k,i;
  int newdim= qh input_dim, newnum= qh num_points;
  signed char *project;
  int size= (qh input_dim+1)*sizeof(*project);
  pointT *newpoints, *coord, *infinity;
  realT paraboloid, maxboloid= 0;
  
  project= (signed char*)qh_memalloc (size);
  memset ((char*)project, 0, size);
  for (k= 0; k<qh input_dim; k++) {   /* skip Delaunay bound */
    if (qh lower_bound[k] == 0 && qh upper_bound[k] == 0) {
      project[k]= -1;
      newdim--;
    }
  }
  if (qh DELAUNAY) {
    project[k]= 1;
    newdim++;
    if (qh ATinfinity)
      newnum++;
  }
  if (newdim != qh hull_dim) {
    fprintf(qh ferr, "qhull internal error (qh_projectinput): dimension after projection %d != hull_dim %d\n", newdim, qh hull_dim);
    qh_errexit(qh_ERRqhull, NULL, NULL);
  }
  if (!(newpoints=(coordT*)malloc(newnum*newdim*sizeof(coordT)))){
    fprintf(qh ferr, "qhull error: insufficient memory to project %d points\n",
           qh num_points);
    qh_errexit(qh_ERRmem, NULL, NULL);
  }
  qh_projectpoints (project, qh input_dim+1, qh first_point,
                    qh num_points, qh input_dim, newpoints, newdim);
  trace1((qh ferr, "qh_projectinput: updating lower and upper_bound\n"));
  qh_projectpoints (project, qh input_dim+1, qh lower_bound,
                    1, qh input_dim+1, qh lower_bound, newdim+1);
  qh_projectpoints (project, qh input_dim+1, qh upper_bound,
                    1, qh input_dim+1, qh upper_bound, newdim+1);
  qh_memfree(project, ((qh input_dim+1)*sizeof(*project)));
  if (qh POINTSmalloc)
    free (qh first_point);
  qh first_point= newpoints;
  qh POINTSmalloc= True;
  if (qh DELAUNAY && qh ATinfinity) {
    coord= qh first_point;
    infinity= qh first_point + qh hull_dim * qh num_points;
    for (k=qh hull_dim-1; k--; )
      infinity[k]= 0.0;
    for (i=qh num_points; i--; ) {
      paraboloid= 0.0;
      for (k=qh hull_dim-1; k--; ) {
        paraboloid += *coord * *coord;
	infinity[k] += *coord;
        coord++;
      }
      *(coord++)= paraboloid;
      maximize_(maxboloid, paraboloid);
    }
    for (k=qh hull_dim-1; k--; )
      *(coord++) /= qh num_points;
    *(coord++)= maxboloid * 1.1;
    qh num_points++;
    trace0((qh ferr, "qh_projectinput: projected points to paraboloid for Delaunay\n"));
  }
} /* projectinput */

  
/*-------------------------------------------------
-projectpoints- project along one or more dimensions
  delete dimension k if project[k] == -1
  add dimension k if project[k] == 1 
  n is size of project
  points, numpoints, dim is old points
  newpoints, newdim is buffer for new points (already allocated)
    newpoints may be points if only adding dimension at end
*/
void qh_projectpoints (signed char *project, int n, realT *points, 
        int numpoints, int dim, realT *newpoints, int newdim) {
  int testdim= dim, oldk=0, newk=0, i,j=0,k;
  realT *newp, *oldp;
  
  for (k= 0; k<n; k++)
    testdim += project[k];
  if (testdim != newdim) {
    fprintf (qh ferr, "qhull internal error (qh_projectpoints): newdim %d should be %d after projection\n",
      newdim, testdim);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  for (j= 0; j<n; j++) {
    if (project[j] == -1)
      oldk++;
    else {
      newp= newpoints+newk++;
      if (project[j] == +1) {
	if (oldk >= dim)
	  continue;
	oldp= points+oldk;
      }else 
	oldp= points+oldk++;
      for (i=numpoints; i--; ) {
        *newp= *oldp;
        newp += newdim;
        oldp += dim;
      }
    }
    if (oldk >= dim)
      break;
  }
  trace1((qh ferr, "qh_projectpoints: projected %d points from dim %d to dim %d\n", 
    numpoints, dim, newdim));
} /* projectpoints */
        

/*-------------------------------------------------
-rand & srand- generate pseudo-random number between 1 and 2^31 -2
  from Park & Miller's minimimal standard random number generator
  Communications of the ACM, 31:1192-1201, 1988.
notes:
  does not use 0 or 2^31 -1
  this is silently enforced by qh_srand()
  copied to rbox.c
  can make 'Rn' much faster by moving qh_rand to qh_distplane
*/
int qh_rand( void) {
#define qh_rand_a 16807
#define qh_rand_m 2147483647
#define qh_rand_q 127773  /* m div a */
#define qh_rand_r 2836    /* m mod a */
  int lo, hi, test;
  int seed = qh rand_seed;

  hi = seed / qh_rand_q;  /* seed div q */
  lo = seed % qh_rand_q;  /* seed mod q */
  test = qh_rand_a * lo - qh_rand_r * hi;
  if (test > 0)
    seed= test;
  else
    seed= test + qh_rand_m;
  qh rand_seed= seed;
  return seed;
} /* rand */

void qh_srand( int seed) {
  if (seed < 1)
    qh rand_seed= 1;
  else if (seed >= qh_rand_m)
    qh rand_seed= qh_rand_m - 1;
  else
    qh rand_seed= seed;
} /* qh_srand */

/*-------------------------------------------------
-randomfactor- return a random factor within qh RANDOMmax of 1.0
  RANDOMa/b definedin global.c
*/
realT qh_randomfactor (void) {
  realT randr;

  randr= qh_RANDOMint;
  return randr * qh RANDOMa + qh RANDOMb;
} /* randomfactor */

/*-------------------------------------------------
-randommatrix- generate a random dimXdim matrix in range (-1,1)
  assumes buffer is dim+1Xdim
returns:
  returns row vector for buffer
  plus row[dim] for scratch
*/
void qh_randommatrix (realT *buffer, int dim, realT **row) {
  int i, k;
  realT **rowi, *coord, realr;

  coord= buffer;
  rowi= row;
  for (i=0; i<dim; i++) {
    *(rowi++)= coord;
    for (k=0; k<dim; k++) {
      realr= qh_RANDOMint;
      *(coord++)= 2.0 * realr/(qh_RANDOMmax+1) - 1.0;
    }
  }
  *rowi= coord;
} /* randommatrix */

        
/*-------------------------------------------------
-rotateinput- rotate input using row matrix
  input points given by qh first_point, num_points, hull_dim
  if qh POINTSmalloc, overwrites input points, else mallocs a new array
  assumes rows[dim] is a scratch buffer
returns:
  sets qh POINTSmalloc
*/
void qh_rotateinput (realT **rows) {
  int size;
  pointT *newpoints;

  if (!qh POINTSmalloc) {
    size= qh num_points*qh hull_dim*sizeof(pointT);
    if (!(newpoints=(coordT*)malloc(size))) {
      fprintf(qh ferr, "qhull error: insufficient memory to rotate %d points\n",
          qh num_points);
      qh_errexit(qh_ERRmem, NULL, NULL);
    }
    memcpy ((char *)newpoints, (char *)qh first_point, size);
    qh first_point= newpoints;
    qh POINTSmalloc= True;
  }
  qh_rotatepoints (qh first_point, qh num_points, qh hull_dim, rows);
}  /* rotateinput */

/*-------------------------------------------------
-rotatepoints- rotate numpoints points by a row matrix
  assumes rows[dim] is a scratch buffer
*/
void qh_rotatepoints (realT *points, int numpoints, int dim, realT **row) {
  realT *point, *rowi, *coord= NULL, sum, *newval;
  int i,j,k;

  for (point= points, j= numpoints; j--; point += dim) {
    newval= row[dim];
    for (i= 0; i<dim; i++) {
      rowi= row[i];
      coord= point;
      for (sum= 0.0, k= dim; k--; )
        sum += *rowi++ * *coord++;
      *(newval++)= sum;
    }
    for (k= dim; k--; )
      *(--coord)= *(--newval);
  }
} /* rotatepoints */  
  

/*-------------------------------------------------
-scaleinput- scale input points using qh low_bound/high_bound
  input points given by qh first_point, num_points, hull_dim
  if qh POINTSmalloc, overwrites input points, else mallocs a new array
returns:
  scales points to low[k], high[k]
  sets qh POINTSmalloc
*/
void qh_scaleinput (void) {
  int size;
  pointT *newpoints;

  if (!qh POINTSmalloc) {
    size= qh num_points*qh hull_dim*sizeof(pointT);
    if (!(newpoints=(coordT*)malloc(size))) {
      fprintf(qh ferr, "qhull error: insufficient memory to scale %d points\n",
          qh num_points);
      qh_errexit(qh_ERRmem, NULL, NULL);
    }
    memcpy ((char *)newpoints, (char *)qh first_point, size);
    qh first_point= newpoints;
    qh POINTSmalloc= True;
  }
  qh_scalepoints (qh first_point, qh num_points, qh hull_dim,
       qh lower_bound, qh upper_bound);
}  /* scaleinput */
  
/*-------------------------------------------------
-scalepoints- scale points to new lowbound and highbound
  retains old bound when newlow= -REALmax or newhigh= +REALmax
  overwrites old points
*/
void qh_scalepoints (pointT *points, int numpoints, int dim,
	realT *newlows, realT *newhighs) {
  int i,k;
  realT shift, scale, *coord, low, high, newlow, newhigh, mincoord, maxcoord;
  boolT nearzero= False;
     
  for (k= 0; k<dim; k++) {
    newhigh= newhighs[k];
    newlow= newlows[k];
    if (newhigh > REALmax/2 && newlow < -REALmax/2)
      continue;
    low= REALmax;
    high= -REALmax;
    for (i= numpoints, coord= points+k; i--; coord += dim) {
      minimize_(low, *coord);
      maximize_(high, *coord);
    }
    if (newhigh > REALmax/2)
      newhigh= high;
    if (newlow < -REALmax/2)
      newlow= low;
    scale= qh_divzero (newhigh - newlow, high - low,
                  qh MINdenom_1, &nearzero);
    if (nearzero) {
      fprintf (qh ferr, "qhull input error: %d'th dimension's new bounds [%2.2g, %2.2g] too wide for\nexisting bounds [%2.2g, %2.2g]\n",
              k, newlow, newhigh, low, high);
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
    shift= (newlow * high - low * newhigh)/(high-low);
    coord= points+k;
    for (i= numpoints; i--; coord += dim)
      *coord= *coord * scale + shift;
    coord= points+k;
    if (newlow < newhigh) {
      mincoord= newlow;
      maxcoord= newhigh;
    }else {
      mincoord= newhigh;
      maxcoord= newlow;
    }
    for (i= numpoints; i--; coord += dim) {
      minimize_(*coord, maxcoord);  /* because of roundoff error */
      maximize_(*coord, mincoord);
    }
    trace0((qh ferr, "qh_scalepoints: scaled %d'th coordinate [%2.2g, %2.2g] to [%.2g, %.2g] for %d points by %2.2g and shifted %2.2g\n",
      k, low, high, newlow, newhigh, numpoints, scale, shift));
  }
} /* scalepoints */    

       
/*-------------------------------------------
-sethalfspace- set coords to dual of halfspace relative to feasible point
  the halfspace is its normal coefficients and offset.
returns:
  false if feasible point is outside of hull (error message reported)
  next value for coords
*/
boolT qh_sethalfspace (int dim, coordT *coords, coordT **nextp, 
         coordT *normal, coordT *offset, coordT *feasible) {
  coordT *normp= normal, *feasiblep= feasible, *coordp= coords;
  realT dist;
  realT r; /*bug fix*/
  int k;
  boolT zerodiv;

  dist= *offset;
  for (k= dim; k--; )
    dist += *(normp++) * *(feasiblep++);
  if (dist > 0)
    goto LABELerroroutside;
  normp= normal;
  if (dist < -qh MINdenom) {
    for (k= dim; k--; )
      *(coordp++)= *(normp++) / -dist;
  }else {
    for (k= dim; k--; ) {
      *(coordp++)= qh_divzero (*(normp++), -dist, qh MINdenom_1, &zerodiv);
      if (zerodiv) 
        goto LABELerroroutside;
    }
  }
  *nextp= coordp;
  if (qh IStracing >= 4) {
    fprintf (qh ferr, "qh_sethalfspace: halfspace at offset %6.2g to point: ", *offset);
    for (k= dim, coordp= coords; k--; )
      fprintf (qh ferr, " %6.2g", r=*coordp++);
    fprintf (qh ferr, "\n");
  }
  return True;
LABELerroroutside:
  feasiblep= feasible;
  normp= normal;
  fprintf(qh ferr, "qhull input error: feasible point is not clearly inside halfspace\nfeasible point: ");
  for (k= dim; k--; )
    fprintf (qh ferr, qh_REAL_1, r=*(feasiblep++));
  fprintf (qh ferr, "\n     halfspace: "); 
  for (k= dim; k--; )
    fprintf (qh ferr, qh_REAL_1, r=*(normp++));
  fprintf (qh ferr, "\n     at offset: ");
  fprintf (qh ferr, qh_REAL_1, *offset);
  fprintf (qh ferr, " and distance: ");
  fprintf (qh ferr, qh_REAL_1, dist);
  fprintf (qh ferr, "\n");
  return False;
} /* sethalfspace */

/*-------------------------------------------------
-sethalfspace_all- generate dual for halfspace intersection with feasible point
     each halfspace is normal coefficients followed by offset 
     the origin is inside the halfspace if the offset is negative
returns:
  unused/untested code: please email barber@tiac.net if this works ok for you
  malloc'd array of count X dim-1 points
  call before qh_init_B or qh_initqhull_globals 
notes:
  If using option 'Fp', also set qh feasible_point. It is a malloc'd array that
  is freed by qh_freebuffers.
*/
coordT *qh_sethalfspace_all (int dim, int count, coordT *halfspaces, pointT *feasible) {
  int i, newdim;
  pointT *newpoints;
  coordT *coordp, *normalp, *offsetp;

  trace0((qh ferr, "qh_sethalfspace_all: compute dual for halfspace intersection\n"));
  newdim= dim - 1;
  if (!(newpoints=(coordT*)malloc(count*newdim*sizeof(coordT)))){
    fprintf(qh ferr, "qhull error: insufficient memory to compute dual of %d halfspaces\n",
          count);
    qh_errexit(qh_ERRmem, NULL, NULL);
  }
  coordp= newpoints;
  normalp= halfspaces;
  for (i= 0; i < count; i++) {
    offsetp= normalp + newdim;
    if (!qh_sethalfspace (newdim, coordp, &coordp, normalp, offsetp, feasible)) {
      fprintf (qh ferr, "The halfspace was at index %d\n", i);
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
    normalp= offsetp + 1;
  }
  return newpoints;
} /* sethalfspace_all */

  
/*-------------------------------------------
-voronoi_center- return Voronoi center for a set of points
  dim is the orginal dimension of the points
notes:
  if non-simplicial, returns center for max simplex of points
  from Bowyer & Woodwark, A Programmer's Geometry, 1983, p. 65
*/
pointT *qh_voronoi_center (int dim, setT *points) {
  pointT *point, **pointp, *point0;
  pointT *center= (pointT*)qh_memalloc (qh center_size);
  setT *simplex;
  int i, j, k, num, size= qh_setsize(points);
  coordT *gmcoord;
  realT *diffp, sum2, *sum2row, *sum2p, det, factor;
  boolT nearzero, infinite;

  if (size == dim+1)
    simplex= points;
  else if (size < dim+1) {
    fprintf (qh ferr, "qhull internal error (qh_voronoi_center):\n  need at least %d points to construct a Voronoi center\n",
	     dim+1);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }else {
    simplex= qh_settemp (dim+1);
    qh_maxsimplex (dim, points, NULL, 0, &simplex);
  }
  num= qh_setsize (simplex);
  point0= SETfirst_(simplex);
  gmcoord= qh gm_matrix;
  for (k=0; k<dim; k++) {
    qh gm_row[k]= gmcoord;
    FOREACHpoint_(simplex) {
      if (point != point0)
        *(gmcoord++)= point[k] - point0[k];
    }
  }
  sum2row= gmcoord;
  for (i=0; i<dim; i++) {
    sum2= 0.0;
    for (k= 0; k<dim; k++) {
      diffp= qh gm_row[k] + i;
      sum2 += *diffp * *diffp;
    }
    *(gmcoord++)= sum2;
  }
  det= qh_determinant (qh gm_row, dim, &nearzero);
  factor= qh_divzero (0.5, det, qh MINdenom, &infinite);
  if (infinite) {
    for (k=dim; k--; )
      center[k]= qh_INFINITE;
    if (qh IStracing)
      qh_printpoints (qh ferr, "qh_voronoi_center: at infinity for ", simplex);
  }else {
    for (i=0; i<dim; i++) {
      gmcoord= qh gm_matrix;
      sum2p= sum2row;
      for (k=0; k<dim; k++) {
	qh gm_row[k]= gmcoord;
	if (k == i) {
	  for (j= dim; j--; )
	    *(gmcoord++)= *sum2p++;
	}else {
	  FOREACHpoint_(simplex) {
	    if (point != point0)
	      *(gmcoord++)= point[k] - point0[k];
	  }
	}
      }
      center[i]= qh_determinant (qh gm_row, dim, &nearzero)*factor + point0[i];
    }
#ifndef qh_NOtrace
    if (qh IStracing >= 3) {
      fprintf (qh ferr, "qh_voronoi_center: det %2.2g factor %2.2g ", det, factor);
      qh_printmatrix (qh ferr, "center:", &center, 1, dim);
      if (qh IStracing >= 5) {
	qh_printpoints (qh ferr, "points", simplex);
	FOREACHpoint_(simplex)
	  fprintf (qh ferr, "p%d dist %.2g, ", qh_pointid (point),
		   qh_pointdist (point, center, dim));
	fprintf (qh ferr, "\n");
      }
    }
#endif
  }
  if (simplex != points)
    qh_settempfree (&simplex);
  return center;
} /* voronoi_center */

