/* qhull - Quickhull algorithm for convex hulls

   qhull() and top-level routines

   see README, qhull.h, unix.c and mac.c

   see qhull_a.h for internal functions
   
   copyright (c) 1993-1995 The Geometry Center        
*/

#include "qhull_a.h" 


/*-------------------------------------------------
-qhull- hull_dim convex hull of num_points starting at first_point
returns:
  returns facet_list, numfacets, etc. 
*/
void qh_qhull (void) {
  setT *maxpoints, *vertices;
  facetT *facet;
  int numpart, i, numoutside;
  realT dist;
  boolT isoutside;

  qh hulltime= (unsigned)clock();
  if (qh DELAUNAY && qh upper_threshold[qh hull_dim-1] > REALmax/2
                  && qh lower_threshold[qh hull_dim-1] < -REALmax/2) {
    for (i= qh_PRINTEND; i--; ) {
      if (qh PRINTout[i] == qh_PRINTgeom && qh DROPdim < 0 
 	  && !qh GOODthreshold && !qh SPLITthresholds)
	break;  /* in this case, don't set upper_threshold */
    }
    if (i < 0) {
      qh upper_threshold[qh hull_dim-1]= 0.0;
      if (!qh GOODthreshold)
	qh SPLITthresholds= True;
    }
  }
  maxpoints= qh_maxmin(qh first_point, qh num_points, qh hull_dim);
  /* qh_maxmin sets DISTround and other precision constants */
  if (qh PRINToptions1st || qh TRACElevel || qh IStracing) {
    if (qh TRACElevel || qh IStracing)
      fprintf (qh ferr, "\nTrace level %d for %s | %s\n", 
         qh IStracing ? qh IStracing : qh TRACElevel, qh rbox_command, qh qhull_command);
    fprintf (qh ferr, "Options selected for qhull %s:\n%s\n", qh_version, qh qhull_options);
  }
  vertices= qh_initialvertices(qh hull_dim, maxpoints, qh first_point, qh num_points); 
  qh_initialhull (vertices);  /* initial qh facet_list */
  qh_partitionall (vertices, qh first_point, qh num_points);
  if (qh PREmerge) {
    qh cos_max= qh premerge_cos;
    qh centrum_radius= qh premerge_centrum;
  }
  if (qh ONLYgood) {
    if (!(qh GOODthreshold || qh GOODpoint
	  || (qh GOODvertex > 0 && !qh MERGING))) {
      fprintf (qh ferr, "qhull input error: 'Qg' (ONLYgood) needs a good threshold ('Pd0D0'), a\n\
good point (QGn or QG-n), or a good vertex without merging (QVn).\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
    if (qh GOODvertex > 0  && !qh MERGING  /* matches qh_partitionall */
	&& !qh_isvertex (qh GOODvertexp, vertices)) {
      facet= qh_findbestnew (qh GOODvertexp, qh facet_list, 
			  &dist, &isoutside, &numpart);
      zadd_(Zdistgood, numpart);
      if (!isoutside) {
        fprintf (qh ferr, "qhull input error: point for QV%d is inside initial simplex\n",
	       qh_pointid(qh GOODvertexp));
        qh_errexit (qh_ERRinput, NULL, NULL);
      }
      if (!qh_addpoint (qh GOODvertexp, facet, False)) {
	qh_settempfree(&vertices);
	qh_settempfree(&maxpoints);
	return;
      }
    }
    qh_findgood (qh facet_list, 0);
  }
  qh_settempfree(&vertices);
  qh_settempfree(&maxpoints);
  qh_resetlists (False /*qh visible_list newvertex_list newfacet_list */);
  qh_buildhull();
  if (!qh STOPpoint && !qh STOPcone) {
    if (qh ZEROall_ok && !qh TESTvneighbors && qh MERGEexact)
      qh_checkzero( qh_ALL);
    if (qh ZEROall_ok && !qh TESTvneighbors) {
      trace2((qh ferr, "qh_qhull: all facets are clearly convex.  Post-merging not needed.\n"));
    }else {
      if (qh MERGEexact || (qh hull_dim > qh_DIMreduceBuild && qh PREmerge))
        qh_postmerge ("First post-merge", qh premerge_centrum, qh premerge_cos, 
             (qh POSTmerge ? False : qh TESTvneighbors));
      else if (!qh POSTmerge && qh TESTvneighbors) 
        qh_postmerge ("For testing vertex neighbors", qh premerge_centrum,
             qh premerge_cos, True); 
      if (qh POSTmerge)
        qh_postmerge ("For post-merging", qh postmerge_centrum, 
             qh postmerge_cos, qh TESTvneighbors);
      if (qh visible_list == qh facet_list) { /* i.e., merging done */
        qh findbestnew= True;
        qh_partitionvisible (/*visible_list, newfacet_list*/ !qh_ALL, &numoutside);
        qh findbestnew= False;
        qh_deletevisible (/*qh visible_list*/);
        qh_resetlists (False /*qh visible_list newvertex_list newfacet_list */);
      }
      if (qh DOcheckmax){
        if (qh REPORTfreq) {
	  qh_buildtracing (NULL, NULL); 
  	  fprintf (qh ferr, "\nTesting all coplanar points.\n");
        }
        qh_check_maxout();
      }
    }
  }
  if (qh_setsize (qhmem.tempstack) != 0) {
    fprintf (qh ferr, "qhull internal error (qh_qhull): temporary sets not empty (%d)\n",
	     qh_setsize (qhmem.tempstack));
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  qh hulltime= (unsigned)clock() - qh hulltime;
  qh QHULLfinished= True;
  trace1((qh ferr, "qh_qhull: algorithm completed\n"));
} /* qhull */

/*-------------------------------------------------
-addpoint-  add point to hull above a facet
  if !facet, locates a facet for the point
    !facet works for Delaunay triangulations
    !facet does not work for lens-shaped hulls
    if point is not outside of the hull, uses qh_partitioncoplanar()
  if checkdist, checks that point is outside of the facet.
    if not outside, uses qh_partitioncoplanar()
    if facet->upperdelaunay, assumes that the point is above facet
  if !checkdist and facet, assumes point is above facet (major damage if below)
returns:
  if unknown point, adds it to qh other_points
  returns False if user requested break
     visible_list, newfacet_list, delvertex_list, NEWfacets may be defined
*/
boolT qh_addpoint (pointT *furthest, facetT *facet, boolT checkdist) {
  int goodvisible, goodhorizon;
  vertexT *vertex;
  facetT *newfacet;
  realT dist, newbalance, pbalance;
  boolT isoutside= False;
  int numpart, numpoints, numnew, firstnew;

  if (qh_pointid (furthest) == -1)
    qh_setappend (&qh other_points, furthest);
  if (!facet 
  || (checkdist && !facet->upperdelaunay)) { /* else missed by qh_findbest */
    if (!facet)
      facet= qh_nonupper( qh facet_list);
    facet= qh_findbest (furthest, facet, !qh_ALL, False, &dist, &isoutside, &numpart);
    zzadd_(Zpartition, numpart);
    if (!isoutside) {
      zinc_(Znotmax);  /* last point of outsideset is no longer furthest. */
      qh_partitioncoplanar (furthest, facet, &dist);
      return True;
    }
  }
  qh_buildtracing (furthest, facet);
  if (qh STOPpoint < 0 && qh furthest_id == -qh STOPpoint-1)
    return False;
  qh_findhorizon (furthest, facet, &goodvisible, &goodhorizon); 
  if (qh ONLYgood && !(goodvisible+goodhorizon)) {
    zinc_(Znotgood);  
    /* last point of outsideset is no longer furthest.  This is ok
       since all points of the outside are likely to be bad */
    qh_resetlists (False /*qh visible_list newvertex_list newfacet_list */);
    return True;
  }
  zzinc_(Zprocessed);
  firstnew= qh facet_id;
  vertex= qh_makenewfacets (furthest /*visible_list, attaches if !ONLYgood */);
  qh_makenewplanes (/* newfacet_list */);
  numnew= qh facet_id - firstnew;
  newbalance= numnew - (realT) (qh num_facets-qh num_visible)
                         * qh hull_dim/qh num_vertices;
  wadd_(Wnewbalance, newbalance);
  wadd_(Wnewbalance2, newbalance * newbalance);
  if (qh ONLYgood && !qh_findgood (qh newfacet_list, goodhorizon)) {
    FORALLnew_facets 
      qh_delfacet (newfacet);
    qh_delvertex (vertex);
    qh_resetlists (True /*qh visible_list newvertex_list newfacet_list */);
    zinc_(Znotgoodnew);
    return True;
  }
  if (qh ONLYgood)
    qh_attachnewfacets(/*visible_list*/);
  qh_matchnewfacets();
  qh_updatevertices();
  if (qh STOPcone && qh furthest_id == qh STOPcone-1)
    return False;  /* visible_list etc. still defined */
  if (qh PREmerge || qh MERGEexact) {
    qh_premerge (vertex, qh premerge_centrum, qh premerge_cos);
    if (zzval_(Ztotmerge) > qh_USEfindbestnew)
      qh findbestnew= True;
    else {
      FORALLnew_facets {
	if (!newfacet->simplicial) {
	  qh findbestnew= True;  /* qh_findbest can not be used */
	  break;
	}
      }
    }
  }else if (qh BESToutside)
    qh findbestnew= True;
  qh_partitionvisible (/*visible_list, newfacet_list*/ !qh_ALL, &numpoints);
  qh findbestnew= False;
  qh findbest_notsharp= False;
  zinc_(Zpbalance);
  pbalance= numpoints - (realT) qh hull_dim /* assumes all points extreme */
                * (qh num_points - qh num_vertices)/qh num_vertices;
  wadd_(Wpbalance, pbalance);
  wadd_(Wpbalance2, pbalance * pbalance);
  qh_deletevisible (/*qh visible_list*/);
  zmax_(Zmaxvertex, qh num_vertices);
  qh NEWfacets= False;
  if (qh IStracing >= 4)
    qh_printfacetlist (qh newfacet_list, NULL, True);
  if (qh CHECKfrequently) {
    if (qh num_facets < 50)
      qh_checkpolygon (qh facet_list);
    else
      qh_checkpolygon (qh newfacet_list);
  }
  if (qh STOPpoint > 0 && qh furthest_id == qh STOPpoint-1)
    return False; 
  qh_resetlists (True /*qh visible_list newvertex_list newfacet_list */);
  trace2((qh ferr, "qh_addpoint: added p%d new facets %d new balance %2.2g point balance %2.2g\n",
    qh_pointid (furthest), numnew, newbalance, pbalance));
  return True;
} /* addpoint */

/*-------------------------------------------------
-buildhull- constructs a hull by adding outside points one at a time
  may be called multiple times
  checks facet and vertex lists for 'visible', 'newfacet', and 'newlist'
notes:
  to recover from STOPcone, call qh_deletevisible and qh_resetlists
*/
void qh_buildhull(void) {
  facetT *facet;
  pointT *furthest;
  vertexT *vertex;
  int id;
  
  trace1((qh ferr, "qh_buildhull: start build hull\n"));
  FORALLfacets {
    if (facet->visible || facet->newfacet) {
      fprintf (qh ferr, "qhull internal error (qh_buildhull): visible or new facet f%d in facet list\n",
                   facet->id);    
      qh_errexit (qh_ERRqhull, facet, NULL);
    }
  }
  FORALLvertices {
    if (vertex->newlist) {
      fprintf (qh ferr, "qhull internal error (qh_buildhull): new vertex f%d in vertex list\n",
                   vertex->id);
      qh_errprint ("ERRONEOUS", NULL, NULL, NULL, vertex);
      qh_errexit (qh_ERRqhull, NULL, NULL);
    }
    id= qh_pointid (vertex->point);
    if ((qh STOPpoint>0 && id == qh STOPpoint-1) ||
	(qh STOPpoint<0 && id == -qh STOPpoint-1) ||
	(qh STOPcone>0 && id == qh STOPcone-1)) {
      trace1((qh ferr,"qh_buildhull: stop point or cone P%d in initial hull\n", id));
      return;
    }
  }
  qh facet_next= qh facet_list;      /* advance facet when processed */
  while ((furthest= qh_nextfurthest (&facet))) {
    qh num_outside--;  /* if ONLYmax, furthest may not be outside */
    if (!qh_addpoint (furthest, facet, qh ONLYmax))
      break;
  }
  if (qh num_outside && !furthest) {
    fprintf (qh ferr, "qhull internal error (qh_buildhull): %d outside points were never processed.\n", qh num_outside);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  trace1((qh ferr, "qh_buildhull: completed the hull construction\n"));
} /* buildhull */
  

/*-------------------------------------------
-buildtracing- for tracing execution of buildhull
  tracks progress with qh lastreport
  updates qh furthest_id (-3 if furthest is NULL)
  also resets visit_id, vertext_visit on wrap around
  if !furthest, prints basic message
  see also qh_tracemerging()
*/
void qh_buildtracing (pointT *furthest, facetT *facet) {
  realT dist= 0;
  float cpu;
  int total, furthestid;
  time_t timedata;
  struct tm *tp;
  vertexT *vertex;

  if (!furthest) {
    time (&timedata);
    tp= localtime (&timedata);
    cpu= (unsigned)clock() - qh hulltime;
    cpu /= qh_SECticks;
    total= zzval_(Ztotmerge) - zzval_(Zcyclehorizon) + zzval_(Zcyclefacettot);
    fprintf (qh ferr, "\n\
At %02d:%02d:%02d & %2.5g CPU secs, qhull has created %d facets and merged %d.\n\
 The current hull contains %d facets and %d vertices.  Last point was p%d\n",
      tp->tm_hour, tp->tm_min, tp->tm_sec, cpu, qh facet_id -1,
      total, qh num_facets, qh num_vertices, qh furthest_id);
    return;
  }
  furthestid= qh_pointid (furthest);
  if (qh TRACEpoint == furthestid) {
    qh IStracing= qh TRACElevel;
    qhmem.IStracing= qh TRACElevel;
  }
  if (qh REPORTfreq && (qh facet_id-1 > qh lastreport+qh REPORTfreq)) {
    qh lastreport= qh facet_id-1;
    time (&timedata);
    tp= localtime (&timedata);
    cpu= (unsigned)clock() - qh hulltime;
    cpu /= qh_SECticks;
    total= zzval_(Ztotmerge) - zzval_(Zcyclehorizon) + zzval_(Zcyclefacettot);
    zinc_(Zdistio);
    qh_distplane (furthest, facet, &dist);
    fprintf (qh ferr, "\n\
At %02d:%02d:%02d & %2.5g CPU secs, qhull has created %d facets and merged %d.\n\
 The current hull contains %d facets and %d vertices.  There are %d\n\
 outside points.  Next is point p%d (v%d), %2.2g above f%d.\n",
      tp->tm_hour, tp->tm_min, tp->tm_sec, cpu, qh facet_id -1,
      total, qh num_facets, qh num_vertices, qh num_outside+1,
      furthestid, qh vertex_id, dist, getid_(facet));
  }else if (qh IStracing >=1) {
    cpu= (unsigned)clock() - qh hulltime;
    cpu /= qh_SECticks;
    qh_distplane (furthest, facet, &dist);
    fprintf (qh ferr, "qh_buildhull: add p%d (v%d) to hull of %d facets (%2.2g above f%d) and %d outside at %4.4g CPU secs.  Previous was p%d.\n",
      furthestid, qh vertex_id, qh num_facets, dist,
      getid_(facet), qh num_outside+1, cpu, qh furthest_id);
  }
  if (qh visit_id > (unsigned) INT_MAX) {
    qh visit_id= 0;
    FORALLfacets
      facet->visitid= qh visit_id;
  }
  if (qh vertex_visit > (unsigned) INT_MAX) {
    qh vertex_visit= 0;
    FORALLvertices
      vertex->visitid= qh vertex_visit;
  }
  qh furthest_id= furthestid;
} /* buildtracing */

/*-------------------------------------------
-errexit2- return exitcode to system after an error
  assumes exitcode non-zero
  for two facets, see qh_errexit() in user.c
*/
void qh_errexit2(int exitcode, facetT *facet, facetT *otherfacet) {
  
  qh_errprint("ERRONEOUS", facet, otherfacet, NULL, NULL);
  qh_errexit (exitcode, NULL, NULL);
} /* errexit2 */


/*-------------------------------------------------
-findhorizon- given a visible facet, find the point's horizon and visible facets
returns:
  qh visible_list to all visible facets 
    marks visible facets with ->visible 
    goodvisible counts visible->good
    initializes num_visible
notes:
  similar to delpoint()
*/
void qh_findhorizon(pointT *point, facetT *facet, int *goodvisible, int *goodhorizon) {
  facetT *neighbor, **neighborp, *visible;
  int numhorizon= 0, coplanar= 0;
  realT dist;
  
  trace1((qh ferr,"qh_findhorizon: find horizon for point p%d facet f%d\n",qh_pointid(point),facet->id));
  *goodvisible= *goodhorizon= 0;
  zinc_(Ztotvisible);
  qh_removefacet(facet);  /* visible_list at end of qh facet_list */
  qh_appendfacet(facet);
  qh num_visible= 1;
  if (facet->good)
    (*goodvisible)++;
  qh visible_list= facet;
  facet->visible= True;
  facet->f.replace= NULL;
  if (qh IStracing >=4)
    qh_errprint ("visible", facet, NULL, NULL, NULL);
  qh visit_id++;
  FORALLvisible_facets {
    visible->visitid= qh visit_id;
    FOREACHneighbor_(visible) {
      if (neighbor->visitid == qh visit_id) 
        continue;
      neighbor->visitid= qh visit_id;
      zzinc_(Znumvisibility);
      qh_distplane(point, neighbor, &dist);
      if (dist > qh MINvisible) {
        zinc_(Ztotvisible);
	qh_removefacet(neighbor);  /* append to end of qh visible_list */
	qh_appendfacet(neighbor);
	neighbor->visible= True;
        neighbor->f.replace= NULL;
	qh num_visible++;
	if (neighbor->good)
	  (*goodvisible)++;
        if (qh IStracing >=4)
          qh_errprint ("visible", neighbor, NULL, NULL, NULL);
      }else {
 	if (dist > - qh MAXcoplanar) {
    	  neighbor->coplanar= True;
          zzinc_(Zcoplanarhorizon);
	  coplanar++;
	  if (qh MERGING) {
	    if (dist > 0) {
	      maximize_(qh max_outside, dist);
	      maximize_(qh max_vertex, dist);
#if qh_MAXoutside
	      maximize_(neighbor->maxoutside, dist);
#endif
	    }else
	      minimize_(qh min_vertex, dist);  /* due to merge later */
	  }
      	  trace2((qh ferr, "qh_findhorizon: point p%d is coplanar to horizon f%d, dist=%2.7g < qh MINvisible (%2.7g)\n",
	      qh_pointid(point), neighbor->id, dist, qh MINvisible));
	}else
    	  neighbor->coplanar= False;
    	zinc_(Ztothorizon);
        numhorizon++;
	if (neighbor->good)
	  (*goodhorizon)++;
        if (qh IStracing >=4)
          qh_errprint ("horizon", neighbor, NULL, NULL, NULL);
      }
    }
  }
  if (!numhorizon) {
    fprintf(qh ferr, "qhull precision error (qh_findhorizon): empty horizon\n\
Point p%d was above all facets.\n", qh_pointid(point));
    qh_printfacetlist (qh facet_list, NULL, True);
    qh_errexit(qh_ERRprec, NULL, NULL);
  }
  trace1((qh ferr, "qh_findhorizon: %d horizon facets (good %d), %d visible (good %d), %d coplanar\n", 
       numhorizon, *goodhorizon, qh num_visible, *goodvisible, coplanar));
  if (qh IStracing >= 4 && qh num_facets < 50) 
    qh_printlists ();
} /* findhorizon */


/*------------------------------------------------
-nextfurthest- returns next furthest point for processing
returns:
  NULL if none available
  visible facet for furthest
  removes empty outside sets  
*/
pointT *qh_nextfurthest (facetT **visible) {
  facetT *facet;
  int size, index;
  realT randr;
  pointT *furthest;

  while ((facet= qh facet_next) != qh facet_tail) {
    if (!facet->outsideset) {
      qh facet_next= facet->next;
      continue;
    }
    SETreturnsize_(facet->outsideset, size);
    if (!size) {
      qh_setfree (&facet->outsideset);
      qh facet_next= facet->next;
      continue;
    }
    if (!qh RANDOMoutside && !qh VIRTUALmemory) {
      *visible= facet;
      return ((pointT*)qh_setdellast (facet->outsideset));
    }
    if (qh RANDOMoutside) {
      randr= qh_RANDOMint;
      randr= randr/(qh_RANDOMmax+1);
      index= floor(qh num_outside * randr);
      FORALLfacet_(qh facet_next) {
        if (facet->outsideset) {
          SETreturnsize_(facet->outsideset, size);
          if (!size)
            qh_setfree (&facet->outsideset);
          else if (size > index) {
            *visible= facet;
            return ((pointT*)qh_setdelnth (facet->outsideset, index));
          }else
            index -= size;
        }
      }
      fprintf (qh ferr, "qhull internal error (qh_nextfurthest): num_outside %d incorrect or random %2.2g >= 1.0\n",
              qh num_outside, randr);
      qh_errexit (qh_ERRqhull, NULL, NULL);
    }else { /* VIRTUALmemory */
      facet= qh facet_tail->previous;
      if (!(furthest= (pointT*)qh_setdellast(facet->outsideset))) {
        if (facet->outsideset)
          qh_setfree (&facet->outsideset);
        qh_removefacet (facet);
        qh_prependfacet (facet, &qh facet_list);
        continue;
      }
      *visible= facet;
      return furthest;
    }
  }
  return NULL;
} /* nextfurthest */

/*-------------------------------------------------
-partitionall- partitions all points into the outsidesets of facets
   vertices= set of vertices used by qh facet_list
     does not partition qh GOODpoint
     if ONLYgood && !MERGING, does not partition GOODvertex
   all facets have ->newfacet for qh_findbestnew in qh_partitionpoint
notes:
   faster if qh facet_list sorted by anticipated size of outside set
*/
void qh_partitionall(setT *vertices, pointT *points, int numpoints){
  setT *pointset;
  vertexT *vertex, **vertexp;
  pointT *point, **pointp, *bestpoint;
  int size, point_i, point_n, point_end, remaining, i, id;
  facetT *facet;
  realT bestdist= -REALmax, dist, distoutside;
    
  trace1((qh ferr, "qh_partitionall: partition all points into outside sets\n"));
  pointset= qh_settemp (numpoints);
  pointp= SETaddr_(pointset, pointT);
  for (i=numpoints, point= points; i--; point += qh hull_dim)
    *(pointp++)= point;
  qh_settruncate (pointset, numpoints);
  FOREACHvertex_(vertices) {
    if ((id= qh_pointid(vertex->point)) >= 0)
      SETelem_(pointset, id)= NULL;
  }
  id= qh_pointid (qh GOODpointp);
  if (id >=0 && qh STOPcone-1 != id && -qh STOPpoint-1 != id)
    SETelem_(pointset, id)= NULL;
  if (qh GOODvertexp && qh ONLYgood && !qh MERGING) { /* matches qhull()*/
    if ((id= qh_pointid(qh GOODvertexp)) >= 0)
      SETelem_(pointset, id)= NULL;
  }
  if (!qh BESToutside) {  /* matches conditional for qh_partitionpoint below */
    if (qh MERGING)
      distoutside= qh_DISToutside; /* defined in user.h */
    else
      distoutside= qh MINoutside;
    zval_(Ztotpartition)= qh num_points - qh hull_dim - 1; /*misses GOOD... */
    remaining= qh num_facets;
    point_end= numpoints;
    FORALLfacets {
      size= point_end/(remaining--) + 100;
      facet->outsideset= qh_setnew (size);
      bestpoint= NULL;
      point_end= 0;
      FOREACHpoint_i_(pointset) {
        if (point) {
          zzinc_(Zpartitionall);
          qh_distplane (point, facet, &dist);
          if (dist < distoutside)
            SETelem_(pointset, point_end++)= point;
          else {
	    qh num_outside++;
            if (!bestpoint) {
              bestpoint= point;
              bestdist= dist;
            }else if (dist > bestdist) {
              qh_setappend (&facet->outsideset, bestpoint);
              bestpoint= point;
              bestdist= dist;
            }else 
              qh_setappend (&facet->outsideset, point);
          }
        }
      }
      if (bestpoint) {
        qh_setappend (&facet->outsideset, bestpoint);
#if !qh_COMPUTEfurthest
	facet->furthestdist= bestdist;
#endif
      }else
        qh_setfree (&facet->outsideset);
      qh_settruncate (pointset, point_end);
    }
  }
  if (qh BESToutside || qh MERGING || qh KEEPcoplanar || qh KEEPinside) {
    qh findbestnew= True;
    FOREACHpoint_i_(pointset) { 
      if (point)
        qh_partitionpoint(point, qh facet_list);
    }
    qh findbestnew= False;
  }
  zzadd_(Zpartitionall, zzval_(Zpartition));
  zzval_(Zpartition)= 0;
  qh_settempfree(&pointset);
  if (qh IStracing >= 4)
    qh_printfacetlist (qh facet_list, NULL, True);
} /* partitionall */


/*-------------------------------------------------
-partitioncoplanar- partition coplanar point to a facet
  if dist NULL, searches for bestfacet, and does nothing if inside
      if qh findbestnew, searches new facets instead of findbest
returns:
  max_ouside updated
  facet->maxoutside's updated at end by qh_check_maxout
  if KEEPcoplanar or KEEPinside
    point assigned to best coplanarset
*/
void qh_partitioncoplanar (pointT *point, facetT *facet, realT *dist) {
  facetT *bestfacet;
  pointT *oldfurthest;
  realT bestdist, dist2;
  int numpart= 0;
  boolT isoutside, istrace= False;

  if (!dist) {
    if (qh findbestnew)
      bestfacet= qh_findbestnew (point, facet, 
			  &bestdist, NULL, &numpart);
    else
      bestfacet= qh_findbest (point, facet, qh_ALL, False, 
                          &bestdist, &isoutside, &numpart);
    zinc_(Ztotpartcoplanar);
    zzadd_(Zpartcoplanar, numpart);
    if (qh KEEPnearinside) {
      if (bestdist < -qh NEARinside) { 
        zinc_(Zcoplanarinside);
        return;
      }
    }else if (!qh KEEPinside && bestdist < -qh MAXcoplanar) {
      zinc_(Zcoplanarinside);
      return;
    }
  }else {
    bestfacet= facet;
    bestdist= *dist;
  }
  if (qh KEEPcoplanar + qh KEEPinside + qh KEEPnearinside) {
    oldfurthest= (pointT*)qh_setlast (bestfacet->coplanarset);
    if (oldfurthest) {
      zinc_(Zcomputefurthest);
      qh_distplane (oldfurthest, bestfacet, &dist2);
    }
    if (!oldfurthest || dist2 < bestdist) {
      qh_setappend(&bestfacet->coplanarset, point);
      if (bestdist > qh max_outside) {
	qh max_outside= bestdist;
	if (bestdist > qh TRACEdist)
	  istrace= True;
      }
    }else
      qh_setappend2ndlast(&bestfacet->coplanarset, point);
  }else { /* !KEEPcoplanar && !KEEPinside */
    if (bestdist > qh max_outside) {
      qh max_outside= bestdist;
      if (bestdist > qh TRACEdist) 
	istrace= True;
    }
  }
  if (istrace) {
    fprintf (qh ferr, "qh_partitioncoplanar: ====== p%d increases max_outside to %2.2g of f%d last p%d\n",
		   qh_pointid(point), bestdist, bestfacet->id, qh furthest_id);
    qh_errprint ("DISTANT", bestfacet, NULL, NULL, NULL);
  }
  trace4((qh ferr, "qh_partitioncoplanar: point p%d is coplanar with facet f%d (or inside) dist %2.2g\n",
	  qh_pointid(point), bestfacet->id, bestdist));
} /* partitioncoplanar */


/*-------------------------------------------------
-partitionpoint- assigns point to a visible facet
  uses qh_findbest with qh BESToutside and newfacets
  qh findbestnew if search new facets instead of findbest()
notes:
  after qh_distplane, this and qh_findbest are most expensive in 3-d
*/
void qh_partitionpoint (pointT *point, facetT *facet) {
  realT bestdist;
  pointT *oldfurthest;
  boolT isoutside;
  facetT *bestfacet;
  int numpart;
#if qh_COMPUTEfurthest
  realT dist;
#endif

  if (qh findbestnew)
    bestfacet= qh_findbestnew (point, facet,
			  &bestdist, &isoutside, &numpart);
  else
    bestfacet= qh_findbest (point, facet, qh BESToutside, True,
			  &bestdist, &isoutside, &numpart);
  zinc_(Ztotpartition);
  zzadd_(Zpartition, numpart);
  if (isoutside) {
    if (!bestfacet->outsideset 
    || !(oldfurthest= (pointT*)qh_setlast (bestfacet->outsideset))) {
      qh_setappend(&(bestfacet->outsideset), point);
      if (!bestfacet->newfacet) {
        qh_removefacet (bestfacet);  /* make sure it's after qh facet_next */
        qh_appendfacet (bestfacet);
      }
#if !qh_COMPUTEfurthest
      bestfacet->furthestdist= bestdist;
#endif
    }else {
#if qh_COMPUTEfurthest
      zinc_(Zcomputefurthest);
      qh_distplane (oldfurthest, bestfacet, &dist);
      if (dist < bestdist) 
	qh_setappend(&(bestfacet->outsideset), point);
      else
	qh_setappend2ndlast(&(bestfacet->outsideset), point);
#else
      if (bestfacet->furthestdist < bestdist) {
	qh_setappend(&(bestfacet->outsideset), point);
	bestfacet->furthestdist= bestdist;
      }else
	qh_setappend2ndlast(&(bestfacet->outsideset), point);
#endif
    }
    qh num_outside++;
    trace4((qh ferr, "qh_partitionpoint: point p%d is outside facet f%d\n",
	  qh_pointid(point), bestfacet->id));
  }else if (bestdist >= -qh MAXcoplanar) {
    zzinc_(Zcoplanarpart);
    if (qh KEEPcoplanar + qh KEEPnearinside || bestdist > qh max_outside) 
      qh_partitioncoplanar (point, bestfacet, &bestdist);
  }else if (qh KEEPnearinside && bestdist > -qh NEARinside) {
    zinc_(Zpartnear);
    qh_partitioncoplanar (point, bestfacet, &bestdist);
  }else {
    zinc_(Zpartinside);
    trace4((qh ferr, "qh_partitionpoint: point p%d is inside all facets, closest to f%d dist %2.2g\n",
	  qh_pointid(point), bestfacet->id, bestdist));
    if (qh KEEPinside)	  
      qh_partitioncoplanar (point, bestfacet, &bestdist);
  }
} /* partitionpoint */

/*-------------------------------------------------
-partitionvisible- partitions points in visible_list to newfacet_list
  1st neighbor (if any) points to a horizon facet or a new facet
  repartitions coplanar points if allpoints (not used)
  qh findbestnew if search new facets instead of findbest()
  qh findbest_notsharp should be clear
*/
void qh_partitionvisible(/*visible_list*/ boolT allpoints, int *numoutside) {
  facetT *visible, *newfacet;
  pointT *point, **pointp;
  int coplanar=0, size, count;
  vertexT *vertex, **vertexp;
  
  if (qh ONLYmax)
    maximize_(qh MINoutside, qh max_vertex);
  *numoutside= 0;
  FORALLvisible_facets {
    if (!visible->outsideset && !visible->coplanarset)
      continue;
    newfacet= visible->f.replace;
    count= 0;
    while (newfacet && newfacet->visible) {
      newfacet= newfacet->f.replace;
      if (count++ > qh facet_id)
	qh_infiniteloop (visible);
    }
    if (!newfacet)
      newfacet= qh newfacet_list;
    if (visible->outsideset) {
      size= qh_setsize (visible->outsideset);
      *numoutside += size;
      qh num_outside -= size;
      FOREACHpoint_(visible->outsideset) 
        qh_partitionpoint (point, newfacet);
    }
    if (visible->coplanarset && (qh KEEPcoplanar + qh KEEPinside + qh KEEPnearinside)) {
      size= qh_setsize (visible->coplanarset);
      coplanar += size;
      FOREACHpoint_(visible->coplanarset) {
        if (allpoints)
          qh_partitionpoint (point, newfacet);
        else
          qh_partitioncoplanar (point, newfacet, NULL);
      }
    }
  }
  FOREACHvertex_(qh del_vertices) {
    if (vertex->point) {
      if (allpoints)
        qh_partitionpoint (vertex->point, qh newfacet_list);
      else
        qh_partitioncoplanar (vertex->point, qh newfacet_list, NULL);
    }
  }
  trace1((qh ferr,"qh_partitionvisible: partitioned %d points from outsidesets and %d points from coplanarsets\n", *numoutside, coplanar));
} /* partitionvisible */



/*----------------------------------------
-printsummary- prints the summary about the computation
  not in io.c so that user_eg.c can prevent io.c from loading
*/
void qh_printsummary(FILE *fp) {
  realT ratio, dist;
  float cpu;
  int size, id, total, numvertices, numcoplanars= 0;
  facetT *facet;

  size= qh num_points + qh_setsize (qh other_points);
  numvertices= qh num_vertices - qh_setsize (qh del_vertices);
  id= qh_pointid (qh GOODpointp);
  if (qh DELAUNAY) 
    numcoplanars= size - numvertices;
  else {
    FORALLfacets {
      if (facet->coplanarset)
	numcoplanars += qh_setsize( facet->coplanarset);
    }
  }
  if (id >=0 && qh STOPcone-1 != id && -qh STOPpoint-1 != id)
    size--;
  if (qh VORONOI)
    fprintf (fp, "\nVoronoi diagram by the convex");
  else if (qh DELAUNAY)
    fprintf (fp, "\nDelaunay triangulation by the convex");
  else if (qh HALFspace)
    fprintf (fp, "\nHalfspace intersection by the convex");
  else
    fprintf (fp, "\nConvex");
  fprintf(fp, " hull of %d points in %d-d:\n\n", size, qh hull_dim);
  fprintf(fp, "  Number of vertices%s: %d\n",
	  qh HALFspace ? " (halfspaces)" : "", numvertices);
  if (numcoplanars) 
    fprintf(fp, "  Number of %s points: %d\n", 
	    qh DELAUNAY ? "similar" : "coplanar", numcoplanars); 
  fprintf(fp, "  Number of facets%s: %d\n", 
	  qh HALFspace ? " (intersections)" : "",
	  qh num_facets - qh num_visible);
  if (qh num_good)
    fprintf(fp, "  Number of good facets: %d\n", qh num_good);
  fprintf(fp, "\nStatistics for: %s | %s", 
                      qh rbox_command, qh qhull_command);
  if (qh ROTATErandom > 0)
    fprintf(fp, " QR%d\n\n", qh ROTATErandom);
  else
    fprintf(fp, "\n\n");
  fprintf(fp, "  Number of points processed: %d\n", zzval_(Zprocessed));
  fprintf(fp, "  Number of hyperplanes created: %d\n", zzval_(Zsetplane));
  fprintf(fp, "  Number of distance tests for qhull: %d\n", zzval_(Zpartition)+
  zzval_(Zpartitionall)+zzval_(Znumvisibility)+zzval_(Zpartcoplanar));
#if 0  /* NOTE: must print before printstatistics() */
  {realT stddev, ave;
  fprintf(fp, "  average new facet balance: %2.2g\n",
	  wval_(Wnewbalance)/zval_(Zprocessed));
  stddev= qh_stddev (zval_(Zprocessed), wval_(Wnewbalance), 
                                 wval_(Wnewbalance2), &ave);
  fprintf(fp, "  new facet standard deviation: %2.2g\n", stddev);
  fprintf(fp, "  average partition balance: %2.2g\n",
	  wval_(Wpbalance)/zval_(Zpbalance));
  stddev= qh_stddev (zval_(Zpbalance), wval_(Wpbalance), 
                                 wval_(Wpbalance2), &ave);
  fprintf(fp, "  partition standard deviation: %2.2g\n", stddev);
  }
#endif
  if (qh MERGING) {
    total= zzval_(Ztotmerge) - zzval_(Zcyclehorizon) + zzval_(Zcyclefacettot);
    fprintf(fp,"  Number of merged facets: %d\n", total);
    fprintf(fp,"  Number of distance tests for merging: %d\n",zzval_(Zbestdist)+
          zzval_(Zcentrumtests)+zzval_(Zdistconvex)+zzval_(Zdistcheck)+
          zzval_(Zdistzero));
  }
  if (!qh RANDOMoutside && qh QHULLfinished) {
    cpu= qh hulltime;
    cpu /= qh_SECticks;
    wval_(Wcpu)= cpu;
    fprintf (fp, "  CPU seconds to compute hull (after input): %2.4g\n", cpu);
  }
  if (qh totarea != 0.0) {
    fprintf(fp, "  %s facet area:   %2.8g\n", 
	    zzval_(Ztotmerge) ? "Approximate" : "Total", qh totarea);
    fprintf(fp, "  %s volume:       %2.8g\n", 
	    zzval_(Ztotmerge) ? "Approximate" : "Total", qh totvol);
  }
  if (!qh FORCEoutput && qh max_outside > qh DISTround) {
    dist= qh max_outside + qh DISTround;   /* agrees with qh_check_points */
    /* 1 DISTround to actual point */
    fprintf(fp, "  Maximum distance of %spoint above facet: %2.2g", 
	    (qh QHULLfinished ? "" : "merged "), dist);
    ratio= dist/(qh ONEmerge+ qh DISTround);
    if (qh MERGING && ratio > 0.05 && (qh ONEmerge > qh MINoutside))
      fprintf (fp, " (%.1fx)\n", ratio);
    else
      fprintf (fp, "\n");
  }
  if (!qh FORCEoutput && qh min_vertex < -qh DISTround) {
    dist= qh min_vertex - qh DISTround;   /* agrees with qh_check_points */
    /* 1 DISTround to actual point */
    fprintf(fp, "  Maximum distance of %svertex below facet: %2.2g", 
	    (qh QHULLfinished ? "" : "merged "), dist);
    ratio= -dist/(qh ONEmerge+qh DISTround);
    if (qh MERGING && ratio > 0.05) 
      fprintf (fp, " (%.1fx)\n", ratio);
    else
      fprintf (fp, "\n");
  }
  fprintf(fp, "\n");
} /* printsummary */


