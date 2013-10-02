/* io.c - Input/Output routines of qhull application

   see README and io.h
   
   see user.c for qh_errprint and qh_printfacetlist

   unix.c calls qh_readpoints and qh_produce_output
   
   unix.c and user.c are the only callers of io.c functions
   This allows the user to avoid loading io.o from qhull.a

   copyright (c) 1993-1995 The Geometry Center        
*/

#include "qhull_a.h"

static int qh_compare_facetarea(const void *p1, const void *p2);
static int qh_compare_facetmerge(const void *p1, const void *p2);
static int qh_compare_facetvisit(const void *p1, const void *p2);
int qh_compare_vertexpoint(const void *p1, const void *p2); /* not used */

/*========= -functions in alphabetical order after produce_output  ========= */

/*-------------------------------------------
-produce_output- prints out the result of qhull in desired format
*/
void qh_produce_output(void) {
  int i, tempsize= qh_setsize (qhmem.tempstack), d_1;

  if (qh VORONOI) {
    qh_clearcenters (qh_ASvoronoi);
    qh_vertexneighbors();
  }
  if (qh GETarea)
    qh_getarea(qh facet_list);
  qh_findgood_all (qh facet_list); 
  if (qh KEEParea || qh KEEPmerge || qh KEEPminArea < REALmax/2)
    qh_markkeep (qh facet_list);
  if (qh PRINTsummary)
    qh_printsummary(qh ferr);
  else if (qh PRINTout[0] == qh_PRINTnone)
    qh_printsummary(qh fout);
  for (i= 0; i< qh_PRINTEND; i++)
    qh_printfacets (qh fout, qh PRINTout[i], qh facet_list, NULL, !qh_ALL);
  qh_allstatistics();
  if (qh PRINTprecision && !qh MERGING)
    qh_printstats (qh ferr, qhstat precision, NULL);
  if (qh PRINTstatistics) {
    qh_collectstatistics();
    qh_printstatistics(qh ferr, "");
    qh_memstatistics (qh ferr);
    d_1= sizeof(setT) + (qh hull_dim - 1) * SETelemsize;
    fprintf(qh ferr, "\
    size in bytes: hashentry %d merge %d ridge %d vertex %d facet %d\n\
         normal %d ridge vertices %d facet vertices or neighbors %d\n",
	    sizeof(hashentryT), sizeof(mergeT), sizeof(ridgeT),
	    sizeof(vertexT), sizeof(facetT),
	    qh normal_size, d_1, d_1 + SETelemsize);
  }
  if (qh_setsize (qhmem.tempstack) != tempsize) {
    fprintf (qh ferr, "qhull internal error (qh_produce_output): temporary sets not empty (%d)\n",
	     qh_setsize (qhmem.tempstack));
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
} /* produce_output */


/*-------------------------------------------------
-compare_vertexpoint- used by qsort() to order vertices by point id 
 them
*/
int qh_compare_vertexpoint(const void *p1, const void *p2) {
  vertexT *a= *((vertexT **)p1), *b= *((vertexT **)p2);
 
  return ((qh_pointid(a->point) > qh_pointid(b->point)?1:-1));
} /* compare_vertexpoint */

/*-------------------------------------------------
-compare_facetarea- used by qsort() to order facets by area
*/
static int qh_compare_facetarea(const void *p1, const void *p2) {
  facetT *a= *((facetT **)p1), *b= *((facetT **)p2);

  if (!a->isarea)
    return -1;
  if (!b->isarea)
    return 1; 
  if (a->f.area > b->f.area)
    return 1;
  else if (a->f.area == b->f.area)
    return 0;
  return -1;
} /* compare_facetarea */

/*-------------------------------------------------
-compare_facetmerge- used by qsort() to order facets by number of merges
*/
static int qh_compare_facetmerge(const void *p1, const void *p2) {
  facetT *a= *((facetT **)p1), *b= *((facetT **)p2);
 
  return (a->nummerge - b->nummerge);
} /* compare_facetvisit */

/*-------------------------------------------------
-compare_facetvisit- used by qsort() to order facets by visit id or id
*/
static int qh_compare_facetvisit(const void *p1, const void *p2) {
  facetT *a= *((facetT **)p1), *b= *((facetT **)p2);
  int i,j;

  if (!(i= a->visitid))
    i= - a->id;
  if (!(j= b->visitid))
    j= - b->id;
  return (i - j);
} /* compare_facetvisit */

/*-------------------------------------------------
-countfacets- count good facets for printing and set visitid
returns:
  each facet with ->visitid indicating 1-relative position
      ->visitid==0 indicates not good
  numfacets, numsimplicial, total neighbors, numridges, coplanars
  if NEWfacets, does not count visible facets (matches qh_printafacet)
*/
void qh_countfacets (facetT *facetlist, setT *facets, boolT printall,
    int *numfacetsp, int *numsimplicialp, int *totneighborsp, int *numridgesp, int *numcoplanarsp) {
  facetT *facet, **facetp;
  int numfacets= 0, numsimplicial= 0, numridges= 0, totneighbors= 0, numcoplanars= 0;

  FORALLfacet_(facetlist) {
    if ((facet->visible && qh NEWfacets)
    || (!printall && qh_skipfacet(facet)))
      facet->visitid= 0;
    else {
      facet->visitid= ++numfacets;
      totneighbors += qh_setsize (facet->neighbors);
      if (facet->simplicial) 
        numsimplicial++;
      else
        numridges += qh_setsize (facet->ridges);
      if (facet->coplanarset)
        numcoplanars += qh_setsize (facet->coplanarset);
    }
  }
  FOREACHfacet_(facets) {
    if ((facet->visible && qh NEWfacets)
    || (!printall && qh_skipfacet(facet)))
      facet->visitid= 0;
    else {
      facet->visitid= ++numfacets;
      totneighbors += qh_setsize (facet->neighbors);
      if (facet->simplicial)
        numsimplicial++;
      else
        numridges += qh_setsize (facet->ridges);
      if (facet->coplanarset)
        numcoplanars += qh_setsize (facet->coplanarset);
    }
  }
  qh visit_id += numfacets+1;
  *numfacetsp= numfacets;
  *numsimplicialp= numsimplicial;
  *totneighborsp= totneighbors;
  *numridgesp= numridges;
  *numcoplanarsp= numcoplanars;
} /* countfacets */

/*---------------------------------------
-dfacet- print facet by id, for debugging
*/
void dfacet (int id) {
  facetT *facet;

  FORALLfacets {
    if (facet->id == id) {
      qh_printfacet (qh fout, facet);
      break;
    }
  }
} /* dfacet */


/*---------------------------------------
-dvertex- print vertex by id, for debugging
*/
void dvertex (int id) {
  vertexT *vertex;

  FORALLvertices {
    if (vertex->id == id) {
      qh_printvertex (qh fout, vertex);
      break;
    }
  }
} /* dvertex */


/*----------------------------------------
-facet2point- return two temporary projected points for a 2-d facet
  may be non-simplicial, returns mindist
*/
void qh_facet2point(facetT *facet, pointT **point0, pointT **point1, realT *mindist) {
  vertexT *vertex0, *vertex1;
  realT dist;
  
  if (facet->toporient ^ qh_ORIENTclock) {
    vertex0= SETfirst_(facet->vertices);
    vertex1= SETsecond_(facet->vertices);
  }else {
    vertex1= SETfirst_(facet->vertices);
    vertex0= SETsecond_(facet->vertices);
  }
  zadd_(Zdistio, 2);
  qh_distplane(vertex0->point, facet, &dist);
  *mindist= dist;
  *point0= qh_projectpoint(vertex0->point, facet, dist);
  qh_distplane(vertex1->point, facet, &dist);
  minimize_(*mindist, dist);		
  *point1= qh_projectpoint(vertex1->point, facet, dist);
} /* facet2point */


/*-------------------------------------------
-facetvertices- returns temporary set of vertices in a set or list of facets
  optimized for allfacets of facet_list
*/
setT *qh_facetvertices (facetT *facetlist, setT *facets, boolT allfacets) {
  setT *vertices;
  facetT *facet, **facetp;
  vertexT *vertex, **vertexp;

  if (facetlist == qh facet_list && allfacets && !facets) {
    vertices= qh_settemp (qh num_vertices);
    FORALLvertices 
      qh_setappend (&vertices, vertex);
  }else {
    vertices= qh_settemp (qh TEMPsize);
    qh vertex_visit++;
    FORALLfacet_(facetlist) {
      if (!allfacets && qh_skipfacet (facet))
        continue;
      FOREACHvertex_(facet->vertices) {
        if (vertex->visitid != qh vertex_visit) {
          vertex->visitid= qh vertex_visit;
          qh_setappend (&vertices, vertex);
        }
      }
    }
  }
  FOREACHfacet_(facets) {
    if (!allfacets && qh_skipfacet (facet))
      continue;
    FOREACHvertex_(facet->vertices) {
      if (vertex->visitid != qh vertex_visit) {
        vertex->visitid= qh vertex_visit;
        qh_setappend (&vertices, vertex);
      }
    }
  }
  return vertices;
} /* facetvertices */

/*-----------------------------------------
-markkeep- mark facets that meet qh KEEParea, qh KEEPmerge, and qh KEEPminArea
  clears ->good
  recomputes qh num_good
*/
void qh_markkeep (facetT *facetlist) {
  facetT *facet, **facetp;
  setT *facets= qh_settemp (qh num_facets);
  int size, count;

  trace2((qh ferr, "qh_markkeep: only keep %d largest and/or %d most merged facets and/or min area %.2g\n",
          qh KEEParea, qh KEEPmerge, qh KEEPminArea));
  FORALLfacet_(facetlist) {
    if (!facet->visible && facet->good)
      qh_setappend (&facets, facet);
  }
  size= qh_setsize (facets);
  if (qh KEEParea) {
    qsort (SETaddr_(facets, facetT), size,
             sizeof (facetT *), qh_compare_facetarea);
    if ((count= size - qh KEEParea) > 0) {
      FOREACHfacet_(facets) {
        facet->good= False;
        if (--count == 0)
          break;
      }
    }
  }
  if (qh KEEPmerge) {
    qsort (SETaddr_(facets, facetT), size,
             sizeof (facetT *), qh_compare_facetmerge);
    if ((count= size - qh KEEPmerge) > 0) {
      FOREACHfacet_(facets) {
        facet->good= False;
        if (--count == 0)
          break;
      }
    }
  }
  if (qh KEEPminArea < REALmax/2) {
    FOREACHfacet_(facets) {
      if (!facet->isarea || facet->f.area < qh KEEPminArea)
	facet->good= False;
    }
  }
  qh_settempfree (&facets);
  count= 0;
  FORALLfacet_(facetlist) {
    if (facet->good)
      count++;
  }
  qh num_good= count;
} /* markkeep */

/*-----------------------------------------
-order_vertexneighbors- order neighbors for a 2-d or 3-d vertex by adjacency
  does not orient the neighbors
*/
void qh_order_vertexneighbors(vertexT *vertex) {
  setT *newset;
  facetT *facet, *neighbor, **neighborp;

  trace4((qh ferr, "qh_order_vertexneighbors: order neighbors of v%d for 3-d\n", vertex->id));
  newset= qh_settemp (qh_setsize (vertex->neighbors));
  facet= (facetT*)qh_setdellast (vertex->neighbors);
  qh_setappend (&newset, facet);
  while (qh_setsize (vertex->neighbors)) {
    FOREACHneighbor_(vertex) {
      if (qh_setin (facet->neighbors, neighbor)) {
        qh_setdel(vertex->neighbors, neighbor);
        qh_setappend (&newset, neighbor);
        facet= neighbor;
        break;
      }
    }
    if (!neighbor) {
      fprintf (qh ferr, "qhull internal error (qh_order_vertexneighbors): no neighbor of v%d for f%d\n",
        vertex->id, facet->id);
      qh_errexit (qh_ERRqhull, facet, NULL);
    }
  }
  qh_setfree (&vertex->neighbors);
  qh_settemppop ();
  vertex->neighbors= newset;
} /* order_vertexneighbors */

/*-----------------------------------------
-printafacet- print facet in given output format (see qh PRINTout)
  nop if skipfacet() unless printall
    nop if visible facet and NEWfacets and format != PRINTfacets
    must match qh_countfacets
  preserves qh visit_id
  facet->normal may be null if PREmerge/MERGEexact and STOPcone before merge
*/
void qh_printafacet(FILE *fp, int format, facetT *facet, boolT printall) {
  realT color[4], mindist, offset, dist;
  boolT zerodiv;
  coordT *point, *normp, *coordp, **pointp, *feasiblep;
  int k;
  vertexT *vertex, **vertexp;
  facetT *neighbor, **neighborp;

  if (!printall && qh_skipfacet (facet))
    return;
  if (facet->visible && qh NEWfacets && format != qh_PRINTfacets)
    return;
  qh printoutnum++;
  switch (format) {
  case qh_PRINTarea:
    if (facet->isarea) {
      fprintf (fp, qh_REAL_1, facet->f.area);
      fprintf (fp, "\n");
    }
    break;
  case qh_PRINTcoplanars:
    fprintf (fp, "%d", qh_setsize (facet->coplanarset));
    FOREACHpoint_(facet->coplanarset)
      fprintf (fp, " %d", qh_pointid (point));
    fprintf (fp, "\n");
    break;
  case qh_PRINTcentrums:
    qh_printcenter (fp, format, NULL, facet);
    break;
  case qh_PRINTfacets:
    qh_printfacet (fp, facet);
    break;
  case qh_PRINTfacets_xridge:
    qh_printfacetheader (fp, facet);
    break;
  case qh_PRINTgeom:  /* either 2 , 3, or 4-d by qh_printbegin */
    if (!facet->normal)
      break;
    for (k= qh hull_dim; k--; ) {
      color[k]= (facet->normal[k]+1.0)/2.0;
      maximize_(color[k], -1.0);
      minimize_(color[k], +1.0);
    }
    qh_projectdim3 (color, color);
    if (qh PRINTdim != qh hull_dim)
      qh_normalize2 (color, 3, True, NULL, NULL);
    if (qh hull_dim <= 2)
      qh_printfacet2geom (fp, facet, color);
    else if (qh hull_dim == 3) {
      if (facet->simplicial)
        qh_printfacet3geom_simplicial (fp, facet, color);
      else
        qh_printfacet3geom_nonsimplicial (fp, facet, color);
    }else {
      if (facet->simplicial)
        qh_printfacet4geom_simplicial (fp, facet, color);
      else
        qh_printfacet4geom_nonsimplicial (fp, facet, color);
    }
    break;
  case qh_PRINTids:
    fprintf (fp, "%d\n", facet->id);
    break;
  case qh_PRINToff:
  case qh_PRINTincidences:
    if (qh hull_dim == 3) 
      qh_printfacet3vertex (fp, facet, format);
    else if (facet->simplicial || qh hull_dim == 2 || format == qh_PRINToff)
      qh_printfacetNvertex_simplicial (fp, facet, format);
    else
      qh_printfacetNvertex_nonsimplicial (fp, facet, qh printoutvar++);
    break;
  case qh_PRINTinner:
    mindist= REALmax;
    FOREACHvertex_(facet->vertices) {
      qh_distplane (vertex->point, facet, &dist);
      minimize_(mindist, dist);
    }
    offset= facet->offset - mindist - qh DISTround;
    goto LABELprintnorm;
  case qh_PRINTmerges:
    fprintf (fp, "%d\n", facet->nummerge);
    break;
  case qh_PRINTnormals:
    offset= facet->offset;
    goto LABELprintnorm;
  case qh_PRINTouter:
#if qh_MAXoutside
    offset= -facet->maxoutside;
#endif
    if (!qh MERGING || !qh_MAXoutside || qh SKIPcheckmax)
      offset= -fmax_(qh max_outside, qh DISTround);
    offset -= qh DISTround;      /* agrees with qh_check_points */
    /* 1 DISTround to actual point */
    offset += facet->offset;
  LABELprintnorm:
    if (!facet->normal) {
      fprintf (fp, "no normal for facet f%d\n", facet->id);
      break;
    }
    if (qh CDDoutput) 
      fprintf (fp, qh_REAL_1, -offset);
    for (k=0; k<qh hull_dim; k++) 
      fprintf (fp, qh_REAL_1, facet->normal[k]);
    if (!qh CDDoutput) 
      fprintf (fp, qh_REAL_1, offset);
    fprintf (fp, "\n");
    break;
  case qh_PRINTmathematica:  /* either 2 or 3-d by qh_printbegin */
    if (qh hull_dim == 2)
      qh_printfacet2math (fp, facet, qh printoutvar++);
    else 
      qh_printfacet3math (fp, facet, qh printoutvar++);
    break;
  case qh_PRINTneighbors:
    fprintf (fp, "%d", qh_setsize (facet->neighbors));
    FOREACHneighbor_(facet)
      fprintf (fp, " %d", 
	       neighbor->visitid ? neighbor->visitid - 1: - neighbor->id);
    fprintf (fp, "\n");
    break;
  case qh_PRINTpointintersect:
    if (!qh feasible_point) {
      fprintf (fp, "qhull input error (qh_printafacet): option 'Fp' needs qh feasible_point\n");
      qh_errexit( qh_ERRinput, NULL, NULL);
    }
    if (facet->offset > 0)
      goto LABELprintinfinite;
    point= coordp= (coordT*)qh_memalloc (qh normal_size);
    normp= facet->normal;
    feasiblep= qh feasible_point;
    if (facet->offset < -qh MINdenom) {
      for (k= qh hull_dim; k--; )
        *(coordp++)= (*(normp++) / - facet->offset) + *(feasiblep++);
    }else {
      for (k= qh hull_dim; k--; ) {
        *(coordp++)= qh_divzero (*(normp++), facet->offset, qh MINdenom_1,
				 &zerodiv) + *(feasiblep++);
        if (zerodiv) {
          qh_memfree (point, qh normal_size);
          goto LABELprintinfinite;
        }
      }
    }
    qh_printpoint (fp, NULL, point);
    qh_memfree (point, qh normal_size);
    break;
  LABELprintinfinite:
    for (k= qh hull_dim; k--; )
      fprintf (fp, qh_REAL_1, qh_INFINITE);
    fprintf (fp, "\n");   
    break;
  case qh_PRINTpointnearest:
    FOREACHpoint_(facet->coplanarset) {
      int id, id2;
      vertex= qh_nearvertex (facet, point, &dist);
      id= qh_pointid (vertex->point);
      id2= qh_pointid (point);
      fprintf (fp, "%d %d %d " qh_REAL_1 "\n", id, id2, facet->id, dist);
    }
    break;
  case qh_PRINTpoints:  /* VORONOI only by qh_printbegin */
    if (qh CDDoutput)
      fprintf (fp, "1 ");
    qh_printcenter (fp, format, NULL, facet);
    break;
  case qh_PRINTvertices:
    fprintf (fp, "%d", qh_setsize (facet->vertices));
    FOREACHvertex_(facet->vertices)
      fprintf (fp, " %d", qh_pointid (vertex->point));
    fprintf (fp, "\n");
    break;
  }
} /* printafacet */

/*-----------------------------------------
-printbegin- prints header for all output formats
  checks for valid format
  uses qh visit_id for 3/4off
  changes qh interior_point if printing centrums
*/
void qh_printbegin (FILE *fp, int format, facetT *facetlist, setT *facets, boolT printall) {
  int numfacets, numsimplicial, numridges, totneighbors, numcoplanars;
  int i, num;
  facetT *facet, **facetp;
  vertexT *vertex, **vertexp;
  setT *vertices;
  pointT *point, **pointp, *pointtemp;

  qh printoutnum= 0;
  qh_countfacets (facetlist, facets, printall, &numfacets, &numsimplicial, 
      &totneighbors, &numridges, &numcoplanars);
  switch (format) {
  case qh_PRINTnone:
    break;
  case qh_PRINTarea:
    fprintf (fp, "%d\n", numfacets);
    break;
  case qh_PRINTcoplanars:
    fprintf (fp, "%d\n", numfacets);
    break;
  case qh_PRINTcentrums:
    if (qh CENTERtype == qh_ASnone)
      qh_clearcenters (qh_AScentrum);
    fprintf (fp, "%d\n%d\n", qh hull_dim, numfacets);
    break;
  case qh_PRINTfacets:
  case qh_PRINTfacets_xridge:
    if (facetlist)
      qh_printvertexlist (fp, "Vertices and facets:\n", facetlist, facets, printall);
    break;
  case qh_PRINTgeom: 
    if (qh hull_dim > 4)  /* qh_initqhull_globals also checks */
      goto LABELnoformat;
    if (qh VORONOI && qh hull_dim > 3)  /* PRINTdim == DROPdim == hull_dim-1 */
      goto LABELnoformat;
    if (qh hull_dim == 2 && (qh PRINTridges || qh DOintersections))
      fprintf (qh ferr, "qhull warning: output for ridges and intersections not implemented in 2-d\n");
    if (qh hull_dim == 4 && (qh PRINTinner || qh PRINTouter ||
			     (qh PRINTdim == 4 && qh PRINTcentrums)))
      fprintf (qh ferr, "qhull warning: output for outer/inner planes and centrums not implemented in 4-d\n");
    if (qh PRINTdim == 4 && (qh PRINTspheres))
      fprintf (qh ferr, "qhull warning: output for vertices not implemented in 4-d\n");
    if (qh PRINTdim == 4 && qh DOintersections && qh PRINTnoplanes)
      fprintf (qh ferr, "qhull warning: 'Gnh' generates no output in 4-d\n");
    if (qh PRINTdim == 2) {
      fprintf(fp, "{appearance {linewidth 3} LIST # %s | %s\n",
	      qh rbox_command, qh qhull_command);
    }else if (qh PRINTdim == 3) {
      fprintf(fp, "{appearance {+edge -evert linewidth 2} LIST # %s | %s\n",
	      qh rbox_command, qh qhull_command);
    }else if (qh PRINTdim == 4) {
      qh visit_id++;
      num= 0;
      FORALLfacet_(facetlist)    /* get number of ridges to be printed */
        qh_printend4geom (NULL, facet, &num, printall);
      FOREACHfacet_(facets)
        qh_printend4geom (NULL, facet, &num, printall);
      qh ridgeoutnum= num;
      qh printoutvar= 0;  /* counts number of ridges in output */
      fprintf (fp, "LIST # %s | %s\n", qh rbox_command, qh qhull_command);
    }
    if (qh PRINTdots) {
      qh printoutnum++;
      num= qh num_points + qh_setsize (qh other_points);
      if (qh DELAUNAY && qh ATinfinity)
	num--;
      if (qh PRINTdim == 4)
        fprintf (fp, "4VECT %d %d 1\n", num, num);
      else
	fprintf (fp, "VECT %d %d 1\n", num, num);
      for (i= num; i--; ) {
        if (i % 20 == 0)
          fprintf (fp, "\n");
	fprintf (fp, "1 ");
      }
      fprintf (fp, "# 1 point per line\n1 ");
      for (i= num-1; i--; ) {
        if (i % 20 == 0)
          fprintf (fp, "\n");
	fprintf (fp, "0 ");
      }
      fprintf (fp, "# 1 color for all\n");
      FORALLpoints {
        if (!qh DELAUNAY || !qh ATinfinity || qh_pointid(point) != qh num_points-1) {
	  if (qh PRINTdim == 4)
	    qh_printpoint (fp, NULL, point);
	  else
	    qh_printpoint3 (fp, point);
	}
      }
      FOREACHpoint_(qh other_points) {
	if (qh PRINTdim == 4)
	  qh_printpoint (fp, NULL, point);
	else
	  qh_printpoint3 (fp, point);
      }
      fprintf (fp, "0 1 1 1  # color of points\n");
    }
    if (qh PRINTdim == 4  && !qh PRINTnoplanes)
      /* 4dview loads up multiple 4OFF objects slowly */
      fprintf(fp, "4OFF %d %d 1\n", 3*qh ridgeoutnum, qh ridgeoutnum);
    qh PRINTcradius= 2 * qh DISTround;  /* include test DISTround */
    if (qh PREmerge) {
      maximize_(qh PRINTcradius, qh premerge_centrum + qh DISTround);
    }else if (qh POSTmerge)
      maximize_(qh PRINTcradius, qh postmerge_centrum + qh DISTround);
    qh PRINTradius= qh PRINTcradius;
    if (qh PRINTspheres + qh PRINTcoplanar)
      maximize_(qh PRINTradius, qh maxmaxcoord * qh_MINradius);
    if (qh premerge_cos < REALmax/2) {
      maximize_(qh PRINTradius, (1- qh premerge_cos) * qh maxmaxcoord);
    }else if (!qh PREmerge && qh POSTmerge && qh postmerge_cos < REALmax/2) {
      maximize_(qh PRINTradius, (1- qh postmerge_cos) * qh maxmaxcoord);
    }
    maximize_(qh PRINTradius, qh MINvisible); 
    if (qh PRINTdim != 4 &&
	(qh PRINTcoplanar || qh PRINTspheres || qh PRINTcentrums)) {
      vertices= qh_facetvertices (facetlist, facets, printall);
      if (qh PRINTspheres && qh PRINTdim <= 3)
         qh_printspheres (fp, vertices, qh PRINTradius);
      if (qh PRINTcoplanar || qh PRINTcentrums) {
        qh firstcentrum= True;
        if (qh PRINTcoplanar&& !qh PRINTspheres) {
          FOREACHvertex_(vertices) 
            qh_printpointvect2 (fp, vertex->point, NULL,
				qh interior_point, qh PRINTradius);
	}
        FORALLfacet_(facetlist) {
	  if (!printall && qh_skipfacet(facet))
	    continue;
	  if (!facet->normal)
	    continue;
          if (qh PRINTcentrums && qh PRINTdim <= 3)
            qh_printcentrum (fp, facet, qh PRINTcradius);
          FOREACHpoint_(facet->coplanarset)
            qh_printpointvect2 (fp, point, facet->normal, NULL, qh PRINTradius);
          FOREACHpoint_(facet->outsideset)
            qh_printpointvect2 (fp, point, facet->normal, NULL, qh PRINTradius);
        }
        FOREACHfacet_(facets) {
	  if (!printall && qh_skipfacet(facet))
	    continue;
	  if (!facet->normal)
	    continue;
          if (qh PRINTcentrums && qh PRINTdim <= 3)
            qh_printcentrum (fp, facet, qh PRINTcradius);
          FOREACHpoint_(facet->coplanarset)
            qh_printpointvect2 (fp, point, facet->normal, NULL, qh PRINTradius);
          FOREACHpoint_(facet->outsideset)
            qh_printpointvect2 (fp, point, facet->normal, NULL, qh PRINTradius);
        }
      }
      qh_settempfree (&vertices);
    }
    qh visit_id++; /* for printing hyperplane intersections */
    break;
  case qh_PRINTids:
    fprintf (fp, "%d\n", numfacets);
    break;
  case qh_PRINTincidences:
    if (qh VORONOI)
      fprintf (qh ferr, "qhull warning: writing Delaunay.  Use 'p' or 'o' for Voronoi centers\n");
    qh printoutvar= qh vertex_id;  /* centrum id for non-simplicial facets */
    if (qh hull_dim <= 3)
      fprintf(fp, "%d\n", numfacets);
    else
      fprintf(fp, "%d\n", numsimplicial+numridges);
    break;
  case qh_PRINTinner:
  case qh_PRINTnormals:
  case qh_PRINTouter:
    if (qh CDDoutput)
      fprintf (fp, "%s | %s\nbegin\n    %d %d real\n", qh rbox_command, 
              qh qhull_command, numfacets, qh hull_dim+1);
    else
      fprintf (fp, "%d\n%d\n", qh hull_dim+1, numfacets);
    break;
  case qh_PRINTmathematica:  
    if (qh hull_dim > 3)  /* qh_initbuffers also checks */
      goto LABELnoformat;
    if (qh VORONOI)
      fprintf (qh ferr, "qhull warning: output is the Delaunay triangulation\n");
    fprintf(fp, "{\n");
    qh printoutvar= 0;   /* counts number of facets for notfirst */
    break;
  case qh_PRINTmerges:
    fprintf (fp, "%d\n", numfacets);
    break;
  case qh_PRINTpointintersect:
    fprintf (fp, "%d\n%d\n", qh hull_dim, numfacets);
    break;
  case qh_PRINTneighbors:
    fprintf (fp, "%d\n", numfacets);
    break;
  case qh_PRINToff:
    if (qh VORONOI)
      goto LABELnoformat;
    fprintf (fp, "%d\n%d %d %d\n", qh hull_dim,
      qh num_points+qh_setsize (qh other_points), numfacets, totneighbors/2);
    FORALLpoints
      qh_printpoint (qh fout, NULL, point);
    FOREACHpoint_(qh other_points)
      qh_printpoint (qh fout, NULL, point);
    break;
  case qh_PRINTpointnearest:
    fprintf (fp, "%d\n", numcoplanars);
    break;
  case qh_PRINTpoints:
    if (!qh VORONOI)
      goto LABELnoformat;
    if (qh CDDoutput)
      fprintf (fp, "%s | %s\nbegin\n%d %d real\n", qh rbox_command,
             qh qhull_command, numfacets, qh hull_dim);
    else
      fprintf (fp, "%d\n%d\n", qh hull_dim-1, numfacets);
    break;
  case qh_PRINTvertices:
    fprintf (fp, "%d\n", numfacets);
    break;
  case qh_PRINTsummary:
  default:
  LABELnoformat:
    fprintf (qh ferr, "qhull internal error (qh_printbegin): can not use this format for dimension %d\n",
         qh hull_dim);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
} /* printbegin */

/*-----------------------------------------
-printcenter- print facet center as either centrum or Voronoi center
  nop if qh CENTERtype neither CENTERvoronoi nor CENTERcentrum
  if upper envelope of Delaunay triangulation and point at-infinity
    prints qh_INFINITE instead;
  defines facet->center if needed
  string may be NULL or include %d for id.  Don't include other '%' codes.
  if format=PRINTgeom, adds a 0 if otherwise 2-d
*/
void qh_printcenter (FILE *fp, int format, char *string, facetT *facet) {
  int k, num;
  
  if (qh CENTERtype != qh_ASvoronoi && qh CENTERtype != qh_AScentrum)
    return;
  if (string)
    fprintf (fp, string, facet->id);
  if (qh CENTERtype == qh_ASvoronoi) {
    num= qh hull_dim-1;
    if (!facet->normal || !facet->upperdelaunay || !qh ATinfinity) {
      if (!facet->center)
        facet->center= qh_facetcenter (facet->vertices);
      for (k=0; k<num; k++)
        fprintf (fp, qh_REAL_1, facet->center[k]);
    }else {
      for (k=0; k<num; k++)
        fprintf (fp, qh_REAL_1, qh_INFINITE);
    }
  }else /* qh CENTERtype == qh_AScentrum */ {
    num= qh hull_dim;
    if (!facet->center)
      facet->center= qh_getcentrum (facet);
    for (k=0; k<num; k++)
      fprintf (fp, qh_REAL_1, facet->center[k]);
  }
  if (format == qh_PRINTgeom && num == 2)
    fprintf (fp, " 0\n");
  else
    fprintf (fp, "\n");
} /* printcenter */

/*-----------------------------------------
-printcentrum- print centrum for a facet in OOGL format 2-d or 3-d only
  defines ->center if needed
*/
void qh_printcentrum (FILE *fp, facetT *facet, realT radius) {
  pointT *centrum, *projpt;
  boolT tempcentrum= False;
  realT xaxis[4], yaxis[4], normal[4], dist;
  realT green[3]={0, 1, 0};
  vertexT *apex;
  int k;
  
  if (qh CENTERtype == qh_AScentrum) {
    if (!facet->center)
      facet->center= qh_getcentrum (facet);
    centrum= facet->center;
  }else {
    centrum= qh_getcentrum (facet);
    tempcentrum= True;
  }
  fprintf (fp, "{appearance {-normal -edge normscale 0} ");
  if (qh firstcentrum) {
    qh firstcentrum= False;
    fprintf (fp, "{INST geom { define centrum CQUAD  # f%d\n\
-0.3 -0.3 0.0001     0 0 1 1\n\
 0.3 -0.3 0.0001     0 0 1 1\n\
 0.3  0.3 0.0001     0 0 1 1\n\
-0.3  0.3 0.0001     0 0 1 1 } transform { \n", facet->id);
  }else
    fprintf (fp, "{INST geom { : centrum } transform { # f%d\n", facet->id);
  apex= SETfirst_(facet->vertices);
  qh_distplane(apex->point, facet, &dist);
  projpt= qh_projectpoint(apex->point, facet, dist);
  for (k= qh hull_dim; k--; ) {
    xaxis[k]= projpt[k] - centrum[k];
    normal[k]= facet->normal[k];
  }
  if (qh hull_dim == 2) {
    xaxis[2]= 0;
    normal[2]= 0;
  }else if (qh hull_dim == 4) {
    qh_projectdim3 (xaxis, xaxis);
    qh_projectdim3 (normal, normal);
    qh_normalize2 (normal, qh PRINTdim, True, NULL, NULL);
  }
  qh_crossproduct (3, xaxis, normal, yaxis);
  fprintf (fp, "%8.4g %8.4g %8.4g 0\n", xaxis[0], xaxis[1], xaxis[2]);
  fprintf (fp, "%8.4g %8.4g %8.4g 0\n", yaxis[0], yaxis[1], yaxis[2]);
  fprintf (fp, "%8.4g %8.4g %8.4g 0\n", normal[0], normal[1], normal[2]);
  qh_printpoint3 (fp, centrum);
  fprintf (fp, "1 }}}\n"); 
  qh_memfree (projpt, qh normal_size);
  qh_printpointvect (fp, centrum, facet->normal, NULL, radius, green);
  if (tempcentrum)
    qh_memfree (centrum, qh normal_size);
} /* printcentrum */
  
/*-----------------------------------------
-printend- prints trailer for all output formats
*/
void qh_printend (FILE *fp, int format, facetT *facetlist, setT *facets, boolT printall) {
  int num;
  facetT *facet, **facetp;

  if (!qh printoutnum)
    fprintf (qh ferr, "qhull warning: no facets printed\n");
  switch (format) {
  case qh_PRINTgeom:
    if (qh hull_dim == 4 && qh DROPdim < 0  && !qh PRINTnoplanes) {
      qh visit_id++;
      num= 0;
      FORALLfacet_(facetlist)
        qh_printend4geom (fp, facet,&num, printall);
      FOREACHfacet_(facets) 
        qh_printend4geom (fp, facet, &num, printall);
      if (num != qh ridgeoutnum || qh printoutvar != qh ridgeoutnum) {
	fprintf (qh ferr, "qhull internal error (qh_printend): number of ridges %d != number printed %d and at end %d\n", qh ridgeoutnum, qh printoutvar, num);
	qh_errexit (qh_ERRqhull, NULL, NULL);
      }
    }else
      fprintf(fp, "}\n");
    break;
  case qh_PRINTinner:
  case qh_PRINTnormals:
  case qh_PRINTouter:
    if (qh CDDoutput) 
      fprintf (fp, "end\n");
    break;
  case qh_PRINTmathematica:
    fprintf(fp, "}\n");
    break;
  case qh_PRINTpoints:
    if (qh CDDoutput)
      fprintf (fp, "end\n");
    break;
  }
} /* printend */

/*------------------------------------------
-printend4geom- helper function for printbegin/printend
  just counts printed ridges if fp=NULL
  uses visitid
  must agree with printfacet4geom...
*/
void qh_printend4geom (FILE *fp, facetT *facet, int *nump, boolT printall) {
  realT color[3];
  int i, num= *nump;
  facetT *neighbor, **neighborp;
  ridgeT *ridge, **ridgep;
  
  if (!printall && qh_skipfacet(facet))
    return;
  if (qh PRINTnoplanes || (facet->visible && qh NEWfacets))
    return;
  if (!facet->normal)
    return;
  if (fp) {
    for (i=0; i<3; i++) {
      color[i]= (facet->normal[i]+1.0)/2.0;
      maximize_(color[i], -1.0);
      minimize_(color[i], +1.0);
    }
  }
  facet->visitid= qh visit_id;
  if (facet->simplicial) {
    FOREACHneighbor_(facet) {
      if (neighbor->visitid != qh visit_id) {
	if (fp)
          fprintf (fp, "3 %d %d %d %8.4g %8.4g %8.4g 1 # f%d f%d\n",
		 3*num, 3*num+1, 3*num+2, color[0], color[1], color[2],
		 facet->id, neighbor->id);
	num++;
      }
    }
  }else {
    FOREACHridge_(facet->ridges) {
      neighbor= otherfacet_(ridge, facet);
      if (neighbor->visitid != qh visit_id) {
	if (fp)
          fprintf (fp, "3 %d %d %d %8.4g %8.4g %8.4g 1 #r%d f%d f%d\n",
		 3*num, 3*num+1, 3*num+2, color[0], color[1], color[2],
		 ridge->id, facet->id, neighbor->id);
	num++;
      }
    }
  }
  *nump= num;
} /* printend4geom */

/*-----------------------------------------
-printfacet- prints all fields of a facet to fp
  ridges printed in neighbor order
*/
void qh_printfacet(FILE *fp, facetT *facet) {

  qh_printfacetheader (fp, facet);
  if (facet->ridges)
    qh_printfacetridges (fp, facet);
} /* printfacet */


/*----------------------------------------
-printfacet2geom- print facet as part of a 2-d VECT for Geomview
*/
void qh_printfacet2geom(FILE *fp, facetT *facet, realT color[3]) {
  pointT *point0, *point1;
  realT mindist, maxdist;
  int k;

  qh_facet2point (facet, &point0, &point1, &mindist);
  /*
    assume precise calculations in io.c with roundoff covered by qh_GEOMepsilon
    mindist is calculated within io.c.  maxoutside is calculated elsewhere
    so a DISTround error may have occured.
  */
  if (qh MERGING) {
    mindist -= qh PRINTradius;
#if qh_MAXoutside
    maxdist= facet->maxoutside+qh PRINTradius+ qh DISTround;
#else
    maxdist= qh max_outside+qh PRINTradius+ qh DISTround;
#endif
    if (qh PRINTcoplanar || qh PRINTspheres)
      maxdist += qh maxmaxcoord * qh_GEOMepsilon;
  }else 
    mindist= maxdist= 0;
  if (qh PRINTouter || (!qh PRINTnoplanes && !qh PRINTinner))
    qh_printfacet2geom_points(fp, point0, point1, facet, maxdist, color);
  if (qh PRINTinner || (!qh PRINTnoplanes && !qh PRINTouter &&
                maxdist - mindist > 2 * qh maxmaxcoord * qh_GEOMepsilon)) {
    for(k= 3; k--; )
      color[k]= 1.0 - color[k];
    qh_printfacet2geom_points(fp, point0, point1, facet, mindist, color);
  }
  qh_memfree (point1, qh normal_size);
  qh_memfree (point0, qh normal_size); 
} /* printfacet2geom */

/*-----------------------------------------
-printfacet2geom_points- prints a 2-d facet as a VECT with 2 points at
   some offset.   The points are on the facet's plane.
*/
void qh_printfacet2geom_points(FILE *fp, pointT *point1, pointT *point2,
			       facetT *facet, realT offset, realT color[3]) {
  pointT *p1= point1, *p2= point2;

  fprintf(fp, "VECT 1 2 1 2 1 # f%d\n", facet->id);
  if (offset != 0.0) {
    p1= qh_projectpoint (p1, facet, -offset);
    p2= qh_projectpoint (p2, facet, -offset);
  }
  fprintf(fp, "%8.4g %8.4g %8.4g\n%8.4g %8.4g %8.4g\n",
           p1[0], p1[1], 0.0, p2[0], p2[1], 0.0);
  if (offset != 0.0) {
    qh_memfree (p1, qh normal_size);
    qh_memfree (p2, qh normal_size);
  }
  fprintf(fp, "%8.4g %8.4g %8.4g 1.0\n", color[0], color[1], color[2]);
} /* printfacet2geom_points */



/*----------------------------------------
-printfacet2math- print 2-d Mathematica output for a facet
  may be non-simplicial
*/
void qh_printfacet2math(FILE *fp, facetT *facet, int notfirst) {
  pointT *point0, *point1;
  realT mindist;
  
  qh_facet2point (facet, &point0, &point1, &mindist);
  if (notfirst)
    fprintf(fp, ",");
  fprintf(fp, "Line[{{%10.8g, %10.8g}, {%10.8g, %10.8g}}]\n",
	  point0[0], point0[1], point1[0], point1[1]);
  qh_memfree (point1, qh normal_size);
  qh_memfree (point0, qh normal_size);
} /* printfacet2math */


/*-----------------------------------------
-printfacet3geom_nonsimplicial- print Geomview OFF for a 3-d nonsimplicial facet.
  if DOintersections, prints ridges to unvisited neighbors (qh visit_id) 
  uses visitid for intersections and ridges
*/
void qh_printfacet3geom_nonsimplicial(FILE *fp, facetT *facet, realT color[3]) {
  ridgeT *ridge, **ridgep;
  setT *projectedpoints, *vertices;
  vertexT *vertex, **vertexp, *vertexA, *vertexB;
  pointT *projpt, *point, **pointp;
  facetT *neighbor;
  realT dist, mindist= REALmax, maxdist;
  int cntvertices, k;
  realT black[3]={0, 0, 0}, green[3]={0, 1, 0};
 
  vertices= qh_facet3vertex (facet); /* oriented */
  cntvertices= qh_setsize(vertices);
  projectedpoints= qh_settemp(cntvertices);
  FOREACHvertex_(vertices) {
    zinc_(Zdistio);
    qh_distplane(vertex->point, facet, &dist);
    projpt= qh_projectpoint(vertex->point, facet, dist);
    qh_setappend (&projectedpoints, projpt);
    minimize_(mindist, dist);
  }
  /*
    assume precise calculations in io.c with roundoff covered by qh_GEOMepsilon
    mindist is calculated within io.c.  maxoutside is calculated elsewhere
    so a DISTround error may have occured.
  */
  if (qh MERGING) {
    mindist -= qh PRINTradius;
#if qh_MAXoutside
    maxdist= facet->maxoutside+qh PRINTradius+ qh DISTround;
#else
    maxdist= qh max_outside+qh PRINTradius+ qh DISTround;
#endif
    if (qh PRINTcoplanar || qh PRINTspheres)
      maxdist += qh maxmaxcoord * qh_GEOMepsilon;
  }else
    mindist= maxdist= 0;
  if (qh PRINTouter || (!qh PRINTnoplanes && !qh PRINTinner))
    qh_printfacet3geom_points(fp, projectedpoints, facet, maxdist, color);
  if (qh PRINTinner || (!qh PRINTnoplanes && !qh PRINTouter &&
                maxdist - mindist > 2 * qh maxmaxcoord * qh_GEOMepsilon)) {
    for (k=3; k--; )
      color[k]= 1.0 - color[k];
    qh_printfacet3geom_points(fp, projectedpoints, facet, mindist, color);
  }
  FOREACHpoint_(projectedpoints)
    qh_memfree (point, qh normal_size);
  qh_settempfree(&projectedpoints);
  qh_settempfree(&vertices);
  if ((qh DOintersections || qh PRINTridges)
  && (!facet->visible || !qh NEWfacets)) {
    facet->visitid= qh visit_id;
    FOREACHridge_(facet->ridges) {
      neighbor= otherfacet_(ridge, facet);
      if (neighbor->visitid != qh visit_id) {
        if (qh DOintersections)
          qh_printhyperplaneintersection(fp, facet, neighbor, ridge->vertices, black);
        if (qh PRINTridges) {
          vertexA= SETfirst_(ridge->vertices);
          vertexB= SETsecond_(ridge->vertices);
          qh_printline3geom (fp, vertexA->point, vertexB->point, green);
        }
      }
    }
  }
} /* printfacet3geom_nonsimplicial */

/*-----------------------------------------
-printfacet3geom_points- prints a 3-d facet as OFF Geomview object. 
  Facet is determined as a list of points
*/
void qh_printfacet3geom_points(FILE *fp, setT *points, facetT *facet, realT offset, realT color[3]) {
  int k, n= qh_setsize(points), i;
  pointT *point, **pointp;
  setT *printpoints;

  fprintf(fp, "{ OFF %d 1 1 # f%d\n", n, facet->id);
  if (offset != 0.0) {
    printpoints= qh_settemp (n);
    FOREACHpoint_(points) 
      qh_setappend (&printpoints, qh_projectpoint(point, facet, -offset));
  }else
    printpoints= points;
  FOREACHpoint_(printpoints) {
    for (k=0; k<qh hull_dim; k++) {
      if (k == qh DROPdim)
        fprintf(fp, "0 ");
      else
        fprintf(fp, "%8.4g ", point[k]);
    }
    if (printpoints != points)
      qh_memfree (point, qh normal_size);
    fprintf (fp, "\n");
  }
  if (printpoints != points)
    qh_settempfree (&printpoints);
  fprintf(fp, "%d ", n);
  for(i= 0; i < n; i++)
    fprintf(fp, "%d ", i);
  fprintf(fp, "%8.4g %8.4g %8.4g 1.0 }\n", color[0], color[1], color[2]);
} /* printfacet3geom_points */


/*-----------------------------------------
-printfacet3geom_simplicial- print Geomview OFF for a 3-d simplicial facet.
  may flip color
  uses visitid for intersections and ridges
*/
void qh_printfacet3geom_simplicial(FILE *fp, facetT *facet, realT color[3]) {
  setT *points, *vertices;
  vertexT *vertex, **vertexp, *vertexA, *vertexB;
  facetT *neighbor, **neighborp;
  realT maxdist, mindist;
  realT black[3]={0, 0, 0}, green[3]={0, 1, 0};
  int k;

  vertices= qh_facet3vertex (facet);
  points= qh_settemp (qh TEMPsize);
  FOREACHvertex_(vertices)
    qh_setappend(&points, vertex->point);
  /*
    assume precise calculations in io.c with roundoff covered by qh_GEOMepsilon
    mindist may be off by qh DISTround.  Maxoutside is calculated elsewhere
    so a DISTround error may have occured.
  */
  if (qh MERGING) {
#if qh_MAXoutside
    maxdist= facet->maxoutside+qh PRINTradius+ qh DISTround;
#else
    maxdist= qh max_outside+qh PRINTradius+ qh DISTround;
#endif
    if (qh PRINTcoplanar || qh PRINTspheres)
      maxdist += qh maxmaxcoord * qh_GEOMepsilon; 
  }else 
    maxdist= 0;
  if (qh MERGING)
    mindist= -qh PRINTradius - qh DISTround;
  else
    mindist= 0;
  if (qh PRINTouter || (!qh PRINTnoplanes && !qh PRINTinner))
    qh_printfacet3geom_points(fp, points, facet, maxdist, color);
  if (qh PRINTinner || (!qh PRINTnoplanes && !qh PRINTouter &&
              maxdist - mindist > 2 * qh maxmaxcoord * qh_GEOMepsilon)) {
    for (k= 3; k--; )
      color[k]= 1.0 - color[k];
    qh_printfacet3geom_points(fp, points, facet, mindist, color);
  }
  qh_settempfree(&points);
  qh_settempfree(&vertices);
  if ((qh DOintersections || qh PRINTridges)
  && (!facet->visible || !qh NEWfacets)) {
    facet->visitid= qh visit_id;
    FOREACHneighbor_(facet) {
      if (neighbor->visitid != qh visit_id) {
	vertices= qh_setnew_delnthsorted (facet->vertices, qh hull_dim,
	                  SETindex_(facet->neighbors, neighbor), 0);
        if (qh DOintersections)
	   qh_printhyperplaneintersection(fp, facet, neighbor, vertices, black); 
        if (qh PRINTridges) {
          vertexA= SETfirst_(vertices);
          vertexB= SETsecond_(vertices);
          qh_printline3geom (fp, vertexA->point, vertexB->point, green);
        }
	qh_setfree(&vertices);
      }
    }
  }
} /* printfacet3geom_simplicial */

/*-----------------------------------------
-printfacet3vertex- print vertex->point for a 3-d facet (may be non-simplicial)
  prints number of vertices first if format == qh_PRINToff
*/
void qh_printfacet3vertex(FILE *fp, facetT *facet, int format) {
  vertexT *vertex, **vertexp;
  setT *vertices;

  vertices= qh_facet3vertex (facet);
  if (format == qh_PRINToff)
    fprintf (fp, "%d ", qh_setsize (vertices));
  FOREACHvertex_(vertices) 
    fprintf (fp, "%d ", qh_pointid(vertex->point));
  fprintf (fp, "\n");
  qh_settempfree(&vertices);
} /* printfacet3vertex */


/*----------------------------------------
-printfacet3math- print 3-d Mathematica output for a facet
  may be non-simplicial
*/
void qh_printfacet3math (FILE *fp, facetT *facet, int notfirst) {
  vertexT *vertex, **vertexp;
  setT *points, *vertices;
  pointT *point, **pointp;
  boolT firstpoint= True;
  realT dist;
  
  if (notfirst)
    fprintf(fp, ",\n");
  vertices= qh_facet3vertex (facet);
  points= qh_settemp (qh_setsize (vertices));
  FOREACHvertex_(vertices) {
    zinc_(Zdistio);
    qh_distplane(vertex->point, facet, &dist);
    point= qh_projectpoint(vertex->point, facet, dist);
    qh_setappend (&points, point);
  }
  fprintf(fp, "Polygon[{");
  FOREACHpoint_(points) {
    if (firstpoint)
      firstpoint= False;
    else
      fprintf(fp, ",\n");
    fprintf(fp, "{%10.8g, %10.8g, %10.8g}", point[0], point[1], point[2]);
  }
  FOREACHpoint_(points)
    qh_memfree (point, qh normal_size);
  qh_settempfree(&points);
  qh_settempfree(&vertices);
  fprintf(fp, "}]");
} /* printfacet3math */


/*-----------------------------------------
-printfacet4geom_nonsimplicial- print Geomview 4OFF file for a 4d nonsimplicial facet
  prints all ridges to unvisited neighbors (qh visit_id)
    must agree with printend4geom()
  prints in OFF format if DROPdim
*/
void qh_printfacet4geom_nonsimplicial(FILE *fp, facetT *facet, realT color[3]) {
  facetT *neighbor;
  ridgeT *ridge, **ridgep;
  vertexT *vertex, **vertexp;
  pointT *point;
  int k;
  realT dist;
  
  facet->visitid= qh visit_id;
  if (qh PRINTnoplanes || (facet->visible && qh NEWfacets))
    return;
  FOREACHridge_(facet->ridges) {
    neighbor= otherfacet_(ridge, facet);
    if (neighbor->visitid == qh visit_id) 
      continue;
    if (qh DOintersections)
      qh_printhyperplaneintersection(fp, facet, neighbor, ridge->vertices, color);
    else {
      if (qh DROPdim >= 0) 
	fprintf(fp, "OFF 3 1 1 # f%d\n", facet->id);
      else {
	qh printoutvar++;
	fprintf (fp, "# r%d between f%d f%d\n", ridge->id, facet->id, neighbor->id);
      }
      FOREACHvertex_(ridge->vertices) {
	zinc_(Zdistio);
	qh_distplane(vertex->point,facet, &dist);
	point=qh_projectpoint(vertex->point,facet, dist);
	for(k= 0; k < qh hull_dim; k++) {
	  if (k != qh DROPdim)
  	    fprintf(fp, "%8.4g ", point[k]);
  	}
	fprintf (fp, "\n");
	qh_memfree (point, qh normal_size);
      }
      if (qh DROPdim >= 0)
        fprintf(fp, "3 0 1 2 %8.4g %8.4g %8.4g\n", color[0], color[1], color[2]);
    }
  }
} /* printfacet4geom_nonsimplicial */


/*-----------------------------------------
-printfacet4geom_simplicial- print Geomview 4OFF file for a 4d simplicial facet
  prints triangles for unvisited neighbors (qh visit_id)
    must agree with printend4geom()
*/
void qh_printfacet4geom_simplicial(FILE *fp, facetT *facet, realT color[3]) {
  setT *vertices;
  facetT *neighbor, **neighborp;
  vertexT *vertex, **vertexp;
  int k;
  
  facet->visitid= qh visit_id;
  if (qh PRINTnoplanes || (facet->visible && qh NEWfacets))
    return;
  FOREACHneighbor_(facet) {
    if (neighbor->visitid == qh visit_id)
      continue;
    vertices= qh_setnew_delnthsorted (facet->vertices, qh hull_dim,
	                  SETindex_(facet->neighbors, neighbor), 0);
    if (qh DOintersections)
      qh_printhyperplaneintersection(fp, facet, neighbor, vertices, color);
    else {
      if (qh DROPdim >= 0) 
	fprintf(fp, "OFF 3 1 1 # ridge between f%d f%d\n",
		facet->id, neighbor->id);
      else {
	qh printoutvar++;
	fprintf (fp, "# ridge between f%d f%d\n", facet->id, neighbor->id);
      }
      FOREACHvertex_(vertices) {
	for(k= 0; k < qh hull_dim; k++) {
	  if (k != qh DROPdim)
  	    fprintf(fp, "%8.4g ", vertex->point[k]);
  	}
	fprintf (fp, "\n");
      }
      if (qh DROPdim >= 0) 
        fprintf(fp, "3 0 1 2 %8.4g %8.4g %8.4g\n", color[0], color[1], color[2]);
    }
    qh_setfree(&vertices);
  }
} /* printfacet4geom_simplicial */


/*-----------------------------------------
-printfacetNvertex_nonsimplicial- print vertices for an N-d non-simplicial facet
  triangulates each ridge to the id
*/
void qh_printfacetNvertex_nonsimplicial(FILE *fp, facetT *facet, int id) {
  vertexT *vertex, **vertexp;
  ridgeT *ridge, **ridgep;

  if (facet->visible && qh NEWfacets)
    return;
  FOREACHridge_(facet->ridges) {
    fprintf(fp, "%d ", id);
    if ((ridge->top == facet) ^ qh_ORIENTclock) {
      FOREACHvertex_(ridge->vertices)
        fprintf(fp, "%d ", qh_pointid(vertex->point));
    }else {
      FOREACHvertexreverse12_(ridge->vertices)
        fprintf(fp, "%d ", qh_pointid(vertex->point));
    }
    fprintf(fp, "\n");
  }
} /* printfacetNvertex_nonsimplicial */


/*-----------------------------------------
-printfacetNvertex_simplicial- print vertices for an N-d simplicial facet
  also prints PRINToff format for non-simplicial facets
*/
void qh_printfacetNvertex_simplicial(FILE *fp, facetT *facet, int format) {
  vertexT *vertex, **vertexp;

  if (format == qh_PRINToff)
    fprintf (fp, "%d ", qh_setsize (facet->vertices));
  if ((facet->toporient ^ qh_ORIENTclock) || !facet->simplicial) {
    FOREACHvertex_(facet->vertices)
      fprintf(fp, "%d ", qh_pointid(vertex->point));
  }else {
    FOREACHvertexreverse12_(facet->vertices)
      fprintf(fp, "%d ", qh_pointid(vertex->point));
  }
  fprintf(fp, "\n");
} /* printfacetNvertex_simplicial */


/*-----------------------------------------
-printfacetheader- prints header fields of a facet to fp
*/
void qh_printfacetheader(FILE *fp, facetT *facet) {
  pointT *point, **pointp, *furthest;
  facetT *neighbor, **neighborp;
  realT dist;

  if (facet == qh_MERGEridge) {
    fprintf (fp, " MERGEridge\n");
    return;
  }else if (facet == qh_DUPLICATEridge) {
    fprintf (fp, " DUPLICATEridge\n");
    return;
  }else if (!facet) {
    fprintf (fp, " NULLfacet\n");
    return;
  }
  fprintf(fp, "- f%d\n", facet->id);
  fprintf(fp, "    - flags:");
  if (facet->toporient) 
    fprintf(fp, " top");
  else
    fprintf(fp, " bottom");
  if (facet->simplicial)
    fprintf(fp, " simplicial");
  if (facet->upperdelaunay)
    fprintf(fp, " upperDelaunay");
  if (facet->visible)
    fprintf(fp, " visible");
  if (facet->newfacet)
    fprintf(fp, " new");
  if (facet->tested)
    fprintf(fp, " tested");
  if (!facet->good)
    fprintf(fp, " notG");
  if (facet->seen)
    fprintf(fp, " seen");
  if (facet->coplanar)
    fprintf(fp, " coplanar");
  if (facet->mergehorizon)
    fprintf(fp, " mergehorizon");
  if (facet->keepcentrum)
    fprintf(fp, " keepcentrum");
  if (facet->dupridge)
    fprintf(fp, " dupridge");
  if (facet->mergeridge && !facet->mergeridge2)
    fprintf(fp, " mergeridge1");
  if (facet->mergeridge2)
    fprintf(fp, " mergeridge2");
  if (facet->newmerge)
    fprintf(fp, " newmerge");
  if (facet->flipped) 
    fprintf(fp, " flipped");
  if (facet->degenerate)
    fprintf(fp, " degenerate");
  if (facet->redundant)
    fprintf(fp, " redundant");
  fprintf(fp, "\n");
  if (facet->isarea)
    fprintf(fp, "    - area: %2.2g\n", facet->f.area);
  else if (qh NEWfacets && facet->visible && facet->f.replace)
    fprintf(fp, "    - replacement: f%d\n", facet->f.replace->id);
  else if (facet->newfacet) {
    if (facet->f.samecycle && facet->f.samecycle != facet)
      fprintf(fp, "    - shares same visible/horizon as f%d\n", facet->f.samecycle->id);
  }else if (facet->f.newcycle)
    fprintf(fp, "    - was horizon to f%d\n", facet->f.newcycle->id);
  if (facet->nummerge)
    fprintf(fp, "    - merges: %d\n", facet->nummerge);
  qh_printpoint(fp, "    - normal: ", facet->normal);
  fprintf(fp, "    - offset: %10.7g\n", facet->offset);
  if (qh CENTERtype == qh_ASvoronoi || facet->center)
    qh_printcenter (fp, qh_PRINTfacets, "    - center: ", facet);
#if qh_MAXoutside
  if (facet->maxoutside > qh DISTround)
    fprintf(fp, "    - maxoutside: %10.7g\n", facet->maxoutside);
#endif
  if (!SETempty_(facet->outsideset)) {
    furthest= (pointT*)qh_setlast(facet->outsideset);
    if (qh_setsize (facet->outsideset) < 6) {
      fprintf(fp, "    - outside set (furthest p%d):\n", qh_pointid(furthest));
      FOREACHpoint_(facet->outsideset)
	qh_printpoint(fp, "     ", point);
    }else if (qh_setsize (facet->outsideset) < 21) {
      qh_printpoints(fp, "    - outside set:", facet->outsideset);
    }else {
      fprintf(fp, "    - outside set:  %d points.", qh_setsize(facet->outsideset));
      qh_printpoint(fp, "  Furthest", furthest);
    }
#if !qh_COMPUTEfurthest
    fprintf(fp, "    - furthest distance= %2.2g\n", facet->furthestdist);
#endif
  }
  if (!SETempty_(facet->coplanarset)) {
    furthest= (pointT*)qh_setlast(facet->coplanarset);
    if (qh_setsize (facet->coplanarset) < 6) {
      fprintf(fp, "    - coplanar set (furthest p%d):\n", qh_pointid(furthest));
      FOREACHpoint_(facet->coplanarset)
	qh_printpoint(fp, "     ", point);
    }else if (qh_setsize (facet->coplanarset) < 21) {
      qh_printpoints(fp, "    - coplanar set:", facet->coplanarset);
    }else {
      fprintf(fp, "    - coplanar set:  %d points.", qh_setsize(facet->coplanarset));
      qh_printpoint(fp, "  Furthest", furthest);
    }
    zinc_(Zdistio);
    qh_distplane (furthest, facet, &dist);
    fprintf(fp, "      furthest distance= %2.2g\n", dist);
  }
  qh_printvertices (fp, "    - vertices:", facet->vertices);
  fprintf(fp, "    - neighboring facets: ");
  FOREACHneighbor_(facet) {
    if (neighbor == qh_MERGEridge)
      fprintf(fp, " MERGE");
    else if (neighbor == qh_DUPLICATEridge)
      fprintf(fp, " DUP");
    else
      fprintf(fp, " f%d", neighbor->id);
  }
  fprintf(fp, "\n");
} /* printfacetheader */


/*-----------------------------------------
-printfacetridges- prints ridges of a facet to fp
  ridges printed in neighbor order
  assumes the ridges exist
*/
void qh_printfacetridges(FILE *fp, facetT *facet) {
  facetT *neighbor, **neighborp;
  ridgeT *ridge, **ridgep, *firstridge;
  int numridges= 0;


  if (facet->visible && qh NEWfacets) {
    fprintf(fp, "    - ridges (ids may be garbage):");
    FOREACHridge_(facet->ridges)
      fprintf(fp, " r%d", ridge->id);
    fprintf(fp, "\n");
  }else {
    fprintf(fp, "    - ridges:\n");
    FOREACHridge_(facet->ridges)
      ridge->seen= False;
    if (qh hull_dim == 3) {
      ridge= firstridge= SETfirst_(facet->ridges);
      while (ridge && !ridge->seen) {
	ridge->seen= True;
	qh_printridge(fp, ridge);
	numridges++;
	ridge= qh_nextridge3d (ridge, facet, NULL);
	}
    }else {
      FOREACHneighbor_(facet) {
	FOREACHridge_(facet->ridges) {
	  if (otherfacet_(ridge,facet) == neighbor) {
	    ridge->seen= True;
	    qh_printridge(fp, ridge);
	    numridges++;
	  }
	}
      }
    }
    if (numridges != qh_setsize (facet->ridges)) {
      fprintf (fp, "     - all ridges:");
      FOREACHridge_(facet->ridges) 
	fprintf (fp, " r%d", ridge->id);
        fprintf (fp, "\n");
    }
    FOREACHridge_(facet->ridges) {
      if (!ridge->seen) 
	qh_printridge(fp, ridge);
    }
  }
} /* printfacetridges */

/*-----------------------------------------
-printfacets- prints facetlist and/or facet set in output format
*/
void qh_printfacets(FILE *fp, int format, facetT *facetlist, setT *facets, boolT printall) {
  int numfacets, numsimplicial, numridges, totneighbors, numcoplanars;
  facetT *facet, **facetp;
  setT *vertices;
  coordT *center;

  qh old_randomdist= qh RANDOMdist;
  qh RANDOMdist= False;
  if (format == qh_PRINTnone)
    ; /* print nothing */
  else if (format == qh_PRINTaverage) {
    vertices= qh_facetvertices (facetlist, facets, printall);
    center= qh_getcenter (vertices);
    fprintf (fp, "%d 1\n", qh hull_dim);
    qh_printpoint (fp, NULL, center);
    qh_memfree (center, qh normal_size);
    qh_settempfree (&vertices);
  }else if (format == qh_PRINToptions)
    fprintf(fp, "Options selected for qhull %s:\n%s\n", qh_version, qh qhull_options);
  else if (format == qh_PRINTpoints && !qh VORONOI)
    qh_printpoints_out (fp, facetlist, facets, printall);
  else if (format == qh_PRINTqhull)
    fprintf (fp, "%s | %s\n", qh rbox_command, qh qhull_command);
  else if (format == qh_PRINTsize) {
    fprintf (fp, "0\n2 ");
    fprintf (fp, qh_REAL_1, qh totarea);
    fprintf (fp, qh_REAL_1, qh totvol);
    fprintf (fp, "\n");
  }else if (format == qh_PRINTsummary) {
    qh_countfacets (facetlist, facets, printall, &numfacets, &numsimplicial, 
      &totneighbors, &numridges, &numcoplanars);
    vertices= qh_facetvertices (facetlist, facets, printall); 
    fprintf (fp, "7 %d %d %d %d %d %d %d\n2 ", qh hull_dim, 
                qh num_points + qh_setsize (qh other_points),
                qh num_vertices, qh num_facets - qh num_visible,
                qh_setsize (vertices), numfacets, numcoplanars);
    fprintf (fp, qh_REAL_2n, qh max_outside + qh DISTround, 
                  qh min_vertex - qh DISTround);  /* agrees with qh_check_points */
    qh_settempfree (&vertices);
  }else if (format == qh_PRINTvneighbors)
    qh_printvneighbors (fp, facetlist, facets, printall);
  else if (qh VORONOI && format == qh_PRINToff)
    qh_printvoronoi (fp, format, facetlist, facets, printall);
  else if (qh VORONOI && format == qh_PRINTgeom) {
    qh_printbegin (fp, format, facetlist, facets, printall);
    qh_printvoronoi (fp, format, facetlist, facets, printall);
    qh_printend (fp, format, facetlist, facets, printall);
  }else {
    qh_printbegin (fp, format, facetlist, facets, printall);
    FORALLfacet_(facetlist)
      qh_printafacet (fp, format, facet, printall);
    FOREACHfacet_(facets) 
      qh_printafacet (fp, format, facet, printall);
    qh_printend (fp, format, facetlist, facets, printall);
  }
  qh RANDOMdist= qh old_randomdist;
} /* printfacets */


/*-----------------------------------------
-printhelp_degenerate- prints descriptive message
  no message if qh_QUICKhelp
*/
void qh_printhelp_degenerate(FILE *fp) {
  
  if (qh MERGEexact || qh PREmerge) 
    fprintf(fp, "\n\
A Qhull error has occurred.  It should have corrected the above\n\
precision error.  Please report the error to qhull_bug@geom.umn.edu\n");
  else if (!qh_QUICKhelp) {
    fprintf(fp, "\n\
Precision problems were detected during construction of the convex hull.\n\
This occurs because convex hull algorithms assume that calculations are\n\
exact, but floating-point arithmetic has round-off errors.\n\
\n\
To correct for precision problems, use option 'C-0' (2-d to 4-d) or \n\
'Qx' (5-d and higher).  This will merge non-convex facets.  See \"Imprecision\n\
in Qhull\" (qh-impre.html).\n\
\n\
If you do not use options 'C-n', 'A-n' or 'Qx', the output may include\n\
coplanar ridges, concave ridges, and flipped facets.  In 4-d and higher,\n\
Qhull may produce a ridge with four neighbors or two facets with the same \n\
vertices.  Qhull reports these events when they occur.  It stops when a\n\
concave ridge, flipped facet, or duplicate facet occurs.\n");
#if REALfloat
    fprintf (fp, "\
\n\
Qhull is currently using single precision arithmetic.  The following\n\
will probably remove the precision problems:\n\
  - recompile qhull for double precision (#define REALfloat 0 in user.h).\n");
#endif
    fprintf(fp, "\
\n\
If you can not merge facets because you want simplicial facets,\n\
try one or more of the following options.  They can not guarantee an output.\n\
  - use 'Po' to produce output and prevent partitioning for flipped facets\n\
  - use 'V0' to set min. distance to visible facet as 0 instead of roundoff\n\
  - use 'En' to specify a maximum roundoff error less than %2.2g.\n\
  - options 'Qf', 'QR0', and 'QbB' may also help\n",
               qh DISTround);
    if (!qh ATinfinity && qh DELAUNAY)
      fprintf(fp, "\
  - remove 'Qu' to add a point \"at-infinity\"\n");
    fprintf(fp, "\
\n\
To guarantee simplicial output:\n\
  - use exact arithmetic (see \"Imprecision in Qhull\", qh-impre.html)\n\
\n\
To merge non-convex facets\n");
    if (qh hull_dim >= 5) 
      fprintf (fp, "\
  - use option 'Qx'\n");
    else
      fprintf (fp, "\
  - use option 'C-0'\n");
  }
} /* printhelp_degenerate */


/*-----------------------------------------
-printhelp_singular- prints descriptive message
  qh_QUICKhelp just prints the numbers.
*/
void qh_printhelp_singular(FILE *fp) {
  facetT *facet;
  vertexT *vertex, **vertexp;
  realT min, max, *coord, dist;
  int i,k;
  
  fprintf(fp, "\n\
The input to qhull appears to be less than %d dimensional, or a\n\
computation has overflowed.\n\n\
Qhull could not construct a clearly convex simplex from points:\n",
           qh hull_dim);
  qh_printvertexlist (fp, "", qh facet_list, NULL, qh_ALL);
  if (!qh_QUICKhelp)
    fprintf(fp, "\n\
The center point is coplanar with a facet, or a vertex is coplanar\n\
with a neighboring facet.  The maximum round off error for\n\
computing distances is %2.2g.  The center point, facets and distances\n\
to the center point are as follows:\n\n", qh DISTround);
  qh_printpoint (fp, "center point", qh interior_point);
  fprintf (fp, "\n");
  FORALLfacets {
    fprintf (fp, "facet");
    FOREACHvertex_(facet->vertices)
      fprintf (fp, " p%d", qh_pointid(vertex->point));
    zinc_(Zdistio);
    qh_distplane(qh interior_point, facet, &dist);
    fprintf (fp, " distance= %4.2g\n", dist);
  }
  if (!qh_QUICKhelp) {
    if (qh HALFspace) 
      fprintf (fp, "\n\
These points are the dual of the given halfspaces.  They indicate that\n\
the intersection is degenerate.\n");
    fprintf (fp,"\n\
These points either have a maximum or minimum x-coordinate, or\n\
they maximize the determinant for k coordinates.  Trial points\n\
are first selected from points that maximize a coordinate.\n");
    if (qh hull_dim >= qh_INITIALmax)
      fprintf (fp, "\n\
Because of the high dimension, the min x-coordinate and max-coordinate\n\
points are used if the determinant is non-zero.  The option 'Qs' will\n\
do a better, though much slower, job.  Instead of 'Qs', you can change\n\
the points by randomly rotating the input with 'QR0'.\n");
  }
  fprintf (fp, "\nThe min and max coordinates for each dimension are:\n");
  for (k=0; k<qh hull_dim; k++) {
    min= REALmax;
    max= -REALmin;
    for (i=qh num_points, coord= qh first_point+k; i--; coord += qh hull_dim) {
      maximize_(max, *coord);
      minimize_(min, *coord);
    }
    fprintf (fp, "  %d:  %8.4g  %8.4g  difference= %4.4g\n", k, min, max, max-min);
  }
  if (!qh_QUICKhelp) {
    fprintf (fp, "\n\
If the input should be full dimensional, you have several options that\n\
may determine an initial simplex:\n\
  - use 'QbB' to scale the points to the unit cube\n\
  - use 'QR0' to randomly rotate the input for different maximum points\n\
  - use the 'Qs' flag to search all points for the initial simplex\n\
  - use 'En' to specify a maximum roundoff error less than %2.2g.\n\
  - trace execution with 'T3' to see the determinant for each point.\n",
                     qh DISTround);
#if REALfloat
    fprintf (fp, "\
  - recompile qhull for double precision (#define REALfloat 0 in qhull.h).\n");
#endif
    fprintf (fp, "\n\
If the input is lower dimensional:\n\
   - use 'Qbk:0Bk:0' to delete coordinate k from the input.  You should\n\
     pick the coordinate with the least range.  The hull will have the\n\
     correct topology.\n\
   - determine the flat containing the points, rotate the points\n\
     into a coordinate plane, and delete the other coordinates.\n\
   - add one or more points to make the input full dimensional.\n");
  }
} /* printhelp_singular */

/*-----------------------------------------
-printhyperplaneintersection- print Geomview OFF or 4OFF for
   the intersection of two hyperplanes in 3-d or 4-d
*/
void qh_printhyperplaneintersection(FILE *fp, facetT *facet1, facetT *facet2,
		   setT *vertices, realT color[3]) {
  realT d1, d2, costheta, denominator, dist1, dist2, s, t, mindenom, p[4];
  vertexT *vertex, **vertexp;
  int i, k;
  boolT nearzero1, nearzero2;
  
  costheta= qh_getangle(facet1->normal, facet2->normal);
  denominator= 1 - costheta * costheta;
  d1= -facet1->offset;
  d2= -facet2->offset;
  i= qh_setsize(vertices);
  if (qh hull_dim == 3)
    fprintf(fp, "VECT 1 %d 1 %d 1 ", i, i);
  else if (qh hull_dim == 4 && qh DROPdim >= 0)
    fprintf(fp, "OFF 3 1 1 ");
  else
    qh printoutvar++;
  fprintf (fp, "# intersect f%d f%d\n", facet1->id, facet2->id);
  mindenom= 1 / (10.0 * qh maxmaxcoord);
  FOREACHvertex_(vertices) {
    zadd_(Zdistio, 2);
    qh_distplane(vertex->point, facet1, &dist1);
    qh_distplane(vertex->point, facet2, &dist2);
    s= qh_divzero (-dist1 + costheta * dist2, denominator,mindenom,&nearzero1);
    t= qh_divzero (-dist2 + costheta * dist1, denominator,mindenom,&nearzero2);
    if (nearzero1 || nearzero2)
      s= t= 0.0;
    for(k= qh hull_dim; k--; )
      p[k]= vertex->point[k] + facet1->normal[k] * s + facet2->normal[k] * t;
    if (qh PRINTdim <= 3) {
      qh_projectdim3 (p, p);
      fprintf(fp, "%8.4g %8.4g %8.4g # ", p[0], p[1], p[2]);
    }else 
      fprintf(fp, "%8.4g %8.4g %8.4g %8.4g # ", p[0], p[1], p[2], p[3]);
    if (nearzero1+nearzero2)
      fprintf (fp, "p%d (coplanar facets)\n", qh_pointid (vertex->point));
    else
      fprintf (fp, "projected p%d\n", qh_pointid (vertex->point));
  }
  if (qh hull_dim == 3)
    fprintf(fp, "%8.4g %8.4g %8.4g 1.0\n", color[0], color[1], color[2]); 
  else if (qh hull_dim == 4 && qh DROPdim >= 0)  
    fprintf(fp, "3 0 1 2 %8.4g %8.4g %8.4g 1.0\n", color[0], color[1], color[2]);
} /* printhyperplaneintersection */

/*----------------------------------------
-printline3geom- prints a line as a VECT
  0's for DROPdim
  if pointA == pointB, it's a 1 point VECT
*/
void qh_printline3geom (FILE *fp, pointT *pointA, pointT *pointB, realT color[3]) {
  int k;
  realT pA[4], pB[4];

  qh_projectdim3(pointA, pA);
  qh_projectdim3(pointB, pB);
  if ((fabs(pA[0] - pB[0]) > 1e-3) || 
      (fabs(pA[1] - pB[1]) > 1e-3) || 
      (fabs(pA[2] - pB[2]) > 1e-3)) {
    fprintf (fp, "VECT 1 2 1 2 1\n");
    for (k= 0; k < 3; k++)
       fprintf (fp, "%8.4g ", pB[k]);
    fprintf (fp, " # p%d\n", qh_pointid (pointB));
  }else
    fprintf (fp, "VECT 1 1 1 1 1\n");
  for (k=0; k< 3; k++)
    fprintf (fp, "%8.4g ", pA[k]);
  fprintf (fp, " # p%d\n", qh_pointid (pointA));
  fprintf (fp, "%8.4g %8.4g %8.4g 1\n", color[0], color[1], color[2]);
}

/*----------------------------------------
-printneighborhood- print neighborhood of one or two facets
  calls findgood_all if active and !printall
  bumps qh visit_id
*/
void qh_printneighborhood (FILE *fp, int format, facetT *facetA, facetT *facetB, boolT printall) {
  facetT *neighbor, **neighborp, *facet;
  setT *facets;

  if (format == qh_PRINTnone)
    return;
  qh_findgood_all (qh facet_list);
  if (facetA == facetB)
    facetB= NULL;
  facets= qh_settemp (2*(qh_setsize (facetA->neighbors)+1));
  qh visit_id++;
  for (facet= facetA; facet; facet= ((facet == facetA) ? facetB : NULL)) {
    if (facet->visitid != qh visit_id) {
      facet->visitid= qh visit_id;
      qh_setappend (&facets, facet);
    }
    FOREACHneighbor_(facet) {
      if (neighbor->visitid == qh visit_id)
        continue;
      neighbor->visitid= qh visit_id;
      if (printall || !qh_skipfacet (neighbor))
        qh_setappend (&facets, neighbor);
    }
  }
  qh_printfacets (fp, format, NULL, facets, printall);
  qh_settempfree (&facets);
} /* printneighborhood */

/*----------------------------------------
-printpoint- prints the coordinates of a point
  prints point id if string and valid point
  nop if point is NULL
*/
void qh_printpoint(FILE *fp, char *string, pointT *point) {
  int k,id;
  realT r; /*bug fix*/
  
  if (!point)
    return;
  if (string) {
    fputs (string, fp);
    if ((id= qh_pointid(point)) != -1)
      fprintf(fp, " p%d: ", id);
  }
  for(k= qh hull_dim; k--; ) {
    if (string)
      fprintf(fp, " %8.4g", r=*point++);
    else
      fprintf(fp, qh_REAL_1, r=*point++);
  }
  fprintf(fp, "\n");
} /* printpoint */

/*------------------------------------
-printpoint3- prints 2-d , 3-d, or 4-d point as Geomview 3-d coordinates
*/
void qh_printpoint3 (FILE *fp, pointT *point) {
  int k;
  realT p[4];
  
  qh_projectdim3 (point, p);
  for (k=0; k< 3; k++)
    fprintf (fp, "%8.4g ", p[k]);
  fprintf (fp, " # p%d\n", qh_pointid (point));
} /* printpoint3 */

/*----------------------------------------
-printpoints- print pointids for a set of points starting at index 
   see geom.c
*/

/*----------------------------------------
-printpoints_out- prints vertices for facets by their point coordinates
  same format as qhull input
  allows CDDoutput
*/
void qh_printpoints_out (FILE *fp, facetT *facetlist, setT *facets, int printall) {
  int allpoints= qh num_points + qh_setsize (qh other_points);
  int numpoints=0, point_i, point_n;
  setT *vertices, *points;
  facetT *facet, **facetp;
  pointT *point, **pointp;
  vertexT *vertex, **vertexp;
  int id;

  points= qh_settemp (allpoints);
  qh_setzero (points, 0, allpoints);
  vertices= qh_facetvertices (facetlist, facets, printall);
  FOREACHvertex_(vertices) {
    id= qh_pointid (vertex->point);
    if (id >= 0)
      SETelem_(points, id)= vertex->point;
  }
  if (qh KEEPinside || qh KEEPcoplanar || qh KEEPnearinside) {
    FORALLfacet_(facetlist) {
      FOREACHpoint_(facet->coplanarset) {
        id= qh_pointid (point);
        if (id >= 0)
          SETelem_(points, id)= point;
      }
    }
    FOREACHfacet_(facets) {
      FOREACHpoint_(facet->coplanarset) {
        id= qh_pointid (point);
        if (id >= 0)
          SETelem_(points, id)= point;
      }
    }
  }
  qh_settempfree (&vertices);
  FOREACHpoint_i_(points) {
    if (point)
      numpoints++;
  }
  if (qh CDDoutput)
    fprintf (fp, "%s | %s\nbegin\n%d %d real\n", qh rbox_command,
             qh qhull_command, numpoints, qh hull_dim + 1);
  else
    fprintf (fp, "%d\n%d\n", qh hull_dim, numpoints);
  FOREACHpoint_i_(points) {
    if (point) {
      if (qh CDDoutput)
	fprintf (fp, "1 ");
      qh_printpoint (fp, NULL, point);
    }
  }
  if (qh CDDoutput)
    fprintf (fp, "end\n");
  qh_settempfree (&points);
} /* printpoints_out */
  

/*----------------------------------------
-printpointvect2- prints a 2-d, 3-d, or 4-d point as 2 3-d VECT's\
  for an imprecise point, 
*/
void qh_printpointvect2 (FILE *fp, pointT *point, coordT *normal, pointT *center, realT radius) {
  realT red[3]={1, 0, 0}, yellow[3]={1, 1, 0};

  qh_printpointvect (fp, point, normal, center, radius, red);
  qh_printpointvect (fp, point, normal, center, -radius, yellow);
} /* printpointvect2 */

/*----------------------------------------
-printpointvect- prints a 2-d, 3-d, or 4-d point as 3-d VECT's
  relative to normal or to center point
*/
void qh_printpointvect (FILE *fp, pointT *point, coordT *normal, pointT *center, realT radius, realT color[3]) {
  realT diff[4], pointA[4];
  int k;
  
  for (k= qh hull_dim; k--; ) {
    if (center)
      diff[k]= point[k]-center[k];
    else if (normal) 
      diff[k]= normal[k];
    else
      diff[k]= 0;
  }
  if (center)
    qh_normalize2 (diff, qh hull_dim, True, NULL, NULL);
  for (k= qh hull_dim; k--; ) 
    pointA[k]= point[k]+diff[k] * radius;
  qh_printline3geom (fp, point, pointA, color);
} /* printpointvect */  

/*----------------------------------------
-printridge- prints the information in a ridge
*/
void qh_printridge(FILE *fp, ridgeT *ridge) {
  
  fprintf(fp, "     - r%d", ridge->id);
  if (ridge->tested)
    fprintf (fp, " tested");
  if (ridge->nonconvex)
    fprintf (fp, " nonconvex");
  fprintf (fp, "\n");
  qh_printvertices (fp, "           vertices:", ridge->vertices);
  if (ridge->top && ridge->bottom)
    fprintf(fp, "           between f%d and f%d\n",
	    ridge->top->id, ridge->bottom->id);
} /* printridge */

/*-----------------------------------------
-printspheres- prints 3-d vertices as OFF spheres
  inflated octahedron from Stuart Levy earth/mksphere2
*/
void qh_printspheres(FILE *fp, setT *vertices, realT radius) {
  vertexT *vertex, **vertexp;
  pointT *point;

  qh printoutnum++;
  fprintf (fp, "{appearance {-edge -normal normscale 0} {\n\
INST geom {define vsphere OFF\n\
18 32 48\n\
\n\
0 0 1\n\
1 0 0\n\
0 1 0\n\
-1 0 0\n\
0 -1 0\n\
0 0 -1\n\
0.707107 0 0.707107\n\
0 -0.707107 0.707107\n\
0.707107 -0.707107 0\n\
-0.707107 0 0.707107\n\
-0.707107 -0.707107 0\n\
0 0.707107 0.707107\n\
-0.707107 0.707107 0\n\
0.707107 0.707107 0\n\
0.707107 0 -0.707107\n\
0 0.707107 -0.707107\n\
-0.707107 0 -0.707107\n\
0 -0.707107 -0.707107\n\
\n\
3 0 6 11\n\
3 0 7 6	\n\
3 0 9 7	\n\
3 0 11 9\n\
3 1 6 8	\n\
3 1 8 14\n\
3 1 13 6\n\
3 1 14 13\n\
3 2 11 13\n\
3 2 12 11\n\
3 2 13 15\n\
3 2 15 12\n\
3 3 9 12\n\
3 3 10 9\n\
3 3 12 16\n\
3 3 16 10\n\
3 4 7 10\n\
3 4 8 7\n\
3 4 10 17\n\
3 4 17 8\n\
3 5 14 17\n\
3 5 15 14\n\
3 5 16 15\n\
3 5 17 16\n\
3 6 13 11\n\
3 7 8 6\n\
3 9 10 7\n\
3 11 12 9\n\
3 14 8 17\n\
3 15 13 14\n\
3 16 12 15\n\
3 17 10 16\n} transforms { TLIST\n");
  FOREACHvertex_(vertices) {
    point= vertex->point;
    fprintf(fp, "%8.4g 0 0 0 # v%d\n 0 %8.4g 0 0\n0 0 %8.4g 0\n",
      radius, vertex->id, radius, radius);
    qh_printpoint3 (fp, vertex->point);
    fprintf (fp, "1\n");
  }
  fprintf (fp, "}}}\n");
} /* printspheres */


/*----------------------------------------------
-printsummary-
                see qhull.c
*/

/*-------------------------------------------
-printvertex- prints the information in a vertex
*/
void qh_printvertex(FILE *fp, vertexT *vertex) {
  pointT *point;
  int k;
  facetT *neighbor, **neighborp;
  realT r; /*bug fix*/

  if (!vertex) {
    fprintf (fp, "  NULLvertex\n");
    return;
  }
  fprintf(fp, "- p%d (v%d):", qh_pointid(vertex->point), vertex->id);
  point= vertex->point;
  if (point) {
    for(k= qh hull_dim; k--; )
      fprintf(fp, " %5.2g", r=*point++);
  }
  if (vertex->deleted)
    fprintf(fp, " deleted");
  if (vertex->delridge)
    fprintf (fp, " ridgedeleted");
  fprintf(fp, "\n");
  if (vertex->neighbors) {
    fprintf(fp, "  neighbors:");
    FOREACHneighbor_(vertex)
      fprintf(fp, " f%d", neighbor->id);
    fprintf(fp, "\n");
  }
} /* printvertex */


/*-------------------------------------------
-printvertexlist- prints vertices used by a facetlist or facet set
  all facets if printall, else !qh_skipfacet
*/
void qh_printvertexlist (FILE *fp, char* string, facetT *facetlist, 
                         setT *facets, boolT printall) {
  vertexT *vertex, **vertexp;
  setT *vertices;
  
  vertices= qh_facetvertices (facetlist, facets, printall);
  fputs (string, fp);
  FOREACHvertex_(vertices)
    qh_printvertex(fp, vertex);
  qh_settempfree (&vertices);
} /* printvertexlist */


/*-------------------------------------------
-printvertices- prints vertices in a set
*/
void qh_printvertices(FILE *fp, char* string, setT *vertices) {
  vertexT *vertex, **vertexp;
  
  fputs (string, fp);
  FOREACHvertex_(vertices) 
    fprintf (fp, " p%d (v%d)", qh_pointid(vertex->point), vertex->id);
  fprintf(fp, "\n");
} /* printvertices */

/*-------------------------------------------
-printvneighbors- print vertex neighbors of vertices in facetlist and facets
*/
void qh_printvneighbors (FILE *fp, facetT* facetlist, setT *facets, boolT printall) {
  int numfacets, numsimplicial, numridges, totneighbors, numneighbors, numcoplanars;
  setT *vertices, *vertex_points, *coplanar_points;
  int numpoints= qh num_points + qh_setsize (qh other_points);
  vertexT *vertex, **vertexp;
  int vertex_i, vertex_n;
  facetT *facet, **facetp, *neighbor, **neighborp;
  pointT *point, **pointp;

  qh_countfacets (facetlist, facets, printall, &numfacets, &numsimplicial, 
      &totneighbors, &numridges, &numcoplanars);  /* sets facet->visitid */
  fprintf (fp, "%d\n", numpoints);
  qh_vertexneighbors();
  vertices= qh_facetvertices (facetlist, facets, printall);
  vertex_points= qh_settemp (numpoints);
  coplanar_points= qh_settemp (numpoints);
  qh_setzero (vertex_points, 0, numpoints);
  qh_setzero (coplanar_points, 0, numpoints);
  FOREACHvertex_(vertices)
    qh_point_add (vertex_points, vertex->point, vertex);
  FORALLfacet_(facetlist) {
    FOREACHpoint_(facet->coplanarset)
      qh_point_add (coplanar_points, point, facet);
  }
  FOREACHfacet_(facets) {
    FOREACHpoint_(facet->coplanarset)
      qh_point_add (coplanar_points, point, facet);
  }
  FOREACHvertex_i_(vertex_points) {
    if (vertex) { 
      numneighbors= qh_setsize (vertex->neighbors);
      fprintf (fp, "%d", numneighbors);
      if (qh hull_dim == 3)
        qh_order_vertexneighbors (vertex);
      else if (qh hull_dim >= 4)
        qsort (SETaddr_(vertex->neighbors, vertexT), numneighbors,
             sizeof (facetT *), qh_compare_facetvisit);
      FOREACHneighbor_(vertex) 
        fprintf (fp, " %d", 
		 neighbor->visitid ? neighbor->visitid - 1 : - neighbor->id);
      fprintf (fp, "\n");
    }else if ((facet= (facetT*)SETelem_(coplanar_points, vertex_i)))
      fprintf (fp, "1 %d\n",
                 facet->visitid ? facet->visitid - 1 : - facet->id);
    else
      fprintf (fp, "0\n");
  }
  qh_settempfree (&coplanar_points);
  qh_settempfree (&vertex_points);
  qh_settempfree (&vertices);
} /* printvneighbors */

/*----------------------------------------
-printvoronoi- print voronoi diagram in 'o' or 'G' format
  for 'o' format
    prints voronoi centers for each facet and for infinity
    for each vertex, lists ids of printed facets or infinity
    assumes facetlist and facets are disjoint
  for 'G' format
    prints an OFF object
    adds a 0 coordinate to center
    prints infinity but does not list in vertices
notes:
  if 'o', prints a line for each point except "at-infinity"
  if all facets are upperdelaunay, reverses lower and upper hull
*/
void qh_printvoronoi (FILE *fp, int format, facetT *facetlist, setT *facets, boolT printall) {
  int k, numcenters=0, numvertices= 0, numneighbors, numinf, vid=1, vertex_i, vertex_n;
  facetT *facet, **facetp, *neighbor, **neighborp;
  setT *vertices;
  vertexT *vertex;
  boolT islower= False;

  qh printoutnum++;
  qh_vertexneighbors();
  vertices= qh_pointvertex();
  if (qh ATinfinity) 
    SETelem_(vertices, qh num_points-1)= NULL;
  qh visit_id++;
  maximize_(qh visit_id, qh num_facets);
  FORALLfacet_(facetlist) {  /* FIXUP: could merge with below */
    if (printall || !qh_skipfacet (facet)) {
      if (!facet->upperdelaunay)
        islower= True;
    }
  }
  FOREACHfacet_(facets) {
    if (printall || !qh_skipfacet (facet)) {
      if (!facet->upperdelaunay)
        islower= True;
    }
  }
  FORALLfacets {
    if (facet->normal && (facet->upperdelaunay == islower))
      facet->visitid= 0;  /* facetlist or facets may overwrite */
    else
      facet->visitid= qh visit_id;
    facet->seen= False;
  }
  numcenters++;  /* qh_INFINITE */
  FORALLfacet_(facetlist) {
    if (printall || !qh_skipfacet (facet)) {
      facet->visitid= numcenters++;
      facet->seen= True;
    }
  }
  FOREACHfacet_(facets) {
    if (printall || !qh_skipfacet (facet)) {
      facet->visitid= numcenters++;  
      facet->seen= True;
    }
  }
  FOREACHvertex_i_(vertices) {
    if (vertex) {
      numvertices++;
      numneighbors = numinf = 0;
      FOREACHneighbor_(vertex) {
        if (neighbor->seen)
          numneighbors++;
        else if (neighbor->visitid == 0)
	  numinf= 1;
      }
      if (numinf && !numneighbors) {
	SETelem_(vertices, vertex_i)= NULL;
	numvertices--;
      }
    }
  }
  if (format == qh_PRINTgeom) 
    fprintf (fp, "{appearance {+edge -face} OFF %d %d 1 # Voronoi centers and cells\n", 
                numcenters, numvertices);
  else
    fprintf (fp, "%d\n%d %d 1\n", qh hull_dim-1, numcenters, qh_setsize(vertices));
  if (format == qh_PRINTgeom) {
    for (k= qh hull_dim-1; k--; )
      fprintf (fp, qh_REAL_1, 0.0);
    fprintf (fp, " 0 # infinity not used\n");
  }else {
    for (k= qh hull_dim-1; k--; )
      fprintf (fp, qh_REAL_1, qh_INFINITE);
    fprintf (fp, "\n");
  }
  FORALLfacet_(facetlist) {
    if (facet->seen) {
      if (format == qh_PRINTgeom)
        fprintf (fp, "# %d f%d\n", vid++, facet->id);
      qh_printcenter (fp, format, NULL, facet);
    }
  }
  FOREACHfacet_(facets) {
    if (facet->seen) {
      if (format == qh_PRINTgeom)
        fprintf (fp, "# %d f%d\n", vid++, facet->id);
      qh_printcenter (fp, format, NULL, facet);
    }
  }
  FOREACHvertex_i_(vertices) {
    numneighbors= 0;
    numinf=0;
    if (vertex) {
      if (qh hull_dim == 3)
        qh_order_vertexneighbors(vertex);
      else if (qh hull_dim >= 4)
        qsort (SETaddr_(vertex->neighbors, vertexT), 
	     qh_setsize (vertex->neighbors),
	     sizeof (facetT *), qh_compare_facetvisit);
      FOREACHneighbor_(vertex) {
        if (neighbor->seen)
          numneighbors++;
        else if (neighbor->visitid == 0)
	  numinf= 1;
      }
    }
    if (format == qh_PRINTgeom) {
      if (vertex) {
	fprintf (fp, "%d", numneighbors);
	if (vertex) {
	  FOREACHneighbor_(vertex) {
	    if (neighbor->seen)
	      fprintf (fp, " %d", neighbor->visitid);
	  }
	}
	fprintf (fp, " # p%d (v%d)\n", vertex_i, vertex->id);
      }else
	fprintf (fp, " # p%d is coplanar or isolated\n", vertex_i);
    }else {
      if (numinf)
	numneighbors++;
      fprintf (fp, "%d", numneighbors);
      if (vertex) {
        FOREACHneighbor_(vertex) {
	  if (neighbor->seen)
	    fprintf (fp, " %d", neighbor->visitid);
  	  else if (numinf && neighbor->visitid == 0) {
	    numinf= 0;
	    fprintf (fp, " %d", neighbor->visitid);
	  }
	}
      }
      fprintf (fp, "\n");
    }
  }
  if (format == qh_PRINTgeom)
    fprintf (fp, "}\n");
  qh_settempfree (&vertices);
} /* printvoronoi */
  
/*-------------------------------------------
-projectdim3 -- project 2-d 3-d or 4-d point to a 3-d point
  uses DROPdim and hull_dim
  source and destination may be the same
  allocate 4 elements just in case
*/
void qh_projectdim3 (pointT *source, pointT *destination) {
  int i,k;

  
  for (k= 0, i=0; k<qh hull_dim; k++) {
    if (qh hull_dim == 4) {
      if (k != qh DROPdim)
        destination[i++]= source[k];
    }else if (k == qh DROPdim)
      destination[i++]= 0;
    else
      destination[i++]= source[k];
  }
  while (i < 3)
    destination[i++]= 0.0;
} /* projectdim3 */

/*-------------------------------------------
-readfeasible- read feasible point from remainder and qh fin
  checks for qh HALFspace
  assumes dim > 1
returns:
  number of lines read from qh fin
  sets qh FEASIBLEpoint with malloc'd coordinates
  see qh_setfeasible
*/
int qh_readfeasible (int dim, char *remainder) {
  boolT isfirst= True;
  int linecount= 0, tokcount= 0;
  char *s, *t, firstline[qh_MAXfirst+1];
  coordT *coords, value;

  if (!qh HALFspace) {
    fprintf  (qh ferr, "qhull input error: feasible point (dim 1 coords) is only valid for halfspace intersection\n");
    qh_errexit (qh_ERRinput, NULL, NULL);
  }  
  if (qh feasible_string)
    fprintf  (qh ferr, "qhull input warning: feasible point (dim 1 coords) overrides 'Hn,n,n' feasible point for halfspace intersection\n");
  if (!(qh feasible_point= (coordT*)malloc (dim* sizeof(coordT)))) {
    fprintf(qh ferr, "qhull error: insufficient memory for feasible point\n");
    qh_errexit(qh_ERRmem, NULL, NULL);
  }
  coords= qh feasible_point;
  while ((s= (isfirst ?  remainder : fgets(firstline, qh_MAXfirst, qh fin)))) {
    if (isfirst)
      isfirst= False;
    else
      linecount++;
    while (*s) {
      while (isspace(*s))
        s++;
      value= qh_strtod (s, &t);
      if (s == t)
        break;
      s= t;
      *(coords++)= value;
      if (++tokcount == dim) {
        while (isspace (*s))
          s++;
        qh_strtod (s, &t);
        if (s != t) {
          fprintf (qh ferr, "qhull input error: coordinates for feasible point do not finish out the line: %s\n",
               s);
          qh_errexit (qh_ERRinput, NULL, NULL);
        }
        return linecount;
      }
    }
  }
  fprintf (qh ferr, "qhull input error: short input.  Could not read feasible point.\n");
  qh_errexit (qh_ERRinput, NULL, NULL);
  return 0;
} /* readfeasible */

/*-------------------------------------------
-readpoints- read points from qh fin into all_points
    qh fin is lines of coordinates, one per vertex, first line number of points
    if 'rbox D4' gives message
returns:
    number of points, array of point coordinates, dimension, ismalloc True
    if DELAUNAY & !PROJECTinput, projects points to paraboloid
        and clears PROJECTdelaunay
    if HALFspace, reads optional feasible point, reads half-spaces,
        converts to dual.
notes:
    dimension will change in qh_initqhull_globals if PROJECTinput
    uses malloc since qh_mem not initialized
    FIXUP: this routine needs rewriting
*/
coordT *qh_readpoints(int *numpoints, int *dimension, boolT *ismalloc) {
  coordT *points, *coords, *infinity= NULL;
  realT paraboloid, maxboloid= -REALmax, value;
  realT *coordp= NULL, *offsetp= NULL, *normalp= NULL;
  char *s, *t, firstline[qh_MAXfirst+1];
  int diminput=0, numinput=0, dimfeasible= 0, newnum, k, tempi;
  int firsttext=0, firstshort=0, firstlong=0, firstpoint=0;
  int tokcount= 0, linecount=0, maxcount, coordcount=0;
  boolT islong, isfirst= True, wasbegin= False;
  boolT isdelaunay= qh DELAUNAY && !qh PROJECTinput;

  if (qh CDDinput) {
    while ((s= fgets(firstline, qh_MAXfirst, qh fin))) {
      linecount++;
      if (!memcmp (firstline, "begin", 5) || !memcmp (firstline, "BEGIN", 5))
        break;
      if (!*qh rbox_command) {
	strncat(qh rbox_command, s, sizeof (qh rbox_command)-1);
        qh rbox_command[strlen(qh rbox_command)-1]= '\0';
      }
    }
    if (!s) {
      fprintf (qh ferr, "qhull input error: missing \"begin\" for cdd-formated input\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
  }
  while(!numinput && (s= fgets(firstline, qh_MAXfirst, qh fin))) {
    linecount++;
    if (!memcmp (s, "begin", 5) || !memcmp (s, "BEGIN", 5))
      wasbegin= True;
    while (*s) {
      while (isspace(*s))
        s++;
      if (!*s)
        break;
      if (!isdigit(*s)) {
        if (!*qh rbox_command) {
          strncat(qh rbox_command, s, sizeof (qh rbox_command)-1);
          qh rbox_command[strlen(qh rbox_command)-1]= '\0';
	  firsttext= linecount;
        }
        break;
      }
      if (!diminput) 
        diminput= qh_strtol (s, &s);
      else {
        numinput= qh_strtol (s, &s);
        if (numinput == 1 && diminput >= 2 && qh HALFspace) {
          linecount += qh_readfeasible (diminput, s); /* checks if ok */
          dimfeasible= diminput;
          diminput= numinput= 0;
        }else 
          break;
      }
    }
  }
  if (!s) {
    fprintf(qh ferr, "qhull input error: short input file.  Did not find dimension and number of points\n");
    qh_errexit(qh_ERRinput, NULL, NULL);
  }
  if (diminput > numinput) {
    tempi= diminput;	/* exchange dim and n, e.g., for cdd input format */
    diminput= numinput;
    numinput= tempi;
  }
  if (diminput < 2) {
    fprintf(qh ferr,"qhull input error: dimension %d (first number) should be at least 2\n",
	    diminput);
    qh_errexit(qh_ERRinput, NULL, NULL);
  }
  if (!strcmp (qh rbox_command, "./rbox D4")) 
    fprintf (qh ferr, "\n\
This is the qhull test case.  If any errors or core dumps occur,\n\
recompile qhull with 'make new'.  If errors still occur, there is\n\
an incompatibility.  You should try a different compiler.  You can also\n\
change the choices in user.h.  If you discover the source of the problem,\n\
please send mail to qhull_bug@geom.umn.edu.\n\
\n\
Type 'qhull' for a short list of options.\n");
  if (isdelaunay) {
    qh PROJECTdelaunay= False;
    *dimension= diminput+1;
    *numpoints= numinput;
    if (qh ATinfinity)
      (*numpoints)++;
  }else if (qh HALFspace) {
    *dimension= diminput - 1;
    *numpoints= numinput;
    if (diminput < 3) {
      fprintf(qh ferr,"qhull input error: dimension %d (first number, includes offset) should be at least 3 for halfspaces\n",
  	    diminput);
      qh_errexit(qh_ERRinput, NULL, NULL);
    }
    if (dimfeasible) {
      if (dimfeasible != *dimension) {
        fprintf(qh ferr,"qhull input error: dimension %d of feasible point is not one less than dimension %d for halfspaces\n",
          dimfeasible, diminput);
        qh_errexit(qh_ERRinput, NULL, NULL);
      }
    }else 
      qh_setfeasible (*dimension);
  }else if (qh CDDinput) {
    *dimension= diminput-1;
    *numpoints= numinput;
  }else {
    *dimension= diminput;
    *numpoints= numinput;
  }
  qh normal_size= *dimension * sizeof(coordT); /* for tracing with qh_printpoint */
  if (qh HALFspace) {
    qh half_space= coordp= (coordT*) malloc (qh normal_size + sizeof(coordT));
    if (qh CDDinput) {
      offsetp= qh half_space;
      normalp= offsetp + 1;
    }else {
      normalp= qh half_space;
      offsetp= normalp + *dimension;
    }
  } 
  qh maxline= diminput * (qh_REALdigits + 5);
  maximize_(qh maxline, 500);
  qh line= malloc ((qh maxline+1) * sizeof (char));
  *ismalloc= True;  /* use malloc since memory not setup */
  coords= points= qh temp_malloc= 
        (coordT*)malloc(*numpoints*(*dimension)*sizeof(coordT));
  if (!coords || !qh line || (qh HALFspace && !qh half_space)) {
    fprintf(qh ferr, "qhull error: insufficient memory to read %d points\n",
	    numinput);
    qh_errexit(qh_ERRmem, NULL, NULL);
  }
  if (isdelaunay && qh ATinfinity) {
    infinity= points + numinput * (*dimension);
    for (k= diminput; k--; )
      infinity[k]= 0.0;
  }
  maxcount= numinput * diminput;
  paraboloid= 0.0;
  while ((s= (isfirst ?  s : fgets(qh line, qh maxline, qh fin)))) {
    if (!isfirst) {
      linecount++;
      if (*s == 'e' || *s == 'E') {
	if (!memcmp (s, "end", 3) || !memcmp (s, "END", 3)) {
	  if (qh CDDinput )
	    break;
	  else if (wasbegin) 
	    fprintf (qh ferr, "qhull input warning: the input appears to be in cdd format.  If so, use 'Fd'\n");
	}
      }
    }
    islong= False;
    while (*s) {
      while (isspace(*s))
        s++;
      value= qh_strtod (s, &t);
      if (s == t) {
        if (*s && !firsttext) 
          firsttext= linecount;
        if (!islong && !firstshort && coordcount)
          firstshort= linecount;
        break;
      }
      if (!firstpoint)
	firstpoint= linecount;
      s= t;
      if (++tokcount > maxcount)
        continue;
      if (qh HALFspace) {
	if (qh CDDinput && !coordcount) 
	  *(coordp++)= -value; /* offset */
	else
	  *(coordp++)= value;
      }else {
        *(coords++)= value;
        if (qh CDDinput && !coordcount) {
          if (value != 1.0) {
            fprintf (qh ferr, "qhull input error: for cdd format, point at line %d does not start with '1'\n",
                   linecount);
            qh_errexit (qh_ERRinput, NULL, NULL);
          }
          coords--;
        }
      }
      if (isdelaunay) {
	paraboloid += value * value;
	if (qh ATinfinity)
	  infinity[coordcount] += value;
      }
      if (++coordcount == diminput) {
        coordcount= 0;
        if (isdelaunay) {
          *(coords++)= paraboloid;
          maximize_(maxboloid, paraboloid);
          paraboloid= 0.0;
        }else if (qh HALFspace) {
          if (!qh_sethalfspace (*dimension, coords, &coords, normalp, offsetp, qh feasible_point)) {
	    fprintf (qh ferr, "The halfspace was on line %d\n", linecount);
	    if (wasbegin)
	      fprintf (qh ferr, "The input appears to be in cdd format.  If so, you should use option 'Fd'\n");
	    qh_errexit (qh_ERRinput, NULL, NULL);
	  }
          coordp= qh half_space;
        }          
        while (isspace(*s))
          s++;
        if (*s) {
          islong= True;
          if (!firstlong)
            firstlong= linecount;
	}
      }
    }
    if (!islong && !firstshort && coordcount)
      firstshort= linecount;
    if (!isfirst && s - qh line >= qh maxline) {
      fprintf(qh ferr, "qhull input error: line %d contained more than %d characters\n", 
	      linecount, (int) (s - qh line));
      qh_errexit(qh_ERRinput, NULL, NULL);
    }
    isfirst= False;
  }
  if (tokcount != maxcount) {
    newnum= fmin_(numinput, tokcount/diminput);
    fprintf(qh ferr,"\
qhull warning: instead of %d %d-dimensional points, input contains\n\
%d points and %d extra coordinates.  Line %d is the first point,\n\
line %d is the first (if any) comment, line %d is the first short line,\n\
and line %d is the first long line.  Continue with %d points.\n",
       numinput, diminput, tokcount/diminput, tokcount % diminput, firstpoint,
       firsttext, firstshort, firstlong, newnum);
    numinput= newnum;
    if (isdelaunay && qh ATinfinity) {
      for (k= tokcount % diminput; k--; )
	infinity[k] -= *(--coords);
    }else
      coords -= tokcount % diminput;
  }
  if (isdelaunay && qh ATinfinity) {
    for (k= diminput; k--; )
      infinity[k] /= numinput;
    if (coords == infinity)
      coords += diminput;
    else {
      for (k= 0; k < diminput; k++)
	*(coords++)= infinity[k];
    }
    *(coords++)= maxboloid * 1.1;
    *numpoints= numinput+1;
  }else
    *numpoints= numinput;
  free (qh line);
  qh line= NULL;
  if (qh half_space) {
    free (qh half_space);
    qh half_space= NULL;
  }
  qh temp_malloc= NULL;
  trace1((qh ferr,"qh_readpoints: read in %d %d-dimensional points\n",
	  numinput, diminput));
  return(points);
} /* readpoints */


/*-------------------------------------------
-setfeasible- set feasible point from qh feasible_string
returns:
  sets qh FEASIBLEpoint with malloc'd coordinates
notes:
  "n,n,n" already checked by qh_initflags
  see qh_readfeasible
*/
void qh_setfeasible (int dim) {
  int tokcount= 0;
  char *s;
  coordT *coords, value;

  if (!(s= qh feasible_string)) {
    fprintf(qh ferr, "\
qhull input error: halfspace intersection needs a feasible point.\n\
Either prepend the input with 1 point or use 'Hn,n,n'.  See manual.\n");
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  if (!(qh feasible_point= malloc (dim* sizeof(coordT)))) {
    fprintf(qh ferr, "qhull error: insufficient memory for 'Hn,n,n'\n");
    qh_errexit(qh_ERRmem, NULL, NULL);
  }
  coords= qh feasible_point;
  while (*s) {
    value= qh_strtod (s, &s);
    if (++tokcount > dim) {
      fprintf (qh ferr, "qhull input warning: more coordinates for 'H%s' than dimension %d\n",
          qh feasible_string, dim);
      break;
    }
    *(coords++)= value;
    if (*s)
      s++;
  }
  while (++tokcount <= dim)    
    *(coords++)= 0.0;
} /* setfeasible */

/*-------------------------------------------
-skipfacet- returns 'True' if this facet is not to be printed 
  (based on the user provided slice thresholds and 'good' specifications)
*/
boolT qh_skipfacet(facetT *facet) {
  facetT *neighbor, **neighborp;

  if (qh PRINTneighbors) {
    if (facet->good)
      return !qh PRINTgood;
    FOREACHneighbor_(facet) {
      if (neighbor->good)
	return False;
    }
    return True;
  }else if (qh PRINTgood)
    return !facet->good;
  else if (!facet->normal)
    return True;
  return (!qh_inthresholds (facet->normal, NULL));
} /* skipfacet */

