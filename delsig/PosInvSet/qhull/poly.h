/* poly.h -- header file for poly.c and poly2.c

   see README, qhull.h and poly.c

   copyright (c) 1993-1995, The Geometry Center

*/

#ifndef qhDEFpoly
#define qhDEFpoly 1

/*-----------------------------------------------
-constants-

	for calling checkconvex()
-ALGORITHMfault   flag for checkconvex for error during buildhull
-DATAfault        flag for checkconvex for error during initialhull

	set by matchneighbor, used by matchmatch and mark_dupridge
-DUPLICATEridge   flag in facet->neighbors to indicate duplicated ridge
-MERGEridge       flag in facet->neighbors to indicate merged ridge
*/

#define qh_ALGORITHMfault 0
#define qh_DATAfault 1

#define qh_DUPLICATEridge (facetT *) 1
#define qh_MERGEridge (facetT *) 2


/* ============ -structures- ====================
*/

/* ----------------------------------------------
-hashentryT- hash table entry for matching sub-ridges in makecone()
*/

typedef struct hashentryT hashentryT;

struct hashentryT {
  facetT     *facet;        /* facet */
  hashentryT *next;         /* next hash table entry for this bucket */
  unsigned    skipindex;    /* skipped vertex in facet, for orientation */
};

/* =========== -macros- ========================= 
*/

/* ----------------------------------------------
-FOREACH... and FORALL... -- standard for loops
  see qhull.h for notes
*/
#define FORALLfacet_(facetlist) if (facetlist) for(facet=(facetlist);facet && facet->next;facet=facet->next)
#define FORALLnew_facets for(newfacet=qh newfacet_list;newfacet && newfacet->next;newfacet=newfacet->next)
#define FORALLvertex_(vertexlist) for (vertex=(vertexlist);vertex && vertex->next;vertex= vertex->next)
#define FORALLvisible_facets for (visible=qh visible_list; visible && visible->visible; visible= visible->next)
/* FORALLsame - only for newfacets, same_cycle includes facet */
#define FORALLsame_(newfacet) for (same= newfacet->f.samecycle; same != newfacet; same= same->f.samecycle)
/* FORALLsame_cycle_ - only for newfacets, last same is newfacet */
#define FORALLsame_cycle_(newfacet) \
     for (same= newfacet->f.samecycle; \
         same; same= (same == newfacet ?  NULL : same->f.samecycle))

#define FOREACHentry_(entries) FOREACHsetelement_(hashentryT, entries, entry)
#define FOREACHvisible_(facets) FOREACHsetelement_(facetT, facets, visible)
#define FOREACHnewfacet_(facets) FOREACHsetelement_(facetT, facets, newfacet)
#define FOREACHvertexA_(vertices) FOREACHsetelement_(vertexT, vertices, vertexA)
#define FOREACHvertexreverse12_(vertices) FOREACHsetelementreverse12_(vertexT, vertices, vertex)


/* ======= -functions =========== 

see poly.c for definitions

	Facetlist functions
-appendfacet	    appends facet to end of qh facet_list,
-prependfacet	    prepends facet to start of facetlist
-removefacet	    unlinks facet from qh facet_list,
-resetlists	    reset newvertex_list, newfacet_list, visible_list
-initialhull	    construct the initial hull as a simplex of vertices
-nonupper           return first facet without upperdelaunay or flipped
-setvoronoi_all     compute Voronoi centers for all facets
-findfacet          find facet that is furthest below a point 
-findgood           identify good facets for qh ONLYgood
-findgood_all       identify good facets for qh PRINTgood

	Facet functions
-createsimplex	    creates a simplex of facets from a set of vertices

-makenewfacet	    creates a toporient? facet from vertices and apex
-makenewfacets	    make new facets from point, horizon facets, and visible facets
-makenewplanes      make new hyperplanes for facets
-makenew_nonsimplicial make new facets for ridges of visible facets
-makenew_simplicial make new facets for horizon neighbors
-attachnewfacets    attach new facets in qh newfacet_list to the horizon
-makeadjacencies    make adjacencies for non-simplicial facets

	Vertex, ridge, and point functions
-appendvertex	    appends vertex to end of qh vertex_list,
-removevertex	    unlinks vertex from qh vertex_list,
-point              return point for a point id, or NULL if unknown
-pointid            return id for a point, or -1 if not known
-nearvertex         return nearest vertex to point
-vertexintersect    intersects two vertex sets
-vertexintersect_new intersects two vertex sets
-facetintersect	    intersect simplicial facets
-isvertex	    true if point is in the vertex set
-vertexsubset	    returns True if vertexsetA is a subset of vertexsetB
-nextridge3d	    iterate each ridge and vertex for a 3d facet
-facet3vertex	    return oriented vertex set for 3-d facet
-vertexneighhbors   for each vertex in hull, determine facet neighbors
-pointfacet	    return temporary set of facets indexed by point id
-pointvertex	    return temporary set of vertices indexed by point id
-initialvertices    return non-singular set of initial vertices
-updatevertices     update vertex neighbors and delete interior vertices

	Hashtable functions
-addhash            add hash element to linear hash table if not already there
-newhashtable	    allocates a new qh hash_table
-gethash	    return hashvalue for a set with firstindex
-matchduplicates    match duplicate ridges in hashtable
-matchnewfacets	    match newfacets in to their newfacet neighbors
-matchneighbor      try to match subridge of newfacet with a neighbor
-matchvertices	    tests whether a facet and hashentry match at a ridge
-printhashtable	    print hash table

	Allocation and deallocation functions
-newfacet	    creates and allocates space for a facet
-newridge	    creates and allocates space for a ridge
-newvertex	    creates and allocates space for a vertex
-deletevisible	    delete visible facets and vertices
-delfacet	    frees up the memory occupied by a facet
-delridge	    deletes ridge from data structures it belongs to and frees up the
-delvertex	    deletes vertex and its memory
-clearcenters       clear old data from facet->center
	
	Check functions
-check_bestdist	    check that points are not outside their best facet
-check_maxout       updates max_outside, checks all points against bestfacet
-check_output	    performs the checks at the end of qhull algorithm
-check_point        check that point is not outside facet
-check_points	    checks that all points are inside all facets
-checkconvex	    check that each ridge in facetlist is convex
-checkfacet	    checks for consistency errors in facet
-checkflipped	    checks facet orientation to the interior point
-checkflipped_all   checks facet orientation for a facet list
-checkpolygon	    checks the correctness of the structure
-checkvertex        check vertex for consistency
-infiniteloop       report infinite loop error due to facet
-printlists         print out facet list for debugging
*/

/*---------- -prototypes poly.c in alphabetical order -----------*/
void    qh_appendfacet(facetT *facet);
void    qh_appendvertex(vertexT *vertex);
void 	qh_attachnewfacets (void);
boolT   qh_checkflipped (facetT *facet, realT *dist, boolT allerror);
void	qh_delfacet(facetT *facet);
void 	qh_deletevisible(/*qh visible_list, qh horizon_list*/);
setT   *qh_facetintersect (facetT *facetA, facetT *facetB, int *skipAp,int *skipBp, int extra);
unsigned qh_gethash (int hashsize, setT *set, int size, int firstindex, void *skipelem);
facetT *qh_makenewfacet(setT *vertices, boolT toporient, facetT *facet);
void    qh_makenewplanes ( void /* newfacet_list */);
facetT *qh_makenew_nonsimplicial (facetT *visible, vertexT *apex, int *numnew);
facetT *qh_makenew_simplicial (facetT *visible, vertexT *apex, int *numnew);
void    qh_matchneighbor (facetT *newfacet, int newskip, int hashsize,
			  int *hashcount);
void	qh_matchnewfacets (void);
boolT   qh_matchvertices (int firstindex, setT *verticesA, int skipA, 
			  setT *verticesB, int *skipB, boolT *same);
facetT *qh_newfacet(void);
ridgeT *qh_newridge(void);
int     qh_pointid (pointT *point);
void 	qh_removefacet(facetT *facet);
void 	qh_removevertex(vertexT *vertex);
void    qh_updatevertices (void);


/*---------- -prototypes poly2.c in alphabetical order -----------*/

void    qh_addhash (void* newelem, setT *hashtable, int hashsize, unsigned hash);
void 	qh_check_bestdist ();
void    qh_check_maxout (void);
void    qh_check_output (void);
void    qh_check_point (pointT *point, facetT *facet, realT *maxoutside, facetT **errfacet1, facetT **errfacet2);
void   	qh_check_points(void);
void 	qh_checkconvex(facetT *facetlist, int fault);
void    qh_checkfacet(facetT *facet, boolT newmerge, boolT *waserrorp);
void 	qh_checkflipped_all (facetT *facetlist);
void 	qh_checkpolygon(facetT *facetlist);
void    qh_checkvertex (vertexT *vertex);
void 	qh_clearcenters (qh_CENTER type);
void 	qh_createsimplex(setT *vertices);
void 	qh_delridge(ridgeT *ridge);
void    qh_delvertex (vertexT *vertex);
setT   *qh_facet3vertex (facetT *facet);
facetT *qh_findfacet (pointT *point, facetT *facet, 
           realT *dist, boolT *isoutside, int *numpart);
int 	qh_findgood (facetT *facetlist, int goodhorizon);
void 	qh_findgood_all (facetT *facetlist);
void    qh_infiniteloop (facetT *facet);
void 	qh_initialhull(setT *vertices);
setT   *qh_initialvertices(int dim, setT *maxpoints, pointT *points, int numpoints);
vertexT *qh_isvertex (pointT *point, setT *vertices);
vertexT *qh_makenewfacets (pointT *point /*horizon_list, visible_list*/);
void qh_matchduplicates (facetT *atfacet, int atskip, int hashsize, int *hashcount);
vertexT *qh_nearvertex (facetT *facet, pointT *point, realT *bestdistp);
int 	qh_newhashtable(int newsize);
vertexT *qh_newvertex(pointT *point);
ridgeT *qh_nextridge3d (ridgeT *atridge, facetT *facet, vertexT **vertexp);
facetT *qh_nonupper (facetT *facetlist);
pointT *qh_point (int id);
void 	qh_point_add (setT *set, pointT *point, void *elem);
setT   *qh_pointfacet (void /*qh facet_list*/);
setT   *qh_pointvertex (void /*qh facet_list*/);
void 	qh_prependfacet(facetT *facet, facetT **facetlist);
void	qh_printhashtable(FILE *fp);
void    qh_printlists (void);
void    qh_resetlists (boolT stats /*qh newvertex_list newfacet_list visible_list*/);
void    qh_setvoronoi_all (void);
void    qh_vertexintersect(setT **vertexsetA,setT *vertexsetB);
setT   *qh_vertexintersect_new(setT *vertexsetA,setT *vertexsetB);
void    qh_vertexneighbors (void /*qh facet_list*/);
boolT 	qh_vertexsubset(setT *vertexsetA, setT *vertexsetB);


#endif /* qhDEFpoly */






