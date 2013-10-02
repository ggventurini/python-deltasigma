/* merge.h -- header file for merge.c
   
   see README and merge.c
   
   copyright (c) 1993-1995, The Geometry Center

*/

#ifndef qhDEFmerge
#define qhDEFmerge 1


/* ============ -constants- ==============
*/
#define qh_ANGLEredundant 6.0 /* angle for redundant merge in mergeT */
#define qh_ANGLEdegen     5.0 /* angle for degenerate facet in mergeT */
#define qh_ANGLEconcave  1.5  /* [2,4] for angle of concave facets in mergeT,
                                 may be <2 or >4 due to roundoff */

#define qh_MERGEapex     True /* parameter to qh_mergefacet */

/* ============ -structures- ====================
*/

/* ----------------------------------------------
-mergeT- structure used to merge facets
*/

typedef struct mergeT mergeT;
typedef enum {	/* in sort order for facet_mergeset */
  MRGnone= 0,
  MRGcoplanar,		/* centrum coplanar */
  MRGanglecoplanar,	/* angle coplanar */
  			/* could detect half concave ridges */
  MRGconcave,		/* concave ridge */
  MRGflip,		/* flipped facet. facet1 == facet2 */
  MRGridge,		/* duplicate ridge (qh_MERGEridge) */
  		/* degen and redundant go onto degen_mergeset */
  MRGdegen,		/* degenerate facet (not enough neighbors) facet1 == facet2 */
  MRGredundant,		/* redundant facet (vertex subset) */
  			/* merge_degenredundant assumes degen < redundant */
  ENDmrg
} mergeType;

struct mergeT {		/* initialize in qh_appendmergeset */
  realT   angle;        /* angle between normals of facet1 and facet2 */
  facetT *facet1; 	/* will merge facet1 into facet2 */
  facetT *facet2;
  mergeType type;
};


/* =========== -macros- =========================
-FOREACHmerge-  if qh_mergefacet() then must restart since facet_mergeset
                may change.
*/
#define FOREACHmerge_(merges) FOREACHsetelement_(mergeT, merges, merge)

/* ======= -functions and procedures- =========== 

	top-level merge functions
-checkzero         check that facets are clearly convex
-premerge          pre-merge nonconvex facets in newfacet_list for apex
-postmerge         post-merge nonconvex facets as defined by maxcentrum/maxangle
-all_merges        merge all non-convex facets
-merge_nonconvex   merge nonconvex ridge
-flippedmerges	   merge flipped facets into best neighbor
-forcedmerges	     merge across duplicated ridges
-tracemerge        print trace message after merge
-tracemerging      print trace message during post merging

	mergeset functions for identifying merges
-mark_dupridges      add duplicated ridges to facet_mergeset
-getmergeset_initial initial mergeset for facets
-getmergeset	     returns facet_mergeset of facet-neighbor pairs to merged
-degen_redundant_facet      check facet for degen. or redundancy
-degen_redundant_neighbors  append degen/redundant neighbors to facet_mergeset
-test_vneighbors     test vertex neighbors for convexity
-test_appendmerge    facet/neighbor and appends to mergeset if nonconvex
-appendmergeset	     appends an entry to facet_mergeset, angle is optional
-facetdegen	     true if facet already in mergeset as a degenerate

	functions for determining the best merge
-findbest_test	   test neighbor for findbestneighbor()
-findbestneighbor  finds best neighbor (least dist) of a facet for merging

	functions for merging facets
-merge_degenredundant  merge degenerate and redundant facets
-mergefacet	   merges facet1 into facet2
-makeridges	   creates explicit ridges between simplicial facets
-mergesimplex      merge simplicial facet1 into facet2
-mergeneighbors	   merges the neighbors of facet1 into facet2
-mergeridges	   merges the ridge set of facet1 into facet2
-mergevertex_del   delete a vertex because of merging facet1 into facet2
-mergevertex_neighbors merge the vertex neighbors of facet1 to facet2
-mergevertices	   merges the vertex set of facet1 into facet2
-mergevertices2d   merges vertices1 into vertices2 in 2-d case
-mungevertices2d   munges vertices1 into vertices2 in 2-d case  RS951221
-newvertices       register all vertices as new vertices
-updatetested      clear facet2->tested/center and ridge->tested for merge
-willdelete        moves facet to visible list, sets replacement or NULL

	functions for merging a cycle of facets
-basevertices       return temporary set of base vertices for samecycle
-mergecycle_all    merge all facets in ->f.samecycle into facet with ->normal
-mergecycle-       merge a ->f.samecycle into newfacet with ->normal
-mergecycle_neighbors  add neighbors for ->f.samecycle to newfacet
-mergecycle_ridges add ridges for ->f.samecycle to newfacet
-mergecycle_facets finish merge of samecycle into newfacet
-mergecycle_vneighbors create vertex neighbors for newfacet

	functions for renaming a vertex
-reducevertices    reduce vertex sets
-rename_sharedvertex  detect and rename if shared vertex in facet
-redundant_vertex  returns true if detect and rename redundant vertex
-renamevertex	     renames oldvertex as newvertex in ridges 
-renameridgevertex renames oldvertex as newvertex in ridge
-maydropneighbor   drop neighbor relationship if no ridge between facet and neighbor
-remove_extravertices  remove extra vertices in non-simplicial facets
-copynonconvex     copy non-convex flag to all ridges between same neighbors

	functions for identifying vertices for renaming
-find_newvertex    locate new vertex for renaming old vertex
-neighbor_intersections	 return intersection for vertex->neighbors
-vertexridges	     return temporary set of ridges adjacent to a vertex
-vertexridges_facet add adjacent ridges for vertex in facet
-hashridge	       add ridge to hashtable without oldvertex
-hashridge_find	   returns matching ridge in hashtable without oldvertex

	check functions
-checkconnect      check that new facets are connected
-checkridge_boundary  checks that ridges of a facet are boundaryless
*/

/*---------- -prototypes in alphabetical order after pre/postmerge -----------*/

void    qh_premerge (vertexT *apex, realT maxcentrum, realT maxangle);
void    qh_postmerge (char *reason, realT maxcentrum, realT maxangle, 
             boolT vneighbors);
void    qh_all_merges (boolT othermerge, boolT vneighbors);
void    qh_appendmergeset(facetT *facet, facetT *neighbor, mergeType mergetype, realT *angle);
setT   *qh_basevertices( facetT *samecycle);
void    qh_checkconnect (void /* qh new_facets */);
boolT   qh_checkzero (boolT testall);
void    qh_copynonconvex (ridgeT *atridge);
void    qh_degen_redundant_facet (facetT *facet);
void   	qh_degen_redundant_neighbors (facetT *facet, facetT *delfacet);
vertexT *qh_find_newvertex (vertexT *oldvertex, setT *vertices, setT *ridges);
void    qh_findbest_test (boolT testcentrum, facetT *facet, facetT *neighbor,
           facetT **bestfacet, realT *distp, realT *mindistp, realT *maxdistp);
facetT *qh_findbestneighbor(facetT *facet, realT *distp, realT *mindistp, realT *maxdistp);
void 	qh_flippedmerges(facetT *facetlist, boolT *wasmerge);
void 	qh_forcedmerges( boolT *wasmerge);
void	qh_getmergeset(facetT *facetlist);
void 	qh_getmergeset_initial (facetT *facetlist);
void    qh_hashridge (setT *hashtable, int hashsize, ridgeT *ridge, vertexT *oldvertex);
ridgeT *qh_hashridge_find (setT *hashtable, int hashsize, ridgeT *ridge, 
              vertexT *vertex, vertexT *oldvertex, int *hashslot);
void 	qh_makeridges(facetT *facet);
void    qh_mark_dupridges(facetT *facetlist);
void    qh_maydropneighbor (facetT *facet);
int     qh_merge_degenredundant (void);
void    qh_merge_nonconvex( facetT *facet1, facetT *facet2, mergeType mergetype);
void    qh_mergecycle (facetT *samecycle, facetT *newfacet);
void    qh_mergecycle_all (facetT *facetlist, boolT *wasmerge);
void    qh_mergecycle_facets( facetT *samecycle, facetT *newfacet);
void    qh_mergecycle_neighbors(facetT *samecycle, facetT *newfacet);
void    qh_mergecycle_ridges(facetT *samecycle, facetT *newfacet);
void    qh_mergecycle_vneighbors( facetT *samecycle, facetT *newfacet);
void 	qh_mergefacet(facetT *facet1, facetT *facet2, realT *mindist, realT *maxdist, boolT mergeapex);
void 	qh_mergeneighbors(facetT *facet1, facetT *facet2);
void 	qh_mergeridges(facetT *facet1, facetT *facet2);
void    qh_mergesimplex(facetT *facet1, facetT *facet2, boolT mergeapex);
void    qh_mergevertex_del (vertexT *vertex, facetT *facet1, facetT *facet2);
void    qh_mergevertex_neighbors(facetT *facet1, facetT *facet2);
void	qh_mergevertices(setT *vertices1, setT **vertices);
void 	qh_mergevertices2d(setT *vertices1, setT *vertices2);
void    qh_mungevertices2d(facetT *facet1, facetT *facet2); /* RS950202 */
setT   *qh_neighbor_intersections (vertexT *vertex);
void    qh_newvertices (setT *vertices);
boolT   qh_reducevertices (void);
vertexT *qh_redundant_vertex (vertexT *vertex);
boolT   qh_remove_extravertices (facetT *facet);
vertexT *qh_rename_sharedvertex (vertexT *vertex, facetT *facet);
void	qh_renameridgevertex(ridgeT *ridge, vertexT *oldvertex, vertexT *newvertex);
void    qh_renamevertex(vertexT *oldvertex, vertexT *newvertex, setT *ridges,
			facetT *oldfacet, facetT *neighborA);
boolT 	qh_test_appendmerge (facetT *facet, facetT *neighbor);
boolT   qh_test_vneighbors (void /* qh newfacet_list */);
void    qh_tracemerge (facetT *facet1, facetT *facet2);
void    qh_tracemerging (void);
void    qh_updatetested( facetT *facet1, facetT *facet2);
setT   *qh_vertexridges (vertexT *vertex);
void    qh_vertexridges_facet (vertexT *vertex, facetT *facet, setT **ridges);
void    qh_willdelete (facetT *facet, facetT *replace);

#endif /* qhDEFmerge */
