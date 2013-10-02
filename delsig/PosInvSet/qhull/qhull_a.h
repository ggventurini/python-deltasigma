/* qhull_a.h -- all header files for compiling qhull

   see README 
   
   see qhull.h for user-level definitions
   
   see user.h for user-defineable constants
   
   defines internal functions for qhull.c global.c

   copyright (c) 1993-1995, The Geometry Center

   Notes:  grep for ((" and (" to catch fprintf("lkasdjf");
           full parens around (x?y:z)
           
*/

#ifndef qhDEFqhulla
#define qhDEFqhulla



#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#include <string.h>
#include <math.h>
#include <float.h>    /* some compilers will not need float.h */
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include <memory.h>

#include "qhull.h"


#include "mem.h"
#include "set.h"
#include "geom.h"
#include "merge.h"
#include "poly.h"
#include "io.h"
#include "stat.h"


/* ======= -macros- =========== */

/*-----------------------------------------------
-traceN((fp.ferr, "format\n", vars));  calls fprintf if IStracing >= N
  removing tracing reduces code size but doesn't change execution speed
*/
#ifndef qh_NOtrace
#define trace0(args) {if (qh IStracing) fprintf args;}
#define trace1(args) {if (qh IStracing >= 1) fprintf args;}
#define trace2(args) {if (qh IStracing >= 2) fprintf args;}
#define trace3(args) {if (qh IStracing >= 3) fprintf args;}
#define trace4(args) {if (qh IStracing >= 4) fprintf args;}
#define trace5(args) {if (qh IStracing >= 5) fprintf args;}
#else /* qh_NOtrace */
#define trace0(args) {}
#define trace1(args) {}
#define trace2(args) {}
#define trace3(args) {}
#define trace4(args) {}
#define trace5(args) {}
#endif /* qh_NOtrace */

/* ======= -functions =========== 

see corresponding .c file for definitions

	Qhull functions (qhull.c)
-qhull		construct the convex hull of a set of points
-qhull_postmerging  set up for post merging of qhull
-addpoint       add point to hull above a facet
-buildhull	constructs a hull by adding points one at a time
-buildtracing   for tracing execution of buildhull
-errexit2	return exitcode to system after an error for two facets
-findhorizon	find the horizon and visible facets for a point
-nextfurthest   returns next furthest point for processing
-partitionall	partitions all points into the outsidesets of facets
-partitioncoplanar partition coplanar point to a facet
-partitionpoint partitions a point as inside, coplanar or outside a facet
-partitionvisible partitions points in visible_list to newfacet_list

	Global.c internal functions (others in qhull.h)
-freebuffers	free up global memory buffers 
-initbuffers	initialize global memory buffers
-option         append option description to qh qhull_options
-strtod/tol     duplicates strtod/tol
*/

/***** -qhull.c prototypes (alphabetical after qhull) ********************/

void 	qh_qhull (void);
void    qh_qhull_postmerging (void);
boolT   qh_addpoint (pointT *furthest, facetT *facet, boolT checkdist);
void 	qh_buildhull(void);
void    qh_buildtracing (pointT *furthest, facetT *facet);
void 	qh_errexit2(int exitcode, facetT *facet, facetT *otherfacet);
void    qh_findhorizon(pointT *point, facetT *facet, int *goodvisible,int *goodhorizon);
pointT *qh_nextfurthest (facetT **visible);
void 	qh_partitionall(setT *vertices, pointT *points,int npoints);
void    qh_partitioncoplanar (pointT *point, facetT *facet, realT *dist);
void    qh_partitionpoint (pointT *point, facetT *facet);
void 	qh_partitionvisible(boolT allpoints, int *numpoints);
void	qh_printsummary(FILE *fp);

/***** -global.c internal prototypes (alphabetical) ***********************/

void    qh_appendprint (qh_PRINT format);
void 	qh_freebuffers (void);
void    qh_initbuffers (coordT *points, int numpoints, int dim, boolT ismalloc);
void    qh_option (char *option, int *i, realT *r);
int     qh_strtol (const char *s, char **endp);
double  qh_strtod (const char *s, char **endp);

/***** -stat.c internal prototypes (alphabetical) ***********************/

void	qh_allstatA (void);
void	qh_allstatB (void);
void	qh_allstatC (void);
void	qh_allstatD (void);
void	qh_allstatE (void);
void	qh_allstatF (void);
void 	qh_freebuffers (void);
void    qh_initbuffers (coordT *points, int numpoints, int dim, boolT ismalloc);

#endif /* qhDEFqhulla */
