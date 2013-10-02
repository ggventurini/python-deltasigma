/* qhull.h -- user-level header file for using qhull.a library
   
   see README, qhull_a.h
   
   see user.h for all user defineable constants

   copyright (c) 1993-1995, The Geometry Center

   defines qh_qh, global data structure for qhull.  
   
   NOTE: access to qh_qh is via the 'qh' macro.  This allows
   qh_qh to be either a pointer or a structure.  An example
   of using qh is "qh DROPdim" which accesses the DROPdim
   field of qh_qh.  Similarly, access to qh_qhstat is via
   the 'qhstat' macro.
   
   includes function prototypes for qhull.c, geom.c, global.c, io.c, user.c

   use mem.h for mem.c
   use set.h for set.c
   
   see unix.c for an example of using qhull.h
   
   recompile qhull if you change this file
*/

#ifndef qhDEFqhull
#define qhDEFqhull 1

/* =========================== -included files ============== */

#include <setjmp.h>
#include <float.h>
#include <time.h>

#if __MWERKS__
#include  <SIOUX.h>
#include  <Files.h>
#include	<Desk.h>
#endif

#include "user.h"      /* user defineable constants */

/* ============ -types- ==================== 

   Note: could use 'float' for data and 'double' for calculations (realT vs. coordT)
        This requires many type casts, and adjusted error bounds.
        Also C compilers may do expressions in double anyway.

*/


#define coordT realT  /* for stored coordinates and coefficients */
#define pointT coordT                 /* array of hull_dim coordinates */
#define flagT unsigned int            /* Boolean flag */
#define boolT unsigned int            /* Boolean value, nothing much is portable */
#ifdef False
#undef False
#endif
#ifdef True
#undef True
#endif
#define False 0
#define True 1

/* -CENTERtype to distinguish facet->center */
typedef enum {qh_ASnone= 0, qh_ASvoronoi, qh_AScentrum} qh_CENTER;

/* -output formats for printing (qh PRINTout).
notes:
  some of these names are similar to qh names.  The similar names are only
  used in switch statements in qh_printbegin() etc.
*/
typedef enum {qh_PRINTnone= 0, qh_PRINTarea, qh_PRINTaverage, 
  qh_PRINTcoplanars, qh_PRINTcentrums, 
  qh_PRINTfacets, qh_PRINTfacets_xridge, 
  qh_PRINTgeom, qh_PRINTids, qh_PRINTinner, qh_PRINTneighbors, 
  qh_PRINTnormals, qh_PRINTouter, qh_PRINTincidences, 
  qh_PRINTmathematica, qh_PRINTmerges, qh_PRINToff, 
  qh_PRINToptions, qh_PRINTpointintersect, qh_PRINTpointnearest,
  qh_PRINTpoints, qh_PRINTqhull,
  qh_PRINTsize, qh_PRINTsummary, qh_PRINTvertices, qh_PRINTvneighbors,
  qh_PRINTEND} qh_PRINT;

/*---------------------------
       arguments
*/
#define qh_ALL 	    True     /* argument for printall and checkall parameters*/


/*-------------------------------------------
-ERR - qhull exit codes, for indicating errors
*/
#define qh_ERRnone  0    /* no error occurred during qhull */
#define qh_ERRinput 1    /* input inconsistency */
#define qh_ERRsingular 2 /* singular input data */
#define qh_ERRprec  3    /* precision error */
#define qh_ERRmem   4    /* insufficient memory, matches mem.h */
#define qh_ERRqhull 5    /* internal error detected, matches mem.h */

/* ============ -structures- ====================
   each of the following structures is defined by a typedef
   all realT and coordT fields occur at the beginning of a structure
        (otherwise space may be wasted due to alignment)
   define all flags together and pack into 32-bit number
*/

typedef struct vertexT vertexT;
typedef struct ridgeT ridgeT;
typedef struct facetT facetT;
#ifndef DEFsetT
#define DEFsetT 1
typedef struct setT setT;          /* defined in set.h */
#endif

/* ----------------------------------------------
-facetT- specifies a facet. 

  qhull() generates the hull as a list of facets.  
*/

struct facetT {
#if !qh_COMPUTEfurthest
  coordT   furthestdist;/* distance to furthest point of outsideset */
#endif
#if qh_MAXoutside
  coordT   maxoutside;  /* max computed distance of point to facet
  			Before QHULLfinished this is an approximation
  			since maxdist not always set for mergefacet
			Actual outer plane is +DISTround and
			computed outer plane is +2*DISTround */
#endif
  coordT   offset;      /* exact offset of hyperplane from origin */ 
  coordT  *normal;      /* normal of hyperplane, hull_dim coefficients */
  union {               /* in order of testing */
   realT   area;        /* area of facet, only in io.c if  ->isarea */
   facetT *replace;	/*  replacement facet if ->visible and NEWfacets
  			     is NULL only if qh_mergedegen_redundant or interior */
   facetT *samecycle;   /*  cycle of facets from the same visible/horizon intersection,
   			     if ->newfacet */
   facetT *newcycle;    /*  in horizon facet, current samecycle of new facets */ 
  }f;
  coordT  *center;      /*  centrum for convexity, qh CENTERtype == qh_AScentrum */
      			/*  Voronoi center, qh CENTERtype == qh_ASvoronoi */
  facetT  *previous;    /* previous facet in the facet_list */
  facetT  *next;        /* next facet in the facet_list */
  setT    *vertices;    /* vertices for this facet, inverse sorted by id 
                           if simplicial, 1st vertex was apex/furthest */
  setT    *ridges;      /* explicit ridges for nonsimplicial facets.
  			   for simplicial facets, neighbors defines ridge */
  setT    *neighbors;   /* neighbors of the facet.  If simplicial, the kth 
			   neighbor is opposite the kth vertex, and the first
			   neighbor is the horizon facet for the first vertex*/
  setT    *outsideset;  /* set of points outside this facet
		           if non-empty, last point is furthest */
  setT    *coplanarset; /* set of points coplanar with this facet
  			   > min_vertex and <= facet->max_outside
                           a point is assigned to the furthest facet
		           if non-empty, last point is furthest away */
  unsigned visitid;     /* visit_id, for visiting all neighbors,
			   all uses are independent */
  unsigned id;	        /* unique identifier from qh facet_id */
  unsigned nummerge:9;  /* number of merges */
#define qh_MAXnummerge 511 /*     2^9-1 */
  flagT	   newfacet:1;  /* True if facet on qh newfacet_list (new or merged) */
  flagT	   visible:1;   /* True if visible facet (will be deleted) */
  flagT    toporient:1; /* True if created with top orientation
			   after merging, use ridge orientation */
  flagT    simplicial:1;/* True if simplicial facet, ->ridges may be implicit */
  flagT    seen:1;      /* used to perform operations only once, like visitid */
  flagT	   flipped:1;   /* True if facet is flipped */
  flagT    upperdelaunay:1; /* True if facet is upper envelope of Delaunay triangulation */
  		/* flags primarily for output */
  flagT	   good:1;      /* True if a facet marked good for output */
  flagT    isarea:1;    /* True if facet->f.area is defined */
		/* flags for merging */
  flagT    dupridge:1;  /* True if duplicate ridge in facet */
  flagT    mergeridge:1; /* True if facet or neighbor contains a qh_MERGEridge
                            ->normal defined (also defined for mergeridge2) */
  flagT    mergeridge2:1; /* True if neighbor contains a qh_MERGEridge (mark_dupridges */
  flagT    coplanar:1;  /* True if horizon facet is coplanar at last use */
  flagT     mergehorizon:1; /* True if will merge into horizon (->coplanar) */
  flagT	    cycledone:1;/* True if mergecycle_all already done */
  flagT    tested:1;    /* True if facet convexity has been tested (false after merge */
  flagT    keepcentrum:1; /* True if keep old centrum after a merge */
  flagT	   newmerge:1;  /* True if facet is newly merged for reducevertices */
  flagT	   degenerate:1; /* True if facet is degenerate (degen_mergeset) */
  flagT	   redundant:1;  /* True if facet is redundant (degen_mergeset) */
};


/*----------------------------------------------
-ridgeT- specifies a ridge

  a ridge is hull_dim-1 simplex between two neighboring facets.  If the
  facets are non-simplicial, there may be more than one ridge between
  two facets.  E.G. a 4-d hypercube has two triangles between each pair
  of neighboring facets.
*/

struct ridgeT {
  setT    *vertices;    /* vertices belonging to this ridge, inverse sorted by id 
                           NULL if a degen ridge (matchsame) */
  facetT  *top;         /* top facet this ridge is part of */
  facetT  *bottom;      /* bottom facet this ridge is part of */
  unsigned id:24;       /* unique identifier, =>room for 8 flags */
  flagT    seen:1;      /* used to perform operations only once */
  flagT    tested:1;    /* True when ridge is tested for convexity */
  flagT    nonconvex:1; /* True if getmergeset detected a non-convex neighbor 
			   only one ridge between neighbors may have nonconvex */
};

/* ----------------------------------------------
-vertexT- specifies a vertex
*/

struct vertexT {
  vertexT *next;        /* next vertex in vertex_list */
  vertexT *previous;    /* previous vertex in vertex_list */
  pointT  *point;       /* hull_dim coordinates (coordT) */
  setT    *neighbors;   /* neighboring facets of vertex, qh_vertexneighbors()
			   inits in io.c or after first merge */
  unsigned visitid;     /* for use with qh vertex_visit */
  unsigned id:24;       /* unique identifier, =>room for 8 flags */
  flagT    seen:1;      /* used to perform operations only once */
  flagT    seen2:1;     /* another seen flag */
  flagT    delridge:1;  /* vertex was part of a deleted ridge */
  flagT	   deleted:1;   /* true if vertex on qh del_vertices */
  flagT    newlist:1;   /* true if vertex on qh newvertex_list */
};

/* ======= -global variables -qh ============================ 

   all global variables for qhull are in qh, qhmem, and qhstat
   
   qhmem is defined in mem.h and qhstat is defined in stat.h

   access to qh_qh is via the "qh" macro.  See qh_QHpointer in user.h
*/

typedef struct qhT qhT;    
#if qh_QHpointer
#define qh qh_qh->
extern qhT *qh_qh;     /* allocated in global.c */
#else
#define qh qh_qh.
extern qhT qh_qh;
#endif
extern char qh_version[];  /* defined in unix.c */

struct qhT {

  /*-user flags */

  boolT ALLpoints;        /* true 'Qi' if search all points for initial simplex */
  boolT ANGLEmerge;	  /* true 'Qa' if sort potential merges by angle */
  boolT APPROXhull;       /* true 'Wn' if MINoutside set */
  realT MINoutside;       /*   'Wn' min. distance for an outside point */
  boolT ATinfinity;       /* true 'd' if point num_points-1 is "at-infinity"
                             false 'Qu' if upper hull (no point at-infinity) */
  boolT AVOIDold;         /* true 'Q4' if avoid old->new merges */
  boolT BESToutside;      /* true 'Qf' if partition points into best outsideset */
  boolT CDDinput;         /* true 'Pc' if input uses CDD format (1.0/offset first) */
  boolT CDDoutput;        /* true 'PC' if print normals in CDD format (offset first) */
  boolT CHECKfrequently;  /* true 'Tc' if checking frequently */
  realT premerge_cos;     /*   'A-n'   cos_max when pre merging */
  realT postmerge_cos;    /*   'An'    cos_max when post merging */
  boolT DELAUNAY;         /* true 'd' if computing DELAUNAY triangulation */
  boolT DOintersections;  /* true 'Gh' if print hyperplane intersections */
  int   DROPdim;          /* drops dim 'GDn' for 4-d -> 3-d output */
  boolT FORCEoutput;      /* true 'Po' if forcing output despite degeneracies */
  int   GOODpoint;        /* 1+n for 'QGn', good facet if visible/not(-) from point n*/
  pointT *GOODpointp;     /*   the actual point */
  boolT GOODthreshold;    /* true if qh lower_threshold/upper_threshold defined 
  			     false if qh SPLITthreshold */
  int   GOODvertex;       /* 1+n, good facet if vertex for point n */
  pointT *GOODvertexp;     /*   the actual point */
  boolT HALFspace;        /* true 'Hn,n,n' if half-space intersection */
  int   IStracing;        /* trace execution, 0=none, 1=least, 4=most, -1=events */
  int   KEEParea;         /* 'PAn' number of largest facets to keep */
  boolT KEEPcoplanar;     /* true if keeping nearest facet for coplanar points */
  boolT KEEPinside;       /* true if keeping nearest facet for inside points */
  int   KEEPmerge;         /* 'PMn' number of facets to keep with most merges */
  realT KEEPminArea;      /* 'PFn' minimum facet area to keep */
  boolT MERGEexact;	  /* true 'Qx' if exact merges (coplanar, degen, dupridge, flipped) */
  boolT MERGEindependent; /* true 'Q2' if merging independent sets */
  boolT MERGING;          /* true if exact-, pre- or post-merging, with angle and centrum tests */
  realT   premerge_centrum;  /*   'C-n' centrum_radius when pre merging */
  realT   postmerge_centrum; /*   'Cn' centrum_radius when post merging */
  boolT MERGEvertices;	  /* true 'Q3' if merging redundant vertices */
  realT MAXcoplanar;      /* 'Un' max distance below a facet to be coplanar*/
  realT MINvisible;       /* 'Vn' min. distance for a facet to be visible */
  boolT NOnearinside;     /* true 'Q8' if ignore near-inside points when partitioning */
  boolT ONLYgood; 	  /* true 'Qg' if process points with good visible or horizon facets */
  boolT ONLYmax; 	  /* true 'Qm' if only process points that increase max_outside */
  boolT POINTSmalloc;     /* true if qh first_point/num_points allocated */
  boolT POSTmerge;        /* true if merging after buildhull (Cn or An) */
  boolT PREmerge;         /* true if merging during buildhull (C-n or A-n) */
  			/* NOTE: some of these names are similar to qh_PRINT names */
  boolT PRINTcentrums;	  /* true 'Gc' if printing centrums */
  boolT PRINTcoplanar;    /* true 'Gp' if printing coplanar points */
  int	PRINTdim;      	  /* print dimension for Geomview output */
  boolT PRINTdots;        /* true 'Ga' if printing all points as dots */
  boolT PRINTgood;        /* true 'Pg' if printing good facets */
  boolT PRINTinner;	  /* true 'Gi' if printing inner planes */
  boolT PRINTneighbors;	  /* true 'PG' if printing neighbors of good facets */
  boolT PRINTnoplanes;	  /* true 'Gn' if printing no planes */
  boolT PRINToptions1st;  /* true 'FO' if printing options to stderr */
  boolT PRINTouter;	  /* true 'Go' if printing outer planes */
  boolT PRINTprecision;   /* false 'Pp' if not reporting precision problems */
  qh_PRINT PRINTout[qh_PRINTEND]; /* list of output formats to print */
  boolT PRINTridges;      /* true 'Gr' if print ridges */
  boolT PRINTspheres;     /* true 'Gv' if print vertices as spheres */
  boolT PRINTstatistics;  /* true 'Ts' if printing statistics to stderr */
  boolT PRINTsummary;     /* true 's' if printing summary to stderr */
  boolT PROJECTdelaunay;  /* true if DELAUNAY, no readpoints() and
			     need projectinput() for Delaunay */
  int   PROJECTinput;     /* number of projected dimensions 'bn:0Bn:0' */
  boolT QUICKhelp;	  /* true if quick help message for degen input */
  boolT RANDOMdist;       /* true if randomly change distplane and setfacetplane */
  realT RANDOMfactor;     /*    maximum perturbation */
  realT RANDOMa;         /*  qh_randomfactor is randr * RANDOMa + RANDOMb */
  realT RANDOMb;
  boolT RANDOMoutside;    /* true if select a random outside point */
  int	REPORTfreq;       /* buildtracing reports every n facets */
  int   REPORTfreq2;	  /* tracemerging reports every REPORTfreq/2 facets */
  int	ROTATErandom;	  /* 'QRn' seed, 0 time, >= rotate input */
  boolT SCALEinput;       /* true if scaling input, 'Qb' */
  boolT SETroundoff;      /* true 'E' if qh DISTround is predefined */
  boolT SKIPcheckmax;	  /* true 'Q5' if skip qh_check_maxout */
  boolT SKIPconvex;       /* true 'Q6' if skip convexity testing during pre-merge */
  boolT SPLITthresholds;  /* true if upper_/lower_threshold defines a region
                               used only for printing (not for qh ONLYgood) */
  int	STOPcone;         /* 'TCn' 1+n for stopping after cone for point n*/
  int	STOPpoint;        /* 'TVn' 'TV-n' 1+n for stopping after/before(-) 
			                adding point n */
  boolT TESTvneighbors;   /*  true 'Qv' if test vertex neighbors at end */
  int   TRACElevel;       /* 'Tn' conditional IStracing level */
  int   TRACEpoint;       /* 'TPn' start tracing when point n is a vertex */
  realT TRACEdist;        /* 'TWn' start tracing when merge distance too big */
  int   TRACEmerge;       /* 'TMn' start tracing before this merge */
  boolT VERIFYoutput;     /* true 'Tv' if verify output at end of qhull */
  boolT VIRTUALmemory;    /* true 'Q7' if depth-first processing in buildhull */
  boolT VORONOI;	  /* true 'v' if computing Voronoi diagram */

  /* -input constants */
  realT AREAfactor;       /* 1/(hull_dim-1)! for converting det's to area */
  boolT DOcheckmax;       /* true if calling qh_check_maxout (qh_initqhull_globals) */
  char	*feasible_string;  /* feasible point 'Hn,n,n' for half-space intersection */
  coordT *feasible_point;  /*    as coordinates, both malloc'd */
  boolT GETarea;          /* true if need to compute facet areas in io.c */
  boolT KEEPnearinside;   /* true if near-inside points in coplanarset */
  int 	input_dim;	  /* dimension of input, set by initbuffers */
  int 	num_points;       /* number of input points */
  pointT *first_point;    /* array of input points, malloc'd */
  int 	hull_dim;         /* dimension of hull, set by initbuffers */
  char 	qhull_command[256];/* command line that invoked this program */
  char 	rbox_command[256]; /* command line that produced the input points */
  char  qhull_options[512]; /* descriptive list of options */
  int   qhull_optionslen; /*    length of line */
  boolT VERTEXneighbors;  /* true if maintaining vertex neighbors */
  boolT ZEROcentrum;      /* true if 'C-0' or 'C-0 Qx'.  sets ZEROall_ok */
  realT *upper_threshold; /* don't print if facet->normal[k]>=upper_threshold[k]
                             must set either GOODthreshold or SPLITthreshold
  			     if Delaunay, default is 0.0 for upper envelope */
  realT *lower_threshold; /* don't print if facet->normal[k] <=lower_threshold[k] */
  realT *upper_bound;     /* scale point[k] to new upper bound */
  realT *lower_bound;     /* scale point[k] to new lower bound 
  			     project if both upper_ and lower_bound == 0 */

  /* -precision constants, computed in qh_maxmin */
  
  realT ANGLEround;       /* max round off error for angles */
  realT centrum_radius;   /* max centrum radius for convexity (roundoff added) */
  realT cos_max;	  /* max cosine for convexity (roundoff added) */
  realT DISTround;        /* max round off error for distances, 'E' overrides */
  realT maxmaxcoord;      /* max coordinate in any dimension */
  realT MINdenom_1;       /* min. abs. value for 1/x */
  realT MINdenom;         /*    use divzero if denominator < MINdenom */
  realT MINdenom_1_2;     /* min. abs. val for 1/x that allows normalization */
  realT MINdenom_2;       /*    use divzero if denominator < MINdenom_2 */
  realT MINnorm;          /* min. norm for not redoing qh_sethyperplane_det */
  realT *NEARzero;        /* hull_dim array for near zero in gausselim */
  realT NEARinside;       /* keep points for qh_check_maxout if close to facet */
  realT ONEmerge;         /* max distance for merging simplicial facets */
  realT WIDEfacet;        /* size of wide facet for skipping ridge in
			     area computation and locking centrum */
  
  /* -internal constants */

  char qhull[sizeof("qhull")]; /* for checking ownership */
  void *old_stat;         /* pointer to saved qh_qhstat, qh_save_qhull */
  jmp_buf errexit;        /* exit label for qh_errexit, defined by setjmp() */
  char jmpXtra[40];       /* extra bytes in case jmp_buf is defined wrong */
  FILE *fin;              /* pointer to input file, init by qh_meminit */
  FILE *fout;             /* pointer to output file */
  FILE *ferr;             /* pointer to error file */
  pointT *interior_point; /* center point of the initial simplex*/
  int   normal_size;      /* size in bytes for facet normals and point coords*/
  int   center_size;      /* size in bytes for Voronoi centers */
  int   TEMPsize;         /* size for small, temporary sets (in quick mem) */

  /* -list of all facets, from facet_list to facet_tail, see qh_appendfacet */
 
  facetT *facet_list;     /* first facet */
  facetT  *facet_tail;     /* end of facet_list (dummy facet) */
  facetT *facet_next;     /* next facet for buildhull()
    			     all previous facets do not have outside sets*/
  facetT *newfacet_list;  /* list of new facets to end of facet_list */
  facetT *visible_list;   /* list of visible facets preceeding newfacet_list,
                             facet->visible set */
  int       num_visible;  /* current number of visible facets */
  unsigned tracefacet_id;  /* set at init, then can print whenever */
  facetT *tracefacet;     /*   set in newfacet/mergefacet, undone in delfacet*/
  unsigned tracevertex_id;  /* set at buildtracing, can print whenever */
  vertexT *tracevertex;     /*   set in newvertex, undone in delvertex*/
  vertexT *vertex_list;   /* list of all vertices, to vertex_tail */
  vertexT  *vertex_tail;   
  vertexT *newvertex_list; /* list of vertices in newfacet_list, to vertex_tail
                             all vertices have 'newlist' set */
  int 	num_facets;	  /* number of facets in facet_list
			     includes visble faces (num_visible) */
  int 	num_vertices;     /* number of vertices in facet_list */
  int   num_outside;      /* number of points in outsidesets (for tracing) */
  int   num_good;         /* number of good facets (after findgood_all) */
  int 	facet_id;         /* id of next, new facet from newfacet() */
  int 	ridge_id;         /* id of next, new ridge from newridge() */
  unsigned vertex_id;        /* id of next, new vertex from newvertex() */

  /* -variables */
  
  unsigned hulltime;       /* ignore time to set up input and randomize */
  qh_CENTER CENTERtype;       /* current type of facet->center, qh_CENTER */
  int 	furthest_id;      /* pointid of furthest point, for tracing */
  facetT *GOODclosest;    /* closest facet to GOODthreshold in qh_findgood */
  realT max_outside;      /* maximum distance from a point to a facet,
			       before roundoff, not simplicial vertices
			       actual outer plane is +DISTround and
			       computed outer plane is +2*DISTround */
  realT max_vertex;       /* maximum distance (>0) from vertex to a facet,
			       before roundoff, not simplicial vertices */
  realT min_vertex;       /* minimum distance (<0) from vertex to a facet,
			       before roundoff, not simplicial vertices 
			       actual inner plane is -DISTround and
			       computed inner plane is -2*DISTround */
  boolT NEWfacets;        /* true while visible facets invalid due to new or merge
			      from makecone/attachnewfacets to deletevisible */
  boolT findbestnew;	  /* true if partitioning calls qh_findbestnew */
  boolT findbest_notsharp; /* true if new facets are at least 90 degrees */
  boolT NOerrexit;        /* true if qh_errexit is not available */
  realT PRINTcradius;     /* radius for printing centrums */
  realT PRINTradius;      /* radius for printing vertex spheres and points */
  boolT POSTmerging;      /* true when post merging */
  int 	printoutvar;	  /* temporary variable for qh_printbegin, etc. */
  int 	printoutnum;	  /* number of facets printed */
  boolT QHULLfinished;    /* True after qhull() is finished */
  realT totarea;          /* total facet area computed by qh_getarea */
  realT totvol;           /* total volume computed by qh_getarea */
  int 	visit_id;         /* unique id for searching neighborhoods, */
  int 	vertex_visit;     /* unique id for searching vertices */
  boolT ZEROall_ok;       /* True if qh_checkzero always succeeds */
  
  /* -sets */
  setT *facet_mergeset;   /* temporary set of merges to be done */
  setT *degen_mergeset;   /* temporary set of degenerate and redundant merges */
  setT *initial_points;   /* initial simplex for buildhull() */
  setT *hash_table;	  /* hash table for matching ridges in qh_matchfacets 
                             size is setsize() */
  int   num_hashentries;  /* current number of hashentries */
  setT *other_points;     /* additional points (first is qh interior_point) */
  setT *del_vertices;     /* vertices to partition and delete with visible 
                             facets.  Have deleted set for checkfacet */

  /* -buffers */
  coordT *gm_matrix;      /* (dim+1)Xdim matrix for geom.c */
  coordT **gm_row;        /* array of gm_matrix rows */
  char* line;             /* malloc'd input line of maxline+1 chars */
  int maxline;
  coordT *half_space;     /* malloc'd input array for halfspace (qh normal_size+coordT) */
  coordT *temp_malloc;    /* malloc'd input array for points */
  
  /* -statics */
  boolT ERREXITcalled;    /* true during errexit (prevents duplicate calls */
  boolT firstcentrum; 	  /* for qh_printcentrum */
  int  lastreport;        /* for qh_buildtracing */
  int  mergereport;       /* for qh_tracemerging */
  boolT old_randomdist;   /* save in io.c for RANDOMdist */
  int   ridgeoutnum;      /* number of ridges in 4OFF output */
  void *old_qhstat;       /* for saving qh_qhstat in save_qhull() */
  setT *old_tempstack;     /* for saving qhmem.tempstack in save_qhull */
  setT *searchset;        /* set of facets for searching in qh_findbest() */
  int   rand_seed;        /* for qh_rand/qh_srand */
};

/* =========== -macros- ========================= 
-otherfacet_(ridge, facet)   return neighboring facet for a ridge in facet
-getid_(p)		     return id or -1 if NULL
*/

#define otherfacet_(ridge, facet) \
                        (((ridge)->top == (facet)) ? (ridge)->bottom : (ridge)->top)
#define getid_(p)       ((p) ? (p)->id : -1)

/* ---------------------------------------------
-FORALL and FOREACH macros

   These all iterate using a variable of the same name, e.g. FORALLfacets
   and FOREACHfacet_ uses 'facet' declared by 'facetT *facet'.  The macros
   may use auxiliary variables as indicated.

-FORALLfacets                iterate over all facets in facetlist 
-FORALLpoint_(points, num)   iterate over num points (uses 'pointT *pointtemp')
-FORALLvertices              iterate over all vertices in vertex_list

-FOREACHfacet_(facets)	     iterate over facet set (uses 'facetT **facetp')
-FOREACHneighbor_(facet)     iterate over facet->neighbors (uses 'facetT **neighborp')
-FOREACHpoint_(points)       iterate over point set (uses 'pointT **pointp')
-FOREACHridge_(ridges)	     iterate over ridge set (uses 'ridgeT **ridgep')
-FOREACHvertex_(vertice)     iterate over vertex set (uses 'vertexT **vertexp')
-FOREACHadjacent_(vertex)    iterate over adjacent vertices to vertex 
-FOREACHneighbor_(vertex)    iterate over neighboring facets to vertex 

-FOREACHfacet_i_(facets)    iterate over facets by facet_i and facet_n
-FOREACHneighbor_i_(facet)  iterate over facet->neighbors by neighbor_i, neighbor_n
-FOREACHvertex_i_(vertices) iterate over vertices by vertex_i, vertex_n
-FOREACHpoint_i_(points)    iterate over points by point_i, point_n
-FOREACHridge_i_(ridges)    iterate over ridges by ridge_i, ridge_n
-FOREACHneighbor_i_(vertex) iterate over vertex->neighbors by neighbor_i, neighbor_n

 WARNING: nested loops can't use the same variable (define another FOREACH)
 WARNING: strange behavior if don't fully brace when nested (including
        intervening blocks, e.g. FOREACH...{ if () FOREACH...} )
 poly.h defines other FOREACH/FORALL macros
 set.h  defines FOREACHsetelement and contains additional notes
*/
#define FORALLfacets for (facet=qh facet_list;facet && facet->next;facet=facet->next)
#define FORALLpoints FORALLpoint_(qh first_point, qh num_points)
#define FORALLvertices for (vertex=qh vertex_list;vertex && vertex->next;vertex= vertex->next)

#define FORALLpoint_(points, num) for(point= (points), \
      pointtemp= (points)+qh hull_dim*(num); point < pointtemp; point += qh hull_dim)
#define FOREACHfacet_(facets)    FOREACHsetelement_(facetT, facets, facet)
#define FOREACHneighbor_(facet)  FOREACHsetelement_(facetT, facet->neighbors, neighbor)
#define FOREACHpoint_(points)    FOREACHsetelement_(pointT, points, point)
#define FOREACHridge_(ridges)    FOREACHsetelement_(ridgeT, ridges, ridge)
#define FOREACHvertex_(vertices) FOREACHsetelement_(vertexT, vertices,vertex)

#define FOREACHfacet_i_(facets)    FOREACHsetelement_i_(facetT, facets, facet)
#define FOREACHneighbor_i_(facet)  FOREACHsetelement_i_(facetT, facet->neighbors, neighbor)
#define FOREACHpoint_i_(points)    FOREACHsetelement_i_(pointT, points, point)
#define FOREACHridge_i_(ridges)    FOREACHsetelement_i_(ridgeT, ridges, ridge)
#define FOREACHvertex_i_(vertices) FOREACHsetelement_i_(vertexT, vertices,vertex)

/* ======= -functions =========== 

  	see corresponding .c file for definitions

	Qhull functions (see qhull.c and qhull_a.h)
-qhull		construct the convex hull of a set of points
-addpoint       add point to hull (must be above facet)
-printsummary	print summary about the output

	User redefinable functions (see user.c)
-errexit	 	return exitcode to system after an error
-errprint		print erroneous facets, ridge, and vertex
-printfacetlist		print all fields for a list of facets
-user_memsizes          define up to 10 additional quick allocation sizes
  	
	Geometric functions (see geom.c and geom.h for other useful functions)
-gram_schmidt   implements Gram-Schmidt orthogonalization by rows
-projectinput  project input along one or more dimensions + Delaunay projection
-randommatrix   generate a random dimXdim matrix in range (-1,1)
-rotatepoints   rotate numpoints points by a row matrix
-scaleinput    scale input to new lowbound and highbound
-sethalfspace_all generate dual for halfspace intersection with feasible point


	Global init/free functions (see global.c and qhull_a.h)
-freeqhull	     free memory used by qhull
-init_A              called before error handling initialized
-init_B              called after points are defined
-initflags	     set flags and initialized constants from command line
-restore_qhull       restores a saved qhull
-save_qhull          saves qhull for later restoring

	Input/output functions (see io.c and io.h)
-dfacet		print facet by id
-dvertex	print vertex by id
-printsummary	print summary about the output
-produce_output prints out the result of qhull in desired format
-readpoints     read points from input

	Polyhedron functions (see poly.c or poly2.c)
-check_output	check output data structure according to user flags
-check_points	verify that all points are inside the hull
-setvoronoi_all compute Voronoi centers for all facets
-findfacet      find facet that is furthest below a point 
-nearvertex     return nearest vertex to point

-point          return point for a point id, or NULL if unknown
-pointid        return id for a point, or -1 if not known

-facetvertices  returns temporary set of vertices in a set of facets
-pointfacet	return temporary set of facets indexed by point id
-pointvertex	return temporary set of vertices indexed by point id

        Statistics functions (see stat.c)
-printallstatistics print all statistics

	other functions no longer needed by most users
-init_qhull_command  build qhull_command from argc/argv
-initqhull_buffers   initialize global memory buffers
-initqhull_globals   initialize globals
-initqhull_mem	     initialize mem.c for qhull
-initqhull_start     start initialization of qhull
-initthresholds	     set thresholds for printing and scaling from command line
-findbest            find visible facet for a point starting at a facet
-findbestnew         find best newfacet for point

*/

/********* -qhull.c prototypes (duplicated from qhull_a.h) **********************/

void    qh_qhull (void);
boolT   qh_addpoint (pointT *furthest, facetT *facet, boolT checkdist);
void	qh_printsummary(FILE *fp);

/********* -user.c prototypes (alphabetical) **********************/

void 	qh_errexit(int exitcode, facetT *facet, ridgeT *ridge);
void 	qh_errprint(char* string, facetT *atfacet, facetT *otherfacet, ridgeT *atridge, vertexT *atvertex);
void    qh_printfacetlist(facetT *facetlist, setT *facets, boolT printall);
void 	qh_user_memsizes (void);

/***** -geom.c/geom2.c prototypes (duplicated from geom.h) ****************/

facetT *qh_findbest (pointT *point, facetT *facet, boolT bestoutside,
		     boolT newfacets, realT *dist, boolT *isoutside, int *numpart);
facetT *qh_findbestnew (pointT *point, facetT *startfacet,
	   realT *dist, boolT *isoutside, int *numpart);
boolT   qh_gram_schmidt(int dim, realT **rows);
void	qh_printsummary(FILE *fp);
void    qh_projectinput (void);
void    qh_randommatrix (realT *buffer, int dim, realT **row);
void    qh_rotateinput (realT **rows);
void    qh_scaleinput (void);
coordT  *qh_sethalfspace_all (int dim, int count, coordT *halfspaces, pointT *feasible);

/***** -global.c prototypes (alphabetical) ***********************/

void 	qh_freebuffers (void);
void    qh_freeqhull (boolT allmem);
void    qh_init_A (FILE *infile, FILE *outfile, FILE *errfile, int argc, char *argv[]);
void    qh_init_B (coordT *points, int numpoints, int dim, boolT ismalloc);
void 	qh_init_qhull_command (int argc, char *argv[]);
void    qh_initbuffers (coordT *points, int numpoints, int dim, boolT ismalloc);
void 	qh_initflags (char *command);
void 	qh_initqhull_buffers (void);
void 	qh_initqhull_globals (coordT *points, int numpoints, int dim, boolT ismalloc);
void    qh_initqhull_mem (void);
void 	qh_initqhull_start (FILE *infile, FILE *outfile, FILE *errfile);
void 	qh_initthresholds (char *command);
#if qh_QHpointer
void 	qh_restore_qhull (qhT **oldqh);
qhT    *qh_save_qhull (void);
#endif

/***** -io.c prototypes (duplicated from io.h) ***********************/

void    dfacet( int id);
void    dvertex( int id);
void	qh_produce_output(void);
coordT *qh_readpoints(int *numpoints, int *dimension, boolT *ismalloc);


/********* -mem.c prototypes (duplicated from mem.h) **********************/

void qh_meminit (FILE *ferr);
void qh_memfreeshort (int *curlong, int *totlong);

/********* -poly.c/poly2.c prototypes (duplicated from poly.h) **********************/

void    qh_check_output (void);
void    qh_check_points (void);
setT   *qh_facetvertices (facetT *facetlist, setT *facets, boolT allfacets);
facetT *qh_findfacet (pointT *point, facetT *facet, 
           realT *dist, boolT *isoutside, int *numpart);
vertexT *qh_nearvertex (facetT *facet, pointT *point, realT *bestdistp);
pointT *qh_point (int id);
setT   *qh_pointfacet (void /*qh.facet_list*/);
int     qh_pointid (pointT *point);
setT   *qh_pointvertex (void /*qh.facet_list*/);
void    qh_setvoronoi_all (void);

/********* -stat.c prototypes (duplicated from stat.h) **********************/

void    qh_printallstatistics (FILE *fp, char *string);

#endif /* qhDEFqhull */



