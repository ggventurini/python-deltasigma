/* user.h - user redefinable constants

   see README.  see COPYING for copyright information.

   before reading any code, review qhull.h for data structure definitions and 
   the "qh" macro.
   

*/

#ifndef qhDEFuser
#define qhDEFuser 1

/*---------------------------------
-realT -- set the size of floating point numbers

   Select whether to store floating point numbers in single precision (float)
   or double precision (double).
   
   Use 'float' to save about 8% in time and 25% in space.  This is particularly
   help if high-d where convex hulls are space limited.  Using 'float' also
   reduces the printed size of Qhull's output since numbers have 8 digits of 
   precision.
   
   Use 'double' when greater arithmetic precision is needed.  This is needed
   for Delaunay triangulations and Voronoi diagrams when you are not merging 
   facets.

   If 'double' gives insufficient precision, your data probably includes
   degeneracies.  If so you should use facet merging (options 'Qx' or 'C-0 Qc')
   or exact arithmetic (see imprecision section of manual, qh-impre.html).  
   You may also use option 'Po' to force output despite precision errors.

   You may use 'long double', but many format statements need to be changed
   and you may need a 'long double' square root routine.  S. Grundmann
   (sg@eeiwzb.et.tu-dresden.de) has done this.  He reports that the code runs 
   much slower with little gain in precision.    

   WARNING: on some machines,    int f(){realT a= REALmax;return (a == REALmax);}
      returns False.  Use (a > REALmax/2) instead of (a == REALmax).

   REALfloat =   1      all numbers are 'float' type
             =   0      all numbers are 'double' type
*/

#define REALfloat 0

#if (REALfloat == 1)
#define realT float
#define REALmax FLT_MAX
#define REALmin FLT_MIN
#define REALepsilon FLT_EPSILON
#define qh_REALdigits 8   /* maximum number of significant digits */
#define qh_REAL_1 "%6.8g "
#define qh_REAL_2n "%6.8g %6.8g\n"
#define qh_REAL_3n "%6.8g %6.8g %6.8g\n"

#elif (REALfloat == 0)
#define realT double
#define REALmax DBL_MAX
#define REALmin DBL_MIN
#define REALepsilon DBL_EPSILON
#define qh_REALdigits 16    /* maximum number of significant digits */
#define qh_REAL_1 "%6.16g "
#define qh_REAL_2n "%6.16g %6.16g\n"
#define qh_REAL_3n "%6.16g %6.16g %6.16g\n"

#else
#error unknown float option
#endif

/*-----------------------------------
-qh_CPUclock- define the clock() function

   if your system does not use clock() to return CPU ticks, replace
   qh_CPUclock with the corresponding function.  It needs to return
   a 'long int'.
   
   qh_SECticks is the number of clock ticks per second
   
   qh_CPUclock is only used for reporting the total time spent by Qhull

   Set qh_CLOCKtype to
   
     1	   	for CLOCKS_PER_SEC, CLOCKS_PER_SECOND, or microsecond
*/
#define qh_CLOCKtype 1  /* change to the desired number */

#if (qh_CLOCKtype == 1)

#if defined (CLOCKS_PER_SECOND)
#define qh_CPUclock    ((unsigned)clock())  /* return CPU clock */
#define qh_SECticks CLOCKS_PER_SECOND

#elif defined (CLOCKS_PER_SEC)
#define qh_CPUclock    ((unsigned)clock())  /* return CPU clock */
#define qh_SECticks CLOCKS_PER_SEC

#else
#define qh_CPUclock    ((unsigned)clock())  /* return CPU clock */
#define qh_SECticks 1E6
#endif

#else /* qh_CLOCKtype == ? */
#error unknown clock option
#endif

/* qh_RANDOM...

   Random number generators differ between systems.  Most systems provide
   rand() but the period varies.  The period of rand() is not critical
   since qhull does not normally use random numbers.  

   The default generator is Park & Miller's minimal standard random
   number generator [CACM 31:1195 '88].  It is included with Qhull.

   Random numbers are used by rbox to generate point sets.  Random
   numbers are used by qhull to rotate the input ('QRn' option),
   simulate a randomized algorithm ('Qr' option), and to simulate
   roundoff errors ('Rn' option).
   
   If qh_RANDOMmax is wrong, qhull will report a warning and Geomview output
   will likely be invisible.
   
   qh_RANDOMint generates a random integer between 0 and qh_RANDOMmax.  
   
   qh_RANDOMseed sets the random number seed for qh_RANDOMint

   Set qh_RANDOMtype to

     1       for random() with 31 bits (UCB)
     2       for rand() with RAND_MAX or 15 bits (system 5)
     3       for rand() with 31 bits (Sun)
     4       for lrand48() with 31 bits (Solaris)
     5       for qh_rand() with 31 bits (included with Qhull)
*/
#define qh_RANDOMtype 4   /* *** change to the desired number *** */

#if (qh_RANDOMtype == 1)
#define qh_RANDOMmax ((realT)0x7fffffffUL)  /* 31 bits, random()/MAX */
#define qh_RANDOMint random()
#define qh_RANDOMseed_(seed) srandom(seed);

#elif (qh_RANDOMtype == 2)
#ifdef RAND_MAX
#define qh_RANDOMmax ((realT)RAND_MAX)
#else
#define qh_RANDOMmax ((realT)32767)   /* 15 bits (System 5) */
#endif
#define qh_RANDOMint  rand()
#define qh_RANDOMseed_(seed) srand((unsigned)seed);
  
#elif (qh_RANDOMtype == 3)
#define qh_RANDOMmax ((realT)0x7fffffffUL)  /* 31 bits, Sun */
#define qh_RANDOMint  rand()
#define qh_RANDOMseed_(seed) srand((unsigned)seed);

#elif (qh_RANDOMtype == 4)
#define qh_RANDOMmax ((realT)0x7fffffffUL)  /* 31 bits, lrand38()/MAX */
#define qh_RANDOMint lrand48()
#define qh_RANDOMseed_(seed) srand48(seed);

#elif (qh_RANDOMtype == 5)
#define qh_RANDOMmax ((realT)2147483646UL)  /* 31 bits, qh_rand/MAX */
#define qh_RANDOMint qh_rand()
#define qh_RANDOMseed_(seed) qh_srand(seed);
/* unlike rand(), never returns 0 */

#else
#error: unknown random option
#endif

/*------------------------------------------
        constants
*/
#define qh_INFINITE  -10.101 /* on output, indicates Voronoi center at infinity */
#define qh_DEFAULTbox 0.5    /* default box size (Geomview expects 0.5) */
#define qh_ORIENTclock 0     /* 0 for inward pointing normals in Geomview */


/*-----------------------
	performance related constants
	
-HASHfactor     total hash slots / used hash slots
-VERIFYdirect   verify all points against all facets if op count smaller
-INITIALmax     if dim >=, use min/max coordinate points for initial simplex
                  use option 'Qs' to override (much slower)
-INITIALsearch  if INITIALmax, search points up to this dimension
*/

#define qh_HASHfactor 2             /* (int) at worst 50% occupancy for qh hash_table
                                       and normally 25% occupancy */
#define qh_VERIFYdirect 1000000     /* if more tests, use qh_findbest instead */
#define qh_INITIALmax 8             /* if hull_dim >=, build initial simplex
				       from points with non-zero determinants */
#define qh_INITIALsearch 6          /* and search points up to this dimension*/

/*-----------------------
       	Memory constants for calling qh_meminitbuffers in global.c
-MEMalign	memory alignment (see mem.h). If using gcc, best alignment is
              #define qh_MEMalign fmax_(__alignof__(realT),__alignof__(void *))

-MEMbufsize 	memory buffer size
-MEMinitbuf 	initial memory buffer size.  It should hold enough
   		   facets to keep outsidesets in short memory.

*/

#define qh_MEMalign fmax_(sizeof(realT), sizeof(void *))
#define qh_MEMbufsize 0x10000       /* allocate 64K memory buffers */
#define qh_MEMinitbuf 0x20000       /* initially allocate 128K buffer */


/* ======= -global variables -qh ============================ 

   all global variables for qhull are in qh, qhmem, and qhstat
   
   qhmem is defined in mem.h and qhstat is defined in stat.h

   access to qh_qh is via the "qh" macro.  There are two choices

   qh_QHpointer = 1     access globals via a pointer to allocated memory
                        enables qh_saveqhull() and qh_restoreqhull()
			costs about 8% in time and 2% in space

		= 0     qh_qh and qh_qhstat are static data structures
		        only one instance of qhull() can be active at a time
			default value
*/

#define qh_QHpointer 0  /* 1 for dynamic allocation, 0 for global structure */

/*--------------------------------------------------
	Conditional compilation

__MWERKS__      defined by Metrowerks when compiling for the Power Macintosh

-QUICKhelp        1 use abbreviated help messages for degenerate inputs

-COMPUTEfurthest computing furthest saves memory but costs time
                         (about 40% more distance tests for partitioning)
                         
-MAXoutside      keep maxoutside for each facet
		   this takes a realT per facet and slightly slows down qhull
		   it produces better outer planes for geomview output 

-KEEPstatistics   0 removes most of statistic gathering and reporting

To disallow statistice, define qh_KEEPstatistics as 0.  It reduces code size
by about 4%.
                       
To disallow merging, define qh_NOmerge.  This saves about 10% space.

    #define qh_NOmerge
    
To disallow tracing define qh_NOtrace.  This saves about 5% space.

    #define qh_NOtrace
    
see also: qh_NOmem in mem.c, and removing io.o via user_eg.c

*/

#define qh_KEEPstatistics 1   /* 0 to take out statistics */
#define qh_QUICKhelp    0    /* 1 for short help messages */
#define qh_COMPUTEfurthest 0    /* 1 removes facet->furthestdist */
#define qh_MAXoutside 1         /* 0 removes facet->maxoutside */

/* ============ -merge constants- ====================

   These constants effect facet merging.  You probably will not need
   to modify these.  They effect the performance of facet merging.


-BESTcentrum     if > 2*dim+n vertices, findbestneighbor tests centrums (faster)
                 else, findbestneighbor test all vertices (much better merges)

-BESTnonconvex   if > dim+n neighbors, findbestneighbor tests nonconvex ridges.
                 It is needed because findbestneighbor is slow for large facets
		   
-MAXnewmerges    if >n newmerges, merge_nonconvex calls reducevertices_centrums.
                 It isbneeded because postmerge can merge many facets at once

-MAXnewcentrum   if <= dim+n vertices (n approximates the number of merges),
                    reset the centrum in reducevertices_centrum 
                  needed to reduce cost and because centrums may move
		        too much if many vertices in high-d

-COPLANARratio   for 3-d+ merging, qh MINvisible is n*premerge_centrum
                 for 2-d merging, it's premerge_centrum 
		 for non-merging, it's DISTround

-DISToutside	 minimum distance when merging to stop search in 
                 qh_findbestnew or partitionall
                   if too big then O(n^2) behavior for partitioning in cone
		   if very small then important points not processed

-DIMmergeVertex  max dimension for vertex merging (not effective in high-d)

-DIMreduceBuild  max dimension for vertex reduction during build (slow in high-d)

-RATIOnearinside ratio of NEARinside to ONEmerge for retaining inside points for
                 qh_check_maxout.  This is overkill since do not know the correct value.

-USEfindbestnew  cut-off between qh_findbest/qh_findbestnew in #merged facets.

-WIDEcoplanar   if vertex is further than qh WIDEfacet from the hyperplane
                then its ridges are not counted in computing the area, and
		the facet's centrum is frozen. 
		qh WIDEfacet= max(qh MAXoutside,qh_WIDEcoplanar*qh MAXcoplanar,
		                  qh_WIDEcoplanar * qh MINvisible);
*/
#define qh_BESTcentrum 20   /* qh_findbestneighbor tests centrum instead of vertices */
#define qh_BESTcentrum2 2   /*   if BESTcentrum2 * hull_dim + BESTcentrum < #vertices */
#define qh_BESTnonconvex 15    /*findbestneighbor only tests nonconvex if > */
#define qh_DIMmergeVertex 6    /* maximum dimension for vertex merging */
#define qh_DIMreduceBuild 5      /* maximum dimension for reduce vertices during build */
#define qh_MAXnewmerges 2   /*merge_nonconvex calls reducevertices*/
#define qh_MAXnewcentrum 5    /* mergefacet resets centrum if <= */
#define qh_COPLANARratio 3    /* furthest is coplanar if < n*premerge_centrum */
#define qh_DISToutside fmax_(4*qh MINoutside, 2*qh max_outside)
#define qh_RATIOnearinside 5  /* ratio of distance for NEARinside to ONEmerge */
#define qh_USEfindbestnew 50 /* use qh_findbestnew if > n merged facets, or merged cone */
#define qh_WIDEcoplanar 6    /* n*MAXcoplanar or n*MINvisible for a WIDEfacet */

#endif /* qh_DEFuser */

