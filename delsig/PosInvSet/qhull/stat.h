/* stat.h - contains all statistics that are collected for qhull

   see README and stat.c

   copyright (c) 1993-1995, The Geometry Center

   recompile qhull if you change this file

   Integer statistics are Z* while real statistics are W*.  

   qh_KEEPstatistics=0 turns off statistic gathering (except zzdef/zzinc/zzadd)

   define maydebugx to call a routine at every statistic event

*/

#ifndef qhDEFstat
#define qhDEFstat 1
#ifndef qh_KEEPstatistics
#define qh_KEEPstatistics 1
#endif

/*----------------------------------------------------
-statistics, Zxxx for integers, Wxxx for reals

can pick up all statistics by:
    grep '[zw].*_[(][ZW]' *.c >z.x
    remove trailers with query-replace-regexp [,)].*
    remove leaders with  query-replace-regexp [ ^I]+  (
*/

#if qh_KEEPstatistics
enum statistics {     /* alphabetical after Z/W */
    Zacoplanar,
    Wacoplanarmax,
    Wacoplanartot,
    Zangle,
    Wangle,
    Wanglemax,
    Wanglemin,
    Zangletests,
    Wareatot,
    Wareamax,
    Wareamin,
    Zavoidold,
    Wavoidoldmax,
    Wavoidoldtot,
    Zback0,
    Zbestcentrum,
    Zbestdist,
    Zcentrumtests,
    Zcheckpart,
    Zcomputefurthest,
    Zconcave,
    Wconcavemax,
    Wconcavetot,
    Zconcaveridges,
    Zconcaveridge,
    Zcoplanar,
    Wcoplanarmax,
    Wcoplanartot,
    Zcoplanarangle,
    Zcoplanarcentrum,
    Zcoplanarhorizon,
    Zcoplanarinside,
    Zcoplanarpart,
    Zcoplanarridges,
    Wcpu,
    Zcyclefacetmax,
    Zcyclefacettot,
    Zcyclehorizon,
    Wcyclemax,
    Wcycletot,
    Zcyclevertex,
    Zdegen,
    Wdegenmax,
    Wdegentot,
    Zdegenvertex,
    Zdelfacetdup, 
    Zdelridge,
    Zdelvertextot,
    Zdelvertexmax,
    Zdetsimplex,
    Zdistcheck,
    Zdistconvex,
    Zdistgood,
    Zdistio,
    Zdistplane,
    Zdiststat,
    Zdistvertex,
    Zdistzero,
    Zdoc1,
    Zdoc2,
    Zdoc3,
    Zdoc4,
    Zdoc5,
    Zdoc6,
    Zdoc7,
    Zdoc8,
    Zdoc9,
    Zdoc10,
    Zdropdegen,
    Zdropneighbor,
    Zdupflip,
    Zduplicate,
    Wduplicatemax,
    Wduplicatetot,
    Zdupridge,
    Zdupsame,
    Zfindfail,
    Zflipped, 
    Wflippedmax, 
    Wflippedtot, 
    Zflippedfacets,
    Zgauss0,
    Zgoodfacet,
    Zhashlookup,
    Zhashridge,
    Zhashridgetest,
    Zhashtests,
    Zinsidevisible,
    Zintersect,
    Zintersectfail,
    Zintersectmax,
    Zintersectnum,
    Zintersecttot,
    Zmaxneighbors,
    Wmaxout,
    Wmaxoutside,
    Zmaxridges,
    Zmaxvertex,
    Zmaxvertices,
    Zmaxvneighbors,
    Zmemfacets,
    Zmempoints,
    Zmemridges,
    Zmemvertices,
    Zmergeflipdup,
    Zmergehorizon,
    Zmergeinittot,
    Zmergeinitmax,
    Zmergeinittot2,
    Zmergeintohorizon,
    Zmergenew,
    Zmergesettot,
    Zmergesetmax,
    Zmergesettot2,
    Zmergesimplex,
    Zmergevertex,
    Wmindenom,
    Wminvertex,
    Zminnorm,
    Zmultiridge,
    Znearlysingular,
    Zneighbor,
    Wnewbalance,
    Wnewbalance2,
    Znewfacettot,
    Znewfacetmax,
    Znewvertex,
    Wnewvertex,
    Wnewvertexmax,
    Znoarea,
    Znotgood,
    Znotgoodnew,
    Znotmax,
    Znumfacets,
    Znummergemax,
    Znummergetot,
    Znumneighbors,
    Znumridges,
    Znumvertices,
    Znumvisibility,
    Znumvneighbors,
    Zonehorizon,
    Zpartcoplanar,
    Zpartinside,
    Zpartition, 
    Zpartitionall,
    Zpartnear,
    Zpbalance,
    Wpbalance,
    Wpbalance2, 
    Zpostfacets, 
    Zpremergetot,
    Zprocessed,
    Zremvertex,
    Zremvertexdel,
    Zrenameall,
    Zrenamepinch,
    Zrenameshare,
    Zsamevertices,
    Zsearchpoints,
    Zsetplane,
    Ztestvneighbor,
    Ztotcheck,
    Ztothorizon,
    Ztotmerge,
    Ztotpartcoplanar,
    Ztotpartition,
    Ztotridges,
    Ztotvertices,
    Ztotvisible,
    Wvertexmax,
    Wvertexmin,
    Zvertexridge,
    Zvertexridgetot,
    Zvertexridgemax,
    Zvertices,
    Zvisfacettot,
    Zvisfacetmax,
    Zvisvertextot,
    Zvisvertexmax,
    Zwidefacet,
    Zwidevertices,
    ZEND};
#else
enum statistics {     /* for zzdef etc. macros */
  Zback0,
  Zbestdist,
  Zcentrumtests,
  Zconcaveridges,
  Zcoplanarhorizon,
  Zcoplanarpart,
  Zcoplanarridges,
  Zcyclefacettot,
  Zcyclehorizon,
  Zdistcheck,
  Zdistconvex,
  Zdistzero,
  Zdoc1,
  Zdoc2,
  Zdoc3,
  Zflippedfacets,
  Zgauss0,
  Zminnorm,
  Zmultiridge,
  Znearlysingular,
  Znumvisibility,
  Zpartcoplanar,
  Zpartition,
  Zpartitionall,
  Zprocessed,
  Zsetplane,
  Ztotmerge,
    ZEND};
#endif

/* ------------ -ztypes- ---------------------
the type of a statistic sets its initial value.  The type should
be the same as the macro for collecting the statistic
*/
enum ztypes {zdoc,zinc,zadd,zmax,zmin,ZTYPEreal,wadd,wmax,wmin,ZTYPEend};

/*------------ -macros -------------
macros:
  zdef_(type, name, doc, -1)	define a statistic (assumes 'qhstat next= 0;')
  zdef_(type, name, doc, count)	   printed as name/count
  zinc_(name)                   integer statistic is count
  zadd/wadd_(name, value)       integer or real statistic is total value
  zmax/wmax_(name, value)	integer or real statistic is max value
  zmin/wmin_(name, value)	integer or real statistic is min value
  zval/wval_(name)		set or return value of statistic
*/

#define MAYdebugx  /* maydebug() is called frequently to trap an error */
#define zzinc_(id) {MAYdebugx; qhstat stats[id].i++;}
#define zzadd_(id, val) {MAYdebugx; qhstat stats[id].i += (val);}
#define zzval_(id) ((qhstat stats[id]).i)
#define wwval_(id) ((qhstat stats[id]).r)
#define zzdef_(stype,name,string,cnt) qhstat id[qhstat next++]=name; \
   qhstat doc[name]= string; qhstat count[name]= cnt; qhstat type[name]= stype

#if qh_KEEPstatistics
#define zinc_(id) {MAYdebugx; qhstat stats[id].i++;}
#define zadd_(id, val) {MAYdebugx; qhstat stats[id].i += (val);}
#define wadd_(id, val) {MAYdebugx; qhstat stats[id].r += (val);}
#define zmax_(id, val) {MAYdebugx; maximize_(qhstat stats[id].i,(val));}
#define wmax_(id, val) {MAYdebugx; maximize_(qhstat stats[id].r,(val));}
#define zmin_(id, val) {MAYdebugx; minimize_(qhstat stats[id].i,(val));}
#define wmin_(id, val) {MAYdebugx; minimize_(qhstat stats[id].r,(val));}
#define zval_(id) ((qhstat stats[id]).i)
#define wval_(id) ((qhstat stats[id]).r)

#define zdef_(stype,name,string,cnt) qhstat id[qhstat next++]=name; \
   qhstat doc[name]= string; qhstat count[name]= cnt; qhstat type[name]= stype


#else  /* !qh_KEEPstatistics */
#define zinc_(id) {}
#define zadd_(id, val) {}
#define wadd_(id, val) {}
#define zmax_(id, val) {}
#define wmax_(id, val) {}
#define zmin_(id, val) {}
#define wmin_(id, val) {}
#define zval_(id) qhstat tempi
#define wval_(id) qhstat tempr
#define zdef_(type,name,doc,count)
#define ZMAXlevel 1
#endif
  
/* -typedef and extern-  types are defined below */

typedef struct qhstatT qhstatT;     /* global data structure for statistics */
typedef union intrealT intrealT;    /* union of int and realT */

/* -qhstat-

   access to qh_qhstat is via the "qhstat" macro.  There are two choices
   qh_QHpointer = 1     access globals via a pointer
                        enables qh_saveqhull() and qh_restoreqhull()
		= 0     qh_qhstat is a static data structure
		        only one instance of qhull() can be active at a time
			default value
   qh_QHpointer is defined in qhull.h
*/
#if qh_QHpointer
#define qhstat qh_qhstat->
extern qhstatT *qh_qhstat;  /* allocated in stat.c */
#else
#define qhstat qh_qhstat.
extern qhstatT qh_qhstat;  /* allocated in stat.c */
#endif

/*-------------------------------------------
-intrealT-  union of integer and real, used for statistics
*/
union intrealT {
    int i;
    realT r;
};

/*--------------------------------------------
-qhstatT- global data structure for statistics
*/
struct qhstatT {  
  intrealT stats[ZEND];  /* integer and real statistics */
  unsigned char id[ZEND];  /* id's in print order */
  char *doc[ZEND];     /* array of documentation strings */
  short int count[ZEND];   /* -1 if none, else index of count to use */
  char type[ZEND];      /* type, see ztypes above */
  char printed[ZEND];   /* true, if statistic has been printed */
  intrealT init[ZTYPEend];  /* initial values by types, set initstatistics */

  int next;           /* next index for zdef_ */
  int precision;      /* index for precision problems */
  int tempi;
  realT tempr;
};

/* ========== -functions- ===========
   see also qhull.h

-allstatA           define statistics in groups of 20
-qh_allstatistics   reset printed flag for all statistics
-collectstatistics  collect statistics for qh facet_list
-freestatistics     free memory used for statistics
-initstatistics     allocate and initialize statistics
-newstats	    returns True if statistics for zdoc
-nostatistic        true if no statistic to print
-printallstatistics print all statistics
-printstatistics    print statistics to a file
-printstatlevel     print level information for a statistic
-printstats         print statistics for a zdoc group
-stddev		    compute the standard deviation and average from statistics
*/
  
void    qh_allstatA(void);
void    qh_allstatB(void);
void    qh_allstatC(void);
void    qh_allstatD(void);
void    qh_allstatE(void);
void    qh_allstatF(void);
void    qh_allstatG(void);
void    qh_allstatH(void);
void    qh_allstatistics (void);
void    qh_collectstatistics (void);
void	qh_freestatistics (void);
void    qh_initstatistics (void);
boolT 	qh_newstats (int index, int *nextindex);
boolT 	qh_nostatistic (int i);
void    qh_printallstatistics (FILE *fp, char *string);
void    qh_printstatistics (FILE *fp, char *string);
void  	qh_printstatlevel (FILE *fp, int id, int start);
void  	qh_printstats (FILE *fp, int index, int *nextindex);
realT   qh_stddev (int num, realT tot, realT tot2, realT *ave);

#endif   /* qhDEFstat */
