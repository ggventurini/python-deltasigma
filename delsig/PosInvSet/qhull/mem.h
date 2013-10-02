/* mem.h - prototypes for memory management functions

   see README, mem.c and set.h
   
   for error handling, writes message and calls
          qh_errexit (qhmem_ERRmem, NULL, NULL) if insufficient memory
   and
          qh_errexit (qhmem_ERRqhull, NULL, NULL) otherwise
   
   copyright (c) 1993-1995, The Geometry Center
*/

#ifndef qhDEFmem
#define qhDEFmem

/* mem.c implements Quickfit memory allocation for about 20% time
   savings.  If it fails on your machine, try to locate the
   problem, and send the answer to qhull@geom.umn.edu.  If this can
   not be done, define qh_NOmem to use malloc/free instead.

   #define qh_NOmem
*/

/* to avoid bus errors, memory allocation must consider alignment requirements.
   malloc() automatically takes care of alignment.   Since mem.c manages
   its own memory, we need to explicitly specify alignment in
   qh_meminitbuffers().

   A safe choice is sizeof(double).  sizeof(float) may be used if doubles 
   do not occur in data structures and pointers are the same size.  Be careful
   of machines (e.g., DEC Alpha) with large pointers.  If gcc is available, 
   use __alignof__(double) or fmax_(__alignof__(float), __alignof__(void *)).

   see user.h for qhull's alignment
*/

#define qhmem_ERRmem 4    /* matches qh_ERRmem in qhull.h */
#define qhmem_ERRqhull 5  /* matches qh_ERRqhull in qhull.h */

/* On 64-bit machines, a pointer may be larger than an 'int'.  The following
  assumes that a 'long' holds a pointer.  This is true for DEC's Alpha.
*/

typedef unsigned long ptr_intT;  /* for casting a void* to an integer-type */

/*---------------------------------------
-qhmemT - global memory structure for mem.c

   users should ignore qhmem except for writing extensions
   
   qhmem could be swapable like qh and qhstat, but then
   multiple qh's and qhmem's would need to keep in synch.  
   A swapable qhmem would also waste memory buffers.  As long
   as memory operations are atomic, there is no problem with
   multiple qh structures being active at the same time.
   If you need separate address spaces, you can swap the
   contents of qhmem.
*/

typedef struct qhmemT qhmemT;
extern qhmemT qhmem;  /* allocated in mem.c */

struct qhmemT {               /* global memory management variables */
  int      BUFsize;	      /* size of memory allocation buffer */
  int      BUFinit;	      /* initial size of memory allocation buffer */
  int      TABLEsize;         /* actual number of sizes in free list table */
  int      NUMsizes;          /* maximum number of sizes in free list table */
  int      LASTsize;          /* last size in free list table */
  int      ALIGNmask;         /* worst-case alignment, must be 2^n-1 */
  void	 **freelists;          /* free list table, linked by offset 0 */
  int     *sizetable;         /* size of each freelist */
  int     *indextable;        /* size->index table */
  void    *curbuffer;         /* current buffer, linked by offset 0 */
  void    *freemem;           /*   free memory in curbuffer */
  int 	   freesize;          /*   size of free memory in bytes */
  void 	  *tempstack;         /* stack of temporary memory, managed by users */
  FILE    *ferr;              /* file for reporting errors */
  int      IStracing;         /* =5 if tracing memory allocations */
  int cntquick;          /* count of quick allocations */
                         /* remove statistics doesn't effect speed */
  int cntshort;          /* count of short allocations */
  int cntlong;           /* count of long allocations */
  int curlong;           /* current count of inuse, long allocations */
  int freeshort;	      /* count of short memfrees */
  int freelong;	      /* count of long memfrees */
  int totshort;          /* total size of short allocations */
  int totlong;           /* total size of long allocations */
  int maxlong;           /* maximum totlong */
  int cntlarger;         /* count of setlarger's */
  int totlarger;         /* total copied by setlarger */
};


/* ======= -macros =========== 

qh_memalloc_(size, freelistp, object, type)  returns object of size bytes 
	assumes size<=qhmem.LASTsize and void **freelistp is a temp

qh_memfree_(object, size, freelistp) free up quick object
	object may be NULL
	assumes size<=qhmem.LASTsize and void **freelistp is a temp
*/

#ifdef qh_NOmem

#define qh_memalloc_(size, freelistp, object, type) {\
  object= (type*)qh_memalloc (size); }

#define qh_memfree_(object, size, freelistp) {\
  qh_memfree (object, size); }

#else /* !qh_NOmem */

#define qh_memalloc_(size, freelistp, object, type) {\
  freelistp= qhmem.freelists + qhmem.indextable[size];\
  if ((object= (type*)*freelistp)) {\
    qhmem.cntquick++;  \
    *freelistp= *((void **)*freelistp);\
  }else object= (type*)qh_memalloc (size);}

#define qh_memfree_(object, size, freelistp) {\
  if (object) { \
    qhmem .freeshort++;\
    freelistp= qhmem.freelists + qhmem.indextable[size];\
    *((void **)object)= *freelistp;\
    *freelistp= object;}}

#endif /* !qh_NOmem */

/* ======= -functions =========== 

	see mem.c for definitions

	User level functions
-memalloc	allocate memory
-memfree	free memory
-memstatistics  print memory statistics

	Initialization and termination functions
-meminit	initialize memory
-meminitbuffers	initialize memory buffers
-memsize	define a free list for a size
-memsetup	set up memory (activates memalloc/free)
-memfreeshort	free up all memory buffers
*/

/*---------- -prototypes in alphabetical order -----------*/

void *qh_memalloc(int insize);
void qh_memfree (void *object, int size);
void qh_memfreeshort (int *curlong, int *totlong);
void qh_meminit (FILE *ferr);
void qh_meminitbuffers (int tracelevel, int alignment, int numsizes,
			int bufsize, int bufinit);
void qh_memsetup (void);
void qh_memsize(int size);
void qh_memstatistics (FILE *fp);

#endif /* qhDEFmem */
