/* set.c -- implements set manipulations needed for quickhull 

   see README and set.h
   
   copyright (c) 1993-1995 The Geometry Center        
*/

#include <stdio.h>
#include <memory.h>
#include <string.h>
#include "set.h"
#include "mem.h"

#ifndef qhDEFqhull
void    qh_errexit(int exitcode, void *, void *);
#endif

/*----------- internal macros -------------------
-SETsizeaddr_(set) - return pointer to actual size+1 of set (set CANNOT be NULL!!)
    *SETsizeaddr==NULL or e[*SETsizeaddr-1].p==NULL
*/
#define SETsizeaddr_(set) (&((set)->e[(set)->maxsize].i))

/*============ functions in alphabetical order ===================*/
  
/*----------------------------------------
-setaddnth- adds newelem as n'th element of sorted or unsorted set
  setp and newelem must be defined
  set may be a temp set
  nth=0 is first element
  errors if nth is out of bounds
*/
void qh_setaddnth(setT **setp, int nth, void *newelem) {
  int *sizep, oldsize, i;
  void **oldp, **newp;

  if (!*setp || !*(sizep= SETsizeaddr_(*setp))) {
    qh_setlarger(setp);
    sizep= SETsizeaddr_(*setp);
  }
  oldsize= *sizep - 1;
  if (nth < 0 || nth > oldsize) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setaddnth): nth %d is out-of-bounds for set:\n", nth);
    qh_setprint (qhmem.ferr, "", *setp);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  (*sizep)++;
  oldp= SETelemaddr_(*setp, oldsize, void);   /* NULL */
  newp= oldp+1;
  for (i= oldsize-nth+1; i--; )  /* move at least NULL  */
    *(newp--)= *(oldp--);       /* may overwrite *sizep */
  *newp= newelem;
} /* setaddnth */


/*----------------------------------------
-setaddsorted- adds an element to a sorted set
  setp and newelem must be defined
  set may be a temp set
  nop if newelem already in set
*/
void qh_setaddsorted(setT **setp, void *newelem) {
  int newindex=0;
  void *elem, **elemp;

  FOREACHelem_(*setp) {          /* could use binary search instead */
    if (elem < newelem)
      newindex++;
    else if (elem == newelem)
      return;
    else
      break;
  }
  qh_setaddnth(setp, newindex, newelem);
} /* setaddsorted */


/*----------------------------------------
-setappend- appends an element to a set
  set may be a temp set
  *setp and newelem may be NULL
*/
void qh_setappend(setT **setp, void *newelem) {
  int *sizep;
  void **endp;

  if (!newelem)
    return;
  if (!*setp || !*(sizep= SETsizeaddr_(*setp))) {
    qh_setlarger(setp);
    sizep= SETsizeaddr_(*setp);
  }
  *(endp= &((*setp)->e[(*sizep)++ - 1].p))= newelem;
  *(++endp)= NULL;
} /* setappend */

/*----------------------------------------
-setappend_set- appends a set to a set
  *setp and set may be NULL
  setp can not be a temp set
*/
void qh_setappend_set(setT **setp, setT *setA) {
  int *sizep, sizeA, size;
  setT *oldset;

  if (!setA)
    return;
  SETreturnsize_(setA, sizeA);
  if (!*setp)
    *setp= qh_setnew (sizeA);
  sizep= SETsizeaddr_(*setp);
  if (!(size= *sizep))
    size= (*setp)->maxsize;
  else
    size--;
  if (size + sizeA > (*setp)->maxsize) {
    oldset= *setp;
    *setp= qh_setcopy (oldset, sizeA);
    qh_setfree (&oldset);
    sizep= SETsizeaddr_(*setp);
  }
  *sizep= size+sizeA+1;   /* memcpy may overwrite */
  if (sizeA > 0) 
    memcpy((char *)&((*setp)->e[size].p), (char *)&(setA->e[0].p), SETelemsize *(sizeA+1));
} /* setappend_set */


/*----------------------------------------
-setappend2ndlast- makes newelem the next to the last element in set
  set must have at least one element, newelem must be defined
  set may be a temp set
*/
void qh_setappend2ndlast(setT **setp, void *newelem) {
  int *sizep;
  void **endp, **lastp;
  
  if (!*setp || !*(sizep= SETsizeaddr_(*setp))) {
    qh_setlarger(setp);
    sizep= SETsizeaddr_(*setp);
  }
  endp= SETelemaddr_(*setp, (*sizep)++ -1, void); /* NULL */
  lastp= endp-1;
  *(endp++)= *lastp;
  *endp= NULL;    /* may overwrite *sizep */
  *lastp= newelem;
} /* setappend2ndlast */


/*----------------------------------------
-setcheck- check set for validity
*/
void qh_setcheck(setT *set, char *typename, int id) {
  int maxsize, size;
  int waserr= 0;

  if (!set)
    return;
  SETreturnsize_(set, size);
  maxsize= set->maxsize;
  if (size > maxsize || !maxsize) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setcheck): actual size %d of %s%d is greater than max size %d\n",
	     size, typename, id, maxsize);
    waserr= 1;
  }else if (set->e[size].p) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setcheck): %s%d (size %d max %d) is not null terminated.\n",
	     typename, id, maxsize, size-1);
    waserr= 1;
  }
  if (waserr) {
    qh_setprint (qhmem.ferr, "ERRONEOUS", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
} /* setcheck */


/*----------------------------------------
-setcompact- compact NULLs in an unsorted set
  set may be NULL
returns:
  updates set
notes:
  is faster to swap tail of set into holes, like qh_setdel
*/
void qh_setcompact(setT *set) {
  int size;
  void **destp, **elemp, **endp, **firstp;

  if (!set)
    return;
  SETreturnsize_(set, size);
  destp= elemp= firstp= SETaddr_(set, void);
  endp= destp + size;
  while (1) {
    if (!(*destp++ = *elemp++)) {
      destp--;
      if (elemp > endp)
	break;
    }
  }
  qh_settruncate (set, destp-firstp);
} /* setcompact */


/*----------------------------------------
-setcopy- copies a sorted or unsorted set into another
returns:
  new set is actual size of old set plus extra
*/
setT *qh_setcopy(setT *set, int extra) {
  setT *newset;
  int size;

  if (extra < 0)
    extra= 0;
  SETreturnsize_(set, size);
  newset= qh_setnew(size+extra);
  *SETsizeaddr_(newset)= size+1;    /* memcpy may overwrite */
  memcpy((char *)&(newset->e[0].p), (char *)&(set->e[0].p), SETelemsize *(size+1));
  return (newset);
} /* setcopy */


/*----------------------------------------
-setdel- deletes oldelem from unsorted set.
   if found, overwrites newlelem with lastelem
   set may be NULL, oldelem must not be NULL;
returns:
  returns oldelem if it was deleted (use type conversion)
notes:
  only deletes one copy of oldelem in set
*/
void *qh_setdel(setT *set, void *oldelem) {
  void **elemp, **lastp;
  int *sizep;

  if (!set)
    return NULL;
  elemp= SETaddr_(set, void);
  while (*elemp != oldelem && *elemp)
    elemp++;
  if (*elemp) {
    sizep= SETsizeaddr_(set);
    if (!(*sizep)--)         /*  if was a full set */
      *sizep= set->maxsize;  /*     *sizep= (maxsize-1)+ 1 */
    lastp= SETelemaddr_(set, *sizep-1, void);
    *elemp= *lastp;      /* may overwrite itself */
    *lastp= NULL;
    return oldelem;
  }
  return NULL;
} /* setdel */


/*----------------------------------------
-setdellast- return last element of set or NULL (use type conversion)
   delete element from set
   set may be NULL
*/
void *qh_setdellast(setT *set) {
  int setsize;  /* actually, actual_size + 1 */
  int maxsize;
  int *sizep;
  void *returnvalue;
  
  if (!set || !(set->e[0].p))
    return NULL;
  sizep= SETsizeaddr_(set);
  if ((setsize= *sizep)) {
    returnvalue= set->e[setsize - 2].p;
    set->e[setsize - 2].p= NULL;
    (*sizep)--;
  }else {
    maxsize= set->maxsize;
    returnvalue= set->e[maxsize - 1].p;
    set->e[maxsize - 1].p= NULL;
    *sizep= maxsize;
  }
  return returnvalue;
} /* setdellast */


/*----------------------------------------
-setdelnth- deletes nth element from unsorted set 
  errors if nth invalid
  returns the element (use type conversion)
*/
void *qh_setdelnth(setT *set, int nth) {
  void **elemp, **lastp, *elem;
  int *sizep;


  elemp= SETelemaddr_(set, nth, void);
  sizep= SETsizeaddr_(set);
  if (!(*sizep)--)         /*  if was a full set */
    *sizep= set->maxsize;  /*     *sizep= (maxsize-1)+ 1 */
  if (nth < 0 || nth >= *sizep) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setaddnth): nth %d is out-of-bounds for set:\n", nth);
    qh_setprint (qhmem.ferr, "", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  lastp= SETelemaddr_(set, *sizep-1, void);
  elem= *elemp;
  *elemp= *lastp;      /* may overwrite itself */
  *lastp= NULL;
  return elem;
} /* setdelnth */

/*----------------------------------------
-setdelnthsorted- deletes nth element from sorted set
  sort order is undefined
  errors if nth invalid
  returns the element (use type conversion)
  see also: setnew_delnthsorted
*/
void *qh_setdelnthsorted(setT *set, int nth) {
  void **newp, **oldp, *elem;
  int *sizep;

  sizep= SETsizeaddr_(set);
  if (nth < 0 || (*sizep && nth >= *sizep-1) || nth >= set->maxsize) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setaddnth): nth %d is out-of-bounds for set:\n", nth);
    qh_setprint (qhmem.ferr, "", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  newp= SETelemaddr_(set, nth, void);
  elem= *newp;
  oldp= newp+1;
  while ((*(newp++)= *(oldp++)))
    ; /* copy remaining elements and NULL */
  if (!(*sizep)--)         /*  if was a full set */
    *sizep= set->maxsize;  /*     *sizep= (max size-1)+ 1 */
  return elem;
} /* setdelnthsorted */


/*----------------------------------------
-setdelsorted- deletes oldelem from sorted set
  sort order is undefined
  set may be NULL
  returns oldelem if it was deleted
*/
void *qh_setdelsorted(setT *set, void *oldelem) {
  void **newp, **oldp;
  int *sizep;

  if (!set)
    return NULL;
  newp= SETaddr_(set, void);
  while(*newp != oldelem && *newp)
    newp++;
  if (*newp) {
    oldp= newp+1;
    while ((*(newp++)= *(oldp++)))
      ; /* copy remaining elements */
    sizep= SETsizeaddr_(set);
    if (!(*sizep)--)    /*  if was a full set */
      *sizep= set->maxsize;  /*     *sizep= (max size-1)+ 1 */
    return oldelem;
  }
  return NULL;
} /* setdelsorted */


/*----------------------------------------
-setduplicate- duplicate a set of elements
notes:
  use setcopy if retaining old elements
*/
setT *qh_setduplicate (setT *set, int elemsize) {
  void		*elem, **elemp, *newElem;
  setT		*newSet;
  int		size;
  
  if (!(size= qh_setsize (set)))
    return NULL;
  newSet= qh_setnew (size);
  FOREACHelem_(set) {
    newElem= qh_memalloc (elemsize);
    memcpy (newElem, elem, elemsize);
    qh_setappend (&newSet, newElem);
  }
  return newSet;
} /* setduplicate */


/*----------------------------------------
-setequal- returns 1 if two sorted sets are equal, otherwise returns 0
    either set may be NULL
*/
int qh_setequal(setT *setA, setT *setB) {
  void **elemAp, **elemBp;
  int sizeA, sizeB;
  
  SETreturnsize_(setA, sizeA);
  SETreturnsize_(setB, sizeB);
  if (sizeA != sizeB)
    return 0;
  if (!sizeA)
    return 1;
  elemAp= SETaddr_(setA, void);
  elemBp= SETaddr_(setB, void);
  if (!memcmp((char *)elemAp, (char *)elemBp, sizeA*SETelemsize))
    return 1;
  return 0;
} /* setequal */


/*----------------------------------------
-setequal_except- returns 1 if two sorted sets are equal except for 2 elements
  neither set may be NULL
  false if either skip is missing
  if second skip is NULL, 
     can skip any one element
*/
int qh_setequal_except (setT *setA, void *skipelemA, setT *setB, void *skipelemB) {
  void **elemA, **elemB;
  int skip=0;

  elemA= SETaddr_(setA, void);
  elemB= SETaddr_(setB, void);
  while (1) {
    if (*elemA == skipelemA) {
      skip++;
      elemA++;
    }
    if (skipelemB) {
      if (*elemB == skipelemB) {
        skip++;
        elemB++;
      }
    }else if (*elemA != *elemB) {
      skip++;
      if (!(skipelemB= *elemB++))
        return 0;
    }
    if (!*elemA)
      break;
    if (*elemA++ != *elemB++) 
      return 0;
  }
  if (skip != 2 || *elemB)
    return 0;
  return 1;
} /* setequal_except */
  

/*----------------------------------------
-setequal_skip- returns 1 if two sorted sets are equal except for skips
  neither set may be NULL
  false if different size
*/
int qh_setequal_skip (setT *setA, int skipA, setT *setB, int skipB) {
  void **elemA, **elemB, **skipAp, **skipBp;

  elemA= SETaddr_(setA, void);
  elemB= SETaddr_(setB, void);
  skipAp= SETelemaddr_(setA, skipA, void);
  skipBp= SETelemaddr_(setB, skipB, void);
  while (1) {
    if (elemA == skipAp)
      elemA++;
    if (elemB == skipBp)
      elemB++;
    if (!*elemA)
      break;
    if (*elemA++ != *elemB++) 
      return 0;
  }
  if (*elemB)
    return 0;
  return 1;
} /* setequal_skip */
  

/*----------------------------------------
-setfree- frees the space occupied by a sorted or unsorted set
  set may be NULL
*/
void qh_setfree(setT **setp) {
  int size;
  void **freelistp;
  
  if (*setp) {
    size= sizeof(setT) + ((*setp)->maxsize)*SETelemsize; 
    if (size <= qhmem.LASTsize) {
      qh_memfree_(*setp, size, freelistp);
    }else
      qh_memfree (*setp, size);
    *setp= NULL;
  }
} /* setfree */


/*----------------------------------------
-setfree2- frees the space occupied by a set and its elements
  set may be NULL
*/
void qh_setfree2 (setT **setp, int elemsize) {
  void		*elem, **elemp;
  
  FOREACHelem_(*setp)
    qh_memfree (elem, elemsize);
  qh_setfree (setp);
} /* setfree2 */


      
/*----------------------------------------
-setfreelong- frees a set only if it's in long memory
  set may be NULL
*/
void qh_setfreelong(setT **setp) {
  int size;
  
  if (*setp) {
    size= sizeof(setT) + ((*setp)->maxsize)*SETelemsize; 
    if (size > qhmem.LASTsize) {
      qh_memfree (*setp, size);
      *setp= NULL;
    }
  }
} /* setfreelong */


/*----------------------------------------
-setin- returns 1 if setelem is in a set, 0 otherwise
  set may be NULL or unsorted
*/
int qh_setin(setT *set, void *setelem) {
  void *elem, **elemp;

  FOREACHelem_(set) {
    if (elem == setelem)
      return 1;
  }
  return 0;
} /* setin */


/*----------------------------------------
-setindex- returns the index of elem in set.   If none, returns -1
  set may be NULL and may contain nulls.
*/
int qh_setindex(setT *set, void *atelem) {
  void **elem;
  int size, i;

  SETreturnsize_(set, size);
  if (size > set->maxsize)
    return -1;
  elem= SETaddr_(set, void);
  for (i=0; i<size; i++) {
    if (*elem++ == atelem)
      return i;
  }
  return -1;
} /* setindex */


/*----------------------------------------
-setlarger- returns a larger set that contains elements of *setp
  the set is at least twice as large
  updates qhmem.tempstack if needed
*/
void qh_setlarger(setT **oldsetp) {
  int size= 1, *sizep;
  setT *newset, *set, **setp, *oldset;
  void **oldp, **newp;

  if (*oldsetp) {
    oldset= *oldsetp;
    SETreturnsize_(oldset, size);
    qhmem.cntlarger++;
    qhmem.totlarger += size+1;
    newset= qh_setnew(2 * size);
    oldp= SETaddr_(oldset, void);
    newp= SETaddr_(newset, void);
    memcpy((char *)newp, (char *)oldp, (size+1) * SETelemsize);
    sizep= SETsizeaddr_(newset);
    *sizep= size+1;
    FOREACHset_((setT *)qhmem.tempstack) {
      if (set == oldset)
	*(setp-1)= newset;
    }
    qh_setfree(oldsetp);
  }else 
    newset= qh_setnew(3);
  *oldsetp= newset;
} /* setlarger */


/*----------------------------------------
-setlast- return last element of set or NULL (use type conversion)
   set may be NULL
*/
void *qh_setlast(setT *set) {
  int size;

  if (set) {
    size= *SETsizeaddr_(set);
    if (!size) 
      return SETelem_(set, set->maxsize - 1);
    else if (size > 1)
      return SETelem_(set, size - 2);
  }
  return NULL;
} /* setlast */


/*----------------------------------------
-setnew- creates and allocates space for a set
    setsize means the number of elements (NOT including the NULL terminator)
    use qh_settemp/qh_setfreetemp if set is temporary
*/
setT *qh_setnew(int setsize) {
  setT *set;
  int sizereceived, size;
  void **freelistp;

  if (!setsize)
    setsize++;
  size= sizeof(setT) + setsize * SETelemsize;
  if ((unsigned) size <= (unsigned) qhmem.LASTsize) {
    qh_memalloc_(size, freelistp, set, setT);
#ifndef qh_NOmem
    sizereceived= qhmem.sizetable[ qhmem.indextable[size]];
    if (sizereceived > size) 
      setsize += (sizereceived - size)/SETelemsize;
#endif
  }else
    set= (setT*)qh_memalloc (size);
  set->maxsize= setsize;
  set->e[setsize].i= 1;
  set->e[0].p= NULL;
  return (set);
} /* setnew */


/*----------------------------------------
-setnew_delnthsorted- creates a sorted set not containing nth element
  the new set may have prepended undefined entries
  set must be defined
  checks nth
  see also: setdelnthsorted
*/
setT *qh_setnew_delnthsorted(setT *set, int size, int nth, int prepend) {
  setT *newset;
  void **oldp, **newp;
  int tailsize= size - nth -1, newsize;

  if (tailsize < 0) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setaddnth): nth %d is out-of-bounds for set:\n", nth);
    qh_setprint (qhmem.ferr, "", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  newsize= size-1 + prepend;
  newset= qh_setnew(newsize);
  newset->e[newset->maxsize].i= newsize+1;  /* may be overwritten */
  oldp= SETaddr_(set, void);
  newp= SETaddr_(newset, void) + prepend;
  switch (nth) {
  case 0:
    break;
  case 1:
    *(newp++)= *oldp++;
    break;
  case 2:
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    break;
  case 3:
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    break;
  case 4:
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    break;
  default:
    memcpy((char *)newp, (char *)oldp, nth * SETelemsize);
    newp += nth;
    oldp += nth;
    break;
  }
  oldp++;
  switch (tailsize) {
  case 0:
    break;
  case 1:
    *(newp++)= *oldp++;
    break;
  case 2:
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    break;
  case 3:
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    break;
  case 4:
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    *(newp++)= *oldp++;
    break;
  default:
    memcpy((char *)newp, (char *)oldp, tailsize * SETelemsize);
    newp += tailsize;
  }
  *newp= NULL;
  return(newset);
} /* setnew_delnthsorted */


/*----------------------------------------
-setprint- print set elements to fp
notes:
  never errors
*/
void qh_setprint(FILE *fp, char* string, setT *set) {
  int size, k;

  if (!set)
    fprintf (fp, "%s set is null\n", string);
  else {
    SETreturnsize_(set, size);
    fprintf (fp, "%s set=%p maxsize=%d size=%d elems=",
	     string, set, set->maxsize, size);
    if (size > set->maxsize)
      size= set->maxsize+1;
    for (k=0; k<size; k++)
      fprintf(fp, " %p", set->e[k].p);
    fprintf(fp, "\n");
  }
} /* setprint */

/*----------------------------------------
-setreplace- replaces oldelem in set with newelem
   errors if oldelem not in the set
   if newelem is NULL then FOREACH no longer works
*/
void qh_setreplace(setT *set, void *oldelem, void *newelem) {
  void **elemp;
  
  elemp= SETaddr_(set, void);
  while(*elemp != oldelem && *elemp)
    elemp++;
  if (*elemp)
    *elemp= newelem;
  else {
    fprintf (qhmem.ferr, "qhull internal error (qh_setreplace): elem %p not found in set\n",
       oldelem);
    qh_setprint (qhmem.ferr, "", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
} /* setreplace */


/*----------------------------------------
-setsize- returns the size of a set
  same as SETreturnsize_(set)
*/
int qh_setsize(setT *set) {
  int size, *sizep;
  
  if (!set)
    return (0);
  sizep= SETsizeaddr_(set);
  if ((size= *sizep)) {
    size--;
    if (size > set->maxsize) {
      fprintf (qhmem.ferr, "qhull internal error (qh_setsize): current set size %d is greater than maximum size %d\n",
	       size, set->maxsize);
      qh_setprint (qhmem.ferr, "set: ", set);
      qh_errexit (qhmem_ERRqhull, NULL, NULL);
    }
  }else
    size= set->maxsize;
  return size;
} /* setsize */

/*----------------------------------------
-settemp- return a stacked, temporary set
  use settempfree or settempfree_all to release from qhmem.tempstack
  see also qh_setnew
*/
setT *qh_settemp(int setsize) {
  setT *newset;
  
  newset= qh_setnew (setsize);
  qh_setappend ((setT **)&qhmem.tempstack, newset);
  if (qhmem.IStracing >= 5)
    fprintf (qhmem.ferr, "qh_settemp: temp set %p of %d elements, depth %d\n",
       newset, newset->maxsize, qh_setsize ((setT*)qhmem.tempstack));
  return newset;
} /* settemp */

/*----------------------------------------
-settempfree- free temporary set at top of qhmem.tempstack
  nop if NULL
  errors if set not from previous qh_settemp
    locate source by T2 and find mis-matching qh_settemp
*/
void qh_settempfree(setT **set) {
  setT *stackedset;

  if (!*set)
    return;
  stackedset= qh_settemppop ();
  if (stackedset != *set) {
    qh_settemppush(stackedset);
    fprintf (qhmem.ferr, "qhull internal error (qh_settempfree): set %p (size %d) was not last temporary allocated (depth %d, set %p, size %d)\n",
	     *set, qh_setsize(*set), qh_setsize((setT*)qhmem.tempstack)+1,
	     stackedset, qh_setsize(stackedset));
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  qh_setfree (set);
} /* settempfree */

/*----------------------------------------
-settempfree_all- free all temporary sets in qhmem.tempstack
*/
void qh_settempfree_all(void) {
  setT *set, **setp;

  FOREACHset_((setT *)qhmem.tempstack) 
    qh_setfree(&set);
  qh_setfree((setT **)&qhmem.tempstack);
} /* settempfree_all */

/*----------------------------------------
-settemppop- pop and return temporary set from qhmem.tempstack (makes it permanent)
*/
setT *qh_settemppop(void) {
  setT *stackedset;
  
  stackedset= (setT*)qh_setdellast((setT *)qhmem.tempstack);
  if (!stackedset) {
    fprintf (qhmem.ferr, "qhull internal error (qh_settemppop): pop from empty temporary stack\n");
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  if (qhmem.IStracing >= 5)
    fprintf (qhmem.ferr, "qh_settemppop: depth %d temp set %p of %d elements\n",
       qh_setsize((setT*)qhmem.tempstack)+1, stackedset, qh_setsize(stackedset));
  return stackedset;
} /* settemppop */

/*----------------------------------------
-settemppush- push temporary set unto qhmem.tempstack (makes it temporary)
  duplicates settemp() for tracing
*/
void qh_settemppush(setT *set) {
  
  qh_setappend ((setT**)&qhmem.tempstack, set);
  if (qhmem.IStracing >= 5)
    fprintf (qhmem.ferr, "qh_settemppush: depth %d temp set %p of %d elements\n",
    qh_setsize((setT*)qhmem.tempstack), set, qh_setsize(set));
} /* settemppush */

 
/*----------------------------------------
-settruncate- truncate set to size elements
  set must be defined
*/
void qh_settruncate (setT *set, int size) {

  if (size < 0 || size > set->maxsize) {
    fprintf (qhmem.ferr, "qhull internal error (qh_settruncate): size %d out of bounds for set:\n", size);
    qh_setprint (qhmem.ferr, "", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  set->e[set->maxsize].i= size+1;   /* maybe overwritten */
  set->e[size].p= NULL;
} /* settruncate */
    
/*----------------------------------------
-setunique- add element if it isn't already
  returns 1 if it's appended
*/
int qh_setunique (setT **set, void *elem) {

  if (!qh_setin (*set, elem)) {
    qh_setappend (set, elem);
    return 1;
  }
  return 0;
} /* setunique */
    
/*----------------------------------------
-setzero- zero remainder of set and set its size
  set must be defined
*/
void qh_setzero (setT *set, int index, int size) {
  int count;

  if (index < 0 || index >= size || size > set->maxsize) {
    fprintf (qhmem.ferr, "qhull internal error (qh_setzero): index %d or size %d out of bounds for set:\n", index, size);
    qh_setprint (qhmem.ferr, "", set);
    qh_errexit (qhmem_ERRqhull, NULL, NULL);
  }
  set->e[set->maxsize].i=  size+1;  /* may be overwritten */
  count= size - index + 1;   /* +1 for NULL terminator */
  memset ((char *)SETelemaddr_(set, index, void), 0, count * SETelemsize);
} /* setzero */

    
