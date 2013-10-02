/* global.c -- contains all the globals of the qhull application

   see README 
   
   see qhull.h for qh globals and function prototypes

   see qhull_a.h for internal functions

   copyright (c) 1993-1995, The Geometry Center
*/

#include "qhull_a.h"

#if qh_QHpointer
qhT *qh_qh= NULL;
#else
qhT qh_qh= {0};   /* remove "= {0}" if this causes a compiler error.  Also
		     qh_qhstat in stat.c and qhmem in mem.c.  */
#endif

/*-------------------------------------------
-appendprint- append output format to qh PRINTout unless already defined
*/
void qh_appendprint (qh_PRINT format) {
  int i;

  for (i=0; i < qh_PRINTEND; i++) {
    if (qh PRINTout[i] == format)
      break;
    if (!qh PRINTout[i]) {
      qh PRINTout[i]= format;
      break;
    }
  }
} /* appendprint */
     
/*-------------------------------------------
-freebuffers- free up global memory buffers
  must match initbuffers()
*/

void qh_freebuffers (void) {

  trace5((qh ferr, "qh_freebuffers: freeing up global memory buffers\n"));
  qh_memfree (qh NEARzero, qh hull_dim * sizeof(realT));
  qh_memfree (qh lower_threshold, (qh input_dim+1) * sizeof(realT));
  qh_memfree (qh upper_threshold, (qh input_dim+1) * sizeof(realT));
  qh_memfree (qh lower_bound, (qh input_dim+1) * sizeof(realT));
  qh_memfree (qh upper_bound, (qh input_dim+1) * sizeof(realT));
  qh_memfree (qh gm_matrix, (qh hull_dim+1) * qh hull_dim * sizeof(coordT));
  qh_memfree (qh gm_row, (qh hull_dim+1) * sizeof(coordT *));
  qh_setfree (&qh hash_table);
  qh_setfree (&qh other_points);
  qh_setfree (&qh del_vertices);
  qh_setfree (&qh searchset);
  /* qh facet_mergeset is a temp */
  qh NEARzero= qh lower_threshold= qh upper_threshold= NULL;
  qh lower_bound= qh upper_bound= NULL;
  qh gm_matrix= NULL;
  qh gm_row= NULL;
  if (qh line)
    free (qh line);
  if (qh half_space)
    free (qh half_space);
  if (qh feasible_point)
    free (qh feasible_point);
  if (qh feasible_string)
    free (qh feasible_string);
  if (qh temp_malloc)
    free (qh temp_malloc);
  qh line= qh feasible_string= NULL;
  qh half_space= qh feasible_point= qh temp_malloc= NULL;
  if (qh first_point && qh POINTSmalloc) {
    free(qh first_point);
    qh first_point= NULL;
  }
  trace5((qh ferr, "qh_freebuffers: finished\n"));
} /* freebuffers */


/*-------------------------------------------
-freeqhull- free global memory
  if allmem, frees all allocated data structures
  else, frees all long memory
    rest of memory freed by qh_memfreeshort();
*/
void qh_freeqhull (boolT allmem) {
  facetT *facet;
  vertexT *vertex;
  ridgeT *ridge, **ridgep;
  mergeT *merge, **mergep;

  trace1((qh ferr, "qh_freeqhull: free global memory\n"));
  qh NOerrexit= True;  /* no more setjmp */
  if (allmem) {
    qh_clearcenters (qh_ASnone);
    while ((vertex= qh vertex_list)) {
      if (vertex->next)
        qh_delvertex (vertex);
      else {
        qh_memfree (vertex, sizeof(vertexT));
        qh vertex_list= NULL;
      }
    }
  }else if (qh VERTEXneighbors) {
    FORALLvertices 
      qh_setfreelong (&(vertex->neighbors));
  }
  if (allmem) {
    FORALLfacets {
      if (!facet->visible) {
	FOREACHridge_(facet->ridges)
	  ridge->seen= False;
      }
    }
    FORALLfacets {
      FOREACHridge_(facet->ridges)
        ridge->seen ^= True;
    }
    while ((facet= qh facet_list)) {
      FOREACHridge_(facet->ridges) {
        if (ridge->seen) {
          qh_setfree(&(ridge->vertices));
          qh_memfree(ridge, sizeof(ridgeT));
        }else
          ridge->seen= True;
      }
      if (facet->next)
        qh_delfacet (facet);
      else {
        qh_memfree (facet, sizeof(facetT));
        qh facet_list= NULL;
      }
    }
  }else {
    FORALLfacets {
      qh_setfreelong (&(facet->outsideset));
      qh_setfreelong (&(facet->coplanarset));
      if (!facet->simplicial) {
        qh_setfreelong (&(facet->neighbors));
        qh_setfreelong (&(facet->ridges));
        qh_setfreelong (&(facet->vertices));
      }
    }
  }
  qh_setfree (&qh hash_table);
  FOREACHmerge_(qh facet_mergeset)  /* usually empty */
    qh_memfree (merge, sizeof(mergeT));
  qh_freebuffers();
  qh_freestatistics();
  qh_settempfree_all();
#if qh_QHpointer
  free (qh_qh);
  qh_qh= NULL;
#endif
} /* freeqhull */

/*---------------------------------------------
-init_A- called before error handling initialized
  argc/argv is used to initialize qh qhull_command.  argc may be 0
  qh_errexit() may not be used
*/
void qh_init_A (FILE *infile, FILE *outfile, FILE *errfile, int argc, char *argv[]) {
  qh_meminit (errfile);
  qh_initqhull_start (infile, outfile, errfile); 
  qh_init_qhull_command (argc, argv);
} /* init_A */

/*---------------------------------------------
-init_B- called after points are defined
  uses qh qhull_command from qh_init_qhull_command or qh_initflags
  qh_errexit() may be used
*/
void qh_init_B (coordT *points, int numpoints, int dim, boolT ismalloc) { 
  qh_initqhull_globals (points, numpoints, dim, ismalloc);
  qh_initqhull_mem();
  /* mem.c and set.c are initialized */
  qh_initqhull_buffers();
  qh_initthresholds (qh qhull_command);
  if (qh PROJECTinput || (qh DELAUNAY && qh PROJECTdelaunay))
    qh_projectinput();
  if (qh SCALEinput)
    qh_scaleinput();
  if (qh ROTATErandom >= 0) {
    qh_randommatrix (qh gm_matrix, qh hull_dim, qh gm_row);
    if (qh DELAUNAY) {
      int k, lastk= qh hull_dim-1;
      for (k= 0; k < lastk; k++) {
        qh gm_row[k][lastk]= 0.0;
        qh gm_row[lastk][k]= 0.0;
      }
      qh gm_row[lastk][lastk]= 1.0;
    }
    qh_gram_schmidt (qh hull_dim, qh gm_row);
    qh_rotateinput (qh gm_row);
  }
} /* init_B */

/*---------------------------------------------
-init_qhull_command- build qhull_command from argc/argv
*/
void qh_init_qhull_command(int argc, char *argv[]) {
  int i;

  if (argc)
    strcpy (qh qhull_command, argv[0]);
  for (i=1; i<argc; i++) {
    if (strlen (qh qhull_command) + strlen(argv[i]) + 1 < sizeof(qh qhull_command)) {
      strcat (qh qhull_command, " ");
      strcat (qh qhull_command, argv[i]);
    }else {
      fprintf (qh ferr, "qhull input error: more than %d characters in command line\n",
        (int)sizeof(qh qhull_command));
      exit (1);  /* can not use qh_errexit */
    }
  }
} /* init_qhull_command */

/*---------------------------------------------
-initflags- set flags and initialized constants from command line
  ignores first word (e.g., "qhull d"
  see 'prompt' in unix.c for documentation
  see also qh_initthresholds
  strtol/strtod may or may not skip trailing spaces
  sets qh qhull_command to command if needed
*/

void qh_initflags(char *command) {
  int k, i, lastproject;
  char *s= command, *t, *prev_s, *start, key;
  boolT isgeom= False, wasproject;
  realT r;

  if (command != &qh qhull_command[0]) {
    *qh qhull_command= '\0';
    strncat( qh qhull_command, command, sizeof( qh qhull_command));
  }
  while (*s && !isspace(*s))  /* skip program name */
    s++;
  while (*s) {
    while (*s && isspace(*s))
      s++;
    if (*s == '-')
      s++;
    prev_s= s;
    switch (*s++) {
    case 'd':
      qh_option ("delaunay", NULL, NULL);
      qh DELAUNAY= True;
      break;
    case 'f':
      qh_option ("facets", NULL, NULL);
      qh_appendprint (qh_PRINTfacets);
      break;
    case 'i':
      qh_option ("incidence", NULL, NULL);
      qh_appendprint (qh_PRINTincidences);
      break;
    case 'm':
      qh_option ("mathematica", NULL, NULL);
      qh_appendprint (qh_PRINTmathematica);
      break;
    case 'n':
      qh_option ("normals", NULL, NULL);
      qh_appendprint (qh_PRINTnormals);
      break;
    case 'o':
      qh_option ("offFile", NULL, NULL);
      qh_appendprint (qh_PRINToff);
      break;
    case 'p':
      qh_option ("points", NULL, NULL);
      qh_appendprint (qh_PRINTpoints);
      break;
    case 's':
      qh_option ("summary", NULL, NULL);
      qh PRINTsummary= True;
      break;
    case 'v':
      qh_option ("voronoi", NULL, NULL);
      qh VORONOI= True;
      qh DELAUNAY= True;
      break;
    case 'A':
      if (!isdigit(*s) && *s != '.' && *s != '-') 
	fprintf(qh ferr, "qhull warning: no maximum cosine angle given for option A.  Ignored.\n");
      else {
	if (*s == '-') {
	  qh premerge_cos= -qh_strtod (s, &s);
          qh_option ("Angle-premerge-", NULL, &qh premerge_cos);
	  qh PREmerge= True;
	}else {
	  qh postmerge_cos= qh_strtod (s, &s);
          qh_option ("Angle-postmerge", NULL, &qh postmerge_cos);
	  qh POSTmerge= True;
	}
	qh MERGING= True; 
      }
      break;
    case 'C':
      if (!isdigit(*s) && *s != '.' && *s != '-')
	fprintf(qh ferr, "qhull warning: no centrum radius given for option C.  Ignored.\n");
      else {
	if (*s == '-') {
	  qh premerge_centrum= -qh_strtod (s, &s);
          qh_option ("Centrum-premerge-", NULL, &qh premerge_centrum);
	  qh PREmerge= True;
	}else {
	  qh postmerge_centrum= qh_strtod (s, &s);
          qh_option ("Centrum-postmerge", NULL, &qh postmerge_centrum);
	  qh POSTmerge= True;
	}
	qh MERGING= True; 
      }
      break;
    case 'E':
      if (*s == '-')
	fprintf(qh ferr, "qhull warning: negative maximum roundoff given for option A.  Ignored.\n");
      else if (!isdigit(*s))
	fprintf(qh ferr, "qhull warning: no maximum roundoff given for option E.  Ignored.\n");
      else {
	qh DISTround= qh_strtod (s, &s);
        qh_option ("Error-roundoff", NULL, &qh DISTround);
	qh SETroundoff= True;
      }
      break;
    case 'H':
      start= s;
      qh HALFspace= True;
      qh_strtod (s, &t);
      while (t > s)  {
        if (*t && !isspace (*t)) {
	  if (*t == ',')
	    t++;
	  else
	    fprintf (qh ferr, "qhull warning: origin for Halfspace intersection should be 'Hn,n,n,...'\n");
	}
        s= t;
	qh_strtod (s, &t);
      }
      if (start < t) {
        if (!(qh feasible_string= (char*)calloc (t-start+1, 1))) {
          fprintf(qh ferr, "qhull error: insufficient memory for 'Hn,n,n'\n");
          qh_errexit(qh_ERRmem, NULL, NULL);
        }
        strncpy (qh feasible_string, start, t-start);
        qh_option ("Halfspace-about", NULL, NULL);
        qh_option (qh feasible_string, NULL, NULL);
      }else
        qh_option ("Halfspace", NULL, NULL);
      break;
    case 'R':
      if (!isdigit(*s))
	fprintf(qh ferr, "qhull warning: no random perturbation given for option R.  Ignored\n");
      else {
	qh RANDOMfactor= qh_strtod (s, &s);
        qh_option ("Random_perturbation", NULL, &qh RANDOMfactor);
        qh RANDOMdist= True;
      }
      break;
    case 'V':
      if (!isdigit(*s) && *s != '-')
	fprintf(qh ferr, "qhull warning: distance not given for option V.  Ignored\n");
      else {
	qh MINvisible= qh_strtod (s, &s);
        qh_option ("Visible", NULL, &qh MINvisible);
      }
      break;
    case 'U':
      if (!isdigit(*s) && *s != '-')
	fprintf(qh ferr, "qhull warning: distance not given for option U.  Ignored\n");
      else {
	qh MAXcoplanar= qh_strtod (s, &s);
        qh_option ("U-coplanar", NULL, &qh MAXcoplanar);
      }
      break;
    case 'W':
      if (*s == '-')
	fprintf(qh ferr, "qhull warning: negative width for option W.  Ignored.\n");
      else if (!isdigit(*s))
	fprintf(qh ferr, "qhull warning: no hull width given for option W.  Ignored\n");
      else {
	qh MINoutside= qh_strtod (s, &s);
        qh_option ("W-outside", NULL, &qh MINoutside);
        qh APPROXhull= True;
      }
      break;
    /************  sub menus ***************/
#pragma mark option_F
    case 'F':
      while (*s && !isspace(*s)) {
	switch(*s++) {
	case 'a':
	  qh_option ("Farea", NULL, NULL);
	  qh_appendprint (qh_PRINTarea);
	  qh GETarea= True;
	  break;
	case 'A':
	  qh_option ("FArea-total", NULL, NULL);
	  qh GETarea= True;
	  break;
        case 'c':
          qh_option ("Fcoplanars", NULL, NULL);
          qh_appendprint (qh_PRINTcoplanars);
          break;
        case 'C':
          qh_option ("FCentrums", NULL, NULL);
          qh_appendprint (qh_PRINTcentrums);
          break;
	case 'd':
          qh_option ("Fd-cdd-in", NULL, NULL);
	  qh CDDinput= True;
	  break;
	case 'D':
          qh_option ("FD-cdd-out", NULL, NULL);
	  qh CDDoutput= True;
	  break;
	case 'F':
	  qh_option ("FFacets-xridge", NULL, NULL);
          qh_appendprint (qh_PRINTfacets_xridge);
	  break;
        case 'i':
          qh_option ("Finner", NULL, NULL);
          qh_appendprint (qh_PRINTinner);
          break;
        case 'I':
          qh_option ("FIds", NULL, NULL);
          qh_appendprint (qh_PRINTids);
          break;
        case 'm':
          qh_option ("Fmerges", NULL, NULL);
          qh_appendprint (qh_PRINTmerges);
          break;
        case 'n':
          qh_option ("Fneighbors", NULL, NULL);
          qh_appendprint (qh_PRINTneighbors);
          break;
        case 'N':
          qh_option ("FNeighbors-vertex", NULL, NULL);
          qh_appendprint (qh_PRINTvneighbors);
          break;
        case 'o':
          qh_option ("Fouter", NULL, NULL);
          qh_appendprint (qh_PRINTouter);
          break;
	case 'O':
	  if (qh PRINToptions1st) {
	    qh_option ("FOptions", NULL, NULL);
	    qh_appendprint (qh_PRINToptions);
	  }else 
	    qh PRINToptions1st= True;
	  break;
	case 'p':
	  qh_option ("Fpoint-intersect", NULL, NULL);
	  qh_appendprint (qh_PRINTpointintersect);
	  break;
	case 'P':
	  qh_option ("FPoint-nearest", NULL, NULL);
	  qh_appendprint (qh_PRINTpointnearest);
	  break;
	case 'Q':
	  qh_option ("FQhull", NULL, NULL);
	  qh_appendprint (qh_PRINTqhull);
	  break;
        case 's':
          qh_option ("Fsummary", NULL, NULL);
          qh_appendprint (qh_PRINTsummary);
          break;
        case 'S':
          qh_option ("FSize", NULL, NULL);
          qh_appendprint (qh_PRINTsize);
          qh GETarea= True;
          break;
        case 'v':
          qh_option ("Fvertices", NULL, NULL);
          qh_appendprint (qh_PRINTvertices);
          break;
        case 'V':
          qh_option ("FVertex-average", NULL, NULL);
          qh_appendprint (qh_PRINTaverage);
          break;
	default:
	  s--;
	  fprintf (qh ferr, "qhull warning: unknown 'F' output option %c, rest ignored\n", (int)s[-1]);
	  while (*++s && !isspace(*s));
	  break;
	}
      }
      break;
#pragma mark option_G
    case 'G':
      isgeom= True;
      qh_appendprint (qh_PRINTgeom);
      while (*s && !isspace(*s)) {
	switch(*s++) {
        case 'a':
          qh_option ("Gall-points", NULL, NULL);
          qh PRINTdots= True;
          break;
        case 'c':
          qh_option ("Gcentrums", NULL, NULL);
          qh PRINTcentrums= True;
          break;
	case 'h':
          qh_option ("Gintersections", NULL, NULL);
	  qh DOintersections= True;
	  break;
	case 'i':
          qh_option ("Ginner", NULL, NULL);
	  qh PRINTinner= True;
	  break;
	case 'n':
          qh_option ("Gno-planes", NULL, NULL);
	  qh PRINTnoplanes= True;
	  break;
	case 'o':
          qh_option ("Gouter", NULL, NULL);
	  qh PRINTouter= True;
	  break;
	case 'p':
          qh_option ("Gcoplanar", NULL, NULL);
	  qh PRINTcoplanar= True;
	  break;
	case 'r':
          qh_option ("Gridges", NULL, NULL);
	  qh PRINTridges= True;
	  break;
	case 'v':
          qh_option ("Gvertices", NULL, NULL);
	  qh PRINTspheres= True;
	  break;
	case 'D':
	  if (!isdigit (*s))
	    fprintf (qh ferr, "qhull input error: missing dimension for 'GD' option\n");
	  else {
	    if (qh DROPdim >= 0)
	      fprintf (qh ferr, "qhull warning: can only drop one dimension.  Previous 'GD%d' ignored\n",
	           qh DROPdim);
  	    qh DROPdim= qh_strtol (s, &s);
            qh_option ("GDrop-dim", &qh DROPdim, NULL);
          }
	  break;
	default:
	  s--;
	  fprintf (qh ferr, "qhull warning: unknown 'G' print option %c, rest ignored\n", (int)s[0]);
	  while (*++s && !isspace(*s));
	  break;
	}
      }
      break;
#pragma mark option_P
    case 'P':
      while (*s && !isspace(*s)) {
	switch(*s++) {
	case 'd': case 'D':
	  key= s[-1];
	  i= qh_strtol (s, &s);
	  r= 0;
	  if (*s == ':') 
	    r= qh_strtod (++s, &s);
	  if (key == 'd')
  	    qh_option ("Pdrop-facets-dim-less", &i, &r);
  	  else
  	    qh_option ("PDrop-facets-dim-more", &i, &r);
	  break;
        case 'g':
          qh_option ("Pgood-facets", NULL, NULL);
          qh PRINTgood= True;
          break;
        case 'G':
          qh_option ("PGood-facet-neighbors", NULL, NULL);
          qh PRINTneighbors= True;
          break;
        case 'o':
          qh_option ("Poutput-forced", NULL, NULL);
          qh FORCEoutput= True;
          break;
        case 'p':
          qh_option ("Pprecision-ignore", NULL, NULL);
          qh PRINTprecision= False;
          break;
	case 'A':
	  if (!isdigit (*s))
	    fprintf (qh ferr, "qhull input error: missing count for 'PA' option\n");
	  else {
  	    qh KEEParea= qh_strtol (s, &s);
            qh_option ("PArea-keep", &qh KEEParea, NULL);
            qh GETarea= True;
          }
	  break;
	case 'F':
	  if (!isdigit (*s))
	    fprintf (qh ferr, "qhull input error: missing area for 'PF' option\n");
	  else {
  	    qh KEEPminArea= qh_strtod (s, &s);
            qh_option ("PFacet-area-keep", NULL, &qh KEEPminArea);
            qh GETarea= True;
          }
	  break;
	case 'M':
	  if (!isdigit (*s))
	    fprintf (qh ferr, "qhull input error: missing count for 'PM' option\n");
	  else {
  	    qh KEEPmerge= qh_strtol (s, &s);
            qh_option ("PMerge-keep", &qh KEEPmerge, NULL);
          }
	  break;
	default:
	  s--;
	  fprintf (qh ferr, "qhull warning: unknown 'P' print option %c, rest ignored\n", (int)s[-1]);
	  while (*++s && !isspace(*s));
	  break;
	}
      }
      break;
#pragma mark option_Q
    case 'Q':
      lastproject= -1;
      while (*s && !isspace(*s)) {
	switch(*s++) {
	case 'b': case 'B':  /* handled by qh_initthresholds */
	  key= s[-1];
	  if (key == 'b' && *s == 'B') {
	    s++;
	    r= qh_DEFAULTbox;
	    qh_option ("QbBound-unit-box", NULL, &r);
	    break;
	  }
	  k= qh_strtol (s, &s);
	  r= 0.0;
	  wasproject= False;
	  if (*s == ':' && (r= qh_strtod(++s, &s)) == 0.0) {
	    t= s;            /* need true dimension for memory allocation */
	    while (*t && !isspace(*t)) {
	      if (toupper(*t++) == 'B' 
	       && k == qh_strtol (t, &t)
	       && *t++ == ':'
	       && qh_strtod(t, &t) == 0.0) {
	        qh PROJECTinput++;
	        trace2((qh ferr, "qh_initflags: project dimension %d\n", k));
	        qh_option ("Qb-project-dim", &k, NULL);
		wasproject= True;
	        lastproject= k;
	        break;
	      }
	    }
  	  }
	  if (!wasproject) {
	    if (lastproject == k && r == 0.0) 
	      lastproject= -1;  /* doesn't catch all possible sequences */
	    else if (key == 'b') {
	      qh SCALEinput= True;
	      if (r == 0.0)
		r= -qh_DEFAULTbox;
	      qh_option ("Qbound-dim-low", &k, &r);
	    }else {
	      qh SCALEinput= True;
	      if (r == 0.0)
		r= qh_DEFAULTbox;
	      qh_option ("QBound-dim-high", &k, &r);
	    }
	  }
	  break;
	case 'c':
	  qh_option ("Qcoplanar-keep", NULL, NULL);
	  qh KEEPcoplanar= True;
	  break;
	case 'f':
	  qh_option ("Qfurthest-outside", NULL, NULL);
	  qh BESToutside= True;
	  break;
	case 'g':
	  qh_option ("Qgood-facets-only", NULL, NULL);
	  qh ONLYgood= True;
	  break;
	case 'i':
	  qh_option ("Qinside-keep", NULL, NULL);
	  qh KEEPinside= True;
	  break;
	case 'm':
	  qh_option ("Qmax-outside-only", NULL, NULL);
	  qh ONLYmax= True;
	  break;
	case 'r':
	  qh_option ("Qrandom-outside", NULL, NULL);
	  qh RANDOMoutside= True;
	  break;
	case 's':
	  qh_option ("Qsearch-initial-simplex", NULL, NULL);
	  qh ALLpoints= True;
	  break;
	case 'u':
	  qh_option ("QupperHull", NULL, NULL);
	  qh ATinfinity= False;
	  break;
	case 'v':
	  qh_option ("Qvertex-neighbors-convex", NULL, NULL);
	  qh TESTvneighbors= True;
	  break;
	case 'x':
	  qh_option ("Qxact-merge", NULL, NULL);
	  qh MERGEexact= True;
	  qh MERGING= True;
	  break;
	case '1':
	  qh_option ("Q1-no-angle-sort", NULL, NULL);
	  qh ANGLEmerge= False;
	  goto LABELcheckdigit;
	case '2':
	  qh_option ("Q2-no-merge-independent", NULL, NULL);
	  qh MERGEindependent= False;
	LABELcheckdigit:
	  if (isdigit(*s)) 
	    fprintf (qh ferr, "qhull warning: can not follow '1' or '2' with a digit.  '%c' skipped.\n",
	             *s++);
	  break;
	case '3':
	  qh_option ("Q3-no-merge-vertices", NULL, NULL);
	  qh MERGEvertices= False;
	  break;
	case '4':
	  qh_option ("Q4-avoid-old-into-new", NULL, NULL);
	  qh AVOIDold= True;
	  break;
	case '5':
	  qh_option ("Q5-no-check-outer", NULL, NULL);
	  qh SKIPcheckmax= True;
	  break;
	case '6':
	  qh_option ("Q6-no-concave-merge", NULL, NULL);
	  qh SKIPconvex= True;
	  break;
	case '7':
	  qh_option ("Q7-no-breadth-first", NULL, NULL);
	  qh VIRTUALmemory= True;
	  break;
	case '8':
	  qh_option ("Q8-no-near-inside", NULL, NULL);
	  qh NOnearinside= True;
	  break;
	case 'G':
	  i= qh_strtol (s, &t);
	  if (qh GOODpoint) 
	    fprintf (qh ferr, "qhull warning: good point already defined for QGn.  Ignored\n");
          else if (s == t)
	    fprintf (qh ferr, "qhull warning: no good point id given for option QGn.  Ignored\n");
	  else if (i < 0 || *s == '-') { 
 	    qh GOODpoint= i-1;
  	    qh_option ("QGood-if-dont-see-point", &i, NULL);
	  }else {
 	    qh GOODpoint= i+1;
  	    qh_option ("QGood-if-see-point", &i, NULL);
  	  }
 	  s= t;
	  break;
	case 'R':
          if (!isdigit(*s) && *s != '-')
	    fprintf (qh ferr, "qhull warning: missing random seed for option QRn.  Ignored\n");
	  else {
 	    qh ROTATErandom= i= qh_strtol(s, &s);
 	    if (i == -1)
   	      qh_option ("QRandom-seed", NULL, NULL );
   	    else if (i > 0)
   	      qh_option ("QRotate-id", &i, NULL );
          }
	  break;
	case 'V':
	  i= qh_strtol (s, &t);
	  if (qh GOODvertex) 
	    fprintf (qh ferr, "qhull warning: good vertex already defined for QV.  Ignored\n");
          else if (s == t)
	    fprintf (qh ferr, "qhull warning: no good point id given for QV.  Ignored\n");
	  else if (i < 0) {
 	    qh GOODvertex= i - 1;
 	    qh_option ("QV-good-facets-not-point", &i, NULL);
	  }else {
  	    qh_option ("QV-good-facets-point", &i, NULL);
	    qh GOODvertex= i + 1;
          }
 	  s= t;
	  break;
	default:
	  fprintf (qh ferr, "qhull warning: unknown 'Q' qhull option %c, rest ignored\n", (int)s[-1]);
	  while (*++s && !isspace(*s));
	  break;
	}
      }
      break;
#pragma mark option_T
    case 'T':
      while (*s && !isspace(*s)) {
	if (isdigit(*s) || *s == '-')
	  qh IStracing= qh_strtol(s, &s);
	else switch(*s++) {
	case 'c':
          qh_option ("Tcheck-frequently", NULL, NULL);
	  qh CHECKfrequently= True;
	  break;
	case 's':
          qh_option ("Tstatistics", NULL, NULL);
	  qh PRINTstatistics= True;
	  break;
	case 'v':
          qh_option ("Tverify", NULL, NULL);
	  qh VERIFYoutput= True;
	  break;
	case 'z':
          qh_option ("Tz-stdout", NULL, NULL);
	  qh ferr= qh fout;
	  qhmem.ferr= qh fout;
	  break;
	case 'C':
	  if (!isdigit(*s))
	    fprintf (qh ferr, "qhull warning: no point given for trace option C.  Ignored\n");
	  else {
	    i= qh_strtol (s, &s);
	    qh_option ("TCone-stop", &i, NULL);
	    qh STOPcone= i + 1;
          }
	  break;
	case 'F':
	  if (!isdigit(*s))
	    fprintf (qh ferr, "qhull warning: no count of new facets for trace option P.  Ignored\n");
	  else {
	    qh REPORTfreq= qh_strtol (s, &s);
            qh_option ("TFacet-log", &qh REPORTfreq, NULL);
	    qh REPORTfreq2= qh REPORTfreq/2;  /* for tracemerging() */
	  }
	  break;
	case 'P':
	  if (!isdigit(*s))
	    fprintf (qh ferr, "qhull warning: no point given for trace option P.  Ignored\n");
	  else {
	    qh TRACEpoint= qh_strtol (s, &s);
            qh_option ("Trace-point", &qh TRACEpoint, NULL);
          }
	  break;
	case 'M':
	  if (!isdigit(*s))
	    fprintf (qh ferr, "qhull warning: no merge given for trace option M.  Ignored\n");
	  else {
	    qh TRACEmerge= qh_strtol (s, &s);
            qh_option ("Trace-merge", &qh TRACEmerge, NULL);
          }
	  break;
	case 'V':
	  i= qh_strtol (s, &t);
	  if (s == t)
	    fprintf (qh ferr, "qhull warning: no point given for trace option V.  Ignored\n");
	  else if (i < 0) {
	    qh STOPpoint= i - 1;
            qh_option ("TV-stop-before-point", &i, NULL);
	  }else {
	    qh STOPpoint= i + 1;
            qh_option ("TV-stop-after-point", &i, NULL);
          }
          s= t;
	  break;
	case 'W':
	  if (!isdigit(*s))
	    fprintf (qh ferr, "qhull warning: no max width given for trace option D.  Ignored\n");
	  else {
 	    qh TRACEdist= (realT) qh_strtod (s, &s);
            qh_option ("TWide-trace", NULL, &qh TRACEdist);
          }
	  break;
	default:
	  fprintf (qh ferr, "qhull warning: unknown 'T' trace option %c, rest ignored\n", (int)s[-1]);
	  while (*++s && !isspace(*s));
	  break;
	}
      }
      break;
    default:
      fprintf (qh ferr, "qhull warning: unknown flag %c (%x)\n", (int)s[-1],
	       (int)s[-1]);
      break;
    }
    if (s-1 == prev_s && *s && !isspace(*s)) {
      fprintf (qh ferr, "qhull warning: missing space after flag %c (%x); reserved for menu. Skipped.\n",
	       (int)*prev_s, (int)*prev_s);
      while (*s && !isspace(*s))
	s++;
    }
  }
  if (isgeom && !qh FORCEoutput && qh PRINTout[1])
    fprintf (qh ferr, "qhull warning: additional output formats are not compatible with Geomview\n");
  /* set derived values in qh_initqhull_globals */
} /* initflags */


/*-------------------------------------------
-initqhull_buffers- initialize global memory buffers
  must match freebuffers()
*/
void qh_initqhull_buffers (void) {
  int k;

  qh TEMPsize= (qhmem.LASTsize - sizeof (setT))/SETelemsize;
  if (qh TEMPsize <= 0 || qh TEMPsize > qhmem.LASTsize)
    qh TEMPsize= 8;  /* e.g., if qh_NOmem */
  qh other_points= qh_setnew (qh TEMPsize);
  qh del_vertices= qh_setnew (qh TEMPsize);
  qh searchset= qh_setnew (qh TEMPsize);
  qh NEARzero= (realT *)qh_memalloc(qh hull_dim * sizeof(realT));
  qh lower_threshold= (realT *)qh_memalloc((qh input_dim+1) * sizeof(realT));
  qh upper_threshold= (realT *)qh_memalloc((qh input_dim+1) * sizeof(realT));
  qh lower_bound= (realT *)qh_memalloc((qh input_dim+1) * sizeof(realT));
  qh upper_bound= (realT *)qh_memalloc((qh input_dim+1) * sizeof(realT));
  for(k= qh input_dim+1; k--; ) {
    qh lower_threshold[k]= -REALmax;
    qh upper_threshold[k]= REALmax;
    qh lower_bound[k]= -REALmax;
    qh upper_bound[k]= REALmax;
  }
  qh gm_matrix= (coordT *)qh_memalloc((qh hull_dim+1) * qh hull_dim * sizeof(coordT));
  qh gm_row= (coordT **)qh_memalloc((qh hull_dim+1) * sizeof(coordT *));
} /* initqhull_buffers */

/*---------------------------------------------
-initqhull_globals- initialize globals
  ismalloc set if points were malloc'd and qhull should free at end
returns:
  sets qh first_point, num_points, input_dim, hull_dim and others
  modifies hull_dim if ((DELAUNAY and PROJECTdelaunay) or PROJECTinput)
  seeds random number generator (seed=1 if tracing)
  adjust user flags as needed
  also checks hull_dim dependencies and constants
*/
void qh_initqhull_globals (coordT *points, int numpoints, int dim, boolT ismalloc) {
  int seed, pointsneeded, extra= 0, i, randi, k;
  boolT printgeom= False, printmath= False;
  realT randr;
  realT factorial;
  
  time_t timedata;

  trace0((qh ferr, "qh_initqhull_globals: for %s | %s\n", qh rbox_command, 
      qh qhull_command));
  qh POINTSmalloc= ismalloc;
  qh first_point= points;
  qh num_points= numpoints;
  qh hull_dim= qh input_dim= dim;
#ifdef qh_NOmerge
  if (qh MERGING) {
    fprintf (qh ferr, "qhull input error: merging not installed (qh_NOmerge + 'Qx', 'Cn' or 'An')\n");
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
#endif
  if (qh MERGING && !qh POSTmerge && qh premerge_cos > REALmax/2
  && (qh premerge_centrum == 0 || qh premerge_centrum > REALmax/2)) {
    qh ZEROcentrum= True;
    qh ZEROall_ok= True;
  }
  if ((qh KEEParea || qh KEEPminArea < REALmax/2 || qh KEEPmerge || qh DELAUNAY)
   && !(qh PRINTgood || qh PRINTneighbors)) {
    qh PRINTgood= True;
    qh_option ("Pgood", NULL, NULL);
  }
  if (qh DELAUNAY && qh HALFspace) {
    fprintf (qh ferr, "qhull input error: can not use Delaunay ('d') or Voronoi ('v') with halfspace intersection ('H')\n");
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  if (!qh DELAUNAY && !qh ATinfinity) {
    fprintf (qh ferr, "qhull input error: use upper-hull ('Qu') with Delaunay ('d') or Voronoi ('v')\n");
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  qh DOcheckmax= (!qh FORCEoutput && !qh SKIPcheckmax && (qh MERGING || qh APPROXhull));
  qh KEEPnearinside= (qh DOcheckmax && !(qh KEEPinside && qh KEEPcoplanar) 
                          && !qh NOnearinside);
  if (qh MERGING)
    qh CENTERtype= qh_AScentrum;
  else if (qh VORONOI)
    qh CENTERtype= qh_ASvoronoi;
  if (qh TESTvneighbors && !qh MERGING) {
    fprintf(qh ferr, "qhull input error: test vertex neighbors ('Qv') needs a merge option\n");
    qh_errexit (qh_ERRinput, NULL ,NULL);
  }
  if (qh PROJECTinput || (qh DELAUNAY && qh PROJECTdelaunay)) {
    qh hull_dim -= qh PROJECTinput;
    if (qh DELAUNAY) {
      qh hull_dim++;
      extra= 1;
    }
  }
  if (qh hull_dim <= 1) {
    fprintf(qh ferr, "qhull error: dimension %d must be > 1\n", qh hull_dim);
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  for (k= 2, factorial=1.0; k < qh hull_dim; k++)
    factorial *= k;
  qh AREAfactor= 1.0 / factorial;
  trace2((qh ferr, "qh_initqhull_globals: initialize globals.  dim %d numpoints %d malloc? %d projected %d to hull_dim %d\n",
	dim, numpoints, ismalloc, qh PROJECTinput, qh hull_dim));
  qh normal_size= qh hull_dim * sizeof(coordT);
  qh center_size= qh normal_size - sizeof(coordT);
  pointsneeded= qh hull_dim+1;
  if (qh hull_dim > qh_DIMmergeVertex) {
    qh MERGEvertices= False;
    qh_option ("Q3-no-merge-vertices-dim-high", NULL, NULL);
  }
  if (qh GOODpoint)
    pointsneeded++;
  if (qh GOODpoint > 0) 
    qh GOODpointp= qh_point (qh GOODpoint-1);
  else if (qh GOODpoint < 0) 
    qh GOODpointp= qh_point (-qh GOODpoint-1);
  if (qh GOODvertex > 0)
    qh GOODvertexp= qh_point (qh GOODvertex-1);
  else if (qh GOODvertex < 0) 
    qh GOODvertexp= qh_point (-qh GOODvertex-1);
  if ((qh GOODpoint  
       && (qh GOODpointp < qh first_point  /* also catches !GOODpointp */
	   || qh GOODpointp > qh_point (qh num_points-1)))
    || (qh GOODvertex
	&& (qh GOODvertexp < qh first_point  /* also catches !GOODvertexp */
	    || qh GOODvertexp > qh_point (qh num_points-1)))) {
    fprintf (qh ferr, "qhull input error: either QGn or QVn point is > p%d\n",
	     qh num_points-1);
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  if (qh TRACEpoint != -1 || qh TRACEdist < REALmax/2 || qh TRACEmerge) {
    qh TRACElevel= (qh IStracing? qh IStracing : 3);
    qh IStracing= 0;
  }
#ifdef qh_NOtrace
  if (qh IStracing || qh TRACElevel) {
    fprintf (qh ferr, "qhull input error: tracing is not installed (qh_NOtrace in user.h)");
    qh_errexit (qh_ERRqhull, NULL, NULL);
  } 
#endif
  if (qh ROTATErandom == 0 || qh ROTATErandom == -1) {
    seed= time (&timedata);
    if (!qh ROTATErandom)
      qh_option ("QRotate-random", &seed, NULL);
    qh ROTATErandom= seed;
  }else if (qh ROTATErandom > 0)
    seed= qh ROTATErandom;
  else
    seed= 1;
  qh_RANDOMseed_(seed);
  randr= 0.0;
  for (i= 1000; i--; ) {
    randi= qh_RANDOMint;
    randr += randi;
    if (randi > qh_RANDOMmax) {
      fprintf (qh ferr, "\
qhull configuration error (qh_RANDOMmax in user.h):\n\
   random integer %d > qh_RANDOMmax (%.8g)\n",
	       randi, qh_RANDOMmax);
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
  }
  if (randr/1000 < qh_RANDOMmax/10)
    fprintf (qh ferr, "\
qhull configuration warning (qh_RANDOMmax in user.h):\n\
   average of 1000 randoms %.2g much less than expected (%.2g).\n\
   Is qh_RANDOMmax wrong?\n",
	     randr/1000, qh_RANDOMmax/2);
  qh RANDOMa= 2.0 * qh RANDOMfactor/qh_RANDOMmax;
  qh RANDOMb= 1.0 - qh RANDOMfactor;
  if (qh_HASHfactor < 1.1) {
    fprintf(qh ferr, "qhull internal error (qh_initqhull_globals): qh_HASHfactor %d must be at least 1.1.  Qhull uses linear hash probing\n",
      qh_HASHfactor);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  if (numpoints+extra < pointsneeded) {
    fprintf(qh ferr,"qhull input error: not enough points (%d) to construct initial simplex (need %d)\n",
	    numpoints, pointsneeded);
    qh_errexit(qh_ERRinput, NULL, NULL);
  }
  for (i= qh_PRINTEND; i--; ) {
    if (qh PRINTout[i] == qh_PRINTgeom)
      printgeom= True;
    else if (qh PRINTout[i] == qh_PRINTmathematica)
      printmath= True;
    else if (qh PRINTout[i] == qh_PRINTpointintersect && !qh HALFspace) {
      fprintf (qh ferr, "qhull input error: option 'Fp' is used for halfspace intersect ('Hn,n,n')\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }else if (qh PRINTout[i] == qh_PRINTpointnearest && !qh KEEPcoplanar && !qh KEEPinside) 
      fprintf (qh ferr, "qhull input warning: option 'FP' needs option 'Qc' or 'Qi' to record the points\n");
    else if (qh PRINTout[i] == qh_PRINTcentrums && qh VORONOI) {
      fprintf (qh ferr, "qhull input error: option 'Fc' is not available for Voronoi vertices ('v')\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
  }
  if (printmath && (qh hull_dim > 3 || qh VORONOI || qh HALFspace)) {
    fprintf (qh ferr, "qhull input error: Mathematica output is only available for 2-d and 3-d convex hulls and Delaunay triangulations\n");
    qh_errexit (qh_ERRinput, NULL, NULL);
  }
  if (printgeom) {
    if (qh hull_dim > 4) {
      fprintf (qh ferr, "qhull input error: Geomview output is only available for 2-d, 3-d and 4-d\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
    if (qh PRINTnoplanes && !(qh PRINTcoplanar + qh PRINTcentrums
     + qh PRINTdots + qh PRINTspheres + qh DOintersections + qh PRINTridges)) {
      fprintf (qh ferr, "qhull input error: no output specified for Geomview\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
    if (qh VORONOI && (qh hull_dim > 3 || qh DROPdim >= 0)) {
      fprintf (qh ferr, "qhull input error: Geomview output for Voronoi diagrams only for 2-d\n");
      qh_errexit (qh_ERRinput, NULL, NULL);
    }
    if (qh hull_dim == 4 && qh DROPdim == -1 &&
	(qh PRINTcoplanar || qh PRINTspheres || qh PRINTcentrums)) {
      fprintf (qh ferr, "qhull input warning: coplanars, vertices, and centrums output not\n\
available for 4-d output (ignored).  Could use 'GDn' instead.\n");
      qh PRINTcoplanar= qh PRINTspheres= qh PRINTcentrums= False;
    }
  }
  qh PRINTdim= qh hull_dim;
  if (qh DROPdim >=0) {    /* after Geomview checks */
    if (qh DROPdim < qh hull_dim) {
      qh PRINTdim--;
      if (!printgeom || qh hull_dim < 3) 
        fprintf (qh ferr, "qhull input warning: drop dimension 'GD%d' is only available for 3-d/4-d Geomview\n", qh DROPdim);
    }else
      qh DROPdim= -1;
  }else if (qh VORONOI) {
    qh DROPdim= qh hull_dim-1;
    qh PRINTdim= qh hull_dim-1; 
  }
} /* initqhull_globals */
 
/*-----------------------------------------------------
-initqhull_mem- initialize mem.c for qhull
  qh.hull_dim and normal_size determines some of the allocation sizes
  if qh MERGING, includes ridgeT
returns:
  mem.c already for memalloc/memfree (errors if called beforehand)
notes:
  the user can add up to 10 additional sizes for quick allocation (increase numsizes)
  print out memsizes in qh_produceoutput
*/
void qh_initqhull_mem (void) {
  int numsizes;
  int i;

  numsizes= 8+10;
  qh_meminitbuffers (qh IStracing, qh_MEMalign, numsizes, 
                     qh_MEMbufsize,qh_MEMinitbuf);
  qh_memsize(sizeof(vertexT));
  if (qh MERGING) {
    qh_memsize(sizeof(ridgeT));
    qh_memsize(sizeof(mergeT));
  }
  qh_memsize(sizeof(facetT));
  qh_memsize(sizeof(hashentryT));
  i= sizeof(setT) + (qh hull_dim - 1) * SETelemsize;  /* ridge.vertices */
  qh_memsize(i);
  qh_memsize(qh normal_size);        /* normal */
  i += SETelemsize;                 /* facet.vertices, .ridges, .neighbors */
  qh_memsize(i);
  qh_user_memsizes();
  qh_memsetup();
} /* initqhull_mem */

/*-------------------------------------------
-initqhull_start -- start initialization of qhull
  inits statistics
*/ 
void qh_initqhull_start (FILE *infile, FILE *outfile, FILE *errfile) {

  clock(); /* start the clock */
#if qh_QHpointer
  if (!(qh_qh= (qhT *)malloc (sizeof(qhT)))) {
    fprintf (errfile, "qhull error (qh_initqhull_globals): insufficient memory\n");
    exit (qh_ERRmem);  /* no error handler */
  }
  memset((char *)qh_qh, 0, sizeof(qhT));   /* every field is 0, FALSE, NULL */
#else
  memset((char *)&qh_qh, 0, sizeof(qhT));
#endif
  strcat (qh qhull, "qhull");
  qh_initstatistics();
  qh ANGLEmerge= True;
  qh ATinfinity= True;
  qh DROPdim= -1;
  qh ferr= errfile;
  qh fin= infile;
  qh fout= outfile;
  qh furthest_id= -1;
  qh KEEPminArea = REALmax;
  qh lastreport= INT_MIN;
  qh mergereport= INT_MIN;
  qh max_outside= 0.0;
  qh maxmaxcoord= 0.0;
  qh max_vertex= 0.0;
  qh MERGEindependent= True;
  qh mergereport= INT_MAX;
  qh min_vertex= 0.0;
  qh MINdenom_1= fmax_(1.0/REALmax, REALmin);
  qh MINoutside= 0.0;
  qh MINvisible= REALmax;
  qh MAXcoplanar= REALmax;
  qh premerge_centrum= 0.0;
  qh premerge_cos= REALmax;
  qh PRINTprecision= True;
  qh PRINTradius= 0.0;
  qh postmerge_cos= REALmax;
  qh postmerge_centrum= 0.0;
  qh rand_seed= 1;
  qh ROTATErandom= INT_MIN;
  qh MERGEvertices= True;
  qh totarea= 0.0;
  qh TRACEdist= REALmax;
  qh TRACEpoint= -1;
  qh tracefacet_id= -1;  /* recompile to trace an id */
  qh tracevertex_id= -1;
} /* initqhull_start */

/*---------------------------------------------
-initthresholds	set thresholds for printing and scaling from command line
  see 'prompt' in unix.c for documentation
  see also initflags(), 'Qbk' 'QBk' 'Pdk' and 'PDk'
  sets qh GOODthreshold or qh SPLITthreshold if 'Pd0D1' used
*/
void qh_initthresholds(char *command) {
  realT value;
  int index, maxdim, k;
  char *s= command;
  char key;
  
  maxdim= qh input_dim;
  if (qh DELAUNAY && (qh PROJECTdelaunay || qh PROJECTinput))
    maxdim++;
  while (*s) {
    if (*s == '-')
      s++;
    if (*s == 'P') {
      s++;
      while (*s && !isspace(key= *s++)) {
	if (key == 'd' || key == 'D') {
	  if (!isdigit(*s)) {
	    fprintf(qh ferr, "qhull warning: no dimension given for Print option '%c' at: %s.  Ignored\n",
		    key, s-1);
	    continue;
	  }
	  index= qh_strtol (s, &s);
	  if (index >= qh hull_dim) {
	    fprintf(qh ferr, "qhull warning: dimension %d for Print option '%c' is >= %d.  Ignored\n", 
	        index, key, qh hull_dim);
	    continue;
	  }
	  if (*s == ':') {
	    value= qh_strtod(++s, &s);
	    if (fabs((double)value) > 1.0) {
	      fprintf(qh ferr, "qhull warning: value %2.4g for Print option %c is > +1 or < -1.  Ignored\n", 
	              value, key);
	      continue;
	    }
	  }else
	    value= 0.0;
	  if (key == 'd')
	    qh lower_threshold[index]= value;
	  else
	    qh upper_threshold[index]= value;
	}
      }
    }else if (*s == 'Q') {
      s++;
      while (*s && !isspace(key= *s++)) {
	if (key == 'b' && *s == 'B') {
	  s++;
	  for (k=maxdim; k--; ) {
	    qh lower_bound[k]= -qh_DEFAULTbox;
	    qh upper_bound[k]= qh_DEFAULTbox;
	  }
	}else if (key == 'b' || key == 'B') {
	  if (!isdigit(*s)) {
	    fprintf(qh ferr, "qhull warning: no dimension given for Qhull option %c.  Ignored\n",
		    key);
	    continue;
	  }
	  index= qh_strtol (s, &s);
	  if (index >= maxdim) {
	    fprintf(qh ferr, "qhull warning: dimension %d for Qhull option %c is >= %d.  Ignored\n", 
	        index, key, maxdim);
	    continue;
	  }
	  if (*s == ':')
	    value= qh_strtod(++s, &s);
	  else if (key == 'b')
	    value= -qh_DEFAULTbox;
	  else
	    value= qh_DEFAULTbox;
	  if (key == 'b')
	    qh lower_bound[index]= value;
	  else
	    qh upper_bound[index]= value;
	}
      }
    }else {
      while (*s && !isspace (*s))
        s++;
    }
    while (isspace (*s))
      s++;
  }
  for (k= qh hull_dim; k--; ) {
    if (qh lower_threshold[k] > -REALmax/2) {
      qh GOODthreshold= True;
      if (qh upper_threshold[k] < REALmax/2) {
        qh SPLITthresholds= True;
        qh GOODthreshold= False;
        break;
      }
    }else if (qh upper_threshold[k] < REALmax/2)
      qh GOODthreshold= True;
  }
} /* initthresholds */

/*------------------------------------------
-option- add an option description to qh qhull_options
  will be printed with statistics ('Ts') and errors
  strlen(option) < 40
*/
void qh_option (char *option, int *i, realT *r) {
  char buf[200];
  int len, maxlen;

  sprintf (buf, "  %s", option);
  if (i)
    sprintf (buf+strlen(buf), " %d", *i);
  if (r)
    sprintf (buf+strlen(buf), " %2.2g", *r);
  len= strlen(buf);
  qh qhull_optionslen += len;
  maxlen= sizeof (qh qhull_options) - len -1;
  maximize_(maxlen, 0);
  if (qh qhull_optionslen > 80 && maxlen > 0) {
    qh qhull_optionslen= len;
    strncat (qh qhull_options, "\n", maxlen--);
  }    
  strncat (qh qhull_options, buf, maxlen);
} /* option */

#if qh_QHpointer
/*------------------------------------------
-restore_qhull- restores a previously saved qhull
  also restores qh_qhstat and qhmem.tempstack
  errors if current qhull hasn't been saved or freed
  uses qhmem for error reporting
*/
void qh_restore_qhull (qhT **oldqh) {

  if (*oldqh && strcmp ((*oldqh)->qhull, "qhull")) {
    fprintf (qhmem.ferr, "qhull internal error (qh_restore_qhull): %p is not a qhull data structure\n",
                  *oldqh);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  if (qh_qh) {
    fprintf (qhmem.ferr, "qhull internal error (qh_restore_qhull): did not save or free existing qhull\n");
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  if (!*oldqh || !(*oldqh)->old_qhstat) {
    fprintf (qhmem.ferr, "qhull internal error (qh_restore_qhull): did not previously save qhull %p\n",
                  *oldqh);
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  qh_qh= *oldqh;
  *oldqh= NULL;
  qh_qhstat= qh old_qhstat;
  qhmem.tempstack= qh old_tempstack;
  trace1((qh ferr, "qh_restore_qhull: restored qhull from %p\n", *oldqh));
} /* restore_qhull */

/*------------------------------------------
-save_qhull- saves qhull for a later qh_restore_qhull
  also saves qh_qhstat and qhmem.tempstack
returns:
  qhull for a later restore_qhull
  qh_qh=NULL
notes:
  need to initialize qhull or call qh_restore_qhull before continuing
*/
qhT *qh_save_qhull (void) {
  qhT *oldqh;

  if (!qh_qh) {
    fprintf (qhmem.ferr, "qhull internal error (qh_save_qhull): qhull not initialized\n");
    qh_errexit (qh_ERRqhull, NULL, NULL);
  }
  qh old_qhstat= qh_qhstat;
  qh_qhstat= NULL;
  qh old_tempstack= qhmem.tempstack;
  qhmem.tempstack= NULL;
  oldqh= qh_qh;
  qh_qh= NULL;
  trace1((qhmem.ferr, "qh_save_qhull: saved qhull %p\n", oldqh));
  return oldqh;
} /* save_qhull */

#endif

/*-----------------------------------------
-strtol/tod -- internal versions that don't skip trailing spaces
*/
double qh_strtod (const char *s, char **endp) {
  double result;

  result= strtod (s, endp);
  if (s < (*endp) && (*endp)[-1] == ' ')
    (*endp)--;
  return result;
} /* strtod */

int qh_strtol (const char *s, char **endp) {
  int result;

  result= (int) strtol (s, endp, 10);
  if (s< (*endp) && (*endp)[-1] == ' ')
    (*endp)--;
  return result;
} /* strtol */
