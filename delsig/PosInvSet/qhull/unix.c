/* unix.c -- Unix version of qhull

   see README
   
   copyright (c) 1993-1995, The Geometry Center
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "qhull.h"
#if __MWERKS__
#include <SIOUX.h>
#include <Files.h>
#include <console.h>
#include <Desk.h>
#endif

char qh_version[]= "version 2.2 95/12/4"; /* Changes Announce README man make*/
/* also change concise prompt below */
char qh_prompt[]= "\n\
qhull- compute the convex hull of a set of n-d points.\n\
    <http://www.geom.umn.edu/locate/qhull>  %s\n\
input (stdin):\n\
    1st lines= dimension #points (or vice-versa).   rest= point coordinates\n\
    comments ignored (non-numeric to end-of-line).\n\
    for half-space intersection: Use dim plus one and put offsets after\n\
      coefficients.  May be proceeded by a single feasible point\n\
\n\
options:\n\
    d    - Delaunay triangulation by lifting points to a paraboloid\n\
    v    - Voronoi vertices via the Delaunay triangulation\n\
    Hn,n,n,... - halfspace intersection about origin [n,n,n,...]\n\
    Qopts- Qhull control options:\n\
           Qbk:n   - scale point[k] to low bound. 'bB' gives cube.\n\
             QBk:n - upper bound (Bk is %2.2g). 'bk:0Bk:0' removes k-th coord.\n\
           QRn     - random rotation (n=seed, n=0 time, n=-1 time/no rotate)\n\
           Qc      - keep coplanar points with nearest facet\n\
           Qf      - partition point to furthest outside facet\n\
           Qg      - only build good facets (needs 'QGn', 'QVn', or 'PdD')\n\
           Qm      - only process points that would increase max_outside\n\
           Qi      - keep inside points with nearest facet\n\
           Qr      - process random outside points instead of furthest ones\n\
           Qs      - search all points for the initial simplex\n\
           Qu      - for 'd', compute upper hull without point at-infinity\n\
           Qv      - test vertex neighbors for convexity\n\
           Qx      - exact pre-merges (skips coplanar and angle coplanar facets)\n\
           QGn     - good facet if visible from point n, -n for not visible\n\
           QVn     - good facet if it includes point n, -n if not\n\
           Q1	   - sort merges by type instead of angle\n\
           Q2      - merge all non-convex at once instead of independent sets\n\
           Q3      - do not merge redundant vertices\n\
	   Q4      - avoid old->new merges\n\
           Q5      - do not correct outer planes at end of qhull\n\
           Q6      - do not pre-merge concave or coplanar facets\n\
           Q7      - depth-first processing instead of breadth-first\n\
           Q8      - do not process near-inside points\n\
    Topts- Trace options:\n\
           T4      - trace at level n, 4=all, 5=mem/gauss, -1= events\n\
           Tc      - check frequently during execution\n\
           Ts      - print statistics\n\
           Tv      - verify result: structure, convexity, and point inclusion\n\
           Tz      - send all output to stdout\n\
           TFn     - report summary when n or more facets created\n\
           TPn     - turn on tracing when point n added to hull\n\
            TMn    - turn on tracing at merge n\n\
            TWn    - trace merge facets when width > n\n\
           TVn     - stop qhull after adding point n, -n for before (see TCn)\n\
            TCn    - stop qhull after building cone for point n (see TVn)\n\
\n\
Precision options (default detects precision problems without correction):\n\
    Cn   - radius of centrum (roundoff added).  Merge facets if non-convex\n\
     An  - cosine of maximum angle.  Merge facets if cosine > n or non-convex\n\
           C-0 roundoff, A-0.99/C-0.01 pre-merge, A0.99/C0.01 post-merge\n\
    En   - max roundoff error for distance computation\n\
    Rn   - randomly perturb computations by a factor of [1-n,1+n]\n\
    Vn   - min distance above plane for a visible facet (default 3C-n or En)\n\
    Un   - max distance below plane for a coplanar facet (default Vn)\n\
    Wn   - min facet width for outside point (before roundoff, default 2Vn)\n\
\n\
Output formats (may be combined; if none, produces a summary to stdout):\n\
    f    - all fields of all facets\n\
    i    - vertices incident to each facet (centers for Voronoi)\n\
    m    - Mathematica output (2-d and 3-d)\n\
    o    - OFF file format (dim, points and facets; Voronoi cells)\n\
    n    - normals with offsets\n\
    p    - point coordinates (Voronoi centers for option 'v')\n\
    s    - print summary to stderr\n\
    Fopts- additional input/output formats:\n\
    	   Fa      - print area for each facet\n\
           FA      - compute total area and volume for option 's'\n\
           Fc      - print count and coplanar points for each facet\n\
           FC      - print centrum or Voronoi center for each facet\n\
           Fd      - use cdd format for input (homogeneous with offset first)\n\
           FD      - use cdd format for normals (offset first)\n\
           FF      - print facets w/o ridges\n\
           Fi      - print inner planes for each facet\n\
           FI      - print id for each facet\n\
           Fm      - print merge count for each facet (511 max)\n\
           Fn      - print count and neighbors for each facet\n\
           FN      - print count and vertex neighbors for each point\n\
           Fo      - print outer planes (or max_outside) for each facet\n\
           FO      - print all options\n\
           Fp      - print point for each halfspace intersection\n\
           FP      - print nearest vertex for each coplanar point (v,p,f,d)\n\
           FQ      - print Qhull and input command\n\
           Fs      - print summary: #int (6), dimension, #points,\n\
                       tot vertices, tot facets, #vertices in output, #facets\n\
                       #real (2), max outer plane and min vertex\n\
           FS      - print sizes: #int (0), #real(2) tot area, tot vol\n\
           Fv      - print count and vertices for each facet\n\
           FV      - print average vertex (feasible point for 'H')\n\
    Gopts- Geomview output (2-d, 3-d and 4-d; 2-d Voronoi)\n\
           Ga      - all points as dots\n\
            Gp     -  coplanar points and vertices as radii\n\
            Gv     -  vertices as spheres\n\
           Gi      - inner planes only\n\
            Gn     -  no planes\n\
            Go     -  outer planes only\n\
           Gc	   - centrums\n\
           Gh      - hyperplane intersections\n\
           Gr      - ridges\n\
           GDn     - drop dimension n in 3-d and 4-d output\n\
    Popts- Print options:\n\
           PAn     - keep n largest facets by area\n\
           Pdk:n   - drop facet if normal[k] <= n (default 0.0)\n\
           PDk:n   - drop facet if normal[k] >= n\n\
           Pg      - print good facets (needs 'QGn' or 'QVn')\n\
           PFn     - keep facets whose area is at least n\n\
           PG      - print neighbors of good facets\n\
           PMn     - keep n facets with most merges\n\
           Po      - force output.  If error, output neighborhood of facet\n\
           Pp      - do not report precision problems\n\
\n\
    .    - list of all options\n\
    -    - help message for all options\n\
";
/* for opts, don't assign 'e' or 'E' to a flag (already used for exponent),
   save 'Qu' for upperHull, farthest-point diagram */

char qh_prompt2[]= "\n\
qhull- compute the convex hull.  %s\n\
    input (stdin): dimension, #points, point coordinates\n\
options:\n\
    d      - Delaunay triangulation by lifting points to a paraboloid\n\
    v      - Voronoi vertices via the Delaunay triangulation\n\
    H1,1   - Halfspace intersection about [1,1,0,...]\n\
Precision options:\n\
    C-0    - handle precision problems by merging non-convex facets\n\
    A-0.9  - merge facets if cosine of angle > 0.9\n\
    Qx     - handle precision problems in 5-d and higher\n\
    Tv     - verify result: structure, convexity, and point inclusion\n\
Output options (default produces a summary to stdout):\n\
    f      - print all fields of all facets\n\
    G      - Geomview output (2-d, 3-d and 4-d)\n\
    i      - vertices incident to each facet (centers for Voronoi)\n\
    m      - Mathematica output (2-d and 3-d)\n\
    o      - OFF file format (if Voronoi, outputs cells)\n\
    n      - normals with offsets\n\
    p      - point coordinates (centers for Voronoi)\n\
    s      - summary to stderr\n\
    FA     - compute total area and volume\n\
    .      - concise list of all options\n\
    -      - help message for all options\n\
example:\n\
    rbox 100 s | qhull\n\
\n\
";
/* for opts, don't assign 'e' or 'E' to a flag (already used for exponent) */

char qh_prompt3[]= "\n\
All options for qhull %s (upper-case takes a number except for 'F')\n\
\n\
 delaunay       voronoi	       Halfspace      facet_dump     Geomview \n\
 incidences     mathematica    normals        points         off_file\n\
 summary\n\
 Farea          FArea-total    Fcoplanars     FCentrums      Fd-cdd-in\n\
 FD-cdd-out     FFacet-xridge  Finner         FIds           Fmerges\n\
 Fneighbors     FNeigh-vertex  Fouter         FOptions       Fpoint-intersect\n\
 FPoint_near    FQhull         Fsummary       FSize          Fvertices\n\
 FVertex-ave\n\
 Angle_max	Centrum_size\n\
 Error_round    Random_dist    Visible_min    Ucoplanar_max  Wide_outside\n\
 Gvertices      Gpoints        Gall_points    Gno_planes     Ginner\n\
 Gcentrums      Ghyperplanes   Gridges        Gouter         GDrop_dim\n\
 PArea-keep     Pdrop d0:0D0   Pgood          PFacet_area_keep\n\
 PGood_neighbors PMerge-keep   Poutput_forced Pprecision_not\n\
 QbBound 0:0.5  Qcoplanar      Qfurthest      Qgood_only     QGood_point\n\
 Qinside        Qmax_out       Qrandom        QRotate        Qsearch_1st\n\
 QupperHull     QVertex_good   Qvneighbors    Qxact_merge    Q1_no_angle\n\
 Q2_no_independ Q3_no_redundant Q4_no_old     Q5_no_check_out Q6_no_concave\n\
 Q7_depth_first Q8_no_near_in  \n\
 T4_trace       Tcheck_often   Tstatistics    Tverify        Tz_stdout\n\
 TFacet_log\n\
 TPoint_trace   TMerge_trace   TWide_trace    TVertex_stop   TCone_stop\n\
\n\
handle roundoff and precision errors with 'C-0' (2-d to 4-d) or 'Qx' (5-d up)\n\
";

int isatty (int);  /* returns 1 if stdin is a tty
		   if "Undefined symbol" this can be deleted, and in main()
		   also need to change rbox.c */

/*-------------------------------------------------
-main- processes the command line, calls qhull() to do the work, and exits
*/
int main(int argc, char *argv[]) {
  int curlong, totlong, exitcode, numpoints, dim;
  coordT *points;
  boolT ismalloc;

#if __MWERKS__
  char inBuf[BUFSIZ], outBuf[BUFSIZ], errBuf[BUFSIZ];
  SIOUXSettings.showstatusline= false;
  SIOUXSettings.tabspaces= 1;
  SIOUXSettings.rows= 40;
  if (setvbuf (stdin, inBuf, _IOFBF, sizeof(inBuf)) < 0   /* w/o, SIOUX I/O is slow*/
  || setvbuf (stdout, outBuf, _IOFBF, sizeof(outBuf)) < 0
  || (stdout != stderr && setvbuf (stderr, errBuf, _IOFBF, sizeof(errBuf)) < 0)) 
    fprintf (stderr, "qhull internal warning (main): could not change stdio to fully buffered.\n");
  argc= ccommand(&argv);
#endif

  if ((argc == 1) && isatty(0)) {      
    fprintf(stdout, qh_prompt2, qh_version);
    exit(qh_ERRnone);
  }
  if (argc > 1 && *argv[1] == '-' && !*(argv[1]+1)) {
    fprintf(stdout, qh_prompt, qh_version, qh_DEFAULTbox);
    exit(qh_ERRnone);
  }
  if (argc >1 && *argv[1] == '.' && !*(argv[1]+1)) {
    fprintf(stdout, qh_prompt3, qh_version);
    exit(qh_ERRnone);
  }
  qh_init_A (stdin, stdout, stderr, argc, argv);  /* sets qh qhull_command */
  exitcode= setjmp (qh errexit); /* simple statement for CRAY J916 */
  if (!exitcode) {
    qh_initflags (qh qhull_command);
    points= qh_readpoints (&numpoints, &dim, &ismalloc);
    qh_init_B (points, numpoints, dim, ismalloc);
    qh_qhull();
    qh_check_output();
    qh_produce_output();
    if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
      qh_check_points();
    exitcode= qh_ERRnone;
  }
  qh NOerrexit= True;  /* no more setjmp */
  qh_freeqhull(False);
  qh_memfreeshort (&curlong, &totlong);
  if (curlong || totlong) 
    fprintf (stderr, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n",
       totlong, curlong);
  return exitcode;
} /* main */

