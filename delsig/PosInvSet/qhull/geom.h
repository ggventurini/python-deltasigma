/* geom.h -- header file for geometric routines

   see README and geom.c
   
   copyright (c) 1993-1995 The Geometry Center        
*/

#ifndef qhDEFgeom
#define qhDEFgeom 1

/* ============ -macros- ======================== */

/*----------------------------------------------
-fabs_(a)          returns the absolute value of a
-fmax_(a,b)        returns the maximum value of a and b
-fmin_(a,b)        returns the minimum value of a and b
-maximize_(maxval, val)  sets maxval to val if greater
-minimize_(minval, val)  sets minval to val if less
*/
#define fabs_(a) (((a) < 0) ? -(a):(a))
#define fmax_(a,b)  ( (a) < (b) ? (b) : (a) )
#define fmin_(a,b)  ( (a) > (b) ? (b) : (a) )
#define maximize_(maxval, val) {if ((maxval) < (val)) (maxval)= (val);}
#define minimize_(minval, val) {if ((minval) > (val)) (minval)= (val);}

/*-----------------------------------------------
-det2_(a1, a2, 		2-d determinate
       b1, b2)
-det3_(a1, a2, a3,      3-d determinate
       b1, b2, b3,
       c1, c2, c3)
*/
#define det2_(a1,a2,b1,b2) ((a1)*(b2) - (a2)*(b1))
#define det3_(a1,a2,a3,b1,b2,b3,c1,c2,c3) ( (a1)*det2_(b2,b3,c2,c3) \
		- (b1)*det2_(a2,a3,c2,c3) + (c1)*det2_(a2,a3,b2,b3) )  

/*-----------------------------------------------
-dX, dY, dZ- coordinate differences given row pointers rows[]
*/
#define dX(p1,p2)  (*(rows[p1]) - *(rows[p2]))
#define dY(p1,p2)  (*(rows[p1]+1) - *(rows[p2]+1))
#define dZ(p1,p2)  (*(rows[p1]+2) - *(rows[p2]+2))
#define dW(p1,p2)  (*(rows[p1]+3) - *(rows[p2]+3))

/* ======= -functions =========== 

   see geom.c for definitions

      	Geometric functions
-crossproduct   compute the cross product of 2 3-d vectors
-determinant    compute the determinant of a square matrix
-detsimplex     return determinate of a simplex of points
-divzero        divide by a number that's nearly zero
-facetarea_simplex return area for a simplex defined by an apex, base, orient, normal
-gausselim      Gaussian elimination with partial pivoting
-getangle       return cosine of angle (dot product of two qh hull_dim vectors)
-gram_schmidt   implements Gram-Schmidt orthogonalization by rows
-inthresholds   return True if normal within qh lower_/upper_threshold
-maxabsval      return max absolute value of a vector
-maxsimplex	determines maximum simplex for a set of points 
-minabsval      return min absolute value of a dim vector
-mindiff        return index of min abs. difference of two vectors
-normalize      normalize a vector
-normalize2     normalize a vector and report if too small
-pointdist      return distance between two points
-printmatrix    print matrix given by row vectors
-printpoints    print pointids for a set of points starting at index 
-projectpoints  project points along one or more dimensions
-rand/srand     generate random numbers
-randomfactor	return a random factor within qh RANDOMdistmax of 1.0
-randommatrix   generate a random dimXdim matrix in range (-1,1)
-rotatepoints   rotate numpoints points by a row matrix
-scalepoints    scale points to new lowbound and highbound
-sethalfspace_all generate dual for halfspace intersection with feasible point
-sethyperplane_det return hyperplane for oriented simplex, uses determinates
-sethyperplane_gauss return hyperplane for oriented simplex, uses Gaussian elimination
-voronoi_center return Voronoi center for a set of points

      	Qhull's geometric functions
-backnormal     solve for normal x using back substitution over rows U
-distplane      return distance from point to facet (>0 if point is above facet)
-facetarea      return area for a facet
-facetcenter    return Voronoi center for a facet's vertices
-findbest	find visible facet for a point starting at a facet
-findbestnew    find best newfacet for point
-findgooddist   find best good facet visible for point from facet
-getarea        get area of all facets in facetlist, collect statistics
-getcenter      return arithmetic center of a set of vertices
-getcentrum     return centrum for a facet
-getdistance    returns the max and min distance of any vertex from neighbor
-maxmin         return max/min points for each dim., sets max roundoff errors
-orientoutside  make facet outside oriented via qh interior_point
-projectinput   project input using qh DELAUNAY and qh low_bound/high_bound
-projectpoint   project point onto a facet by distance
-rotateinput    rotate input using row matrix
-scaleinput     scale input using qh low_bound/high_bound
-setfacetplane  sets the hyperplane for a facet
-sethalfspace   set coords to dual of halfspace relative to feasible point
*/

/*---------- -prototypes in alphabetical order, infrequent at end -----------*/
void    qh_backnormal (realT **rows, int numrow, int numcol, boolT sign, coordT *normal, boolT *nearzero);
void	qh_distplane (pointT *point, facetT *facet, realT *dist);
facetT *qh_findbest (pointT *point, facetT *facet, boolT bestoutside,
		boolT newfacets, realT *dist, boolT *isoutside, int *numpart);
facetT *qh_findbestnew (pointT *point, facetT *startfacet,
	   realT *dist, boolT *isoutside, int *numpart);
void 	qh_gausselim(realT **rows, int numrow, int numcol, boolT *sign, boolT *nearzero);
realT   qh_getangle(pointT *vect1, pointT *vect2);
pointT *qh_getcenter(setT *vertices);
pointT *qh_getcentrum(facetT *facet);
realT   qh_getdistance(facetT *facet, facetT *neighbor, realT *mindist, realT *maxdist);
void    qh_normalize (coordT *normal, int dim, boolT toporient);
void    qh_normalize2 (coordT *normal, int dim, boolT toporient, 
            realT *minnorm, boolT *ismin);
pointT *qh_projectpoint(pointT *point, facetT *facet, realT dist);

void    qh_setfacetplane(facetT *newfacets);
void 	qh_sethyperplane_det (int dim, coordT **rows, coordT *point0, 
              boolT toporient, coordT *normal, realT *offset, boolT *nearzero);
void 	qh_sethyperplane_gauss (int dim, coordT **rows, pointT *point0, 
	     boolT toporient, coordT *normal, coordT *offset, boolT *nearzero);

  /* ======= infrequently used code in geom2.c =========== */
void    qh_crossproduct (int dim, realT vecA[3], realT vecB[3], realT vecC[3]);
realT 	qh_determinant (realT **rows, int dim, boolT *nearzero);
realT   qh_detsimplex(pointT *apex, setT *points, int dim, boolT *nearzero);
realT   qh_divzero(realT numer, realT denom, realT mindenom1, boolT *zerodiv);
realT   qh_facetarea (facetT *facet);
realT   qh_facetarea_simplex (int dim, coordT *apex, setT *vertices, 
          vertexT *notvertex,  boolT toporient, coordT *normal, realT *offset);
pointT *qh_facetcenter (setT *vertices);
boolT   qh_findbestsharp (pointT *point, facetT **bestfacet, 
           realT *bestdist, int *numpart);
facetT *qh_findgooddist (pointT *point, facetT *facetA, realT *distp, facetT **facetlist);
void    qh_getarea (facetT *facetlist);
boolT   qh_gram_schmidt(int dim, realT **rows);
boolT   qh_inthresholds (coordT *normal, realT *angle);
realT  *qh_maxabsval (realT *normal, int dim);
setT   *qh_maxmin(pointT *points, int numpoints, int dimension);
void    qh_maxsimplex (int dim, setT *maxpoints, pointT *points, int numpoints, setT **simplex);
realT   qh_minabsval (realT *normal, int dim);
int     qh_mindiff (realT *vecA, realT *vecB, int dim);
boolT   qh_orientoutside (facetT *facet);
coordT  qh_pointdist(pointT *point1, pointT *point2, int dim);
void    qh_printmatrix (FILE *fp, char *string, realT **rows, int numrow, int numcol);
void    qh_printpoints (FILE *fp, char *string, setT *points);
void    qh_projectinput (void);
void 	qh_projectpoints (signed char *project, int n, realT *points, 
             int numpoints, int dim, realT *newpoints, int newdim);
int     qh_rand( void);
void    qh_srand( int seed);
realT   qh_randomfactor (void);
void    qh_randommatrix (realT *buffer, int dim, realT **row);
void    qh_rotateinput (realT **rows);
void    qh_rotatepoints (realT *points, int numpoints, int dim, realT **rows);
void    qh_scaleinput (void);
void 	qh_scalepoints (pointT *points, int numpoints, int dim,
  		realT *newlows, realT *newhighs);
boolT   qh_sethalfspace (int dim, coordT *coords, coordT **nextp, 
              coordT *normal, coordT *offset, coordT *feasible);
coordT *qh_sethalfspace_all (int dim, int count, coordT *halfspaces, pointT *feasible);
pointT *qh_voronoi_center (int dim, setT *points);

#endif /* qhDEFgeom */



