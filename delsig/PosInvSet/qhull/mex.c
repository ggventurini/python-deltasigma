/* hull.c - A MEX-file interface to qhull 
Returns vertices, edges, normals and offsets.
*/

#include "qhull_a.h" 
#include "mex.h"
char qh_version[] = "mex interface to qhull V2.2 by R. Schreier 95.12.21";

/* Fill up the return matrices, keeping everything in order */
Matrix* stuff2D( double *v, double *n, double *o){
    Matrix *E = mxCreateFull(2, qh num_vertices, REAL);
    double *e = mxGetPr(E);
    vertexT *vertex = qh vertex_list;
    facetT *facet;
    coordT *coord;
    int i,j;
    qh_vertexneighbors();
    qh visit_id++;
    for( i=0; i< qh num_vertices; ++i ){
        coord = vertex->point;
        for( j=0; j< qh hull_dim; j++)
            *v++ = *coord++;
        vertex->visitid = qh visit_id;

        *e++ = i+1;
        if( i == qh num_vertices-1 )
            *e = 1;
        else
            *e++ = i+2;

        facet = (facetT *)(SETfirst_(vertex->neighbors));
        if( facet->visitid == qh visit_id )
            facet = (facetT *)(SETsecond_(vertex->neighbors));
        facet->visitid = qh visit_id;

        coord = facet->normal;
        for( j=0; j< qh hull_dim; j++)
            *n++ = *coord++;

        *o++ = facet->offset;

        vertex = (vertexT *)(SETfirst_(facet->vertices));
        if( vertex->visitid == qh visit_id )
            vertex = (vertexT *)(SETsecond_(facet->vertices));
    }
    return E;
}

/* Fill up the return matrices for the 3D case */
Matrix *stuff3D( double *v, double *n, double *o){
    Matrix *E;
    double *e;
    facetT *facet, *neighbor;
    ridgeT **ridgep, *ridge;	/* needed by the FOREACHridge_ macro */
    vertexT *vertex;
    coordT *coord;
    int i, vi=1, *vertex_numbers;
    char string[256];

    vertex_numbers = (int*)mxCalloc(qh vertex_id, sizeof(int));
    FORALLvertices{
        for( i=0, coord = vertex->point; i< qh hull_dim; i++)
            *v++ = *coord++;
        if( vertex->id >= qh vertex_id ){
            sprintf(string,"Internal Error. vertex->id=%d, qh vertex_id=%d.", vertex->id, qh vertex_id );
            mexErrMsgTxt(string);
        }
        vertex_numbers[vertex->id] = vi++;
    }

    {
    int num_edges = 0;
    FORALLfacets
	if (facet->simplicial)
	    num_edges += qh hull_dim;
	else
	    num_edges += qh_setsize(facet->ridges);
    E = mxCreateFull(2, num_edges/2 , REAL);
    }

    e = mxGetPr(E);
    qh visit_id++;
    FORALLfacets{
        facet->visitid = qh visit_id;
        if( facet->simplicial )
            qh_makeridges(facet); /* creates explicit ridges for simplicial facets*/
        FOREACHridge_(facet->ridges){
            neighbor = otherfacet_(ridge, facet);
            if (neighbor->visitid != qh visit_id) {
                vertex = (vertexT*)SETfirst_(ridge->vertices);
                *e++ = vertex_numbers[vertex->id];
                vertex = (vertexT*)SETsecond_(ridge->vertices);
                *e++ = vertex_numbers[vertex->id];
            }
        } /* FOREACHridge() */
    }

    FORALLfacets{
        for( i=0, coord = facet->normal; i< qh hull_dim; i++)
            *n++ = *coord++;
        *o++ = facet->offset;
    }
    return E;
} /* stuff3D(...) */

mexPrintPoint( double *x, char *fmt ){
    int i;
    for(i=0; i<qh hull_dim; ++i)
	mexPrintf(fmt,*x++);
}

int note_edge( int **connect, int vi, int vj, int id){
    int *conn;
    if( vi<vj ) 
	conn = &connect[vj][vi];
    else
	conn = &connect[vi][vj];

    if( *conn == 0 ){
	/* New edge. Mark it. */
	*conn = id;
	return 1;
    } else if( *conn == id ) {	
	/* Duplicate edge. Delete it */
	*conn = 0;
	return -1;
    } else 
	/* Edge seen from intersection of another facet pair. Do nothing. */
	return 0;
}

/* Return an int that is unique to the (m,n)/(n,m) pair.  m,n >=0, m!=n. */
int id_number( int m, int n ){
    if( n>m )
	return n*(n-1)/2 +m;
    else
	return m*(m-1)/2 +n;
}

/* Fill up the return matrices for the 4D case */
Matrix *stuff4D( double *v, double *n, double *o){
    Matrix *E;
    double *e;
    int num_edges = 0;
    facetT *facet, *neighbor;
    ridgeT **ridgep, *ridge;	/* needed by the FOREACHridge_ macro */
    vertexT **vertexp, *vertex;	/* needed by the FOREACHvertex_ macro */
    coordT *coord;
    coordT centrum[4], centrum_1, centrum_2, n1[4], n2[4];
    int i, j, vi=1, vj, *vertex_numbers;
    char string[256];
    int **connect, conn_id;

    /* Fill the vertex array */
    vertex_numbers = (int*)mxCalloc(qh vertex_id, sizeof(int));
    FORALLvertices{
        for( i=0, coord = vertex->point; i< qh hull_dim; i++)
            *v++ = *coord++;
        if( vertex->id >= qh vertex_id ){
            sprintf(string,"Internal Error. vertex->id=%d, qh vertex_id=%d.", vertex->id, qh vertex_id );
            mexErrMsgTxt(string);
        }
        vertex_numbers[vertex->id] = vi++;
    }

    /* Fill the normal and offset arrays */
    FORALLfacets{
        for( i=0, coord = facet->normal; i< qh hull_dim; i++)
            *n++ = *coord++;
        *o++ = facet->offset;
    }

    /* Initialize the connection matrix */
    connect = (int**)mxCalloc(qh num_vertices, sizeof(int*));
    connect[0] = (int*)mxCalloc((qh num_vertices-1)*(qh num_vertices)/2, sizeof(int));
    for( i=0; i<qh num_vertices -1; ++i)
	connect[i+1] = connect[i] + i;

    /* Fill the connection matrix and keep track of the number of edges */
    qh visit_id++;
    FORALLfacets{
	facet->visitid = qh visit_id;
	if( facet->simplicial )
	    qh_makeridges(facet); /* creates explicit ridges for simplicial facets*/
	FOREACHridge_(facet->ridges){
	    int vi, vip, vi0=-1;
	    neighbor = otherfacet_(ridge, facet);
	    conn_id = id_number(facet->id, neighbor->id);
	    if (neighbor->visitid == qh visit_id) /* Ridge already seen */
		continue;
	    FOREACHvertex_(ridge->vertices){
		vi = vertex_numbers[vertex->id]-1;
		if( vi0 == -1 )
		    vi0 = vi;
		else
		    num_edges += note_edge(connect,vi,vip,conn_id);
		vip = vi;
	    }
	    num_edges += note_edge(connect,vi,vi0,conn_id);
	} /* FOREACHridge() */
    } /* FORALLfacets */

    /* Use the connection matrix to fill the edge list */
    E = mxCreateFull(2, num_edges, REAL);
    e = mxGetPr(E);
    for(i=1; i<qh num_vertices ; ++i){
	for(j=0; j<i; ++j){
	    if( connect[i][j] ){
                *e++ = i+1;
                *e++ = j+1;
	    }
	}
    }
    return E;
} /* stuff4D(...) */

void cleanup(){
    int curlong, totlong;
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
    if (curlong || totlong)      /* optional */
        fprintf (stderr, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n",
            totlong, curlong);
}

void mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[]){
    Matrix *V, *E, *N, *O;
    double *Vp, *Np, *Op;
    int exitcode, dim;
    char *cmd = "qhull(mex) Qcx C0.001 A0.999", string[256];

    /* check the arguments */
    if( nrhs > 2 )
        mexErrMsgTxt("hull error: only two input arguments are allowed.");
    if( nrhs < 1 )
        mexErrMsgTxt("hull error: no input argument.");
    if( mxIsNumeric(prhs[0]) != 1 )
        mexErrMsgTxt("hull error: The first argument is not numeric.");
    if( nrhs > 1 ){	/* Use arg2 as the argument string for qhull */
	Matrix *str = prhs[1];
	int strlen = mxGetM(str)*mxGetN(str)+1;
	if( mxIsString(str) != 1 )
	    mexErrMsgTxt("hull error: the second argument is not a string.");
	cmd = mxCalloc(strlen,sizeof(char));
	mxGetString(str, cmd, strlen);
	}
    if( nlhs != 4 )
        mexErrMsgTxt("hull error: 4 left-hand arguments are needed.");
    dim = mxGetM( prhs[0] );
    if( dim > 4 )
        mexErrMsgTxt("hull error: Unable to handle dimension >4.");

    qh_meminit (stderr);
    qh_initqhull_start (stdin, stdout, stderr);

    if (!(exitcode= setjmp (qh errexit))) {
        strcat (qh rbox_command, "hull");    /* for qh_printstatistics() */
        strcat (qh qhull_command, cmd);    /* Needed to set parameters */
        qh_initflags (qh qhull_command);
        qh_initqhull_globals (mxGetPr(prhs[0]), mxGetN(prhs[0]), dim, False);
        qh_initqhull_mem();
        qh_initqhull_buffers();

        qh_qhull();

        V = plhs[0] = mxCreateFull(dim, qh num_vertices, REAL);
        N = plhs[2] = mxCreateFull(dim, qh num_facets, REAL);
        O = plhs[3] = mxCreateFull(1, qh num_facets, REAL);

	Vp = mxGetPr(V);
	Np = mxGetPr(N);
	Op = mxGetPr(O);

/* DBG 
qh_printfacetlist(qh facet_list, NULL, True);
*/
        switch(dim){
        case 2:
            E = plhs[1] = stuff2D(Vp,Np,Op);
            break;
        case 3:
            E = plhs[1] = stuff3D(Vp,Np,Op);
            break;
        case 4:
            E = plhs[1] = stuff4D(Vp,Np,Op);
            break;
        default:
            break;
        }

        exitcode= qh_ERRnone;
    }

    qh NOerrexit= True;
    cleanup();
} /* mexFunction */

/*-------------------------------------------
-errexit- return exitcode to system after an error
  assumes exitcode non-zero
  prints useful information
  see qh_errexit2() in qhull.c for 2 facets
*/
void qh_errexit(int exitcode, facetT *facet, ridgeT *ridge) {

    if (qh ERREXITcalled) {
        fprintf (qh ferr, "\nqhull error while processing previous error.  Exit program\n");
        exit(1);
    }
    qh ERREXITcalled= True;
    if (!qh QHULLfinished)
        qh hulltime= clock() - qh hulltime;
    qh_errprint("ERRONEOUS", facet, NULL, ridge, NULL);
    fprintf (qh ferr, "\nWhile executing: %s | %s\n", qh rbox_command, qh qhull_command);
    if (qh furthest_id) {
        fprintf(qh ferr, "Last point added to hull was p%d.", qh furthest_id);
        if (zzval_(Ztotmerge))
            fprintf(qh ferr, "  Last merge was #%d.", zzval_(Ztotmerge));
    }
    if (qh ROTATErandom > 0)
        fprintf(qh ferr, "  Repeat with 'QR%d'.", qh ROTATErandom);
    if (qh furthest_id || qh ROTATErandom > 0)
        fprintf (qh ferr, "\n");
    if (qh FORCEoutput && (qh QHULLfinished || (!facet && !ridge)))
        qh_produce_output();
    else {
        if (exitcode != qh_ERRsingular && zzval_(Zsetplane) > qh hull_dim+1) {
            fprintf (qh ferr, "\nAt error exit:\n");
            qh_printsummary (qh ferr);
            if (qh PRINTstatistics) {
                qh_collectstatistics();
                qh_printstatistics(qh ferr, "at error exit");
                qh_memstatistics (qh ferr);
            }
        }
        if (qh PRINTprecision)
            qh_printstats (qh ferr, qhstat precision, NULL);
    }
    if (!exitcode)
        exitcode= qh_ERRqhull;
    else if (exitcode == qh_ERRsingular)
        qh_printhelp_singular(qh ferr);
    else if (exitcode == qh_ERRprec && !qh PREmerge)
        qh_printhelp_degenerate (qh ferr);
    if (qh NOerrexit) {
        fprintf (qh ferr, "qhull error while ending program.  Exit program\n");
        exit(1);
    }
    qh NOerrexit= True;
    longjmp(qh errexit, exitcode);
} /* errexit */


/*-------------------------------------------
-errprint- prints out the information of the erroneous object
    any parameter may be NULL, also prints neighbors and geomview output
*/
void qh_errprint(char *string, facetT *atfacet, facetT *otherfacet, ridgeT *atridge, vertexT *atvertex) {
    int i;

    if (atfacet) {
        fprintf(qh ferr, "%s FACET:\n", string);
        qh_printfacet(qh ferr, atfacet);
    }
    if (otherfacet) {
        fprintf(qh ferr, "%s OTHER FACET:\n", string);
        qh_printfacet(qh ferr, otherfacet);
    }
    if (atridge) {
        fprintf(qh ferr, "%s RIDGE:\n", string);
        qh_printridge(qh ferr, atridge);
        if (atridge->top && atridge->top != atfacet && atridge->top != otherfacet)
            qh_printfacet(qh ferr, atridge->top);
        if (atridge->bottom
            && atridge->bottom != atfacet && atridge->bottom != otherfacet)
            qh_printfacet(qh ferr, atridge->bottom);
        if (!atfacet)
            atfacet= atridge->top;
        if (!otherfacet)
            otherfacet= otherfacet_(atridge, atfacet);
    }
    if (atvertex) {
        fprintf(qh ferr, "%s VERTEX:\n", string);
        qh_printvertex (qh ferr, atvertex);
    }
    if (qh FORCEoutput && atfacet && !qh QHULLfinished) {
        fprintf(qh ferr, "ERRONEOUS and NEIGHBORING FACETS to output\n");
        for (i= 0; i< qh_PRINTEND; i++)
            qh_printneighborhood (qh fout, qh PRINTout[i], atfacet, otherfacet,
                !qh_ALL);
    }
} /* errprint */


/*-----------------------------------------
-printfacetlist- print all fields for a list and/or set of facets to .ferr
  includes all vertices
  if !printall, only prints good facets
*/
void qh_printfacetlist(facetT *facetlist, setT *facets, boolT printall) {
    facetT *facet, **facetp;

    qh_printbegin (qh ferr, qh_PRINTfacets, facetlist, facets, printall);
    FORALLfacet_(facetlist)
        qh_printafacet(qh ferr, qh_PRINTfacets, facet, printall);
    FOREACHfacet_(facets)
        qh_printafacet(qh ferr, qh_PRINTfacets, facet, printall);
    qh_printend (qh ferr, qh_PRINTfacets, facetlist, facets, printall);
} /* printfacetlist */


/*-----------------------------------------
-user_memsizes- allocate up to 10 additional, quick allocation sizes
*/
void qh_user_memsizes (void) {

    /* qh_memsize (size); */
} /* user_memsizes */
