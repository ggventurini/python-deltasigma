/*************************************************************
  AI-2-MIF   Version 1.0                     DATE:  5/6/94

	BY:  Deron Jackson
    E-mail:  (djackson@mit.edu)
  Platform:  Windows/DOS
  Compiler:  Borland C++ Version 2.0

  This code is written explicitly for converting MATLAB 4.2's
  Adobe Illustrator output into Frame MIF (Maker Interchange
  Format).  This allows MATLAB graphics to be imported into
  Frame.  The resulting importing graphic is then EDITABLE in
  Frame.

  NOTES:
  -----
  1) This code is a quick and somewhat dirty solution
  to the problem of imporing MATLAB graphics.  The code
  should support any MATLAB "line" graphic.  This includes
  any 2D, 3D, mesh and contour plots consisting only of lines
  and text.  However any SURFACE or BITMAPS graphics are not
  supported.

  2) The code contains almost no error checking so
  use it at your own risk.

  3) I am quite certain that this code will NOT WORK with
  anything other than a MATLAB AI file.  I wrote the program
  by deciphering the symbol meanings from a number of MATLAB
  outputs.  I had no references to the TRUE AI format.

  4) Object clipping is rather crudely implemented.  This
  keeps graphs from exceeding their axis edges.  However I have
  not accounted for a few unlikely possibilities.

  5) Matlab plots using NON-SOLID linetypes (dashed, dotted,
  etc.) which have a LARGE number of points will appear
  solid in FRAME.  You can fix this by ungrouping the graphic
  after it is imported and using the Frame SMOOTH command on
  the individual traces.  Frame seems to display the linetype
  correctly when the POLYLINE is smoothed.

  6) You will notice that if you SCALE the graphic once
  imported into Frame the text size will not scale with it.
  This is anoying but I don't know a way around it.  You
  can easily scale the text manually using the "character
  designer" menu.

  7) All colors are mapped to Frame BLACK or NONE.  Color
  support would not be to bad to add but I have no need for
  it.

  Modifications
  Changed handling of text so that the MIF special characters
  /\<> are escaped by a \.		(R. Schreier, 07/29/96)
  
  Mex-ified.				(R. Schreier, 07/30/97)

  Added color and text rotation 	(R. Schreier, 08/20/00)

  Change a leading minus signs into an EnDash 	(R. Schreier, 08/20/00)

*************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "mex.h"

#ifdef MACINTOSH
#include <unix.h>
#endif

#define PPI   72.0
#define TRUE    1
#define FALSE   0

#define debug 0

/* Make these GLOBALS */
FILE *aifile;
FILE *miffile;
int finished=TRUE;
int need_start=FALSE;
int need_finish=FALSE;
int grp=3;
int pen=0;
int points=1;
int rec_bound=FALSE;
int clip_in=FALSE;
int clip_out=FALSE;
float minx,maxx,miny,maxy;
float xoffset=0,yoffset=0;
int color = 0;
char *ColorTable[]={"Black", "Black", "Yellow", "Black", 
                  "Magenta", "Black", "Red", "Black",
                  "Cyan", "Black", "Green", "Black",
                  "Blue", "Black", "Black", "Black",
                  "Black", "Black", "Black", "Black", 
                  "Black", "Black", "Black", "Black"}; /* White and black are both 0 0 0 1??? */

/***********************************************************
CHECK_FINISH()
This procedure finishes a polyline command before adding
anything else to the MIF file.
************************************************************/
void check_finish()
{
   if(finished==FALSE && need_start==FALSE)
   {
      if(points==2)     /* All 2 pt. lines go in group 3 */
	 grp=3;
      else
	 grp=1;     /* Multiple point lines go in group 1 */
      fprintf(miffile,"   <GroupID %d> <NumPoints %d>>\n",grp,points);
      finished=TRUE;
   }
}


/***********************************************************
REC_PROCEDURE()
This procedure records the limits of the current bounding
box.  This is used to CLIP data.
************************************************************/
void rec_procedure(float x,float y, int start) {
   if(start)
   {
      minx=x;
      maxx=x;
      miny=y;
      maxy=y;
   }
   else
   {
      if(x<minx)
	 minx=x;
      if(x>maxx)
	 maxx=x;
      if(y<miny)
	 miny=y;
      if(y>maxy)
	 maxy=y;
   }
}


/***********************************************************
CLIP_END()
This procedure clips the end point of a line to the endge
of the bounding box.
************************************************************/
void clip_end(xs,ys,xe,ye)
float xs,ys;
float *xe,*ye;
{
   float xn,yn;
   if(!rec_bound)
   {
      if(*xe>maxx)
      {
	 if(*xe!=xs)
	    yn=((maxx-xs)/(*xe-xs))*(*ye-ys)+ys;
	 xn=maxx;
      }
      if(*xe<minx)
      {
	 if(*xe!=xs)
	    yn=((xs-minx)/(xs-*xe))*(*ye-ys)+ys;
	 xn=minx;
      }
      if(*ye>maxy)
      {
	 if(*ye!=ys)
	    xn=((maxy-ys)/(*ye-ys))*(*xe-xs)+xs;
	 yn=maxy;
      }
      if(*ye<miny)
      {
	 if(*ye!=ys)
	    xn=((ys-miny)/(ys-*ye))*(*xe-xs)+xs;
	 yn=miny;
      }
      *xe=xn;
      *ye=yn;
   }
}

/***********************************************************
CLIP_START()
This procedure clips the start point of a line to the edge
of the bounding box.
************************************************************/
void clip_start(xs,ys,xe,ye)
float *xs,*ys;
float xe,ye;
{
   float xn,yn;
   if(!rec_bound)
   {
      if(*xs>maxx)
      {
	 if(*xs!=xe)
	    yn=((*xs-maxx)/(*xs-xe))*(ye-*ys)+*ys;
	 xn=maxx;
      }
      if(*xs<minx)
      {
	 if(*xs!=xe)
	    yn=((minx-*xs)/(xe-*xs))*(ye-*ys)+*ys;
	 xn=minx;
      }
      if(*ys>maxy)
      {
	 if(*ys!=ye)
	    xn=((*ys-maxy)/(*ys-ye))*(xe-*xs)+*xs;
	 yn=maxy;
      }
      if(*ys<miny)
      {
	 if(*ys!=ye)
	    xn=((miny-*ys)/(ye-*ys))*(xe-*xs)+*xs;
	 yn=miny;
      }
      *xs=xn;
      *ys=yn;
   }
}


/***********************************************************
IN_BOUNDS()
This procedure checks to see if a point is inside the
bounding box.
************************************************************/
int in_bounds(x,y)
float x,y;
{
   if(x>=minx && x<=maxx && y>=miny && y<=maxy)
	 return(TRUE);
   return(FALSE);
}


/***********************************************************
DRAW_LINE()
This procedure actually writes the MIF POLYLINE command.
************************************************************/
void draw_line(xs,ys,xe,ye)
float xs,ys;
float xe,ye;
{
   if(need_start){
      if(debug){
         printf("POLYLINE\n");
	     printf("  Start Point==> X: %f  Y: %f",xs,ys);
	  }
      fprintf(miffile,"<PolyLine <Pen %d> <Fill 15> <ObColor `%s'>\n",pen,ColorTable[color]);
      fprintf(miffile,"   <Point %f%c %f%c>",(xs+xoffset)/PPI,34,(-ys+yoffset)/PPI,34);
      if(debug && clip_in){
		 printf("  CLIPPED\n");
		 fprintf(miffile," # CLIPPED\n");
      }
      else{
	   if(debug) printf("\n");
	 fprintf(miffile,"\n");
      }
      need_start=FALSE;
      finished=FALSE;
   }
   if(debug) printf("  Next Point==> X: %f  Y: %f",xe,ye);
   fprintf(miffile,"   <Point %f%c %f%c>",(xe+xoffset)/PPI,34,(-ye+yoffset)/PPI,34);
   if(clip_out)
   {
      if(debug) printf("  CLIPPED\n");
      fprintf(miffile," # CLIPPED\n");
   }
   else
   {
      if(debug) printf("\n");
      fprintf(miffile,"\n");
   }
   if(need_finish)
   {
      if(points==2)     /* All 2 pt. lines go in group 3 */
	 grp=3;
      else
	 grp=1;     /* Multiple point lines go in group 1 */
      fprintf(miffile,"   <GroupID %d> <NumPoints %d>>\n",grp,points);
      need_finish=FALSE;
      finished=TRUE;
      need_start=TRUE;
      points=1;
   }
}


/*******************************************/
/*            Main                         */
/*******************************************/
void ai2mif(char* filename){
   char c, aifilename[256],miffilename[256];
   char msg5[256];
   char msg4[256];
   char msg3[256];
   char msg2[256];
   char msg1[256];
   char msg[256];
   char ffamily[256];
   char talign[10];
   float fsize;
   float lwidth=0.5;
   float p_maxx,p_minx,p_maxy,p_miny;
   float lastx,lasty,newx,newy,tempx,tempy;
   float textx,texty;
   float seg1,seg2,seg3,seg4;
   float tempf,tempf2, pi=3.14159265358979;
   int tlength=1;
   float tangle=0;
   int i;
   int skip_line=FALSE;
    

   strcpy(talign,"Left"); /* Default LEFT text alignment */
   strcpy(aifilename,filename);
   strcat(aifilename,".ai");
   strcpy(miffilename,filename);
   strcat(miffilename,".mif");
   aifile = fopen(aifilename,"r");
#ifdef MACINTOSH
_fcreator = 'Fra5';
_ftype = 'TEXT';
#endif
   miffile = fopen(miffilename,"w");
   
   fputs("<MIFFile 4.00> # Generated by AI2MIF\n",miffile);
   fputs("# AI2MIF Written by D. Jackson\n",miffile);
   fputs("#        Modified by R. Schreier\n",miffile);
   fputs("# If you have comments e-mail djackson@mit.edu\n",miffile);
   fputs("<PenWidth 0.500>\n",miffile);
   fputs("<HeadCap Round>\n",miffile);
   fputs("<TailCap Round>\n",miffile);


   /* Read through HEADER.  Pull out the BOUNDING BOX */
   while(!feof(aifile) && strncmp("%%EndComments",msg,13))
   {
      fscanf(aifile,"%s", msg);
      if(!strncmp("%%BoundingBox:",msg,14))
      {
	 fscanf(aifile,"%f",&minx);
	 fscanf(aifile,"%f",&miny);
	 fscanf(aifile,"%f",&maxx);
	 fscanf(aifile,"%f",&maxy);
	 xoffset=-(minx+maxx)/2.0;
	 yoffset=(miny+maxy)/2.0;
	 p_minx=minx;                 /* Hold on to these for later */
	 p_maxx=maxx;
	 p_miny=miny;
	 p_maxy=maxy;
      }
   }

   /* Now LOOP through each string and interpret */
   while(!feof(aifile))
   {
      strcpy(msg5,msg4); /* This is a record of the last 5 fields. msg1 is the most recent */
      strcpy(msg4,msg3);
      strcpy(msg3,msg2);
      strcpy(msg2,msg1);
      strcpy(msg1,msg);

      if (fscanf(aifile, "%s", msg))
      {
   /* Set the font */
	 if(!strncmp("z",msg,1) && !strncmp("/_",msg5,2))
	 {
	    sprintf(ffamily,"%s",msg5+2);  /* Grab Font FAMILY */
	    sscanf(msg4,"%f",&fsize);      /* Grab Font Size */
	    sscanf(msg1,"%f",&tempf);      /* Find alignment */
	    if(tempf>0)
	       strcpy(talign,"Center");
	    else
	       strcpy(talign,"Left");
	    if(debug) printf("FONT CHANGE ==> %s, %f, %s\n",ffamily,fsize,talign);
	 }
   /* This sets the pen width */
	 if(!strncmp("w",msg1,1))
	 {
	    check_finish();
	    sscanf(msg2,"%f",&lwidth);
	    fprintf(miffile,"<PenWidth %f>\n",lwidth);
	    if(debug) printf("LINWIDTH ==> %f\n",lwidth);
	 }
   /* This means a change in linestyle */
	 if(!strncmp("d",msg,1))
	 {
	    check_finish();
	    fputs("<DashedPattern <DashedStyle ",miffile);
	    if(!strncmp("[]",msg2,2))
	    {
	       fputs("Solid>>\n",miffile);
	       if(debug) printf("LINESTYLE ==> Solid\n");
	    }
	    else if(!strncmp("[",msg2,1))
	    {
	       sscanf(msg2,"[%f]",&seg1);
	       fprintf(miffile,"Dashed> <NumSegments 2> <DashSegment %f pt> <DashSegment %f pt>>\n",seg1,seg1);
	       /* Notice that I repeated the dash segment twice.  Frame will BOMB with a FATAL ERROR */
	       /* if you try define only 1 segment.                                                  */
	       if(debug) printf("LINESTYLE ==> Dashed\n");
	    }
	    else if(!strncmp("[",msg3,1))
	    {
	       sscanf(msg2,"%f]",&seg2);
	       sscanf(msg3,"[%f",&seg1);
	       fprintf(miffile,"Dashed> <NumSegments 2> <DashSegment %f pt> <DashSegment %f pt>>\n",seg1,seg2);
	       if(debug) printf("LINESTYLE ==> Dotted\n");
	    }
	    else if(!strncmp("[",msg5,1))
	    {
	       sscanf(msg2,"%f]",&seg4);
	       sscanf(msg3,"%f",&seg3);
	       sscanf(msg4,"%f",&seg2);
	       sscanf(msg5,"[%f",&seg1);
	       fprintf(miffile,"Dashed> <NumSegments 4> <DashSegment %f pt> <DashSegment %f pt> <DashSegment %f pt> <DashSegment %f pt>>\n",seg1,seg2,seg3,seg4);
	       if(debug) printf("LINESTYLE ==> Dash-Dot\n");
	    }
	 }
   /* This comes after the text position */
	 if(!strncmp("]e",msg1,2)) {
	    check_finish();
	    sscanf(msg5,"%f",&tempf2);   /* -sin(theta) */
	    sscanf(msg4,"%f",&tempf);   /* cos(theta) */
		tangle = atan2(-tempf2,tempf)*180/pi;
	    sscanf(msg3,"%f",&textx);
	    sscanf(msg2,"%f",&texty);
	    sscanf(msg,"%d",&tlength);

	    /* Place all text on GROUP 2 */
	    fprintf(miffile,"<TextLine <GroupID 2> <Angle %0.1f>\n",tangle);
	    fprintf(miffile,"   <Font <FFamily `%s'> <FSize %f> <FPlain Yes> <FBold No> <FColor `%s'>>\n",
		    ffamily,fsize,ColorTable[color]);
	    fprintf(miffile,"   <TLOrigin %f%c %f%c> <TLAlignment %s>",
		    (textx+xoffset)/PPI,34,(-texty+yoffset)/PPI,34,talign);

	    c=fgetc(aifile);   /* Toss away a blank space */
	    c=fgetc(aifile);   /* Toss away a ( character */
	    if(debug) printf("TEXT STRING ==> ");
	    for(i=0;i<tlength;i++) {  /* Read text string */
	       c=fgetc(aifile);
	       if( i==0 )
		  if( c=='-' ) {
		      fprintf(miffile," <Char EnDash> <String `");
		      c=fgetc(aifile);
		      --tlength;
		      }
		   else
		      fprintf(miffile," <String `");
	       if(c=='\\') /* Check next char for '()' */
		  {
		  char next_c=fgetc(aifile);
		  if(next_c=='(' || next_c==')')
		     c=next_c;
		  else if(next_c!='\\')
		     ungetc(next_c,aifile);
		  }
	       if(c=='\\' || c=='/' || c=='>' || c=='<') /* Escape \/>< */
	          fputc('\\',miffile);
	       fputc(c,miffile);
	       if(debug) printf("%c",c);
	    }
	    fprintf(miffile,"'>>\n");
	    if(debug) printf("\n  Location==> X: %f  Y: %f\n",textx,texty);
	 }
   /* This is a color definition.  I map everything to BLACK */
	 if(!strncmp("k",msg1,1) || !strncmp("K",msg1,1))
	 {
	    tempf2=0;
	    sscanf(msg5,"%f",&tempf);
	    /* RS 2000.08.17: record CMYK color as a binary number 
	       and decode as follows: 1000 = 8 = cyan, ...*/
	    tempf2 += 8*tempf;
	    sscanf(msg4,"%f",&tempf);
	    tempf2 += 4*tempf;
	    sscanf(msg3,"%f",&tempf);
	    tempf2 += 2*tempf;
	    sscanf(msg2,"%f",&tempf);
	    tempf2 += 1*tempf;
	    color = (int)tempf2;
	    if( color<0 ) color=0;
	    if( color>15 ) color=0;

	    if(tempf2==0)
	       pen=15;
	    else
	       pen=0;
	 }
   /* This moves the marker to start a line */
	 if(!strncmp("m",msg1,1))
	 {
	    check_finish();
	    sscanf(msg3,"%f",&lastx);
	    sscanf(msg2,"%f",&lasty);
	    if(rec_bound)                /* Record bounding box */
	       rec_procedure(lastx,lasty,1);
	    need_start=TRUE;
	    points=1;
	 }
   /* This means draw a line between this point and the marker */
	 if(!strncmp("L",msg1,1))
	 {
	    sscanf(msg3,"%f",&newx);
	    sscanf(msg2,"%f",&newy);
	    tempx=newx;
	    tempy=newy;
	    if(rec_bound)
	       rec_procedure(newx,newy,0);
	    clip_in=!in_bounds(lastx,lasty); /* If not in bounds need to clip */
	    clip_out=!in_bounds(newx,newy);   /* TRUE if not in bounds */
	    if(clip_in && clip_out)
	    {
	       skip_line=TRUE;
	    }
	    else if(clip_in)
	    {
	       clip_start(&lastx,&lasty,newx,newy);
	    }
	    else if(clip_out)
	    {
	       clip_end(lastx,lasty,&newx,&newy);
	       need_finish=TRUE;
	    }
	    if(skip_line==FALSE)
	    {
	       points++;
	       draw_line(lastx,lasty,newx,newy);
	    }
	    else
	       skip_line=FALSE;
	    lastx=tempx;
	    lasty=tempy;
	 }
	 if(!strncmp("q",msg1,1))   /* Start bound box */
	 {
	    rec_bound=TRUE;
	    pen=15;
	 }
	 if(!strncmp("W",msg1,1))   /* End bound box */
	 {
	    rec_bound=FALSE;
	    check_finish();
	    if(debug) printf("BOUNDING BOX: X: %f %f  Y: %f %f\n",minx,maxx,miny,maxy);
	    fprintf(miffile,"# BOUNDING BOX: X: %f %f  Y: %f %f\n",
		    (minx+xoffset)/PPI,(maxx+xoffset)/PPI,(-miny+yoffset)/PPI,(-maxy+yoffset)/PPI);
	 }
	 if(!strncmp("Q",msg1,1))   /* This means reset the bounding box */
	 {
	    maxx=p_maxx;
	    minx=p_minx;
	    maxy=p_maxy;
	    miny=p_miny;
	    check_finish();
	    if(debug) printf("BOUNDING BOX: X: %f %f  Y: %f %f\n",minx,maxx,miny,maxy);
	    fprintf(miffile,"# BOUNDING BOX: X: %f %f  Y: %f %f\n",
		    (minx+xoffset)/PPI,(maxx+xoffset)/PPI,(-miny+yoffset)/PPI,(-maxy+yoffset)/PPI);
	 }
      }
      else
      {
	 printf("\nError reading from AI file.\n");
	 exit(1);
      }
   }

   /* This catches the last grouping since there are no more M's */
   if(finished==FALSE && need_start==FALSE)
   {
      if(points==2)     /* All 2 pt. lines go in group 3 */
	 grp=3;
      else
	 grp=1;     /* Multiple point lines go in group 1 */
      fprintf(miffile,"<GroupID %d> <NumPoints %d>>\n",grp,points);
      finished=TRUE;
   }

   /* Now group the TEXT, SINGLE LINES, and then the whole thing */
   fprintf(miffile,"<Group <ID 2> <GroupID 1>>\n");
   fprintf(miffile,"<Group <ID 3> <GroupID 1>>\n");
   fprintf(miffile,"<Group <ID 1>>\n");

   fclose(aifile);
   fclose(miffile);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    char filename[256];
    if( nrhs < 1 )
		strcpy(filename, "temp");
	else{
		if( !mxIsChar(prhs[0]) )
			mexErrMsgTxt("The first argument must be a string");
		else
			mxGetString(prhs[0],filename,256);
	}
    ai2mif(filename);
}


