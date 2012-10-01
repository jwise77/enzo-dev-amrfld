// $Id: hypre-init.cpp 10 2010-03-18 22:35:19Z jbordner $

/// @file      hypre-init.cpp
/// @author    James Bordner
/// @brief     Creates a hypre-solve input file
///
/// Usage: 
///
///   hypre-init N0 np0 np1 np2 num_levels offset serial [cg|mg|gmres|pcg] [periodic|dirichlet]
///
/// Creates an input file where unigrid base has global size N0**3
/// distributed on a <np0,np1,np2> processor grid, with num_levels total
/// levels.  The point mass may be offset or not ("offset" = [0|1]), and
/// the problem may be run serially or with MPI ("serial" = [0|1]).  Solvers
/// include Conjugate Gradient (cg), FAC / multigrid (mg), GMRES (gmres), or
/// AMG-preconditioned multigrid [NOT IMPLEMENTED] (pmg).  Boundary
/// conditions are dirichlet or periodic.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define DIR_MODE 0755
#define STRLEN   80
#define BOXSIZE 8e9

//----------------------------------------------------------------------
void usage (char ** argv)
{
  // Print usage and exit abnormally

  printf ("\n"
	  "Usage: %s N0 np0 np1 np2 num_levels offset serial [cg|mg|gmres|pcg] [periodic|dirichlet]\n",argv[0]);
  exit(1);
}

void find_limit_centered(
			 double lp3[3], // OUT: Lower point defining grid patch
			 double up3[3], // OUT: Upper point defining grid patch
			 int li3[3], // OUT: Lower index of grid
			 int     N0, //  IN: Global grid size
			 int  n3[3], //  IN: Local grid size
			 int ip3[3], //  IN: rank in MPI grid topology
			 int np3[3], //  IN: size of MPI grid topology
			 int r       //  IN: r0 ^ level
			 )
{



  for (int i=0; i<3; i++) {
    lp3[i] = (-BOXSIZE/2 +  ip3[i]   *BOXSIZE / np3[i])/r;
    up3[i] = (-BOXSIZE/2 + (ip3[i]+1)*BOXSIZE / np3[i])/r;
    li3[i] = N0*(r-1)/2 + ip3[i]*n3[i];
  }
}
//----------------------------------------------------------------------

void find_limit_offset(
		       double lp3[3], // OUT: Lower point defining grid patch
		       double up3[3], // OUT: Upper point defining grid patch
		       int li3[3], // OUT: Lower index of grid
		       int     N0, //  IN: Global grid size
		       int  n3[3], //  IN: Local grid size
		       int ip3[3], //  IN: rank in MPI grid topology
		       int np3[3], //  IN: size of MPI grid topology
		       int r       //  IN: r0 ^ level
		       )
{
  for (int i=0; i<3; i++) {
    lp3[i] = BOXSIZE* ip3[i]   /r/np3[i];
    up3[i] = BOXSIZE*(ip3[i]+1)/r/np3[i];
    li3[i] = ip3[i]*n3[i];
  }
}

//----------------------------------------------------------------------

int main(int argc, char **argv)
{

//-----------------------------------------------------------------------
// Parse command-line arguments
//-----------------------------------------------------------------------

  const int expected_args = 9;
  if (argc != expected_args + 1) {
    printf ("Number of arguments %d is not %d\n",argc-1,expected_args);
    usage(argv);
  }

  int N0, np3[3],num_levels,is_offset,is_serial;
  char * solver, *boundary;

  N0         = atoi(argv[1]);
  np3[0]     = atoi(argv[2]);
  np3[1]     = atoi(argv[3]);
  np3[2]     = atoi(argv[4]);
  num_levels = atoi(argv[5]);
  is_offset  = atoi(argv[6]);
  is_serial  = atoi(argv[7]);
  solver     = argv[8];
  boundary   = argv[9];

  int np = np3[0]*np3[1]*np3[2];

  //-----------------------------------------------------------------------
  // Generate input file
  //-----------------------------------------------------------------------

  // Get input, output, and directory names

  char dir[STRLEN],infile[STRLEN];
  sprintf (dir,    "N%d.P%d%d%d.L%d.O%d.S%d.%s.%s",N0,np3[0],np3[1],np3[2],num_levels,is_offset,is_serial,solver,boundary);
  sprintf (infile,"in.%s",dir);
  FILE * fp = fopen (infile,"w");

  // Create new directory...

  if (mkdir (dir,DIR_MODE) != 0) {
    fprintf (stderr,"Cannot make directory %s.\n",dir);
    exit(1);
  }

  // And change path to the new directory

  char path[STRLEN];
  getcwd(path, STRLEN);
  strcat (path,"/");
  strcat (path,dir);
  chmod (path,DIR_MODE);

  // Copy required files to new directory

  link ("../hypre-solve",".");

  // Compute local grid sizes
  // WARNING: assumes # processors along each dimension evenly divides problem size

  int i,n3[3];

  for (i=0; i<3; i++) {
    n3[i] = N0 / np3[i];  
    if ((n3[i] * np3[i]) != N0) {
      fprintf (stderr,"Problem size %d not evenly divisible by np3[%d] = %d",
	       N0,i,np3[i]);
      exit(1);
    }
  }

  // Loop though processors, generating a grid per processor in each level

  int ip3[3];
  double lp3[3],up3[3];
  int li3[3];
  int id_point[8];
  double point_pos[8][3];
  int level;

  int r0 = 2;  // refinement factor
  int r = 1;   // r0^level

  // Loop over levels

  for (level = 0; level < num_levels; level++, r *= r0) {

    // Loop over processors

    for (ip3[2] = 0; ip3[2] < np3[2]; ip3[2]++) {
      for (ip3[1] = 0; ip3[1] < np3[1]; ip3[1]++) {
	for (ip3[0] = 0; ip3[0] < np3[0]; ip3[0]++) {

	  if (is_offset == false) {
	    find_limit_centered(lp3,up3,li3,N0,n3,ip3,np3,r);
	  } else {
	    find_limit_offset(lp3,up3,li3,N0,n3,ip3,np3,r);
	  }
               
	  int ip = ip3[0] + np3[0]*(ip3[1] + np3[1]*ip3[2]);

	  int id = ip + np*level;

	  int id_parent;
	  if (level == 0) {
	    id_parent = -1;
	  } else {
	    int offset = np*(level-1);
	    int i0 = ((1-is_offset)*r0*np3[0] + 4*ip3[0])/(4*r0);
	    int i1 = ((1-is_offset)*r0*np3[1] + 4*ip3[1])/(4*r0);
	    int i2 = ((1-is_offset)*r0*np3[2] + 4*ip3[2])/(4*r0);
	    id_parent = offset + i0 + np3[0]*(i1 + np3[1]*i2);
	  }

	  fprintf (fp, "grid %d %d %d "
		   "%g %g %g  %g %g %g  "
		   "%d %d %d  %d %d %d\n",
		   id,id_parent,(is_serial) ?  0 : ip,
		   lp3[0],lp3[1],lp3[2], up3[0],up3[1],up3[2],
		   li3[0],li3[1],li3[2],  n3[0], n3[1], n3[2]);

	  // Store id for point.  Requires coarse-to-finer level loop
	  // since it overwrites values.

	  //	  if (level == num_levels-1) {
	    for (int k0=0; k0<2; k0++) {
	      double p0 = (k0-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5 
		+ is_offset*BOXSIZE*pow(0.5,num_levels);
	      for (int k1=0; k1<2; k1++) {
		double p1 = (k1-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5  
		  + is_offset*BOXSIZE*pow(0.5,num_levels);
		for (int k2=0; k2<2; k2++) {
		  double p2 = (k2-0.5)*BOXSIZE/(N0*pow(2.0,level-1))*0.5  
		    + is_offset*BOXSIZE*pow(0.5,num_levels);
		  int k = k0 + 2*(k1 + 2*k2);
		  if (lp3[0] < p0 && p0 < up3[0] &&
		      lp3[1] < p1 && p1 < up3[1] &&
		      lp3[2] < p2 && p2 < up3[2]) {
		    id_point[k] = id;
		    point_pos[k][0] = p0;
		    point_pos[k][1] = p1;
		    point_pos[k][2] = p2;
		  }
		}
	      }
	      //	    }
	}
	}
      }
    }
  }
  // Create problem

  fprintf (fp, "dimension 3\n");
  if (is_offset) {
    fprintf (fp, "domain    3 0e9 0e9 0e9  8e9 8e9 8e9\n");
  } else {
    fprintf (fp, "domain    3 -4e9 -4e9 -4e9  4e9 4e9 4e9\n");
  }
  fprintf (fp, "boundary  %s\n", boundary);
  for (int k=0; k<8; k++) {
    fprintf (fp, "point     1e43   %g %g %g   %d\n",
	     point_pos[k][0],point_pos[k][1],point_pos[k][2],id_point[k]);
  }
  fprintf (fp, "discret constant\n");
  fprintf (fp, "solver_restol 1e-2\n");
  if (strcmp(solver,"mg")==0) {
    fprintf (fp, "solver %s\n",((num_levels==1) ? "pfmg" : "fac") );
    fprintf (fp, "solver_itmax 20\n");
  } else   if (strcmp(solver,"cg")==0) {
    fprintf (fp, "solver bicgstab\n");
    fprintf (fp, "solver_itmax 100\n");
  } else   if (strcmp(solver,"gmres")==0) {
    fprintf (fp, "solver gmres\n");
    fprintf (fp, "solver_itmax 500\n");
  } else   if (strcmp(solver,"pcg")==0) {
    fprintf (fp, "solver bicgstab-boomer\n");
    fprintf (fp, "solver_itmax 100\n");
  } else {
    usage(argv);
  }

  //========================================================================

  // Exit normally

  fclose (fp);
    
  exit(0);
}

