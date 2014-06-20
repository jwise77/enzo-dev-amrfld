/// @file      AMRsolve_hierarchy.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @author    Daniel Reynolds (reynolds@smu.edu)
/// @brief     Implementation of the AMRsolve_Hierarchy class

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <vector>
#include <map>
#include <string>

#include "HYPRE_sstruct_ls.h"

//----------------------------------------------------------------------

const bool trace           = false;
const bool debug_hierarchy = false;
const bool debug_detailed  = false;
const bool geomview        = false;

const int maximum_level    = 1000;

//----------------------------------------------------------------------

#include "AMRsolve_scalar.h"
#include "AMRsolve_constants.h"
#include "AMRsolve_error.h"
#include "AMRsolve_point.h"
#include "AMRsolve_domain.h"
#include "AMRsolve_faces.h"
#include "AMRsolve_mpi.h"
#include "AMRsolve_grid.h"
#include "AMRsolve_level.h"
#include "AMRsolve_hierarchy.h"
#include "AMRsolve_defs.h"

//----------------------------------------------------------------------

#ifdef AMR_SOLVE /* defined if we're interfacing with enzo */

#include "AMRsolve_enzo.h"

// Define what Enzo thinks is a float

#  ifdef CONFIG_BFLOAT_4
#    define ENZO_FLOAT float
#  endif
#  ifdef CONFIG_BFLOAT_8
#    define ENZO_FLOAT double
#  endif

#endif /* ifdef AMR_SOLVE */

//----------------------------------------------------------------------

// #include "debug.h"
//----------------------------------------------------------------------

/// AMRsolve_Hierarchy class constructor.
AMRsolve_Hierarchy::AMRsolve_Hierarchy() throw()
  : dimension_(0)
{
  grids0_.push_back(0);  // Create a NULL sentinel
  levels0_.push_back(0); // Create a NULL sentinel
}
	  
//----------------------------------------------------------------------

/// AMRsolve_Hierarchy class descructor.
AMRsolve_Hierarchy::~AMRsolve_Hierarchy() throw()
{
  deallocate_();
}

//======================================================================

#ifdef AMR_SOLVE

int AMRsolve_Hierarchy::enzo_attach_grav(LevelHierarchyEntry *LevelArray[],
					 int level_coarsest,
					 int level_finest) throw()
{
  _TRACE_;
  set_dim(3);

  // Determine Grid ID's

  // Count grids
  int id = 0;

  // Traverse levels in enzo hierarchy
  std::map<grid *,int> enzo_id;
  if (debug_hierarchy) printf("DEBUG: %s:%d Called enzo_attach_grav (%d %d)\n",
			      __FILE__,__LINE__,level_coarsest,level_finest);

  // Set defaults for level_coarsest and level_finest:
  //    if no levels specified: default to entire hierarchy
  //    if one level specified: default to just that level
  if (level_coarsest == LEVEL_UNKNOWN && level_finest == LEVEL_UNKNOWN) {
    level_coarsest = 0;
    level_finest   = INT_MAX;
  } else if (level_coarsest != LEVEL_UNKNOWN && level_finest == LEVEL_UNKNOWN) {
    level_finest   = level_coarsest;
  }

  if (!(level_coarsest <= level_finest)) {
    ERROR("level_finest must not be less than level_coarsest");
  }

  if (debug_hierarchy) printf("DEBUG: %s:%d Attaching to Enzo (%d %d)\n",
			      __FILE__,__LINE__,level_coarsest,level_finest);

  // Loop over levels
  for (int level = level_coarsest; 
       level <= level_finest && LevelArray[level] != NULL;
       level++) {

    // Traverse grids in enzo level
    for (LevelHierarchyEntry * itEnzoLevelGrid = LevelArray[level];
	 itEnzoLevelGrid != NULL;
	 itEnzoLevelGrid = itEnzoLevelGrid->NextGridThisLevel) {

      grid* enzo_grid = itEnzoLevelGrid->GridData;
      if (enzo_grid == NULL) {
	fprintf(stderr,"GridData is NULL\n");
	return(1);
      }

      // Save grid id's to determine parents
      enzo_id[enzo_grid] = id;

      // Set owning processor for this grid
      int ip = enzo_grid->ProcessorNumber;
      
      // Collect required enzo grid data
      int id_parent;
      if (level == level_coarsest) {
	id_parent = -1;
      } else {
	grid* enzo_parent = itEnzoLevelGrid->GridHierarchyEntry->ParentGrid->GridData;
	if (enzo_parent == NULL) {
	  fprintf(stderr,"ParentGrid->GridData is NULL\n");
	  return(1);
	}
	id_parent = enzo_id[enzo_parent];
      }

      // Compute global grid size at this level
      long long nd[3];
      for (int dim=0; dim<3; dim++) {
	nd[dim] = (long long) ((DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
		       enzo_grid->CellWidth[dim][0] + 0.5);
      }
	
      // Compute active size of this grid, and location information
      Scalar xl[3],xu[3];
      int il[3],n[3];
      for (int dim=0; dim<3; dim++) {
	xl[dim] = enzo_grid->GridLeftEdge[dim];
	xu[dim] = enzo_grid->GridRightEdge[dim];
	n[dim]  = enzo_grid->GridEndIndex[dim] - enzo_grid->GridStartIndex[dim] + 1;
	il[dim] = nint((xl[dim] - DomainLeftEdge[dim]) 
		       / (DomainRightEdge[dim]-DomainLeftEdge[dim]) * nd[dim]);
      }

      // Create a new amr_solve grid
      AMRsolve_Grid* grid = new AMRsolve_Grid(id,id_parent,ip,xl,xu,il,n);
      insert_grid(grid);

      // if this is a local grid, allocate data and set pointers to Enzo data
      if (grid->is_local()) {

	// allocate the u and f arrays on this grid
	grid->allocate();

	// Set pointers to the relevant Enzo data arrays for setting up the 
	// linear system later on.
	Scalar* phi   = enzo_grid->AccessPotentialField();
	Scalar* gmass = enzo_grid->AccessGravitatingMassField();
	if (phi == NULL) {
	  fprintf(stderr,"PotentialField is NULL\n");
	  return(1);
	}
	if (gmass == NULL) {
	  fprintf(stderr,"GravitatingMassField is NULL\n");
	  return(1);
	}
	grid->set_phi(phi);
	grid->set_gmass(gmass);

	// Set Enzo gravity ghost zone information for this grid
	int grav_ghosts = (enzo_grid->GravitatingMassFieldDimension[0] - n[0])/2;
	int enzo_ghosts[][2] = { {grav_ghosts, grav_ghosts}, 
				 {grav_ghosts, grav_ghosts}, 
				 {grav_ghosts, grav_ghosts} };
	grid->set_GravGhosts(enzo_ghosts);
      }
      id++;
    } // grids within level

  } // level

  return(0);
  _TRACE_;
}
#endif

//======================================================================

#ifdef AMR_SOLVE

void AMRsolve_Hierarchy::enzo_attach_fld(LevelHierarchyEntry *LevelArray[],
					 int NumBins,
					 int level_coarsest,
					 int level_finest) throw()
{
  _TRACE_;
  set_dim(3);

  // Determine Grid ID's

  // Count grids
  int id = 0;

  // Traverse levels in enzo hierarchy
  std::map<grid *,int> enzo_id;
  if (debug_hierarchy) printf("DEBUG: %s:%d Called enzo_attach_fld (%d %d)\n",
			      __FILE__,__LINE__,level_coarsest,level_finest);

  // Set defaults for level_coarsest and level_finest:
  //    if no levels specified: default to entire hierarchy
  //    if only one level specified: default to just that level
  if (level_coarsest == LEVEL_UNKNOWN && level_finest == LEVEL_UNKNOWN) {
    level_coarsest = 0;
    level_finest   = INT_MAX;
  } else if (level_coarsest != LEVEL_UNKNOWN && level_finest == LEVEL_UNKNOWN) {
    level_finest   = level_coarsest;
  }

  if (!(level_coarsest <= level_finest)) {
    ERROR("level_finest must not be less than level_coarsest");
  }

  if (debug_hierarchy) printf("DEBUG: %s:%d Attaching to Enzo (%d %d)\n",
			      __FILE__,__LINE__,level_coarsest,level_finest);

  // Loop over levels
  for (int level = level_coarsest; 
       level <= level_finest && LevelArray[level] != NULL;
       level++) {

    // Traverse grids in enzo level
    for (LevelHierarchyEntry* itEnzoLevelGrid = LevelArray[level];
	 itEnzoLevelGrid != NULL;
	 itEnzoLevelGrid = itEnzoLevelGrid->NextGridThisLevel) {

      grid* enzo_grid = itEnzoLevelGrid->GridData;

      // Save grid id's to determine parents
      enzo_id[enzo_grid] = id;

      // Set owning processor for this grid
      int ip = enzo_grid->ProcessorNumber;
      
      // Collect required enzo grid data
      int id_parent;
      if (level == level_coarsest) {
	id_parent = -1;
      } else {
	grid* enzo_parent = itEnzoLevelGrid->GridHierarchyEntry->ParentGrid->GridData;
	id_parent = enzo_id[enzo_parent];
      }

      // Compute global grid size at this level
      long long nd[3];
      for (int dim=0; dim<3; dim++) {
	nd[dim] = (long long) ((DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
		       enzo_grid->CellWidth[dim][0] + 0.5);
      }

      // Compute active size of this grid, and location information
      Scalar xl[3],xu[3];
      int il[3],n[3];
      for (int dim=0; dim<3; dim++) {
	xl[dim] = enzo_grid->GridLeftEdge[dim];
	xu[dim] = enzo_grid->GridRightEdge[dim];
	n[dim]  = enzo_grid->GridEndIndex[dim] - enzo_grid->GridStartIndex[dim] + 1;
	il[dim] = nint((xl[dim] - DomainLeftEdge[dim]) 
		       / (DomainRightEdge[dim]-DomainLeftEdge[dim]) * nd[dim]);
      }

      // Create a new amr_solve grid
      AMRsolve_Grid* grid = new AMRsolve_Grid(id,id_parent,ip,xl,xu,il,n);
      insert_grid(grid);

      // if this is a local grid, allocate data and set pointers to Enzo data
      if (grid->is_local()) {

	// allocate the u and f arrays on this grid
	grid->allocate();

	// Set pointers to the relevant Enzo data arrays for setting up the 
	// linear system later on.
	Scalar *ptr;
	if (NumBins > 0) {
	  ptr = enzo_grid->AccessRadiationFrequency0();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 0 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(0,ptr);
	  ptr = enzo_grid->AccessEmissivity0();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 0 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(0,ptr);
	}
	if (NumBins > 1) {
	  ptr = enzo_grid->AccessRadiationFrequency1();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 1 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(1,ptr);
	  ptr = enzo_grid->AccessEmissivity1();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 1 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(1,ptr);
	}
	if (NumBins > 2) {
	  ptr = enzo_grid->AccessRadiationFrequency2();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 2 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(2,ptr);
	  ptr = enzo_grid->AccessEmissivity2();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 2 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(2,ptr);
	}
	if (NumBins > 3) {
	  ptr = enzo_grid->AccessRadiationFrequency3();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 3 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(3,ptr);
	  ptr = enzo_grid->AccessEmissivity3();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 3 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(3,ptr);
	}
	if (NumBins > 4) {
	  ptr = enzo_grid->AccessRadiationFrequency4();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 4 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(4,ptr);
	  ptr = enzo_grid->AccessEmissivity4();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 4 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(4,ptr);
	}
	if (NumBins > 5) {
	  ptr = enzo_grid->AccessRadiationFrequency5();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 5 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(5,ptr);
	  ptr = enzo_grid->AccessEmissivity5();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 5 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(5,ptr);
	}
	if (NumBins > 6) {
	  ptr = enzo_grid->AccessRadiationFrequency6();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 6 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(6,ptr);
	  ptr = enzo_grid->AccessEmissivity6();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 6 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(6,ptr);
	}
	if (NumBins > 7) {
	  ptr = enzo_grid->AccessRadiationFrequency7();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 7 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(7,ptr);
	  ptr = enzo_grid->AccessEmissivity7();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 7 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(7,ptr);
	}
	if (NumBins > 8) {
	  ptr = enzo_grid->AccessRadiationFrequency8();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 8 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(8,ptr);
	  ptr = enzo_grid->AccessEmissivity8();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 8 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(8,ptr);
	}
	if (NumBins > 9) {
	  ptr = enzo_grid->AccessRadiationFrequency9();
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get E bin 9 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_E(9,ptr);
	  ptr = enzo_grid->AccessEmissivity9();  
	  if (!ptr) {
	    printf("grid %i, proc %i: cannot get eta bin 9 from Enzo!\n",grid->id(),grid->ip());
	    assert(ptr);
	  }
	  grid->set_eta(9,ptr);
	}
	ptr = enzo_grid->AccessKPhHI();       // temporary field for old radiation arrays
	if (!ptr) {
	  printf("grid %i, proc %i: cannot get E0 from Enzo!\n",grid->id(),grid->ip());
	  assert(ptr);
	}
	grid->set_E0(ptr);
	ptr = enzo_grid->AccessPhotoGamma();  // temporary field for opacity arrays
	if (!ptr) {
	  printf("grid %i, proc %i: cannot get PhotoGamma from Enzo!\n",grid->id(),grid->ip());
	  assert(ptr);
	}
	grid->set_kap(ptr);

	// Set Enzo ghost zone information for this grid
	int enzo_ghosts[][2] = { {DEFAULT_GHOST_ZONES, DEFAULT_GHOST_ZONES}, 
				 {DEFAULT_GHOST_ZONES, DEFAULT_GHOST_ZONES}, 
				 {DEFAULT_GHOST_ZONES, DEFAULT_GHOST_ZONES} };
	grid->set_Ghosts(enzo_ghosts);
      }
      id++;
    } // grids within level

  } // level

  _TRACE_;
}
#endif

//----------------------------------------------------------------------

#ifdef AMR_SOLVE
void AMRsolve_Hierarchy::enzo_detach() throw()
{
  _TRACE_;
  deallocate_();
  _TRACE_;
}
#endif

//======================================================================

/// Insert a grid into the hierarchy.
/** Inserts the grid's pointer into the grids0_ vector, at the
  * position given by the grid's id(), and always maintaining a 0
  * sentinel at the end of grids0_ */
void AMRsolve_Hierarchy::insert_grid(AMRsolve_Grid* pgrid) throw()
{
  // Insert grid into the list of all grids in the Hierarchy
  if (unsigned(pgrid->id()) + 1 >= grids0_.size()) {
    grids0_.resize (pgrid->id() + 2);
    grids0_[grids0_.size() - 1] = 0;
  }
  grids0_[pgrid->id()] = pgrid;
}

//======================================================================

/// Initialize the hierarchy after all grids have been inserted
/** After all grids are inserted into the hierarchy, this function
    determines each grid's parent, containing level, children, 
    neighbors, face-zone categories, and global indices. */
void AMRsolve_Hierarchy::initialize(AMRsolve_Domain& domain, AMRsolve_Mpi& mpi,
				    bool is_periodic[3]) throw()
{
  if (debug_hierarchy) printf("AMRsolve_Hierarchy::init_levels()\n");

  init_grid_parents_();
  init_grid_levels_();
  init_indices_(is_periodic);     // DEPENDENCY: Requires init_grid_levels_()
  init_extents_(is_periodic);     // DEPENDENCY: Requires init_grid_levels_()
  init_grid_children_();
  init_grid_neighbors_();
  _TRACE_;
  init_grid_faces_(domain, mpi);  // DEPENDENCY: Requires init_indices_()
  _TRACE_;

  geomview_grids(mpi);
}

//------------------------------------------------------------------------

//======================================================================

/// Write hierarchy grids to geomview files grid-L<level>-P<processor>.vect 
void AMRsolve_Hierarchy::geomview_grids(AMRsolve_Mpi& mpi) throw()
{
  if (geomview) {
    ItHierarchyLevels itl (*this);
    while (AMRsolve_Level * level = itl++) {
      char filename[20];
      sprintf(filename,"grid-L%d-P%d.vect",level->index(),mpi.ip());
      FILE * fp = fopen (filename,"w");
      level->geomview_grid_local (fp);
      fclose (fp);
    }
  }
}

//======================================================================

void AMRsolve_Hierarchy::print() throw()
{
  printf("AMRsolve_Hierarchy\n");
  printf("   levels      = %d\n",int(levels0_.size()-1));
  printf("   grids       = %d\n",int(grids0_.size()-1));
  printf("   dimension   = %d\n",dimension_);
  for (int i=0; i<num_levels(); i++) {
    level(i).print();
  }
}

//----------------------------------------------------------------------

void AMRsolve_Hierarchy::write(FILE* fp) throw()
{
  fprintf(fp,"AMRsolve_Hierarchy\n");
  for (int i=0; i<num_levels(); i++) {
    level(i).write();
  }
}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

/// Determines the parent grid of each grid in the hierarchy.
void AMRsolve_Hierarchy::init_grid_parents_() throw()
{
  ItHierarchyGridsAll itg (*this);
  while (AMRsolve_Grid * g = itg++) {
    AMRsolve_Grid * p = (g->id_parent() >= 0) ? & return_grid(g->id_parent()) : 0;
    this->set_parent(g,p);
  }
}

//------------------------------------------------------------------------

/// Determines the level of each grid in the hierarchy.
void AMRsolve_Hierarchy::init_grid_levels_() throw()
{
  bool done = false;
  int pass_count = 0;

  while (! done && pass_count++ < maximum_level) {

    done = true;

    // Loop through grids, and try to determine their level if unknown
    ItHierarchyGridsAll itg (*this);
    while (AMRsolve_Grid * g = itg++) {

      // If grid's level < 0, then we haven't determined its level yet
      if (g->level() < 0) { 

	// Grids without parents must be in level 0 ...
	if (parent(*g) == 0) {
	  g->set_level(0);
	  insert_in_level_ (0,*g);
	}
	
	// ... grids with parents of known level have level = 1 +
	// parent level ...
	else if (parent(*g)->level() >= 0) {
	  int level = parent(*g)->level() + 1;
	  g->set_level(level);
	  insert_in_level_ (level,*g);
	} 
	
	// ... otherwise a grid's parents is in an unknown level, so
	// we're not done yet
	else {
	  done = false;
	}
      }
    }
  } // while

  if (pass_count >= maximum_level) {
    printf("%s:%d Hmmm, pass_count >= maximum_level\n",__FILE__,__LINE__);
    ItHierarchyGridsAll itg (*this);
    while (AMRsolve_Grid * g = itg++) {
      char filename[20];
      sprintf(filename,"debug-grid-%d",g->id());
      g->write("nofields",filename);
    }
    ERROR("too many passes in init_grid_levels_ loop: dumped grids");
  }
}

//------------------------------------------------------------------------

///  Initialize grid children.
void AMRsolve_Hierarchy::init_grid_children_() throw()
{
  ItHierarchyGridsAll itg (*this);
  while (AMRsolve_Grid * g = itg++) {
    // If a grid has a parent, then the grid is the parent's child
    if (parent(*g) != 0) {
      parent(*g)->set_child(*g);
    }
  }
}


//------------------------------------------------------------------------

///  Find each grid's neighbors.
/**  If a grid is in level 0, then find its neighbors by comparing
     with all other grids in the level.  Otherwise, if a grid is in a
     level > 0, then we can save time by only testing grids that are
     children of the parent, and children of all the parents
     neighbors. */
void AMRsolve_Hierarchy::init_grid_neighbors_() throw()
{

  int i,j;

  // For level == 0, test all pairs
  for (i=0; i<level(0).num_grids(); i++) {

    AMRsolve_Grid* g1 = & level(0).return_grid(i);

    // Starting index of loop was i+1, but need to include possibility
    // of a grid being a neighbor with itself for periodic problems
    for (j=i; j<level(0).num_grids(); j++) {

      AMRsolve_Grid* g2 = & level(0).return_grid(j);

      if (g1->is_adjacent(*g2,period_domain_)) {
	assert_neighbors(*g1,*g2);
      }
    }
  }

  // For levels > 0, only test parents' children, and parents'
  // neighbors' children.  If two grids' parents are not the same or
  // not neighbors, then they necessarily cannot be neighbors.
  int k,j1,j2;
  for (k=1; k<num_levels(); k++) {
    for (i=0; i<level(k).num_grids(); i++) {

      AMRsolve_Grid* g1 = & level(k).return_grid(i);

      // Check parents' children
      _TRACE_;
      for (j=0; j<parent(*g1)->num_children(); j++) {
	AMRsolve_Grid* g2 = & parent(*g1)->child(j);
	if (g1->is_adjacent(*g2,period_domain_) && g1->id() > g2->id()) {
	  assert_neighbors (*g1,*g2);
	}
      }

      // Check parents' neighbors' children
      _TRACE_;
      for (j1=0; j1<parent(*g1)->num_neighbors(); j1++) {
	AMRsolve_Grid* gn = & parent(*g1)->neighbor(j1);
	for (j2=0; j2<gn->num_children(); j2++) {
	  AMRsolve_Grid* g2 = & gn->child(j2);
	  if (g1->is_adjacent(*g2,period_domain_) && g1->id() > g2->id()) {
	    if (debug_hierarchy) printf("DEBUG 3 grids %d and %d are neighbors\n",
			       g1->id(),g2->id());
	    assert_neighbors(*g1,*g2);
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------

///  Find each grid's face categories.
/** After all grid inter-connections are determined, this function
    determines the neighbor structure for each face-zone. Only
    performed for local grids, and grids neighboring local grids. 
    Allocates AMRsolve_Faces for grids that haven't had them allocated
    yet. */
void AMRsolve_Hierarchy::init_grid_faces_(AMRsolve_Domain& domain,
					  AMRsolve_Mpi& mpi) throw()

{
  int axis, face;
  int ig0,in0,il0,iu0;
  int ig1,in1,il1,iu1;
  int j0,j1;            // face axes

  // ------------------------------------------------------------
  // Determine pointers to neighboring grid for each face-zone
  // ------------------------------------------------------------

  ItHierarchyLevels itl (*this);
  AMRsolve_Level* level;

  while ((level = itl++)) {

    ItLevelGridsAll itg (*level);

    // Determine level periodicity
    int period3[3];
    period_index(period3,level->index());

    while (AMRsolve_Grid* grid = itg++) {

      if (debug_hierarchy) grid->print();
      ItGridNeighbors itn (*grid);

      // Set "adjacent" pointers for adjacent grids in same level
      AMRsolve_Grid* adjacent;
      while ((adjacent = itn++)) {
	if (grid->is_local() || adjacent->is_local()) {
	  int count = 0;
	  while (grid->neighbor_shared_face
		 (*adjacent,axis,face,il0,il1,iu0,iu1,period3,count)) {
	    for (ig0=il0; ig0<=iu0; ig0++) {
	      for (ig1=il1; ig1<=iu1; ig1++) {
		grid->faces().adjacent(axis,face,ig0,ig1) = adjacent;
	      }
	    }
	  }
	}
      }

      // Set "adjacent" pointers for adjacent non-parent coarse grids
      if (parent (*grid) != NULL) {

	ItGridNeighbors itpn (*parent(*grid));

	while ((adjacent = itpn++)) {
	  if (grid->is_local() || adjacent->is_local()) {
	    int count=0;
	    while (grid->coarse_shared_face
		   (*adjacent,axis,face,il0,il1,iu0,iu1,period3,count)) {
	      for (ig0=il0; ig0<=iu0; ig0++) {
		for (ig1=il1; ig1<=iu1; ig1++) {
		  grid->faces().adjacent(axis,face,ig0,ig1) = adjacent;
		}
	      }
	    }
	  }
	}

	// Set remaining "adjacent" pointers to parent grid
	if (grid->is_local() || parent(*grid)->is_local()) {
	  int count = 0;

	  while (grid->parent_interior_face
		 (*parent(*grid),axis,face,il0,il1,iu0,iu1,count)) {
	    for (ig0=il0; ig0<=iu0; ig0++) {
	      for (ig1=il1; ig1<=iu1; ig1++) {
		if (grid->faces().adjacent(axis,face,ig0,ig1) == NULL) {
		  grid->faces().adjacent(axis,face,ig0,ig1) = parent (*grid);
		}
	      }
	    }
	  }
	}
      }

      if (debug_detailed) {
	int ig3[3][2];
	grid->indices(ig3);
	for (axis=0; axis<3; axis++) {
	  int j0 = (axis+1)%3;
	  int j1 = (axis+2)%3;
	  int n0 = ig3[j0][1] - ig3[j0][0];
	  int n1 = ig3[j1][1] - ig3[j1][0];
	  for (face=0; face<2; face++) {
	    for (ig0=0; ig0<n0; ig0++) {
	      for (ig1=0; ig1<n1; ig1++) {
		AMRsolve_Grid * neighbor = grid->faces().adjacent(axis,face,ig0,ig1);
		printf("ip=%d grid=%d neighbor=%d (axis=%d face=%d) [%d,%d]\n",
		       pmpi->ip(),grid->id(),neighbor?neighbor->id():-1,
		       axis,face,ig0,ig1);
	      }
	    }
	  }
	}
	fflush(stdout);
      }
    }
  }

  // ------------------------------------------------------------
  // 1. Categorize boundary face-zones
  // ------------------------------------------------------------

  // Domain boundary
  double db3[3][2]; 
  domain.lower(db3[0][0],db3[1][0],db3[2][0]);
  domain.upper(db3[0][1],db3[1][1],db3[2][1]);

  int ih0[3][2];
  this->indices0(ih0);
  while ((level = itl++)) {
    ItLevelGridsLocal itgl (*level);
    while (AMRsolve_Grid* grid = itgl++) {

      // Grid boundary
      double gb3[3][2];
      grid->x_lower(gb3[0][0],gb3[1][0],gb3[2][0]);
      grid->x_upper(gb3[0][1],gb3[1][1],gb3[2][1]);

      // Grid spacing (should be level-independent)
      double h3[3];
      grid->h(h3[0],h3[1],h3[2]);

      for (axis = 0; axis < 3; axis++) {
	for (face = 0; face < 2; face++) {

	  // Label as boundary if on domain boundary AND not periodic

	  bool is_boundary = 
	    fabs(gb3[axis][face] - db3[axis][face]) < 0.5*h3[axis];

	  if (!is_periodic(axis) && is_boundary) {
	    grid->faces().label(axis,face,AMRsolve_Faces::_boundary_);
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // 2. Categorize covered (and temporary adjacent-covered) face-zones
  // ------------------------------------------------------------

  for (int ilevel = num_levels()-1; ilevel>0; ilevel--) {

    AMRsolve_Level* level = &this->level(ilevel);

    ItLevelGridsAll itg (*level);
    while (AMRsolve_Grid* grid = itg++) {
      AMRsolve_Grid* parent = this->parent(*grid);
      if (grid->is_local() || parent->is_local()) {
	int count = 0;
	while (grid->parent_shared_face 
	       (*parent, axis, face, il0,il1,iu0,iu1,count)) {
	  int ig3[3][2];
	  grid->indices(ig3);
	  for (ig0=il0; ig0<=iu0; ig0++) {
	    for (ig1=il1; ig1<=iu1; ig1++) {
	      // Label covered
	      parent->faces().label(axis,face,ig0,ig1) = 
		AMRsolve_Faces::_covered_;
	      // Label adjacent-covered
	      AMRsolve_Grid* neighbor = parent->faces().adjacent(axis,face,ig0,ig1);
	      if (neighbor != NULL && neighbor->level() == grid->level()) {
		int in3[3][2];
		neighbor->indices(in3);
		j0 = (axis+1)%3;
		j1 = (axis+2)%3;
		in0 = ig0 + in3[j0][0] - ig3[j0][0];
		in1 = ig1 + in3[j1][0] - ig3[j1][0];
		if (neighbor->faces().label(axis,1-face,in0,in1) != 
		    AMRsolve_Faces::_covered_) {
		  neighbor->faces().label(axis,1-face,in0,in1) 
		    = AMRsolve_Faces::_adjacent_covered_;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // 3. Categorize fine and neighbor face-zones
  // ------------------------------------------------------------

  for (int ilevel = 0; ilevel<num_levels(); ilevel++) {
    AMRsolve_Level* level = &this->level(ilevel);

    // Determine level periodicity
    int period3[3];
    period_index (period3,level->index());

    ItLevelGridsAll itg (*level);
    while (AMRsolve_Grid* grid = itg++) {
      //      int ig3[3][2];
      //      grid->indices(ig3);
      ItGridNeighbors itn (*grid);
      while (AMRsolve_Grid* neighbor = itn++) {
	int count = 0;
	while (grid->neighbor_shared_face 
	       (*neighbor, axis, face, il0,il1,iu0,iu1,period3,count)) {
	  for (ig0=il0; ig0<=iu0; ig0++) {
	    for (ig1=il1; ig1<=iu1; ig1++) {
	      AMRsolve_Faces::Label& fz = 
		grid->faces().label(axis,face,ig0,ig1);
	      if (fz == AMRsolve_Faces::_adjacent_covered_) {
		fz = AMRsolve_Faces::_fine_;
	      } else if (fz != AMRsolve_Faces::_covered_) {
		fz = AMRsolve_Faces::_neighbor_;
	      }
	    }
	  }
	}
      }	  
    }
  }

  // ------------------------------------------------------------
  // 4. Categorize coarse face-zones
  // ------------------------------------------------------------

  while (AMRsolve_Level* level = itl++) {
    ItLevelGridsAll itga(*level);
    while (AMRsolve_Grid* grid = itga++) {
      int ig3[3][2];
      grid->indices(ig3);
      for (axis=0; axis<3; axis++) {
	int j0 = (axis+1)%3;
	int j1 = (axis+2)%3;
	int n0 = ig3[j0][1]-ig3[j0][0];
	int n1 = ig3[j1][1]-ig3[j1][0];
	for (face=0; face<2; face++) {
	  for (ig0=0; ig0<n0; ig0++) {
	    for (ig1=0; ig1<n1; ig1++) {
	      AMRsolve_Faces::Label& fz = 
		grid->faces().label(axis,face,ig0,ig1);
	      if (fz==AMRsolve_Faces::_unknown_) 
		fz = AMRsolve_Faces::_coarse_;
	    }
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // Dump out any requested geomview files
  // ------------------------------------------------------------

  if (geomview) {
    char filename[80];
    for (int ilevel=0; ilevel<this->num_levels(); ilevel++) {
      for (AMRsolve_Faces::Label type=AMRsolve_Faces::_first_; 
	   type<=AMRsolve_Faces::_last_; 
	   type = AMRsolve_Faces::Label(type+1)) {
	sprintf(filename,"facezones-L%d-P%d-%s.quad",
		ilevel,mpi.ip(),AMRsolve_Faces::LabelName[type]);
	FILE* fp = fopen(filename,"w");
	this->level(ilevel).geomview_face_types(fp,&type,1);
	fclose(fp);
      }
    }
  }
}

//======================================================================

void AMRsolve_Hierarchy::init_indices_(bool is_periodic[3]) throw()
{
  _TRACE_;
  // Determine problem size the hard way
  int i,lower[3],upper[3];
  for (i=0; i<3; i++) {
    lower[i] = INT_MAX;
    upper[i] = INT_MIN;
  }

  _TRACE_;

  ItLevelGridsAll itg (level(0));

  while (AMRsolve_Grid*grid = itg++) {
    _TRACE_;
    for (i=0; i<3; i++) {
      lower[i] = MIN(grid->index_lower(i),lower[i]);
      upper[i] = MAX(grid->index_upper(i),upper[i]);
    }
    _TRACE_;
  }

  _TRACE_;

  for (i=0; i<3; i++) {
    il0_[i] = lower[i];
    n0_[i] = upper[i] - lower[i] + 1;
  }

  _TRACE_;

  for (i=0; i<3; i++) {
    period_index_[i] = is_periodic[i] ? n0_[i] : 0.0;
  }

  _TRACE_;
}

//======================================================================

void AMRsolve_Hierarchy::init_extents_(bool is_periodic[3]) throw()
{
  _TRACE_;

  // Determine hierarchy extents from grid extents
  ItLevelGridsAll itg (level(0));

  int i;
  for (i=0; i<3; i++) {
    xl_[i] = SCALAR_MAX;
    xu_[i] = -SCALAR_MAX;
  }

  while (AMRsolve_Grid *grid = itg++) {

    // Get grid extents
    Scalar gl[3],gu[3];
    grid->x_lower(gl[0],gl[1],gl[2]);
    grid->x_upper(gu[0],gu[1],gu[2]);

    // Adjust hierarchy extents
    for (i=0; i<3; i++) {
      xl_[i] = MIN(xl_[i],gl[i]);
      xu_[i] = MAX(xu_[i],gu[i]);
    }
  }

  for (i=0; i<3; i++) {
    period_domain_[i] = is_periodic[i] ? xu_[i] - xl_[i] : 0.0;
  }

  _TRACE_;
}

//----------------------------------------------------------------------

void AMRsolve_Hierarchy::insert_in_level_(int level, AMRsolve_Grid& grid) throw()
{
  // Resize levels0_[] if needed
  if (unsigned(level + 1) >= levels0_.size()) {
    levels0_.resize (level + 2);
    levels0_[levels0_.size() - 1] = 0;
  }
  if (levels0_.at(level) == 0) {
    if (debug_hierarchy) printf("DEBUG: creating new Level at %d\n",level);
    levels0_[level] = new AMRsolve_Level(level);
  }
  levels0_[level]->insert_grid (grid);
}

//----------------------------------------------------------------------

void AMRsolve_Hierarchy::deallocate_() throw()
{
  // Delete levels
  for (unsigned i=0; i<levels0_.size(); i++) {
    if (levels0_[i] != NULL) {
      // AMRsolve_Level objects deleted here
      // AMRsolve_Grid objects are deleted below
      // AMRsolve_Level objects do not delete containing AMRsolve_Grid objects
      delete levels0_[i];
      levels0_[i] = 0;
    }
  }
  levels0_.resize(0);
  levels0_.push_back(0);

  // Delete grids list
  for (unsigned i=0; i<grids0_.size(); i++) {
    if (grids0_[i] != NULL) {
      // AMRsolve_Grid objects deleted here
      delete grids0_[i];
      grids0_[i] = 0;
    }
  }
  grids0_.resize(0);
  grids0_.push_back(0);
}

//----------------------------------------------------------------------
