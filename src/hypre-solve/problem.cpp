// $Id: problem.cpp 19 2010-03-19 23:04:28Z jbordner $

/// @file      problem.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @brief     Implementation of the Problem class

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <map>
#include <string>
#include <vector>

#include "HYPRE_sstruct_ls.h"

#include "newgrav-hypre-solve.h"


//----------------------------------------------------------------------

#define debug false
#define trace false

//----------------------------------------------------------------------

#include "newgrav-scalar.h"
#include "newgrav-constants.h"
#include "newgrav-error.h"
#include "newgrav-parameters.h"
#include "newgrav-point.h"
#include "newgrav-faces.h"
#include "newgrav-mpi.h"
#include "newgrav-domain.h"
#include "newgrav-grid.h"
#include "newgrav-level.h"
#include "newgrav-hierarchy.h"
#include "newgrav-problem.h"

//======================================================================

Problem::Problem () throw ()
{
  //
}

//----------------------------------------------------------------------

Problem::~Problem () throw ()
{
  deallocate_ ();
  //
}

//----------------------------------------------------------------------

Problem::Problem (const Problem & p) throw ()
{
  domain_    = p.domain_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
}

//----------------------------------------------------------------------

Problem & Problem::operator = (const Problem & p) throw ()
{
  domain_    = p.domain_;
  points_    = p.points_;
  hierarchy_ = p.hierarchy_;
  return *this;
}

//----------------------------------------------------------------------

void Problem::read (std::string filename) throw ()
{

  bool boundary_defined = false;

  parameters_.read(filename);
  ItParameters itp (parameters_);

  while (itp++) {

    std::string key   = itp.key();
    std::string value = itp.value();

    if (key == "dimension") {

      int d = atoi(value.c_str());
      hierarchy_.set_dim(d);
      Point::set_dim(d);

      // Domain

    } else if (key == "domain") {

      domain_.input (value);

      // Grid ...

    } else if (key == "grid") {

      hierarchy_.insert_grid(new Grid(value));

    } else if (key == "point") {

      points_.push_back(new Point(value));      

    } else if (key == "boundary") {
      boundary_defined = true;
    }
  }

  if (! boundary_defined ) {
    fprintf (stderr,"Input parameter 'boundary' must be defined\n");
    MPI_Finalize();
    exit(1);
  }
}

//----------------------------------------------------------------------

void Problem::print () throw ()
{
  int i;
  domain_.print();
  hierarchy_.print();
  for (i=0; i<num_points(); i++)  point(i).print();
}

//----------------------------------------------------------------------

void Problem::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;
  int i;
  domain_.write(fp);
  hierarchy_.write(fp);
  for (i=0; i<num_points(); i++)  point(i).write(fp);
}

//----------------------------------------------------------------------

void Problem::deallocate_ () throw ()
{
  for (unsigned i=0; i<points_.size(); i++) {
    delete points_[i];
    points_[i] = 0;
  }
  points_.resize(0);
}
