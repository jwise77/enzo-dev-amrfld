// $Id: test_Parameters.cpp 2148 2011-03-25 01:29:44Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Parameters.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 21 16:04:03 PST 2008
/// @brief    Program implementing unit tests for the Parameters class
//----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <string>

#include "parameters.hpp"
#include "test.hpp"

//----------------------------------------------------------------------

/// @def      CLOSE
/// @brief    Local definition for testing whether two scalars are close
#define MACH_EPS 1e-15
#define CLOSE(a,b) ( (fabs(a)-fabs(b))/(fabs(a)) < 2*MACH_EPS)

//----------------------------------------------------------------------

void generate_input()
{
  FILE * fp = fopen ("test.in","w");

  // Groups
  // 
  //  Integer::
  //  List
  //  Logical::
  //  Logical_expr::var_logical_1
  //  Scalar::
  //  Scalar::const_scalar_1
  //  Scalar::const_scalar_2
  //  Scalar_expr::var_scalar_1
  //  Scalar_expr::var_scalar_2
  //  String

  fprintf (fp, "Logical {\n");
  fprintf (fp, "  test_true  = true;\n");
  fprintf (fp, "  test_false = false;\n");
  fprintf (fp, "}\n");
  fprintf (fp, "   \n");
  fprintf (fp, "Integer {\n");
  fprintf (fp, "  test_1 = 1;\n");
  fprintf (fp, "  test_37 = 37\n");
  fprintf (fp, "}\n");

  fprintf (fp, "Scalar {\n");
  fprintf (fp, "  test_1_5 = 1.5;\n");
  fprintf (fp, "  test_37_25 = 37.25;\n");
  fprintf (fp, "}\n");

  fprintf (fp, "Scalar {\n");
  fprintf (fp, "  const_scalar_1 {\n");
  fprintf (fp, "    num1 = 24.5 + 6.125;\n");
  fprintf (fp, "    num2 = 24.5 - 6.125;\n");
  fprintf (fp, "    num3 = 24.5 * 6.125;\n");
  fprintf (fp, "    num4 = 24.5 / 6.125;\n");
  fprintf (fp, "  };\n");
  fprintf (fp, "  const_scalar_2 {\n");
  fprintf (fp, "    num1 = 24.5 + 6.125*2.0;\n");
  fprintf (fp, "    num2 = 24.5*3.0 - 6.125;\n");
  fprintf (fp, "    num3 = (24.5 + 6.125*2.0) - (24.5*3.0 - 6.125);\n");
  fprintf (fp, "  }\n");
  fprintf (fp, "}\n");

  fprintf (fp, "String {\n");
  fprintf (fp, "  str1 = \"testing\";\n");
  fprintf (fp, "  str2 = \"one\";\n");
  fprintf (fp, "}\n");

  fprintf (fp, "Scalar_expr {\n");
  fprintf (fp, "  var_scalar_1 {\n");
  fprintf (fp, "    num1 = x;\n");
  fprintf (fp, "    num2 = x + 3.0;\n");
  fprintf (fp, "    num3 = x+y+z+t;\n");
  fprintf (fp, "  }\n");
  fprintf (fp, "}\n");

  fprintf (fp, " Scalar_expr var_scalar_2 {\n");
  fprintf (fp, "  num1 = sin(x);\n");
  fprintf (fp, "  num2 = atan(y/3.0+3.0*t);\n");
  fprintf (fp, "}\n");

  fprintf (fp, " Logical_expr {\n");
  fprintf (fp, "  var_logical_1 {\n");
  fprintf (fp, "    num1 = x < y;\n");
  fprintf (fp, "    num2 = x + y >= t + 3.0;\n");
  fprintf (fp, "    num3 = x == y;\n");
  fprintf (fp, "  }\n");
  fprintf (fp, "}\n");

  fprintf (fp, " List {\n");
  fprintf (fp, "  num1 = [1.0, true, 37, \"string\", x-y+2.0*z, x+y+t > 0.0 ];\n");
  fprintf (fp, "}\n");
  fclose(fp);

}

//======================================================================

int main()
{

  unit_init (0,1);

  unit_class("Parameters");

  //----------------------------------------------------------------------
  // test parameter
  //----------------------------------------------------------------------

  Parameters * parameters = new Parameters();

  // Generate test.in to make sure it exists

  generate_input();

  // Read

  unit_func("read");
  parameters->read ( "test.in" );

  // set_current_group()

  unit_func("set_current_group");
  parameters->set_current_group("Group");
  unit_assert(parameters->current_group() == "Group");

  // set_current_subgroup()
  unit_func("set_current_subgroup");
  parameters->set_current_subgroup("subgroup_1");
  unit_assert(parameters->current_subgroup() == "subgroup_1");
  parameters->set_current_subgroup("subgroup_2");
  unit_assert(parameters->current_subgroup() == "subgroup_2");

  // set_current_group() without set_current_subgroup() clears subgroup to ""

  unit_func("set_current_group");
  parameters->set_current_group("Group2");
  unit_assert(parameters->current_group() == "Group2");
  unit_assert(parameters->current_subgroup() == "");

  // value_logical()

  unit_func("value_logical");

  parameters->set_current_group("Logical");
  
  unit_assert (parameters->value_logical("test_true")  == true);
  unit_assert (parameters->value_logical("test_false") == false);
  unit_assert (parameters->value_logical("none",true) == true);
  unit_assert (parameters->value_logical("none",false) == false);

  // bool l,ld;

  // parameters->value("test_true",parameter_logical,&l);
  // unit_assert (l == true);

  // parameters->value("test_false",parameter_logical,&l);
  // unit_assert (l == false);

  // ld = true;
  // parameters->value("none",parameter_logical,&l,&ld);
  // unit_assert (l == ld);

  // ld = false;
  // parameters->value("none",parameter_logical,&l,&ld);
  // unit_assert (l == ld);

  // set_logical()

  unit_func("set_logical");

  parameters->set_logical("test_true",false);
  unit_assert (parameters->value_logical("test_true") == false);
  parameters->set_logical("none",true);
  unit_assert (parameters->value_logical("none") == true);
  parameters->set_logical("none_l1",true);
  unit_assert (parameters->value_logical("none_l1") == true);
  parameters->set_logical("none_l2",false);
  unit_assert (parameters->value_logical("none_l2") == false);

  // value_integer()

  unit_func("value_integer");

  parameters->set_current_group("Integer");
  
  unit_assert (parameters->value_integer("test_1")  == 1);
  unit_assert (parameters->value_integer("test_37") == 37);
  unit_assert (parameters->value_integer("none",58) == 58);

  // int i,id;
  // parameters->value("test_1",parameter_integer,&i);
  // unit_assert (i == 1);
  // parameters->value("test_37",parameter_integer,&i);
  // unit_assert (i == 37);
  // id = 58;
  // parameters->value("none",parameter_integer,&i,&id);
  // unit_assert (i == id);

  // set_integer()

  unit_func("set_integer");

  parameters->set_integer("test_1",2);
  unit_assert (parameters->value_integer("test_1") == 2);
  parameters->set_integer("none",3);
  unit_assert (parameters->value_integer("none") == 3);
  parameters->set_integer("none2",4);
  unit_assert (parameters->value_integer("none2") == 4);


  // value_scalar()
  
  unit_func("value_scalar");

  parameters->set_current_group("Scalar");
  
  unit_assert (parameters->value_scalar("test_1_5")  == 1.5);
  unit_assert (parameters->value_scalar("test_37_25") == 37.25);
  unit_assert (parameters->value_scalar("none",58.75) == 58.75);

  // double d,dd;
  
  // parameters->value("test_1_5",parameter_scalar,&d);
  // unit_assert (d  == 1.5);
  // parameters->value("test_37_25",parameter_scalar,&d);
  // unit_assert (d == 37.25);
  // dd = 58.75;
  // parameters->value("none",parameter_scalar,&d,&dd);
  // unit_assert (d == dd);

  // set_scalar()

  unit_func("set_scalar");

  parameters->set_scalar("test_1_5",27.0);
  unit_assert (parameters->value_scalar("test_1_5") == 27.0);
  parameters->set_scalar("none_s",1.5);
  unit_assert (parameters->value_scalar("none_s") == 1.5);

  // Constant scalar expressions
  // subgroups

  unit_func("value_scalar");

  parameters->set_current_group("Scalar");

  parameters->set_current_subgroup("const_scalar_1");

  unit_assert(parameters->value_scalar("num1") == 30.625);
  unit_assert(parameters->value_scalar("num2") == 18.375);
  unit_assert(parameters->value_scalar("num3") == 150.0625);
  unit_assert(parameters->value_scalar("num4") == 4.0000000000);

  parameters->set_current_subgroup("const_scalar_2");

  unit_assert(parameters->value_scalar("num1") == 36.750);
  unit_assert(parameters->value_scalar("num2") == 67.375);
  unit_assert(parameters->value_scalar("num3") == -30.625);

  // Strings

  unit_func("value_string");

  parameters->set_current_group("String");
  unit_assert(strcmp(parameters->value_string("str1"),"testing")==0);
  unit_assert(strcmp(parameters->value_string("str2","blah"),"one")==0);
  unit_assert(strcmp(parameters->value_string("none","blah"),"blah")==0);

  // const char *s, *sd = "blah";
  // parameters->value("str1",parameter_string,&s);
  // unit_assert(strcmp(s,"testing")==0);
  // parameters->value("str2",parameter_string,&s,&sd);
  // unit_assert(strcmp(s,"one")==0);
  // parameters->value("none",parameter_string,&s,&sd);
  // unit_assert(strcmp(s,"blah")==0);

  // set_string()

  unit_func("set_string");

  parameters->set_string("str1","yahoo");
  unit_assert (strcmp(parameters->value_string("str1"),"yahoo")==0);
  parameters->set_string("none_str","hello");
  unit_assert (strcmp(parameters->value_string("none_str"),"hello")==0);

  // Variable scalar expressions

  unit_func("evaluate_scalar");

  double x[] = { 1, 2, 3};
  double y[] = {5 , 4, 3};
  double z[] = {8, 9, 10};
  double t[] = {-1, 2, -7};
  double values_scalar[] = {0,0,0};
  double deflts_scalar[] = {-1,-2,-3};

  parameters->set_current_group("Scalar_expr");

  parameters->set_current_subgroup("var_scalar_1");

  parameters->evaluate_scalar("num1",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0]==x[0]);
  unit_assert (values_scalar[1]==x[1]);
  unit_assert (values_scalar[2]==x[2]);

  
  parameters->evaluate_scalar("num2",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0]==x[0]+3.0);
  unit_assert (values_scalar[1]==x[1]+3.0);
  unit_assert (values_scalar[2]==x[2]+3.0);

  parameters->evaluate_scalar("num3",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0]==x[0]+y[0]+z[0]+t[0]);
  unit_assert (values_scalar[1]==x[1]+y[1]+z[1]+t[1]);
  unit_assert (values_scalar[2]==x[2]+y[2]+z[2]+t[2]);

  parameters->set_current_subgroup("var_scalar_2");

  parameters->evaluate_scalar("num1",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (CLOSE(values_scalar[0],sin(x[0])));
  unit_assert (CLOSE(values_scalar[1],sin(x[1])));
  unit_assert (CLOSE(values_scalar[2],sin(x[2])));

  parameters->evaluate_scalar("num2",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (CLOSE(values_scalar[0],atan(y[0]/3.0+3*t[0])));
  unit_assert (CLOSE(values_scalar[1],atan(y[1]/3.0+3*t[1])));
  unit_assert (CLOSE(values_scalar[2],atan(y[2]/3.0+3*t[2])));

  // Logical expressions

  unit_func("evaluate_logical");

  bool values_logical[] = {false, false, false};
  bool deflts_logical[] = {true, false,true};

  parameters->set_current_group("Logical_expr");
  parameters->set_current_subgroup("var_logical_1");

  parameters->evaluate_logical("num1",3,values_logical,deflts_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0] < y[0]));
  unit_assert (values_logical[1] == (x[1] < y[1]));
  unit_assert (values_logical[2] == (x[2] < y[2]));

  parameters->evaluate_logical("num2",3,values_logical,deflts_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0] + y[0] >= t[0] + 3.0));
  unit_assert (values_logical[1] == (x[1] + y[1] >= t[1] + 3.0));
  unit_assert (values_logical[2] == (x[2] + y[2] >= t[2] + 3.0));

  parameters->evaluate_logical("num3",3,values_logical,deflts_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0] == y[0]));
  unit_assert (values_logical[1] == (x[1] == y[1]));
  unit_assert (values_logical[2] == (x[2] == y[2]));

  // Lists

  parameters->set_current_group("List");

  unit_func("list_length");
  unit_assert(parameters->list_length("num1") == 6);

  unit_func("list_value_scalar");
  unit_assert(parameters->list_value_scalar(0,"num1") == 1.0);

  unit_func("list_value_logical");
  unit_assert(parameters->list_value_logical(1,"num1") == true);

  unit_func("list_value_integer");
  unit_assert(parameters->list_value_integer(2,"num1") == 37);

  unit_func("list_value_string");
  unit_assert(strcmp(parameters->list_value_string(3,"num1"),"string")==0);

  unit_func("list_evaluate_scalar");
  parameters->list_evaluate_scalar(4,"num1",3,values_scalar,deflts_scalar,x,y,z,t);
  unit_assert (values_scalar[0] == (x[0]-y[0]+2.0*z[0]));
  unit_assert (values_scalar[1] == (x[1]-y[1]+2.0*z[1]));
  unit_assert (values_scalar[2] == (x[2]-y[2]+2.0*z[2]));

  unit_func("list_evaluate_logical");
  parameters->list_evaluate_logical(5,"num1",3,values_logical,deflts_logical,x,y,z,t);
  unit_assert (values_logical[0] == (x[0]+y[0]+t[0] > 0 ));
  unit_assert (values_logical[1] == (x[1]+y[1]+t[1] > 0 ));
  unit_assert (values_logical[2] == (x[2]+y[2]+t[2] > 0 ));

  // group_count(), group(i)

  unit_func("group_count");

  int num_groups = parameters->group_count();
  // Groups:
  //  Integer
  //  List
  //  Logical
  //  Logical_expr
  //  Scalar
  //  Scalar_expr
  //  String
  unit_assert (num_groups == 7); 

  unit_func("subgroup_count");
  std::map <std::string,int> num_subgroups;
  num_subgroups["Integer"]      = 1; // ""
  num_subgroups["List"]         = 1; // ""
  num_subgroups["Logical"]      = 1; // "var_logical_1"
  num_subgroups["Logical_expr"] = 1; // ""
  num_subgroups["Scalar"]       = 3; // "","const_scalar_1","const_scalar_2"
  num_subgroups["Scalar_expr"]  = 2; // "var_scalar_1","var_scalar_2"
  num_subgroups["String"]       = 1; // ""

  for (int i=0; i<num_groups; i++) {
    std::string group = parameters->group(i);
    parameters->set_current_group(group);
    unit_assert (num_subgroups[group] == parameters->subgroup_count());
  }
  
  // Write

  unit_func("write");
  parameters->write ( "test.out" );
  unit_assert(unit_incomplete);

  unit_finalize();

}

