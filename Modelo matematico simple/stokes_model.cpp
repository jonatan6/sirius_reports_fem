#include <iostream>
#include <algorithm>
#include <math.h>
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/elem.h"
using namespace libMesh;
void assemble_stokes (EquationSystems & es,
                      const std::string & system_name);
int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");
  libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc or --enable-laspack");
  Mesh mesh(init.comm());
  MeshTools::Generation::build_square (mesh,
                                       15, 15,
                                       0., 1.,
                                       0., 1.,
                                       QUAD9);
  mesh.print_info();
  EquationSystems equation_systems (mesh);
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Stokes");
  system.add_variable ("u", SECOND);
  system.add_variable ("v", SECOND);
  system.add_variable ("p", FIRST);
  system.attach_assemble_function (assemble_stokes);
  equation_systems.init ();
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 250;
  equation_systems.parameters.set<Real>        ("linear solver tolerance") = TOLERANCE;
  equation_systems.print_info();
  equation_systems.get_system("Stokes").solve();
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems ("out.e",
                                            equation_systems);
#endif 
  return 0;
}
void assemble_stokes (EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to (system_name, "Stokes");
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem> ("Stokes");
  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int p_var = system.variable_number ("p");
  FEType fe_vel_type = system.variable_type(u_var);
  FEType fe_pres_type = system.variable_type(p_var);
  UniquePtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
  UniquePtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  const std::vector<Real> & JxW = fe_vel->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe_vel->get_dphi();
  const std::vector<std::vector<Real> > & psi = fe_pres->get_phi();
  const DofMap & dof_map = system.get_dof_map();
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke);
  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_p;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);
      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_p_dofs; j++)
              Kup(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](0);
          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_p_dofs; j++)
              Kvp(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](1);
          for (unsigned int i=0; i<n_p_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              Kpu(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](0);
          for (unsigned int i=0; i<n_p_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              Kpv(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](1);
        } 
      {
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor_ptr(s) == libmesh_nullptr)
            {
              UniquePtr<const Elem> side (elem->build_side_ptr(s));
              for (unsigned int ns=0; ns<side->n_nodes(); ns++)
                {
                  const Real yf = side->point(ns)(1);
                  const Real penalty = 1.e10;
                  const Real u_value = (yf > .99) ? 1. : 0.;
                  const Real v_value = 0.;
                  for (unsigned int n=0; n<elem->n_nodes(); n++)
                    if (elem->node_id(n) == side->node_id(ns))
                      {
                        Kuu(n,n) += penalty;
                        Kvv(n,n) += penalty;
                        Fu(n) += penalty*u_value;
                        Fv(n) += penalty*v_value;
                      }
                } 
            } 
      } 
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    } 
}
