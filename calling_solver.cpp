#include <iostream>
#include <algorithm>
#include <math.h>

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/fe.h"

#include "libmesh/quadrature_gauss.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"

#include "libmesh/dof_map.h"

using namespace libMesh;

void assemble_poisson(EquationSystems & es,
                      const std::string & system_name);

Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.);

int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  Mesh mesh(init.comm());

  MeshTools::Generation::build_square (mesh,
                                       15, 15,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD9);

  mesh.print_info();

  EquationSystems equation_systems (mesh);

  equation_systems.add_system<LinearImplicitSystem> ("Poisson");

  equation_systems.get_system("Poisson").add_variable("u", SECOND);

  equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);

  equation_systems.init();

  equation_systems.print_info();

  equation_systems.get_system("Poisson").solve();

#if defined(LIBMESH_HAVE_VTK) && !defined(LIBMESH_ENABLE_PARMESH)

  VTKIO (mesh).write_equation_systems ("out.pvtu", equation_systems);


  return 0;
}



void assemble_poisson(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{

  libmesh_assert_equal_to (system_name, "Poisson");

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");

  const DofMap & dof_map = system.get_dof_map();

  FEType fe_type = dof_map.variable_type(0);

  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

  QGauss qrule (dim, FIFTH);

  fe->attach_quadrature_rule (&qrule);

  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  QGauss qface(dim-1, FIFTH);

  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();

  const std::vector<Point> & q_point = fe->get_xyz();

  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el ; ++el)
    {
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);

      fe->reinit (elem);


      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          for (unsigned int i=0; i<phi.size(); i++)
            for (unsigned int j=0; j<phi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              }

          {
            const Real x = q_point[qp](0);
            const Real y = q_point[qp](1);
            const Real eps = 1.e-3;


            const Real fxy = -(exact_solution(x, y-eps) +
                               exact_solution(x, y+eps) +
                               exact_solution(x-eps, y) +
                               exact_solution(x+eps, y) -
                               4.*exact_solution(x, y))/eps/eps;

            for (unsigned int i=0; i<phi.size(); i++)
              Fe(i) += JxW[qp]*fxy*phi[i][qp];
          }
        }

      {

        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              const std::vector<std::vector<Real> > & phi_face = fe_face->get_phi();

              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              const std::vector<Point> & qface_point = fe_face->get_xyz();

              fe_face->reinit(elem, side);

              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  const Real xf = qface_point[qp](0);
                  const Real yf = qface_point[qp](1);

                  const Real penalty = 1.e10;

                  const Real value = exact_solution(xf, yf);

                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

                  for (unsigned int i=0; i<phi_face.size(); i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                }
            }
      }


      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

}
