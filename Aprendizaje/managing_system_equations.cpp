#include <iostream>
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"

using namespace libMesh;



int main (int argc, char ** argv)
{
  LibMeshInit init (argc, argv);

  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  libMesh::out << "Running " << argv[0];
  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];
  libMesh::out << std::endl << std::endl;

  Mesh mesh(init.comm());

  MeshTools::Generation::build_square (mesh, 5, 5);

  EquationSystems equation_systems (mesh);

  equation_systems.parameters.set<bool> ("test") = true;

  equation_systems.parameters.set<Real> ("dummy") = 42.;

  equation_systems.parameters.set<Real> ("nobody") = 0.;

  equation_systems.add_system<TransientLinearImplicitSystem> ("Simple System");

  equation_systems.get_system("Simple System").add_variable("u", FIRST);

  equation_systems.add_system<ExplicitSystem> ("Complex System");

  equation_systems.get_system("Complex System").add_variable("c", FIRST);
  equation_systems.get_system("Complex System").add_variable("T", FIRST);
  equation_systems.get_system("Complex System").add_variable("dv", SECOND, MONOMIAL);

  equation_systems.init();

  mesh.print_info();

  equation_systems.print_info();

  if (argc > 1)
    if (argv[1][0] != '-')
      {
        libMesh::out << "<<< Writing system to file " << argv[1]
                     << std::endl;

        equation_systems.write (argv[1], WRITE);

        equation_systems.clear ();

        libMesh::out << ">>> Reading system from file " << argv[1]
                     << std::endl << std::endl;

        equation_systems.read (argv[1], READ);

        equation_systems.print_info();
      }

  return 0;
}
