#include <iostream>
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

using namespace libMesh;

int main (int argc, char ** argv)
{

  // En este programa se implementa una cuadricula de 5x5 (se explica en el informe) y se carga la informacion de una malla en un archivo xdr

  LibMeshInit init (argc, argv);

  if (argc < 4)
    {
      libmesh_error_msg("Usage: " << argv[0] << " -d 2 in.mesh [-o out.mesh]");
    }

  const unsigned int dim = std::atoi(argv[2]);

  libmesh_example_requires(dim <= LIBMESH_DIM, "2D/3D support");

  Mesh mesh(init.comm());

  std::string input_filename = argv[3];
#ifndef LIBMESH_HAVE_XDR
  libmesh_example_requires(input_filename.rfind(".xdr") >=
                           input_filename.size(), "XDR support");
#endif

  mesh.read (argv[3]);

  mesh.print_info();

  if (argc >= 6 && std::string("-o") == argv[4])
    {
      std::string output_filename = argv[5];
#ifndef LIBMESH_HAVE_XDR
      libmesh_example_requires(output_filename.rfind(".xdr") >=
                               output_filename.size(), "XDR support");
#endif

      mesh.write (argv[5]);
    }
  return 0;
}
