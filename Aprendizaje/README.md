# USO DE LIBMESH

CONFIGURACION DE GPU

Actualmente un requisito fundamental de los frameworks y aplicaciones en general, es poseer mecanismos que permitan el aprovechamiento de GPU y acelerar el desempeño computacional. Libmesh se apoya en la librería PETSC que permite ejecutar solvers para algebra lineal, además de ello, PETSC permite la paralización con GPU mediante solvers que usan formas matriciales, lo cual explota al máximo las capacidades del hardware. Para este propósito se pueden aprovechar dos tecnologías: CUDA y OpenCL.

En las instrucciones a continuacion, todos los binarios que están disponibles para sistemas Windows son autoinstalables, es decir que basta con seguir los pasos del asistente para llevar a cabo la correcta instalación de los paquetes. En Linux, a menos que se mencione lo contrario, se debe contar con el paquete _&#39;build-essentials&#39;_ (disponible en casi todos los gestores de paquetes para sistemas Linux), para realizar la compilación ___./configure__, __make__, __make install__. Solo algunos de los paquetes (como los de CUDA), están disponibles en los gestores de paquetes como cabeceras genéricas.

## CUDA

Es una arquitectura de cálculo paralelo de NVIDIA que aprovecha la gran potencia de la GPU (unidad de procesamiento gráfico) para proporcionar un incremento extraordinario del rendimiento del sistema. Los sistemas informáticos están pasando de realizar el &quot;procesamiento central&quot; en la CPU a realizar &quot;procesamiento&quot; repartido entre la CPU y la GPU. Para posibilitar este nuevo paradigma computacional, NVIDIA ha inventado la arquitectura de cálculo paralelo CUDA, que ahora se incluye en las GPUs GeForce, ION Quadro y Tesla GPUs, lo cual representa una base instalada considerable para los desarrolladores de aplicaciones.

Con base en lo anterior, sólo los dispositivos NVIDIA anteriormente mencionados pueden realizar aceleración de algoritmos a través de CUDA. Para ello, se debe contar con:

- --La suite de desarrollo CUDA que se puede obtener del siguiente vinculo: [https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads). Hay versiones de CUDA para sistemas operativos Windows, Mac OSX y Linux. Minimo versión 4.0
- --Thrust, el conjunto de algoritmos paralelos que reensamblan las librerías estándar de C++ que se pueden obtener en: [http://thrust.github.io/](http://thrust.github.io/)
- --La librería para algebra lineal dispersa y grafos computacionales, Cusp. Puede obtenerse en: [http://cusplibrary.github.io/](http://cusplibrary.github.io/)

Se deben desactivar los controladores gráficos Nouveau para la tarjeta de video y remplazarlos por sus correspondientes versiones de NVIDIA, actualmente 369. Este paso es importante pues el controlador Nouveau (desarrollado de manera libre por la Xorg fundation) no soporta todos los componentes necesarios para la suite CUDA. Dentro de la carpeta config/examples/ hay un archivo .py que al ejecutarse configura en PETSC las banderas necesarias para linkear las librerías respectivas de CUDA.

## OpenCL

Es un framework que consta de una interfaz de programación de aplicaciones y de un lenguaje de programación. Juntos permiten crear aplicaciones con paralelismo a nivel de datos y de tareas que pueden ejecutarse tanto en unidades centrales de procesamiento como unidades de procesamiento gráfico. El lenguaje está basado en C99, eliminando cierta funcionalidad y extendiéndolo con operaciones vectoriales. Los siguientes requerimientos deben cumplirse para el uso de OpenCL en PETSc:

- --ViennaCL: Se trata de una librería de algebra lineal orientada al cómputo multinúcleo en CPU, GPU y MIC, basado en CUDA, OpenCL y OpenMP. Es capaz de funcionar según los requerimientos instalados y provee una capa de abstracción que permite el uso simultaneo de las tecnologías mencionadas. Puede obtenerse de: [http://viennacl.sourceforge.net/viennacl-download.html](http://viennacl.sourceforge.net/viennacl-download.html), hay binarios para sistemas Windows, y tarballs compilables para linux
- --Los headers requeridos para cada una de las arquitecturas pueden obtenerse para:
  - NVIDIA: [https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads) (incluido en la suite CUDA)
  - AMD: [http://developer.amd.com/tools-and-sdks/opencl-zone/amd-accelerated-parallel-processing-app-sdk/](http://developer.amd.com/tools-and-sdks/opencl-zone/amd-accelerated-parallel-processing-app-sdk/) (incluido en el kit SDK para AMD)
  - Intel: [https://software.intel.com/en-us/intel-inde](https://software.intel.com/en-us/intel-inde) (Incluido en la suite INDE)

Las cabeceras genéricas en sistemas Linux puede obtenerse bajo el nombre _&#39;opencl-headers&#39;_ desde la mayoría de gestores de paquetes.

- --Tener instalado los drivers dedicados para cada modelo de GPU (ya mencionados previamente).
- --En los archivos de instalación hay un archivo .py que configura las banderas necesarias para linkear las librerías.

PETSc usa estructuras de datos vectoriales para abstraer operaciones paralelas out-of-box, es decir que operaciones vectoriales básicas con las estructuras de proporcionadas por el framework son realizadas en GPU sin necesidad de implementarlas manualmente. Sin embargo, es posible usar los paquetes mencionados anteriormente para especificar en detalle la manera en que deben ejecutarse las operaciones paralelas deseadas.

Además, es posible combinar los paquetes de MPI y métodos de memoria compartida para consolidar una estructura de desarrollo más compleja y con altos niveles de paralelismo. Cabe mencionar que la optimización mediante paralelismo es difícil de realizar; las operaciones pueden distribuirse, pero ello no asegura la reducción de tiempos de ejecución con respecto a implementaciones secuenciales que hagan uso de GPU. La filosofía principal consta en poder compensar con tiempo de computo, el costo que implica transferir los datos desde la memoria host al device y viceversa; una vez esta barrera es superada, la eficiencia del algoritmo paralelo supera a la de la implementación secuencial.
