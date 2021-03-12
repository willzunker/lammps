
#include "lammpsplugin.h"

#include "comm.h"
#include "error.h"
#include "pointers.h"
#include "version.h"

#include <cstring>

namespace LAMMPS_NS {
  class Hello : protected Pointers {
    public:
      Hello(class LAMMPS *lmp) : Pointers(lmp) {};
      void command(int, char **);
  };
}

using namespace LAMMPS_NS;

void Hello::command(int argc, char **argv)
{
   if (argc != 1) error->all(FLERR,"Illegal hello command");
   if (comm->me == 0)
     utils::logmesg(lmp,fmt::format("Hello, {}!\n",argv[0]));
}

static void hellocreator(LAMMPS *lmp, int argc, char **argv)
{
    Hello hello(lmp);
    hello.command(argc,argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style   = "command";
  plugin.name    = "hello";
  plugin.info    = "Hello world command v1.0";
  plugin.author  = "Axel Kohlmeyer (akohlmey@gmail.com)";
  plugin.creator.v3 = (lammpsplugin_factory3 *) &hellocreator;
  plugin.handle  = handle;
  (*register_plugin)(&plugin,lmp);
}
