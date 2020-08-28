
import sys,os,unittest
from lammps import lammps

has_mpi=False
has_mpi4py=False
try:
    from mpi4py import __version__ as mpi4py_version
    # tested to work with mpi4py versions 2 and 3
    has_mpi4py = mpi4py_version.split('.')[0] in ['2','3']
except:
    pass

try:
    lmp = lammps()
    has_mpi = lmp.has_mpi_support
    lmp.close()
except:
    pass

class PythonOpen(unittest.TestCase):

    def setUp(self):
        self.machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            self.machine=os.environ['LAMMPS_MACHINE_NAME']

    def testNoArgs(self):
        """Create LAMMPS instance without any arguments"""

        lmp=lammps(name=self.machine)
        self.assertIsNot(lmp.lmp,None)
        self.assertEqual(lmp.opened,1)
        self.assertEqual(has_mpi4py,lmp.has_mpi4py)
        self.assertEqual(has_mpi,lmp.has_mpi_support)
        lmp.close()
        self.assertIsNone(lmp.lmp,None)
        self.assertEqual(lmp.opened,0)

    def testWithArgs(self):
        """Create LAMMPS instance with a few arguments"""
        lmp=lammps(name=self.machine,
                   cmdargs=['-nocite','-sf','opt','-log','none'])
        self.assertIsNot(lmp.lmp,None)
        self.assertEqual(lmp.opened,1)

    @unittest.skipIf(not (has_mpi and has_mpi4py),"Skipping MPI test since LAMMPS is not parallel or mpi4py is not found")
    def testWithMPI(self):
        from mpi4py import MPI
        mycomm=MPI.Comm.Split(MPI.COMM_WORLD, 0, 1)
        lmp=lammps(name=self.machine,comm=mycomm)
        self.assertIsNot(lmp.lmp,None)
        self.assertEqual(lmp.opened,1)
        lmp.close()

if __name__ == "__main__":
    unittest.main()
