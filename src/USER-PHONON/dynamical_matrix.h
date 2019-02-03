//
// Created by charlie sievers on 6/21/18.
//

#ifdef COMMAND_CLASS

CommandStyle(dynamical_matrix,DynamicalMatrix)

#else

#ifndef LMP_DYNAMICAL_MATRIX_H
#define LMP_DYNAMICAL_MATRIX_H

#include "pointers.h"

namespace LAMMPS_NS {

    class DynamicalMatrix : protected Pointers {
    public:
        DynamicalMatrix(class LAMMPS *);
        virtual ~DynamicalMatrix();
        void command(int, char **);
        void setup();

    protected:
        int eflag,vflag;            // flags for energy/virial computation
        int external_force_clear;   // clear forces locally or externally


        int triclinic;              // 0 if domain is orthog, 1 if triclinic
        int pairflag;

        int pair_compute_flag;            // 0 if pair->compute is skipped
        int kspace_compute_flag;          // 0 if kspace->compute is skipped

        int nvec;                   // local atomic dof = length of xvec

        void energy_force(int);
        void force_clear();
        virtual void openfile(const char* filename);

    private:
        void options(int, char **);
        void calculateMatrix();
        void create_groupmap();
        void writeMatrix();
        void convert_units(const char *style);
        void displace_atom(int local_idx, int direction, int magnitude);

        double conversion;
        double conv_energy;
        double conv_distance;
        double conv_mass;
        double del;
        int igroup,groupbit;
        int scaleflag;
        int me;
        bigint dynlen;
        int *groupmap;
        double **dynmat;
        double **final_dynmat;

        int compressed;            // 1 if dump file is written compressed, 0 no
        int binaryflag;            // 1 if dump file is written binary, 0 no
        int file_opened;           // 1 if openfile method has been called, 0 no
        int file_flag;             // 1 custom file name, 0 dynmat.dat

        FILE *fp;
    };
}


#endif //LMP_DYNAMICAL_MATRIX_H
#endif