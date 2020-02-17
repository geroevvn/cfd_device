/*
 * Fvm_tvd_implicit.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef FVM_TVD_IMPLICIT_H_
#define FVM_TVD_IMPLICIT_H_

#include "../mesh/Grid.h"
#include "../mesh/Mesh.h"
#include "../iterators/MeshIterator.h"
#include "../iterators/FilterIterator.h"
#include "../iterators/BndIterator.h"
#include "Method.h"
#include "MatrixSolver/MatrixSolver.h"
#include <vector>


using namespace std;

class FVM_TVD_IMPLICIT: public Method
{
public:
		FVM_TVD_IMPLICIT();
		~FVM_TVD_IMPLICIT();

		virtual void init(char* xmlFileName);
		virtual void run();
		virtual void done();

private:
		static const int PLUS_JACOBIAN = 0;
		static const int MINUS_JACOBIAN = 1;

		double			TMAX;
		int				STEP_MAX;
		int 			FILE_STEP_SAVE;
		int 			LOG_STEP_SAVE;
		double			TAU;
		double			CFL;

		double* Flux;
		double* Flux1;
		double* Flux2;
		double** temp_mat;

private:
        Grid* grid;
        Mesh* msh;

        vector<string> bndInletNames;
        vector<string> bndOutletNames;
        vector<string> bndWallNames;


        int check_bnd_cond();
        void clear5(double**);
        void clear_vec(double*);
        double** allocate_mem();
        void free_mem(double**);
        void matrix_A(double**, double**, double*, double**, int);
        void matrix_Ap(double**, double**, double*, double**);
        void matrix_Am(double**, double**, double*, double**);
        void eigen_values(double*, double, double, double, double, const Point&);
        void left_eigen_vecs(double**, double, double, double, double, double, const Point&);
        void right_eigen_vecs(double**, double, double, double, double, double, const Point&);
        void calc_F(double*, const CellFluidDynamicsProps&);
        void calc_H(double*, const CellFluidDynamicsProps&);
        void calc_G(double*, const CellFluidDynamicsProps&);
        void flux_Lax_Friedrichs(double*, const CellFluidDynamicsProps&, const CellFluidDynamicsProps&, const Point&);
        void save(int);

        inline double _max(double, double);
};

#endif /* FVM_TVD_IMPLICIT_H_ */
