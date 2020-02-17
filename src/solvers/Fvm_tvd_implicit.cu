/*
 * Fvm_tvd_implicit.cpp
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#include "Fvm_tvd_implicit.h"
#include "../tinyxml/tinyxml.h"
#include "../mesh_properties/CellFluidDynamicsProps.h"

#include <string.h>
#include <vector>
#include <algorithm>

#include <cuda.h>
#include <cuda_runtime.h>
#include <unistd.h>

#include <ctime>

#define INLET_TYPE 4
#define OUTLET_TYPE 1
#define WALL_TYPE 2


#define TETRA_FACES_CNT 4
#define INLET_BND_SIZE 50

__constant__ float* TAU_dev;
__constant__ float inlet_ro_dev[INLET_BND_SIZE];
__constant__ float inlet_ru_dev[INLET_BND_SIZE];
__constant__ float inlet_rv_dev[INLET_BND_SIZE];
__constant__ float inlet_rw_dev[INLET_BND_SIZE];
__constant__ float inlet_rE_dev[INLET_BND_SIZE];
__constant__ float inlet_P_dev[INLET_BND_SIZE];
__constant__ float inlet_gamma_dev[INLET_BND_SIZE];



double FVM_TVD_IMPLICIT::_max(double x, double y)
{
	return (x > y) ? x : y;
}


FVM_TVD_IMPLICIT::FVM_TVD_IMPLICIT()
{
	grid = 0;

	Flux = new double[5];
	Flux1 = new double[5];
	Flux2 = new double[5];
}

FVM_TVD_IMPLICIT::~FVM_TVD_IMPLICIT()
{
	//cout << "FVM_TVD_IMPLICIT" << endl;
	delete [] Flux;
	delete [] Flux1;
	delete [] Flux2;

	if(grid != 0)
	{
		delete grid;
	}
}



void FVM_TVD_IMPLICIT::init(char* xmlFileName)
{
	TiXmlDocument doc( xmlFileName );
	bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
	if (!loadOkay)
	{
		Logger::Instance()->logging()->error("Failed to open file : \"%s\"", xmlFileName);
		Logger::Instance()->EXIT(doc.ErrorId());
	}

	double ro, u, v, w, P, gamma;

	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild( "task" );


	node0 = task->FirstChild("mesh");
	const char* fileType = task->FirstChild("mesh")->FirstChild("fileType")->ToElement()->Attribute("value");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");

	if(fileType == 0)
	{
		Logger::Instance()->logging()->error("Filetype of Mesh error");
		Logger::Instance()->EXIT(-1);
	}

	if(fName == 0)
	{
		Logger::Instance()->logging()->error("Filename of Mesh error");
		Logger::Instance()->EXIT(-1);
	}

	grid = new Grid(fileType);
	grid->read(fName);

	msh = grid->get_mesh();

	int steadyVal = 1;
	node0 = task->FirstChild("control");
	//node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_STEP_SAVE);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &LOG_STEP_SAVE);

	/*
	const char * flxStr = node0->FirstChild("FLUX")->ToElement()->Attribute("value");
	if (strcmp(flxStr, "GODUNOV") == 0) {
		FLUX = FLUX_GODUNOV;
	}
	else if (strcmp(flxStr, "LAX") == 0) {
		FLUX = FLUX_LAX;
	}
	else {
		FLUX = FLUX_GODUNOV;
	}

	if (steadyVal == 0) {
		STEADY = false;
	} else {
		STEADY = true;
		node1 = node0->FirstChild("CFL");
		node1->FirstChild("start")->ToElement()->Attribute("value", &CFL);
		node1->FirstChild("scale")->ToElement()->Attribute("value", &scaleCFL);
		node1->FirstChild("max")->ToElement()->Attribute("value", &maxCFL);
		node1->FirstChild("step")->ToElement()->Attribute("value", &stepCFL);
		node1->FirstChild("max_limited_cells")->ToElement()->Attribute("value", &maxLimCells);
	}


	int smUsing = 1;
	node0 = task->FirstChild("smoothing");
	node0->FirstChild("using")->ToElement()->Attribute("value", &smUsing);
	node0->FirstChild("coefficient")->ToElement()->Attribute("value", &SMOOTHING_PAR);
	SMOOTHING = (smUsing == 1);


	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute( "min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute( "max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute( "max", &limitUmax);


	node0 = task->FirstChild("materials");
	node0->ToElement()->Attribute("count", &matCount);;
	materials = new Material[matCount];
	TiXmlNode* matNode = node0->FirstChild("material");
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		matNode->ToElement()->Attribute("id", &mat.id);
		node1 = matNode->FirstChild("name");
		el = node1->ToElement();
		mat.name = el->GetText();
		node1 = matNode->FirstChild("parameters");
		node1->FirstChild( "M"  )->ToElement()->Attribute( "value", &mat.M  );
		node1->FirstChild( "Cp" )->ToElement()->Attribute( "value", &mat.Cp );
		node1->FirstChild( "K"  )->ToElement()->Attribute( "value", &mat.K  );
		node1->FirstChild( "ML" )->ToElement()->Attribute( "value", &mat.ML );
		matNode = matNode->NextSibling("material");
	}
	*/

	node0 = task->FirstChild("regions");
	int regCount;
	node0->ToElement()->Attribute("count", &regCount);

	TiXmlNode* regNode = node0->FirstChild("region");
	for (int i = 0; i < regCount; i++)
	{
		node1 = regNode->FirstChild("parameters");

		node1->FirstChild( "ro" )->ToElement()->Attribute( "value", &ro );
		node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &u );
		node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &v );
		node1->FirstChild( "Vz" )->ToElement()->Attribute( "value", &w );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &P );
		node1->FirstChild( "Gamma"  )->ToElement()->Attribute( "value", &gamma );

		for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
		{
			it->cellFDP.ro = ro;
			it->cellFDP.ru = u * ro;
			it->cellFDP.rv = v * ro;
			it->cellFDP.rw = w * ro;
			it->cellFDP.gamma = gamma;
			it->cellFDP.P = P;

			it->cellFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, gamma);
		}

		regNode = regNode->NextSibling("region");
	}


	node0 = task->FirstChild("boundaries");
	double bCount;
	node0->ToElement()->Attribute("count", &bCount);
	TiXmlNode* bNode = node0->FirstChild("boundCond");

	for (int i = 0; i < bCount; i++)
	{
		const char * name = bNode->FirstChild("name")->ToElement()->GetText();
		const char * str = bNode->FirstChild("type")->ToElement()->GetText();

		if (strcmp(str, "BOUND_WALL") == 0)
		{
			bndWallNames.push_back(name);
			for (Mesh::FaceIterator it = msh->beginBndFace(name), ite = msh->endBndFace(name); it != ite; ++it)
			{
				it->bnd_type = Face::BND_TYPE_WALL;
			}
		}
		else if (strcmp(str, "BOUND_OUTLET") == 0)
		{
			bndOutletNames.push_back(name);
			for (Mesh::FaceIterator it = msh->beginBndFace(name), ite = msh->endBndFace(name); it != ite; ++it)
			{
				it->bnd_type = Face::BND_TYPE_OUTLET;
			}
		}
		else if (strcmp(str, "BOUND_INLET") == 0)
		{
			bndInletNames.push_back(name);

			node1 = bNode->FirstChild("parameters");

			node1->FirstChild( "ro" )->ToElement()->Attribute( "value", &ro );
			node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &u );
			node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &v );
			node1->FirstChild( "Vz" )->ToElement()->Attribute( "value", &w );
			node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &P );
			node1->FirstChild( "Gamma"  )->ToElement()->Attribute( "value", &gamma );

			for (Mesh::FaceIterator it = msh->beginBndFace(name), ite = msh->endBndFace(name); it != ite; ++it)
			{
				it->bnd_type = Face::BND_TYPE_INLET;
				it->faceFDP.ro = ro;
				it->faceFDP.ru = u * ro;
				it->faceFDP.rv = v * ro;
				it->faceFDP.rw = w * ro;
				it->faceFDP.gamma = gamma;
				it->faceFDP.P = P;

				it->faceFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, gamma);
			}
		}
		else
		{
			Logger::Instance()->logging()->error("Unsupported boundary condition type \"%s\"", str);
			Logger::Instance()->EXIT(1);
		}

		bNode = bNode->NextSibling("boundCond");
	}

	bool check = check_bnd_cond();

	if(!check)
	{
		Logger::Instance()->logging()->error("Boundary names from \"%s\" != boundary names from \"%s\"", xmlFileName, fName);
		Logger::Instance()->EXIT(1);
	}

	save(0);
}

int FVM_TVD_IMPLICIT::check_bnd_cond()
{
	vector<string> v1, v2;

	for(int i = 0; i < bndInletNames.size(); i++)
	{
		v1.push_back(bndInletNames[i]);
	}

	for(int i = 0; i < bndOutletNames.size(); i++)
	{
		v1.push_back(bndOutletNames[i]);
	}

	for(int i = 0; i < bndWallNames.size(); i++)
	{
		v1.push_back(bndWallNames[i]);
	}

	for( map<string, vector<Face*> >::iterator it = msh->bnd_faces.begin(); it != msh->bnd_faces.end(); ++it)
	{
		v2.push_back(it->first);
	}

	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	return ( v1.size() == v2.size() && std::equal(v1.begin(), v1.end(), v2.begin()) );
}


void FVM_TVD_IMPLICIT::done()
{
	//free_mem(temp_mat);

	//delete [] Flux;
	//delete [] Flux1;
	//delete [] Flux2;
}


void FVM_TVD_IMPLICIT::flux_Lax_Friedrichs(double* Flux, const CellFluidDynamicsProps& cfdp1, const CellFluidDynamicsProps& cfdp2, const Point& n)
{
	double v_n1 = (cfdp1.ru * n.x + cfdp1.rv * n.y + cfdp1.rw * n.z) / cfdp1.ro;
	double v_n2 = (cfdp2.ru * n.x + cfdp2.rv * n.y + cfdp2.rw * n.z) / cfdp2.ro;

	Flux1[0] = cfdp1.ro * v_n1;
	Flux1[1] = cfdp1.ru * v_n1 + cfdp1.P * n.x;
	Flux1[2] = cfdp1.rv * v_n1 + cfdp1.P * n.y;
	Flux1[3] = cfdp1.rw * v_n1 + cfdp1.P * n.z;
	Flux1[4] = ( cfdp1.rE + cfdp1.P ) * v_n1;

	Flux2[0] = cfdp2.ro * v_n2;
	Flux2[1] = cfdp2.ru * v_n2 + cfdp2.P * n.x;
	Flux2[2] = cfdp2.rv * v_n2 + cfdp2.P * n.y;
	Flux2[3] = cfdp2.rw * v_n2 + cfdp2.P * n.z;
	Flux2[4] = ( cfdp2.rE + cfdp2.P ) * v_n2;

	double eigen_val1 = sqrt(cfdp1.gamma * cfdp1.P / cfdp1.ro) + abs( v_n1 );
	double eigen_val2 = sqrt(cfdp2.gamma * cfdp2.P / cfdp2.ro) + abs( v_n2 );
	double alpha = _max(eigen_val1, eigen_val2);

	Flux[0] = 0.5 * ( Flux1[0] + Flux2[0] - alpha * (cfdp2.ro - cfdp1.ro) );
	Flux[1] = 0.5 * ( Flux1[1] + Flux2[1] - alpha * (cfdp2.ru - cfdp1.ru) );
	Flux[2] = 0.5 * ( Flux1[2] + Flux2[2] - alpha * (cfdp2.rv - cfdp1.rv) );
	Flux[3] = 0.5 * ( Flux1[3] + Flux2[3] - alpha * (cfdp2.rw - cfdp1.rw) );
	Flux[4] = 0.5 * ( Flux1[4] + Flux2[4] - alpha * (cfdp2.rE - cfdp1.rE) );
}


#define TAUU 1E-5


__global__ void calc_fluxes(float* ro, float* ru, float* rv, float* rw, float* rE, float* P, float* gamma, float* n_x, float* n_y, float* n_z, float* S, int* inds_cell, float* fluxes, int nc, float t)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	//printf("%d \n", tid);

	float flux1[5];
	float flux2[5];

	while(tid < nc)
	{
		fluxes[5 * tid + 0] = 0;
		fluxes[5 * tid + 1] = 0;
		fluxes[5 * tid + 2] = 0;
		fluxes[5 * tid + 3] = 0;
		fluxes[5 * tid + 4] = 0;

		for(int i = 0; i < TETRA_FACES_CNT; i++)
		{
			int ind = inds_cell[TETRA_FACES_CNT * tid + i];

			if( ind >= 0 )
			{
				float v_n = (ru[tid] * n_x[TETRA_FACES_CNT * tid + i] + rv[tid] * n_x[TETRA_FACES_CNT * tid + i] + rw[tid] * n_z[TETRA_FACES_CNT * tid + i]) / ro[tid];
				float v_n_i = (ru[ind] * n_x[TETRA_FACES_CNT * tid + i] + rv[ind] * n_x[TETRA_FACES_CNT * tid + i] + rw[ind] * n_z[TETRA_FACES_CNT * tid + i]) / ro[ind];

				flux1[0] = ro[tid] * v_n;
				flux1[1] = ru[tid] * v_n + P[tid] * n_x[TETRA_FACES_CNT * tid + i];
				flux1[2] = rv[tid] * v_n + P[tid] * n_y[TETRA_FACES_CNT * tid + i];
				flux1[3] = rw[tid] * v_n + P[tid] * n_z[TETRA_FACES_CNT * tid + i];
				flux1[4] = ( rE[tid] + P[tid] ) * v_n;

				flux2[0] = ro[ind] * v_n_i;
				flux2[1] = ru[ind] * v_n_i + P[ind] * n_x[TETRA_FACES_CNT * tid + i];
				flux2[2] = rv[ind] * v_n_i + P[ind] * n_y[TETRA_FACES_CNT * tid + i];
				flux2[3] = rw[ind] * v_n_i + P[ind] * n_z[TETRA_FACES_CNT * tid + i];
				flux2[4] = ( rE[ind] + P[ind] ) * v_n_i;

				float eigen_val1 = sqrtf(gamma[tid] * P[tid] / ro[tid]) + fabsf( v_n );
				float eigen_val2 = sqrtf(gamma[ind] * P[ind] / ro[ind]) + fabsf( v_n_i );
				float alpha = fmaxf(eigen_val1, eigen_val2);

				fluxes[5 * tid + 0] = fluxes[5 * tid + 0] + 0.5 * ( flux1[0] + flux2[0] - alpha * (ro[ind] - ro[tid]) ) * S[TETRA_FACES_CNT * tid + i];
				fluxes[5 * tid + 1] = fluxes[5 * tid + 1] + 0.5 * ( flux1[1] + flux2[1] - alpha * (ru[ind] - ru[tid]) ) * S[TETRA_FACES_CNT * tid + i];
				fluxes[5 * tid + 2] = fluxes[5 * tid + 2] + 0.5 * ( flux1[2] + flux2[2] - alpha * (rv[ind] - rv[tid]) ) * S[TETRA_FACES_CNT * tid + i];
				fluxes[5 * tid + 3] = fluxes[5 * tid + 3] + 0.5 * ( flux1[3] + flux2[3] - alpha * (rw[ind] - rw[tid]) ) * S[TETRA_FACES_CNT * tid + i];
				fluxes[5 * tid + 4] = fluxes[5 * tid + 4] + 0.5 * ( flux1[4] + flux2[4] - alpha * (rE[ind] - rE[tid]) ) * S[TETRA_FACES_CNT * tid + i];
			}
			else
			{
				if( ind == - OUTLET_TYPE )
				{
					float v_n = (ru[tid] * n_x[TETRA_FACES_CNT * tid + i] + rv[tid] * n_x[TETRA_FACES_CNT * tid + i] + rw[tid] * n_z[TETRA_FACES_CNT * tid + i]) / ro[tid];

					flux1[0] = ro[tid] * v_n;
					flux1[1] = ru[tid] * v_n + P[tid] * n_x[TETRA_FACES_CNT * tid + i];
					flux1[2] = rv[tid] * v_n + P[tid] * n_y[TETRA_FACES_CNT * tid + i];
					flux1[3] = rw[tid] * v_n + P[tid] * n_z[TETRA_FACES_CNT * tid + i];
					flux1[4] = ( rE[tid] + P[tid] ) * v_n;

					fluxes[5 * tid + 0] = fluxes[5 * tid + 0] + flux1[0] * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 1] = fluxes[5 * tid + 1] + flux1[1] * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 2] = fluxes[5 * tid + 2] + flux1[2] * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 3] = fluxes[5 * tid + 3] + flux1[3] * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 4] = fluxes[5 * tid + 4] + flux1[4] * S[TETRA_FACES_CNT * tid + i];
				}
				else if( ind == - WALL_TYPE )
				{
					float v_n = (ru[tid] * n_x[TETRA_FACES_CNT * tid + i] + rv[tid] * n_x[TETRA_FACES_CNT * tid + i] + rw[tid] * n_z[TETRA_FACES_CNT * tid + i]) / ro[tid];

					float ru2 = ru[tid] - 2 * v_n * ro[tid] * n_x[TETRA_FACES_CNT * tid + i];
					float rv2 = rv[tid] - 2 * v_n * ro[tid] * n_y[TETRA_FACES_CNT * tid + i];
					float rw2 = rw[tid] - 2 * v_n * ro[tid] * n_z[TETRA_FACES_CNT * tid + i];
					float rE2 =	P[tid] / (gamma[tid] - 1) + 0.5 * (ru2*ru2 + rv2*rv2 + rw2*rw2) / ro[tid];

					float v_n_i = (ru2 * n_x[TETRA_FACES_CNT * tid + i] + rv2 * n_x[TETRA_FACES_CNT * tid + i] + rw2 * n_z[TETRA_FACES_CNT * tid + i]) / ro[tid];

					flux1[0] = ro[tid] * v_n;
					flux1[1] = ru[tid] * v_n + P[tid] * n_x[TETRA_FACES_CNT * tid + i];
					flux1[2] = rv[tid] * v_n + P[tid] * n_y[TETRA_FACES_CNT * tid + i];
					flux1[3] = rw[tid] * v_n + P[tid] * n_z[TETRA_FACES_CNT * tid + i];
					flux1[4] = ( rE[tid] + P[tid] ) * v_n;

					flux2[0] = ro[tid] * v_n_i;
					flux2[1] = ru2 * v_n_i + P[tid] * n_x[TETRA_FACES_CNT * tid + i];
					flux2[2] = rv2 * v_n_i + P[tid] * n_y[TETRA_FACES_CNT * tid + i];
					flux2[3] = rw2 * v_n_i + P[tid] * n_z[TETRA_FACES_CNT * tid + i];
					flux2[4] = ( rE2 + P[tid] ) * v_n_i;

					float alpha = sqrtf(gamma[tid] * P[tid] / ro[tid]) + fmaxf(fabsf( v_n ), fabsf( v_n_i ));

					fluxes[5 * tid + 0] = fluxes[5 * tid + 0] + 0.5 * ( flux1[0] + flux2[0] ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 1] = fluxes[5 * tid + 1] + 0.5 * ( flux1[1] + flux2[1] - alpha * (ru2 - ru[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 2] = fluxes[5 * tid + 2] + 0.5 * ( flux1[2] + flux2[2] - alpha * (rv2 - rv[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 3] = fluxes[5 * tid + 3] + 0.5 * ( flux1[3] + flux2[3] - alpha * (rw2 - rw[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 4] = fluxes[5 * tid + 4] + 0.5 * ( flux1[4] + flux2[4] - alpha * (rE2 - rE[tid]) ) * S[TETRA_FACES_CNT * tid + i];
				}
				else
				{
					int inlet_ind = -ind - INLET_TYPE; // -4 -> 0, -5 -> 1, -6 -> 2 ...
					/*
					float ro2 = inlet_ro_dev[inlet_ind];
					float ru2 = inlet_ru_dev[inlet_ind];
					float rv2 = inlet_rv_dev[inlet_ind];
					float rw2 = inlet_rw_dev[inlet_ind];
					float rE2 = inlet_rE_dev[inlet_ind];
					float P2 = inlet_P_dev[inlet_ind];
					float gamma2 = inlet_gamma_dev[inlet_ind];
					*/
					/*
					float ro2 = 0.43567;
					float ru2 = ro2 * 219.521;
					float rv2 = ro2 * 8.85526;
					float rw2 = ro2 * 0;
					float P2 = 28263.7;
					float gamma2 = 1.4;
					float rE2 = P2 / (gamma2 - 1) + 0.5 * (ru2*ru2 + rv2*rv2 + rw2*rw2) / ro2;
					*/

					float ro2 = 1.4;
					float ru2 = ro2 * 0;
					float rv2 = ro2 * 0;
					float rw2 = ro2 * 3;
					float P2 = 1;
					float gamma2 = 1.4;
					float rE2 = P2 / (gamma2 - 1) + 0.5 * (ru2*ru2 + rv2*rv2 + rw2*rw2) / ro2;


					float v_n = (ru[tid] * n_x[TETRA_FACES_CNT * tid + i] + rv[tid] * n_x[TETRA_FACES_CNT * tid + i] + rw[tid] * n_z[TETRA_FACES_CNT * tid + i]) / ro[tid];
					float v_n_i = (ru2 * n_x[TETRA_FACES_CNT * tid + i] + rv2 * n_x[TETRA_FACES_CNT * tid + i] + rw2 * n_z[TETRA_FACES_CNT * tid + i]) / ro2;

					flux1[0] = ro[tid] * v_n;
					flux1[1] = ru[tid] * v_n + P[tid] * n_x[TETRA_FACES_CNT * tid + i];
					flux1[2] = rv[tid] * v_n + P[tid] * n_y[TETRA_FACES_CNT * tid + i];
					flux1[3] = rw[tid] * v_n + P[tid] * n_z[TETRA_FACES_CNT * tid + i];
					flux1[4] = ( rE[tid] + P[tid] ) * v_n;

					flux2[0] = ro2 * v_n_i;
					flux2[1] = ru2 * v_n_i + P2 * n_x[TETRA_FACES_CNT * tid + i];
					flux2[2] = rv2 * v_n_i + P2 * n_y[TETRA_FACES_CNT * tid + i];
					flux2[3] = rw2 * v_n_i + P2 * n_z[TETRA_FACES_CNT * tid + i];
					flux2[4] = ( rE2 + P2 ) * v_n_i;

					float eigen_val1 = sqrtf(gamma[tid] * P[tid] / ro[tid]) + fabsf( v_n );
					float eigen_val2 = sqrtf(gamma2 * P2 / ro2) + fabsf( v_n_i );
					float alpha = fmaxf(eigen_val1, eigen_val2);

					fluxes[5 * tid + 0] = fluxes[5 * tid + 0] + 0.5 * ( flux1[0] + flux2[0] - alpha * (ro2 - ro[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 1] = fluxes[5 * tid + 1] + 0.5 * ( flux1[1] + flux2[1] - alpha * (ru2 - ru[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 2] = fluxes[5 * tid + 2] + 0.5 * ( flux1[2] + flux2[2] - alpha * (rv2 - rv[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 3] = fluxes[5 * tid + 3] + 0.5 * ( flux1[3] + flux2[3] - alpha * (rw2 - rw[tid]) ) * S[TETRA_FACES_CNT * tid + i];
					fluxes[5 * tid + 4] = fluxes[5 * tid + 4] + 0.5 * ( flux1[4] + flux2[4] - alpha * (rE2 - rE[tid]) ) * S[TETRA_FACES_CNT * tid + i];
				}
			}
		}

		//printf("%f %f %f %f %f %f\n",fluxes[5 * tid + 0], fluxes[5 * tid + 1], fluxes[5 * tid + 2], fluxes[5 * tid + 3], fluxes[5 * tid + 4], t);

		tid += blockDim.x * gridDim.x;
	}
}

__device__ float calc_P_dev(float ro, float rE, float ru, float rv, float rw, float gamma)
{
	return (rE - 0.5 *(ru*ru + rv*rv + rw*rw) / ro) * (gamma - 1);
}


__global__ void calc_new_values(float* ro, float* ru, float* rv, float* rw, float* rE, float* P, float* gamma, float* V, float* fluxes, int nc)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	while(tid < nc)
	{
		ro[tid] = ro[tid] - TAUU * fluxes[5 * tid + 0] / V[tid];
		ru[tid] = ru[tid] - TAUU * fluxes[5 * tid + 1] / V[tid];
		rv[tid] = rv[tid] - TAUU * fluxes[5 * tid + 2] / V[tid];
		rw[tid] = rw[tid] - TAUU * fluxes[5 * tid + 3] / V[tid];
		rE[tid] = rE[tid] - TAUU * fluxes[5 * tid + 4] / V[tid];

		P[tid] = calc_P_dev(ro[tid], rE[tid], ru[tid], rv[tid], rw[tid], gamma[tid]);
		//printf("%f %f %f %f %f %f\n",fluxes[5 * tid + 0], fluxes[5 * tid + 1], fluxes[5 * tid + 2], fluxes[5 * tid + 3], fluxes[5 * tid + 4], P[tid]);
		tid += blockDim.x * gridDim.x;
	}
}

void FVM_TVD_IMPLICIT::run()
{
	Logger::Instance()->logging()->info("TMAX = %e STEP_MAX = %d", TMAX, STEP_MAX);

	unsigned int nc = msh->cells.size();

	double t = 0;
	int step = 0;

	int num_bytes_cells = nc * sizeof(float);

	float* buffer = new float[ nc ];

	float* V_dev = 0;
	float* fluxes_dev = 0;
	float* ro_dev = 0;
	float* ru_dev = 0;
	float* rv_dev = 0;
	float* rw_dev = 0;
	float* rE_dev = 0;
	float* P_dev = 0;
	float* gamma_dev = 0;


	cudaMalloc( (void**)&V_dev, num_bytes_cells );
	cudaMalloc( (void**)&fluxes_dev, 5 * nc * sizeof(float) );
    cudaMalloc( (void**)&ro_dev, num_bytes_cells );
    cudaMalloc( (void**)&ru_dev, num_bytes_cells );
    cudaMalloc( (void**)&rv_dev, num_bytes_cells );
    cudaMalloc( (void**)&rw_dev, num_bytes_cells );
    cudaMalloc( (void**)&rE_dev, num_bytes_cells );
    cudaMalloc( (void**)&P_dev, num_bytes_cells );
    cudaMalloc( (void**)&gamma_dev, num_bytes_cells );


    cudaMemcpyToSymbol(TAU_dev, &TAU, sizeof(float), cudaMemcpyHostToDevice);

    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->V;
	}
	cudaMemcpy(V_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);


	////////////////////////

    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.ro;
	}
    cudaMemcpy(ro_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);


    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.ru;
	}
    cudaMemcpy(ru_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);


    for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.rv;
	}
	cudaMemcpy(rv_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);


	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.rw;
	}
	cudaMemcpy(rw_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.rE;
	}
	cudaMemcpy(rE_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.P;
	}
	cudaMemcpy(P_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);


	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		buffer[it->index] = it->cellFDP.gamma;
	}
	cudaMemcpy(gamma_dev, buffer, num_bytes_cells, cudaMemcpyHostToDevice);


	//////////////////


	int num_bytes_faces = TETRA_FACES_CNT * nc * sizeof(float);

	int* buffer_ind = new int[TETRA_FACES_CNT * nc];
	float* buffer_faces = new float[TETRA_FACES_CNT * nc];

	float* n_x_dev = 0;
	float* n_y_dev = 0;
	float* n_z_dev = 0;
	float* S_dev = 0;

	int* inds_cell_dev = 0;

	cudaMalloc( (void**)&n_x_dev, num_bytes_faces );
	cudaMalloc( (void**)&n_y_dev, num_bytes_faces );
	cudaMalloc( (void**)&n_z_dev, num_bytes_faces );
	cudaMalloc( (void**)&S_dev, num_bytes_faces );
	cudaMalloc( (void**)&inds_cell_dev, TETRA_FACES_CNT * nc * sizeof(int) );

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		for(int i = 0; i < TETRA_FACES_CNT; i++)
		{
			int ic = it->f[i]->in_cell;

			if(it->f[i]->c[ic] == msh->cells[it->index])
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = -it->f[i]->n.x;
			}
			else
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = it->f[i]->n.x;
			}
		}
	}
	cudaMemcpy(n_x_dev, buffer_faces, num_bytes_faces, cudaMemcpyHostToDevice);

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		for(int i = 0; i < TETRA_FACES_CNT; i++)
		{
			int ic = it->f[i]->in_cell;

			if(it->f[i]->c[ic] == msh->cells[it->index])
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = -it->f[i]->n.y;
			}
			else
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = it->f[i]->n.y;
			}
		}
	}
	cudaMemcpy(n_y_dev, buffer_faces, num_bytes_faces, cudaMemcpyHostToDevice);

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		for(int i = 0; i < TETRA_FACES_CNT; i++)
		{
			int ic = it->f[i]->in_cell;
			if(it->f[i]->c[ic] == msh->cells[it->index])
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = -it->f[i]->n.z;
			}
			else
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = it->f[i]->n.z;
			}
		}
	}
	cudaMemcpy(n_z_dev, buffer_faces, num_bytes_faces, cudaMemcpyHostToDevice);

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		for(int i = 0; i < TETRA_FACES_CNT; i++)
		{
			buffer_faces[TETRA_FACES_CNT * it->index + i] = it->f[i]->S;
		}
	}
	cudaMemcpy(S_dev, buffer_faces, num_bytes_faces, cudaMemcpyHostToDevice);

	for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		for(int i = 0; i < TETRA_FACES_CNT; i++)
		{
			Cell* other_cell = (msh->cells[it->index] == it->f[i]->c[0]) ? it->f[i]->c[1] : it->f[i]->c[0];

			if(other_cell != 0)
			{
				buffer_faces[TETRA_FACES_CNT * it->index + i] = other_cell->index;
			}
			else
			{
				switch(it->f[i]->bnd_type)
				{
					case Face::BND_TYPE_WALL:
						buffer_ind[TETRA_FACES_CNT * it->index + i] = - WALL_TYPE;
						break;

					case Face::BND_TYPE_OUTLET:
						buffer_ind[TETRA_FACES_CNT * it->index + i] = - OUTLET_TYPE;
						break;

					case Face::BND_TYPE_INLET:
					{
						buffer_ind[TETRA_FACES_CNT * it->index + i] = - INLET_TYPE;

						float _ro = it->f[i]->faceFDP.ro;
						float _ru = it->f[i]->faceFDP.ru;
						float _rv = it->f[i]->faceFDP.rv;
						float _rw = it->f[i]->faceFDP.rw;
						float _rE = it->f[i]->faceFDP.rE;
						float _P = it->f[i]->faceFDP.P;
						float _gamma = it->f[i]->faceFDP.gamma;

						cudaMemcpyToSymbol(inlet_ro_dev, &_ro, sizeof(float), cudaMemcpyHostToDevice);
						cudaMemcpyToSymbol(inlet_ru_dev, &_ru, sizeof(float), cudaMemcpyHostToDevice);
						cudaMemcpyToSymbol(inlet_rv_dev, &_rv, sizeof(float), cudaMemcpyHostToDevice);
						cudaMemcpyToSymbol(inlet_rw_dev, &_rw, sizeof(float), cudaMemcpyHostToDevice);
						cudaMemcpyToSymbol(inlet_rE_dev, &_rE, sizeof(float), cudaMemcpyHostToDevice);
						cudaMemcpyToSymbol(inlet_P_dev, &_P, sizeof(float), cudaMemcpyHostToDevice);
						cudaMemcpyToSymbol(inlet_gamma_dev, &_gamma, sizeof(float), cudaMemcpyHostToDevice);

						break;
					}
					default:
						Logger::Instance()->logging()->error("There is boundary face that has not BND_TYPE(WALL, OUTLET, INLET)");
				}
			}
		}
	}
	cudaMemcpy(inds_cell_dev, buffer_ind, TETRA_FACES_CNT * nc * sizeof(int), cudaMemcpyHostToDevice);

	free(buffer_ind);
	free(buffer_faces);


	Logger::Instance()->logging()->info("complete...");


	dim3 threads = dim3(2, 1, 1);
	//dim3 blocks = dim3(nc / threads.x + threads.x, 1, 1);
	dim3 blocks = dim3(10, 1, 1);

	Logger::Instance()->logging()->info("Solving the equation (FVM_TVD_IMPLICIT)");
	float* _buffer = new float[ 5*nc ];
	while(t < TMAX && step < STEP_MAX)
	{
		long time_start, time_end;
		time_start = clock();

		t += TAU;
		step++;


		calc_fluxes<<<blocks, threads>>> (ro_dev, ru_dev, rv_dev, rw_dev, rE_dev, P_dev, gamma_dev, n_x_dev, n_y_dev, n_z_dev, S_dev, inds_cell_dev, fluxes_dev, nc, t);
		//sleep(1); printf("%f\n",t);

		calc_new_values<<<blocks, threads>>>(ro_dev, ru_dev, rv_dev, rw_dev, rE_dev, P_dev, gamma_dev, V_dev, fluxes_dev, nc);

		cudaMemcpy(_buffer, fluxes_dev, 5 * nc * sizeof(float), cudaMemcpyDeviceToHost);
		//printf("%f %f %f %f %f\n",_buffer[5 * 4 + 0], _buffer[5 * 4 + 1], _buffer[5 * 4 + 2], _buffer[5 * 4 + 3], _buffer[5 * 4 + 4]);

		/*
		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndWallNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndWallNames); it != ite; ++it)
		{
			c1 = it->c[0]->index;

            it->faceFDP.ro = it->c[0]->cellFDP.ro;
			it->faceFDP.rE = it->c[0]->cellFDP.rE;
			it->faceFDP.P = it->c[0]->cellFDP.P;
			it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

			double rvel_n = it->c[0]->cellFDP.ru * it->n.x + it->c[0]->cellFDP.rv * it->n.y + it->c[0]->cellFDP.rw * it->n.z;


			it->faceFDP.ru = it->c[0]->cellFDP.ru - 2 * rvel_n * it->n.x;
			it->faceFDP.rv = it->c[0]->cellFDP.rv - 2 * rvel_n * it->n.y;
			it->faceFDP.rw = it->c[0]->cellFDP.rw - 2 * rvel_n * it->n.z;


			flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;

			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
		}


		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
		{
			c1 = it->c[0]->index;

			flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;

			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
		}


		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
		{
			c1 = it->c[0]->index;

			it->faceFDP.ro = it->c[0]->cellFDP.ro;
			it->faceFDP.ru = it->c[0]->cellFDP.ru;
			it->faceFDP.rv = it->c[0]->cellFDP.rv;
			it->faceFDP.rw = it->c[0]->cellFDP.rw;
			it->faceFDP.rE = it->c[0]->cellFDP.rE;
			it->faceFDP.P = it->c[0]->cellFDP.P;
			it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

			flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;
			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
		}



		for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite; ++it)
		{
			oc = it->out_cell;
			ic = it->in_cell;
			c1 = it->c[oc]->index;
			c2 = it->c[ic]->index;

			flux_Lax_Friedrichs(Flux, it->c[oc]->cellFDP, it->c[ic]->cellFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;
				right5[ c2 ][i] += Flux[i] * it->S;
			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->c[1]->cellFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);


			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);
			matrix_A(A_minus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::MINUS_JACOBIAN);


			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j]  *= it->S;
					A_minus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
			solverMtx->addMatrElement(c1, c2, A_minus);

			pc.x = -it->n.x;
			pc.y = -it->n.y;
			pc.z = -it->n.z;

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);
			matrix_A(A_minus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::MINUS_JACOBIAN);


			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j]  *= it->S;
					A_minus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c2, c2, A_plus);
			solverMtx->addMatrElement(c2, c1, A_minus);
		}


		for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
		{
			c1 = it->index;
			double V_tau = it->V / TAU;

			for(int i = 0; i < 5; i++)
			{
				mtx5[i][i] = V_tau;
			}

			solverMtx->addMatrElement(c1, c1, mtx5);
			solverMtx->setRightElement(c1, right5[c1]);
		}


		solveErr = solverMtx->solve(eps, max_iter);
		*/

		time_end = clock();

		if(step % FILE_STEP_SAVE == 0)
		{
			cudaMemcpy(buffer, ro_dev, num_bytes_cells, cudaMemcpyDeviceToHost);
			for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
			{
				it->cellFDP.ro = buffer[it->index];
			}

			cudaMemcpy(buffer, ru_dev, num_bytes_cells, cudaMemcpyDeviceToHost);
			for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
			{
				it->cellFDP.ru = buffer[it->index];
			}

			cudaMemcpy(buffer, rv_dev, num_bytes_cells, cudaMemcpyDeviceToHost);
			for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
			{
				it->cellFDP.rv = buffer[it->index];
			}

			cudaMemcpy(buffer, rw_dev, num_bytes_cells, cudaMemcpyDeviceToHost);
			for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
			{
				it->cellFDP.rw = buffer[it->index];;
			}

			cudaMemcpy(buffer, P_dev, num_bytes_cells, cudaMemcpyDeviceToHost);
			for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
			{
				it->cellFDP.P = buffer[it->index];
			}

			save(step);
		}

		if(step % LOG_STEP_SAVE == 0)
		{
			Logger::Instance()->logging()->info("step : %d\ttime step : %.16f\ttime: %d ticks", step, t, time_end - time_start);
		}

	}

    save(step);
    Logger::Instance()->logging()->info("complete...");

}




void FVM_TVD_IMPLICIT::save(int step)
{
    FILE *out;
    char c[20];

    sprintf(c, "res_%d.vtk", step);
    out = fopen(c, "w");
    fprintf(out, "# vtk DataFile Version 3.0\n");
    //The header can be used to describe the data
    fprintf(out, "GASDIN data file\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(out, "POINTS %d double\n", msh->pCount);
    for (int i = 0; i < msh->pCount; i++)
    {
        fprintf(out, "%f %f %f\n", msh->points[i].x, msh->points[i].y, msh->points[i].z);
    }

    int cellCount = msh->cells.size();

    /*
    cellSize + cellCount :
    cellSize + one number for each cell - count of points in this cell
    */
    fprintf(out, "CELLS %d %d\n", cellCount, msh->cnt_of_points + cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        fprintf(out, "%d", msh->cells[i]->pCount);

        for (int k = 0; k < msh->cells[i]->pCount; k++)
        {
            fprintf(out, " %d", msh->cells[i]->p[k]->index);
        }

        fprintf(out, "\n");
    }

    fprintf(out, "CELL_TYPES %d\n", cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        switch (msh->cells[i]->type)
        {
			case Mesh::TYPE_TETRAHEDRON:
			{
				fprintf(out, "10\n"); //10 - VTK_TETRA
				break;
			}
			case Mesh::TYPE_WEDGE:
			{
				fprintf(out, "13\n"); //13 - VTK_WEDGE
				break;
			}
			case Mesh::TYPE_HEXAHEDRON:
			{
				fprintf(out, "12\n"); //12 - VTK_HEXAHEDRON
				break;
			}
        }
    }


	fprintf(out, "CELL_DATA %d\nSCALARS Density double 1\nLOOKUP_TABLE default\n", cellCount);
	for (int i = 0; i < cellCount; i++)
	{
	  fprintf(out, "%25.16f\n", msh->cells[i]->cellFDP.ro);
	}

	fprintf(out, "SCALARS Pressure double 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < cellCount; i++)
	{
	   fprintf(out, "%25.16f\n", msh->cells[i]->cellFDP.P);
	}

	fprintf(out, "SCALARS Mach double 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < cellCount; i++)
	{
        double ro = msh->cells[i]->cellFDP.ro;
        double u = msh->cells[i]->cellFDP.ru / ro;
        double v = msh->cells[i]->cellFDP.rv / ro;
        double w = msh->cells[i]->cellFDP.rw / ro;
        double c_2 =  msh->cells[i]->cellFDP.gamma * msh->cells[i]->cellFDP.P / ro;

	   fprintf(out, "%25.16f\n", sqrt( (u*u + v*v + w*w) / c_2 ) );
	}

    fprintf(out, "VECTORS Velocity double \n");
    for (int i = 0; i < cellCount; i++)
    {
	   double ro = msh->cells[i]->cellFDP.ro;
	   fprintf(out, "%f %f %f\n", msh->cells[i]->cellFDP.ru/ro, msh->cells[i]->cellFDP.rv/ro, msh->cells[i]->cellFDP.rw/ro);
    }

    fclose(out);
}
