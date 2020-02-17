
#include <cuda.h>
#include <cuda_runtime.h>
#include "global.h"
#include "solvers/Solver.h"
#include <ctime>
#include <stdio.h>
#include <iostream>




int main(int argc, char** argv)
{
	double te;

	Logger::Instance()->open_log_file("task.log");
    te = clock();

	Method* m = Solver::initMethod("task.xml");
    Solver::runMethod(m);
    Solver::destroyMethod(m);


    cout << endl << "time of execution : " << (clock() - te) / CLOCKS_PER_SEC /60 /60 << " hours"<< endl;
    Logger::Instance()->close_log_file();


    return 0;
}
