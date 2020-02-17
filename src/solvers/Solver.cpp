/*
 * Solver.cpp
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#include "Solver.h"
#include "Method.h"
#include "Fvm_tvd_implicit.h"
#include "../tinyxml/tinyxml.h"
#include <string.h>



Method* Solver::initMethod(char* fileName)
{
	Method * m;
	int num_of_method;


    TiXmlDocument doc( fileName );
    bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
    if (!loadOkay)
    {
        Logger::Instance()->logging()->error("Failed to open file : \"%s\"", fileName);
        //cout << doc.ErrorDesc();
        Logger::Instance()->EXIT(doc.ErrorId());
    }

    TiXmlNode* task = 0;
    TiXmlElement* el = 0;
    TiXmlNode* node0 = 0;
    TiXmlNode* node1 = 0;
    task = doc.FirstChild( "task" );

    const char* methodName = task->ToElement()->Attribute("method");

    if (strcmp("FVM_TVD_IMPLICIT", methodName) == 0)
    {
        m = new FVM_TVD_IMPLICIT();
    }
    else
    {
        Logger::Instance()->logging()->error("Unsupported method : \"%s\"", methodName);
        Logger::Instance()->EXIT(-1);
    }

    m->init(fileName);


	return m;
}

void Solver::runMethod(Method* m)
{
    m->run();
}

void Solver::destroyMethod(Method* m)
{
	m->done();
	delete m;
}
