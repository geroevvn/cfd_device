#include "Mesh.h"
#include "../mesh_properties/BndFaceTemperature.h"
#include "../mesh_properties/CellTemperature.h"
#include "../mesh_properties/CellFluidDynamicsProps.h"
#include "../iterators/MeshIterator.h"
#include "../iterators/PointIterator.h"
#include "../iterators/FilterIterator.h"
#include "../iterators/BndIterator.h"


///////////////////////

inline Point Point::operator+(const Point& b) const
{
    Point vec;

    vec.x = this->x + b.x;
    vec.y = this->y + b.y;
    vec.z = this->z + b.z;

    return vec;
}

inline Point Point::operator/(double lyambda) const
{
    Point vec;

    vec.x = this->x / lyambda;
    vec.y = this->y / lyambda;
    vec.z = this->z / lyambda;

    return vec;
}


double Point::scalar_product(Point* vec1, Point* vec2)
{
    return vec1->x*vec2->x + vec1->y*vec2->y + vec1->z*vec2->z;
}

void Point::normalize(Point* p1)
{
    double norm = sqrt(p1->x * p1->x + p1->y * p1->y + p1->z * p1->z);

    p1->x /= norm;
    p1->y /= norm;
    p1->z /= norm;
}

Edge::Edge(Point* p1, Point* p2)
{
    p[0] = p1;
    p[1] = p2;
}

double Edge::getlength()
{
    return sqrt(pow(p[0]->x - p[1]->x, 2) + pow(p[0]->y - p[1]->y, 2) + pow(p[0]->z - p[1]->z, 2));
}

double Edge::getlength(Point* p1, Point* p2)
{
    return sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2) + pow(p1->z - p2->z, 2));
}

Point Edge::getVector(Point* p1, Point* p2)
{
    Point vec;

    vec.x = p2->x - p1->x;
    vec.y = p2->y - p1->y;
    vec.z = p2->z - p1->z;

    return vec;
}

////////////////////////

Face::~Face()
{
    delete [] p;
    delete [] e;
}

Face::Face(const int type, Point* p1, Point* p2, Point* p3)
{
	this->bnd_type = 0;
    this->type = type;
    pCount = 3;
    eCount = 3;
    isUsed = false;
    c[0] = 0;
    c[1] = 0;

    p = new Point*[pCount];
    e = new Edge*[eCount];

    p[0] = p1;
    p[1] = p2;
    p[2] = p3;

    Point vec1 = Edge::getVector(p[0],p[1]);
    Point vec2 = Edge::getVector(p[0],p[2]);

    n.x = vec1.y*vec2.z - vec1.z*vec2.y;
    n.y = vec1.z*vec2.x - vec1.x*vec2.z;
    n.z = vec1.x*vec2.y - vec1.y*vec2.x;

    Point::normalize(&n);
}

void Face::initFace(Cell* c, Edge* e1, Edge* e2, Edge* e3)
{
  if( !isUsed ) { // инициализируем в первый раз

    e[0] = e1;
    e[1] = e2;
    e[2] = e3;

    this->c[0] = c;
 }
 else  {
    this->c[1] = c;
 }
  isUsed = true;
}


Face::Face(const int type, Point* p1, Point* p2, Point* p3, Point* p4)
{
	this->bnd_type = 0;
    this->type = type;
    pCount = 4;
    eCount = 4;
    isUsed = false;
    c[0] = 0;
    c[1] = 0;

    p = new Point*[pCount];
    e = new Edge*[eCount];

    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    p[3] = p4;

    Point vec1 = Edge::getVector(p[0],p[1]);
    Point vec2 = Edge::getVector(p[0],p[2]);

    n.x = vec1.y*vec2.z - vec1.z*vec2.y;
    n.y = vec1.z*vec2.x - vec1.x*vec2.z;
    n.z = vec1.x*vec2.y - vec1.y*vec2.x;

    Point::normalize(&n);
}

void Face::initFace(Cell* c, Edge* e1, Edge* e2, Edge* e3, Edge* e4)
{
  if( !isUsed ) { // инициализируем в первый раз

    e[0] = e1;
    e[1] = e2;
    e[2] = e3;
    e[3] = e4;

    this->c[0] = c;
 }
 else  {
    this->c[1] = c;
 }
  isUsed = true;
}

// Вычисление площади грани
void Face::area()
{
    S = 0;
     switch( type )
     {
      case Mesh::TYPE_TRIANGLE: // треугольник( код формата unv - 91)
      {
            double a = Edge::getlength(p[0], p[1]);
            double b = Edge::getlength(p[0], p[2]);
            double c = Edge::getlength(p[1], p[2]);

            double pp = (a + b + c)/2;

            S = sqrt(pp*(pp-a)*(pp-b)*(pp-c));
        break;
      }
      case Mesh::TYPE_QUADRILATERAL: // четырехугольник ( код формата unv - 94)
      {
            double a = Edge::getlength(p[0], p[1]);
            double b = Edge::getlength(p[1], p[2]);
            double diag = Edge::getlength(p[0], p[2]);
            double c = Edge::getlength(p[2], p[3]);
            double d = Edge::getlength(p[3], p[0]);

            double p1 = (a + b + diag)/2;
            double p2 = (c + d + diag)/2;

            S = sqrt(p1*(p1-a)*(p1-b)*(p1-diag)) + sqrt(p2*(p2-c)*(p2-d)*(p2-diag));
        break;
      }
    }
}

// Вычисление расстояния между центрами
void Face::calc_h() {

    Point center1 = c[0]->center;

    if(c[1] != 0) {

        Point center2 = c[1]->center;
        h = Edge::getlength(&center1, &center2);
    }
    else {

        Point vec1 = Edge::getVector(p[0],p[1]);
        Point vec2 = Edge::getVector(p[0],p[2]);

        Point ort_vec;
        ort_vec.x = vec1.y*vec2.z - vec1.z*vec2.y;
        ort_vec.y = vec1.z*vec2.x - vec1.x*vec2.z;
        ort_vec.z = vec1.x*vec2.y - vec1.y*vec2.x;

        Point::normalize(&ort_vec);

        Point vec_res = Edge::getVector(p[0], &(c[0]->center));

        h = 2 * abs(Point::scalar_product(&vec_res, &ort_vec));
    }
}

///////////////////


Cell::~Cell()
{
    delete [] p;
    delete [] e;
    delete [] f;
}

Cell::Cell(const int type, Point* p1, Point* p2, Point* p3, Point* p4)
{
    this->type = type;
    pCount = 4;
    eCount = 6;
    fCount = 4;

    p = new Point*[pCount];
    e = new Edge*[eCount];
    f = new Face*[fCount];

    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    p[3] = p4;

    center = (*p[0] + *p[1] + *p[2] + *p[3]) / 4;
}

Cell::Cell(const int type, Point* p1, Point* p2, Point* p3, Point* p4, Point* p5, Point* p6)
{
    this->type = type;
    pCount = 6;
    eCount = 9;
    fCount = 5;

    p = new Point*[pCount];
    e = new Edge*[eCount];
    f = new Face*[fCount];

    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    p[3] = p4;
    p[4] = p5;
    p[5] = p6;

    center = (*p[0] + *p[1] + *p[2] + *p[3] + *p[4] + *p[5]) / 6;
}

Cell::Cell(const int type, Point* p1, Point* p2, Point* p3, Point* p4, Point* p5, Point* p6, Point* p7, Point* p8)
{
    this->type = type;
    pCount = 8;
    eCount = 12;
    fCount = 6;

    p = new Point*[pCount];
    e = new Edge*[eCount];
    f = new Face*[fCount];

    p[0] = p1;
    p[1] = p2;
    p[2] = p3;
    p[3] = p4;
    p[4] = p5;
    p[5] = p6;
    p[6] = p7;
    p[7] = p8;

    center = (*p[0] + *p[1] + *p[2] + *p[3] + *p[4] + *p[5] + *p[6] + *p[7]) / 8;
}

// Вычисление объема клетки
void Cell::volume()
{
    switch( type ) {

     case Mesh::TYPE_TETRAHEDRON:
        {

        Point vec_a = Edge::getVector(p[0], p[1]);
        Point vec_b = Edge::getVector(p[0], p[2]);
        Point vec_c = Edge::getVector(p[0], p[3]);

        V = abs( vec_a.x*(vec_b.y*vec_c.z - vec_b.z*vec_c.y) - vec_a.y*(vec_b.x*vec_c.z - vec_b.z*vec_c.x) + vec_a.z*(vec_b.x*vec_c.y - vec_b.y*vec_c.x) ) / 6;

        break;
        }

     case Mesh::TYPE_WEDGE:
        {
         Point vec_10 = Edge::getVector(p[1], p[0]);
         Point vec_12 = Edge::getVector(p[1], p[2]);
         Point vec_13 = Edge::getVector(p[1], p[3]);

         V = abs( vec_10.x*(vec_12.y*vec_13.z - vec_12.z*vec_13.y) - vec_10.y*(vec_12.x*vec_13.z - vec_12.z*vec_13.x) + vec_10.z*(vec_12.x*vec_13.y - vec_12.y*vec_13.x));

         Point vec_14 = Edge::getVector(p[1], p[4]);

         V += abs( vec_12.x*(vec_13.y*vec_14.z - vec_13.z*vec_14.y) - vec_12.y*(vec_13.x*vec_14.z - vec_13.z*vec_14.x) + vec_12.z*(vec_13.x*vec_14.y - vec_13.y*vec_14.x));

         Point vec_23 = Edge::getVector(p[2], p[3]);
         Point vec_25 = Edge::getVector(p[2], p[5]);
         Point vec_24 = Edge::getVector(p[2], p[4]);

         V += abs( vec_23.x*(vec_24.y*vec_25.z - vec_24.z*vec_25.y) - vec_23.y*(vec_24.x*vec_25.z - vec_24.z*vec_25.x) + vec_23.z*(vec_24.x*vec_25.y - vec_24.y*vec_25.x));

         V /= 6;

         break;
        }

     case Mesh::TYPE_HEXAHEDRON:
        {
        Point vec_01 = Edge::getVector(p[0], p[1]);
        Point vec_03 = Edge::getVector(p[0], p[3]);
        Point vec_04 = Edge::getVector(p[0], p[4]);

        V = abs( vec_01.x*(vec_04.y*vec_03.z - vec_04.z*vec_03.y) - vec_01.y*(vec_04.x*vec_03.z - vec_04.z*vec_03.x) + vec_01.z*(vec_04.x*vec_03.y - vec_04.y*vec_03.x) );

        Point vec_12 = Edge::getVector(p[1], p[2]);
        Point vec_13 = Edge::getVector(p[1], p[3]);
        Point vec_14 = Edge::getVector(p[1], p[4]);

        V += abs( vec_12.x*(vec_13.y*vec_14.z - vec_13.z*vec_14.y) - vec_12.y*(vec_13.x*vec_14.z - vec_13.z*vec_14.x) + vec_12.z*(vec_13.x*vec_14.y - vec_13.y*vec_14.x) );

        Point vec_15 = Edge::getVector(p[1], p[5]);

        V += abs( vec_12.x*(vec_14.y*vec_15.z - vec_14.z*vec_15.y) - vec_12.y*(vec_14.x*vec_15.z - vec_14.z*vec_15.x) + vec_12.z*(vec_14.x*vec_15.y - vec_14.y*vec_15.x) );

        Point vec_76 = Edge::getVector(p[7], p[6]);
        Point vec_73 = Edge::getVector(p[7], p[3]);
        Point vec_74 = Edge::getVector(p[7], p[4]);

        V += abs( vec_76.x*(vec_74.y*vec_73.z - vec_74.z*vec_73.y) - vec_76.y*(vec_74.x*vec_73.z - vec_74.z*vec_73.x) + vec_76.z*(vec_74.x*vec_73.y - vec_74.y*vec_73.x) );

        Point vec_62 = Edge::getVector(p[6], p[2]);
        Point vec_63 = Edge::getVector(p[6], p[3]);
        Point vec_64 = Edge::getVector(p[6], p[4]);

        V += abs( vec_62.x*(vec_63.y*vec_64.z - vec_63.z*vec_64.y) - vec_62.y*(vec_63.x*vec_64.z - vec_63.z*vec_64.x) + vec_62.z*(vec_63.x*vec_64.y - vec_63.y*vec_64.x) );

        Point vec_65 = Edge::getVector(p[6], p[5]);

        V += abs( vec_62.x*(vec_64.y*vec_65.z - vec_64.z*vec_65.y) - vec_62.y*(vec_64.x*vec_65.z - vec_64.z*vec_65.x) + vec_62.z*(vec_64.x*vec_65.y - vec_64.y*vec_65.x) );

        V /= 6;

        break;
        }
    }
}

////////////////////

bool IsBoundaryFace::operator()(Face* f) { Cell** c = f->getCells(); return (c[1] == NULL); }
bool IsInnerFace::operator()(Face* f) { Cell** c = f->getCells(); return (c[1] != NULL); }

////////////////////

void Mesh::createPoints(Point* p, unsigned int pCount)
{
        this->pCount = pCount;
        points = new Point[pCount];
        for(int i = 0; i < pCount; i++)
        {
            points[i] = p[i];
            points[i].index = i;
        }
}

void Mesh::createPoints(const vector<Point>& p)
{
		//cout << "ASDASD";
        this->pCount = p.size();
        points = new Point[pCount];

        for(int i = 0; i < pCount; i++)
        {
            points[i] = p[i];
            points[i].index = i;
        }
}

Mesh::Mesh(){}

Mesh::~Mesh()
{
        delete [] points;

        for(vector<Edge*>::iterator it = edges.begin(); it != edges.end(); ++it)
            delete (*it);

         for(vector<Face*>::iterator it = faces.begin(); it != faces.end(); ++it)
            delete (*it);

         for(vector<Cell*>::iterator it = cells.begin(); it != cells.end(); ++it)
            delete (*it);

        // cout << "MEsh" << endl;
}


Mesh::CellIterator Mesh::beginCell() { return CellIterator(cells.begin()); }
Mesh::CellIterator Mesh::endCell() { return CellIterator(cells.end()); }

Mesh::FaceIterator Mesh::beginFace() { return FaceIterator(faces.begin()); }
Mesh::FaceIterator Mesh::endFace() { return FaceIterator(faces.end()); }

Mesh::EdgeIterator Mesh::beginEdge() { return EdgeIterator(edges.begin()); }
Mesh::EdgeIterator Mesh::endEdge() { return EdgeIterator(edges.end()); }

PointIterator Mesh::beginPoint() { return PointIterator(points); }
PointIterator Mesh::endPoint() { return PointIterator(&points[pCount]); }

IsBoundaryFace predicateBoundary;
Mesh::BoundaryFaceIterator Mesh::beginBoundaryFace()
{
    FaceIterator bf = beginFace();
    return BoundaryFaceIterator(predicateBoundary, bf, endFace(), bf);
}

Mesh::BoundaryFaceIterator Mesh::endBoundaryFace()
{
    FaceIterator ef = endFace();
    return BoundaryFaceIterator(predicateBoundary, beginFace(), ef, ef);
}

/*
IsInnerFace predicateInner;
Mesh::InnerFaceIterator Mesh::beginInnerFace()
{
    FaceIterator bf = beginFace();
    return InnerFaceIterator(predicateInner, bf, endFace(), bf);
}

Mesh::InnerFaceIterator Mesh::endInnerFace()
{
    FaceIterator ef = endFace();
    return InnerFaceIterator(predicateInner, beginFace(), ef, ef);
}
*/
////
Mesh::FaceIterator Mesh::beginInnerFace(){ return FaceIterator(inner_faces.begin()); }
Mesh::FaceIterator Mesh::endInnerFace(){ return FaceIterator(inner_faces.end()); }
////
Mesh::FaceIterator Mesh::beginBndFace(string str) { return FaceIterator(bnd_faces[str].begin()); }
Mesh::FaceIterator Mesh::endBndFace(string str) { return FaceIterator(bnd_faces[str].end()); }

Mesh::BndFaceIterator Mesh::beginBndFace(map<std::string, vector<Face*>> *m, vector<string> *str)
{
    auto it = str->begin();
    auto it1 = str->end();
    if(it1 == it)
    {
    	vector<Face*> f(0);
    	return BndFaceIterator(m, it, str->end(), f.begin());
    }

    return BndFaceIterator(m, it, str->end(), bnd_faces[*it].begin());
}
Mesh::BndFaceIterator Mesh::endBndFace(map<std::string, vector<Face*>> *m, vector<string> *str)
{
    auto it = str->end();
    auto it1 = str->begin();
    if(it == it1)
    {
    	vector<Face*> f(0);
    	return BndFaceIterator(m, it, it, f.end());
    }

    return BndFaceIterator(m, it, it, bnd_faces[*(it - 1)].end());
}

void Mesh::iterateCells(iterateCellsFunc f)
{
    for (CellIterator it = beginCell(), ite = endCell(); it != ite; ++it)
    {
        f(&(*it));
    }
}

void Mesh::iterateFaces(iterateFacesFunc f)
{
    for (FaceIterator it = beginFace(), ite = endFace(); it != ite; ++it)
    {
        f(&(*it));
    }
}

void Mesh::iterateEdges(iterateEdgesFunc f)
{
    for (EdgeIterator it = beginEdge(), ite = endEdge(); it != ite; ++it)
    {
        f(&(*it));
    }
}

void Mesh::iteratePoints(iteratePointsFunc f)
{
    for (PointIterator it = beginPoint(), ite = endPoint(); it != ite; ++it)
    {
        f(&(*it));
    }
}

///////////////////////////////

/*int pointC = 0;
void funcCells(Cell *c)
{
    pointC++;
}

void funcPoints(Point *p)
{
    std::cout << p->x << std::endl;
}
*/
