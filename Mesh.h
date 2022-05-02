#pragma once

#include <vector>

using namespace std;

struct Vertex
{
	float* coords, * normals; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vertices
	vector< int > triList;
	vector <int> quadList;
	vector< int > edgeList;

	Vertex(int i, float* c) : idx(i), coords(c) {};
};

struct Edge
{
	int idx; //edges[idx]
	vector<int> quadIdx; //an edge consists of 2 quads, we need them for calculating
						 //edge point on catmull-clark algorithm
	int midPointIdx; //mid point of the edge, will be needed for catmull-clark
	int v1i, v2i; //endpnts
	float length;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) {};
};

struct Triangle
{
	int idx; //tris[idx]
	int midPointIdx;
	int v1i, v2i, v3i;
	float normals[3];
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {
	};
};

struct Quad {
	int idx;
	int facePointIdx = -1;
	int v1i, v2i, v3i, v4i;
	Quad(int id, int v1, int v2, int v3, int v4) : idx(id), v1i(v1), v2i(v2), v3i(v3), v4i(v4) {};
};

class Mesh
{
private:
	void addTriangle(int v1, int v2, int v3);
	void addSqrt3Triangle(int v1, int v2, int v3, vector<Triangle*>& sqrt3Triangles, vector<Edge*>& sqrt3Edges);

	void removeSqrt3Triangle(int triIdx, vector<Triangle*>& sqrt3Triangles);
	void replaceSqrt3Triangle(int triIdx, int v1, int v2, int v3, vector<Triangle*>& sqrt3Triangles);

	void setTriNormal(int triIdx);
	void setSqrt3TriNormal(int triIdx, vector<Triangle*>& sqrt3Triangles);

	void addEdge(int v1, int v2, int quadIdx=-1);
	void addSqrt3Edge(int v1, int v2, vector<Edge*>& sqrt3Edges);
	void replaceSqrt3Edge(int edgeIdx, int v1, int v2, vector<Edge*> sqrt3Edges);
	void removeSqrt3Edge(int edgeIdx, vector<Edge*>& sqrt3Edges);

	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);
	int getEdgeFromVerts(int v1, int v2);
	void addQuad(int v1, int v2, int v3, int v4);
	int isCommonNeighbourExist(int v1, int v2);
	float* calculateMidPoint(float* c1, float* c2,float c1_coeff=1.0f,float c2_coeff=1.0f);
	float calculateEdgeLength(int v1i, int v2i); //v1i and v2i are directly taken from edge struct

	void copyEdges(const vector<Edge*> edgeV);
	void copyQuads(const vector<Quad*> quadV);
	void copyTriangles(const vector<Triangle*> triangleV);

	void removeTriFromVert(int v, int tri);
	/*
		Catmull Clark subdivision algorithm will require updating the Quads and Edges of vertices.
		We will use mesh->quads and mesh->edges to calculate face points and edge points.
		But also we need to get use of verts->edgeList and verts->quadList so we cannot simply delete
		our mesh structure's current edges and quads lists. To rewrite them we will use seperate lists on 
		our CatmullClark function and pass these lists by ref to functions listed below.
	*/
	void addCatmullClarkEdge(int v1, int v2, vector<Edge*>& cmEdges, int quadIdx = -1); //this function will add edges to given vector
	void addCatmullClarkQuad(int v1, int v2, int v3, int v4, vector<Quad*>& cmQuads, vector<Edge*>& cmEdges, vector<int>& verts_oldEdgeList);
	int getCatmullClarkEdgeFromVerts(int v1, int v2, vector<Edge*>& cmEdges, vector<int>& verts_oldEdgeList);



public:
	vector<Vertex*> verts;
	vector<Triangle*> tris;
	vector<Edge*> edges;
	vector<Quad*> quads;
	Mesh() {};
	void createCube(float side);
	void createCubeQuad(float sideLength);
	void loadOffTriMesh(const char* name);
	void loadOffQuadMesh(const char* name);
	void convertQuadbySplit();
	void CatmullClarkSubdiv();
	void Sqrt3SubDiv(bool checkNormals=false);


};
