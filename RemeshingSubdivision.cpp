#include "Mesh.h"
#include <stdio.h>
#include <cstdlib>
#include <time.h>
int main(int, char** argv)
{
	bool quadratedMesh = false; // if mesh is a quad mes OR converted to a quad mesh, set true
	bool triangulatedMesh = false;
	bool drawEdges = false;

	Mesh* triMesh = new Mesh();
	Mesh* quadMesh = new Mesh();

	quadMesh->loadOffTriMesh("horse0.off"); //horse.off is a triangular mesh originally
	quadMesh->convertQuadbySplit();
	quadMesh->CatmullClarkSubdiv();
	quadratedMesh = true;


	triMesh->loadOffTriMesh("horse0.off");
	triMesh->Sqrt3SubDiv();
	triangulatedMesh = true;
	
	
	//Saving the new mesh
	if (quadratedMesh)
	{
		FILE* quadratedMesh;
		fopen_s(&quadratedMesh,"result_quad.off", "w+");
		fprintf(quadratedMesh, "OFF\n%d %d %d\n", quadMesh->verts.size(), quadMesh->quads.size(), quadMesh->edges.size());
		for (int i = 0; i < quadMesh->verts.size(); i++) {
			fprintf(quadratedMesh, "%0.6g %0.6g %0.6g\n", quadMesh->verts[i]->coords[0], quadMesh->verts[i]->coords[1], quadMesh->verts[i]->coords[2]);
		}
		for (int i = 0; i < quadMesh->quads.size(); i++) {
			fprintf(quadratedMesh, "4  %d %d %d %d\n", quadMesh->quads[i]->v1i, quadMesh->quads[i]->v2i, quadMesh->quads[i]->v3i, quadMesh->quads[i]->v4i);
		}
		fclose(quadratedMesh);

	}

	if(triangulatedMesh)
	{
		FILE* triangleMesh;
		fopen_s(&triangleMesh,"result_tri.off", "w+");
		fprintf(triangleMesh, "OFF\n%d %d %d\n", triMesh->verts.size(), triMesh->tris.size(), triMesh->edges.size());
		for (int i = 0; i < triMesh->verts.size(); i++) {
			fprintf(triangleMesh, "%0.6g %0.6g %0.6g\n", triMesh->verts[i]->coords[0], triMesh->verts[i]->coords[1], triMesh->verts[i]->coords[2]);
		}
		for (int i = 0; i < triMesh->tris.size(); i++) {
			fprintf(triangleMesh, "3  %d %d %d\n", triMesh->tris[i]->v1i, triMesh->tris[i]->v2i, triMesh->tris[i]->v3i);
		}
		fclose(triangleMesh);
	}
	printf("\n\nSubdivisions completed.\n\n");
	
}