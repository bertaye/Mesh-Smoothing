#include "Mesh.h"
#include "Painter.h"
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCone.h>
float getAreaFromThreePoints(float* p1, float* p2, float* p3);
int main(int, char** argv)
{
	//Coin3D Initializations
	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer* viewer = new SoWinExaminerViewer(window);
	SoSeparator* root = new SoSeparator;
	root->ref();

	Mesh* mesh = new Mesh();
	Mesh* baseMesh = new Mesh();
	Painter* painter = new Painter();

	mesh->createCube(5.0f);


	//mesh->loadOffTriMesh("horse0.off");
	//baseMesh->loadOffTriMesh("horse0.off");


	mesh->convertQuadbySplit();
	//mesh->Sqrt3SubDiv();
	//mesh->Sqrt3SubDiv();
	//mesh->loadOffQuadMesh("quadrated_result.off");

	//printf("mesh converted to quads\n\n\n\n\n\n\n\n\n\n\n\n\n\n");

	//mesh->createCubeQuad(10.0f);

	mesh->CatmullClarkSubdiv();
	mesh->CatmullClarkSubdiv();
	//mesh->CatmullClarkSubdiv();
	//mesh->CatmullClarkSubdiv();
	
    ////NOT WORKING... BECAUSE ITS  O(N1*N2) :(	root->addChild(painter->compareAndPaint(baseMesh, mesh));
	
	bool drawQuad = true; // if mesh is a quad mes OR converted to a quad mesh, set true!!!
	bool drawEdges = false;
	
	root->addChild(painter->getShapeSep(mesh, drawEdges, drawQuad));
	//root->addChild(painter->paintNeighbours(mesh, 3,true));
	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	
	if (drawQuad)
	{
		FILE* quadratedMesh;
		quadratedMesh = fopen("result.off", "w+");
		fprintf(quadratedMesh, "OFF\n%d %d %d\n", mesh->verts.size(), mesh->quads.size(), mesh->edges.size());
		for (int i = 0; i < mesh->verts.size(); i++) {
			fprintf(quadratedMesh, "%0.6g %0.6g %0.6g\n", mesh->verts[i]->coords[0], mesh->verts[i]->coords[1], mesh->verts[i]->coords[2]);
		}
		for (int i = 0; i < mesh->quads.size(); i++) {
			fprintf(quadratedMesh, "4  %d %d %d %d\n", mesh->quads[i]->v1i, mesh->quads[i]->v2i, mesh->quads[i]->v3i, mesh->quads[i]->v4i);
		}
		fclose(quadratedMesh);

	}
	else {
		FILE* triangleMesh;
		triangleMesh = fopen("result.off", "w+");
		fprintf(triangleMesh, "OFF\n%d %d %d\n", mesh->verts.size(), mesh->tris.size(), mesh->edges.size());
		for (int i = 0; i < mesh->verts.size(); i++) {
			fprintf(triangleMesh, "%0.6g %0.6g %0.6g\n", mesh->verts[i]->coords[0], mesh->verts[i]->coords[1], mesh->verts[i]->coords[2]);
		}
		for (int i = 0; i < mesh->tris.size(); i++) {
			fprintf(triangleMesh, "3  %d %d %d\n", mesh->tris[i]->v1i, mesh->tris[i]->v2i, mesh->tris[i]->v3i);
		}
		fclose(triangleMesh);
	}

	float totalSurfaceArea=0.0f, totalSurfaceArea_base=0.0f;
	if (!drawQuad) {

		for (int i = 0; i < mesh->tris.size(); i++) {
			if (i < baseMesh->tris.size()) {

				totalSurfaceArea_base += getAreaFromThreePoints(baseMesh->verts[baseMesh->tris[i]->v1i]->coords,
					baseMesh->verts[baseMesh->tris[i]->v2i]->coords,
					baseMesh->verts[baseMesh->tris[i]->v3i]->coords);
			}

			totalSurfaceArea += getAreaFromThreePoints(mesh->verts[mesh->tris[i]->v1i]->coords,
				mesh->verts[mesh->tris[i]->v2i]->coords,
				mesh->verts[mesh->tris[i]->v3i]->coords);

		}
		FILE* comparedResults;
		comparedResults = fopen("comparedResults.txt", "w+");
		fprintf(comparedResults, "\t\t\tBASE MESH \t\t\t DIVIDED MESH\n");
		fprintf(comparedResults, "Face Count: \t\t %d \t\t\t\t %d\n", baseMesh->tris.size(), mesh->quads.size());
		fprintf(comparedResults, "Total Surface Area:\t%.9g\t\t\t%.9g\n", totalSurfaceArea_base, totalSurfaceArea);

		fclose(comparedResults);


	}
	else {
		for (int i = 0; i < mesh->quads.size(); i++) {
			if (i < baseMesh->tris.size()) {

				totalSurfaceArea_base += getAreaFromThreePoints(baseMesh->verts[baseMesh->tris[i]->v1i]->coords,
					baseMesh->verts[baseMesh->tris[i]->v2i]->coords,
					baseMesh->verts[baseMesh->tris[i]->v3i]->coords);
			}

			totalSurfaceArea += getAreaFromThreePoints(mesh->verts[mesh->quads[i]->v1i]->coords,
				mesh->verts[mesh->quads[i]->v2i]->coords,
				mesh->verts[mesh->quads[i]->v3i]->coords);

			totalSurfaceArea += getAreaFromThreePoints(mesh->verts[mesh->quads[i]->v2i]->coords,
				mesh->verts[mesh->quads[i]->v3i]->coords,
				mesh->verts[mesh->quads[i]->v4i]->coords);

		}
		FILE* comparedResults;
		comparedResults = fopen("comparedResults.txt", "w+");
		fprintf(comparedResults, "\t\t\tBASE MESH \t\t\t DIVIDED MESH\n");
		fprintf(comparedResults, "Face Count: \t\t %d \t\t\t\t %d\n", baseMesh->tris.size(), mesh->quads.size());
		fprintf(comparedResults, "Total Surface Area:\t%.9g\t\t\t%.9g\n", totalSurfaceArea_base, totalSurfaceArea);

		fclose(comparedResults);
	}
	
	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
}

float getAreaFromThreePoints(float* p1, float* p2, float* p3) {
	float U[3], V[3], crossP[3];
	U[0] = p2[0] - p1[0];
	U[1] = p2[1] - p1[1];
	U[2] = p2[2] - p1[2];

	V[0] = p3[0] - p1[0];
	V[1] = p3[1] - p1[1];
	V[2] = p3[2] - p1[2];

	crossP[0] = U[1] * V[2] - U[2] * V[1];
	crossP[1] = U[2] * V[0] - U[0] * V[2];
	crossP[2] = U[0] * V[1] - U[1] * V[0];

	return sqrt(pow(crossP[0], 2) + pow(crossP[1], 2) + pow(crossP[2], 2))/2.0f;

}