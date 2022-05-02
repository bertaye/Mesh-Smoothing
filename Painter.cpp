#include "Painter.h"

SoSeparator* Painter::getShapeSep(Mesh* mesh, bool drawThickEdges, bool drawQuad)
{
	
	SoSeparator* res = new SoSeparator();
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(1, 1, 1); 

	bool youWantToPaintEachVertexDifferently = false;



	res->addChild(mat);

	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints);
	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->verts.size(); c++)
		coords->point.set1Value(c, mesh->verts[c]->coords[0], mesh->verts[c]->coords[1], mesh->verts[c]->coords[2]);
	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	if(drawQuad)
		for (int c = 0; c < mesh->quads.size(); c++)
		{
			//printf("mesh->quads->[%d] : %d\n", c, mesh->quads[c]->v1i);
			//printf("mesh->quads->[%d] : %d\n", c, mesh->quads[c]->v2i);
			//printf("mesh->quads->[%d] : %d\n", c, mesh->quads[c]->v3i);
			//printf("mesh->quads->[%d] : %d\n", c, mesh->quads[c]->v4i);
			//printf("*****\n");


			faceSet->coordIndex.set1Value(c*5, mesh->quads[c]->v1i);
			faceSet->coordIndex.set1Value(c*5 + 1, mesh->quads[c]->v2i);
			faceSet->coordIndex.set1Value(c*5 + 2, mesh->quads[c]->v3i);
			faceSet->coordIndex.set1Value(c*5 + 3, mesh->quads[c]->v4i);
			faceSet->coordIndex.set1Value(c*5 + 4, -1);
		}
	if(!drawQuad)
		for (int c = 0; c < mesh->tris.size(); c++)
		{
			faceSet->coordIndex.set1Value(c * 4, mesh->tris[c]->v1i);
			faceSet->coordIndex.set1Value(c * 4 + 1, mesh->tris[c]->v2i);
			faceSet->coordIndex.set1Value(c * 4 + 2, mesh->tris[c]->v3i);
			faceSet->coordIndex.set1Value(c * 4 + 3, -1);
		}

	if (drawThickEdges)
	{
		SoSeparator* thickEdgeSep = new SoSeparator;
		SoMaterial* ma = new SoMaterial;
		ma->diffuseColor.set1Value(0, 0.0f, 0.0f, 1.0f);
		thickEdgeSep->addChild(ma);
		SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 1.0f;	thickEdgeSep->addChild(sty);

		SoIndexedLineSet* ils = new SoIndexedLineSet;
		SoCoordinate3* co = new SoCoordinate3;

		//assumes no edge in sedges is removed
		for (unsigned int se = 0; se < mesh->edges.size(); se++)
		{
			SbVec3f end1 = mesh->verts[mesh->edges[se]->v1i]->coords + SbVec3f(0.0f, 0.0f, 0.0f),
				end2 = mesh->verts[mesh->edges[se]->v2i]->coords + SbVec3f(0.0f, 0.0f, 0.0f);
			co->point.set1Value(2 * se, end1);
			co->point.set1Value(2 * se + 1, end2);
		}

		for (unsigned int ci = 0; ci < mesh->edges.size(); ci++)
		{
			ils->coordIndex.set1Value(3 * ci, 2 * ci);	ils->coordIndex.set1Value(3 * ci + 1, 2 * ci + 1);
			ils->coordIndex.set1Value(3 * ci + 2, -1); //end this edge with -1
		}
		thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);
		res->addChild(thickEdgeSep);
	}

	res->addChild(coords);
	res->addChild(faceSet);

	return res;
}
SoSeparator* Painter::DrawLines(Mesh* mesh, int source, int targ, int N, int* prev) {
	SoSeparator* lineSep = new SoSeparator;
	Painter* paint = new Painter();
	//material
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 1.0f, 0.0f, 0.0f);
	lineSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 8.0f;	lineSep->addChild(sty);

	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;
	vector<int> path = getPath(prev, source, targ, N);
	for (unsigned int se = 0; se < path.size(); se++)
	{
		SbVec3f end1 = mesh->verts[path[se]]->coords + SbVec3f(0.0f, 0.0f, 0.0f);

		co->point.set1Value(se, end1);
	}

	for (unsigned int ci = 0; ci < path.size()-1; ci++)
	{
		ils->coordIndex.set1Value(3 * ci, ci);	ils->coordIndex.set1Value(3 * ci + 1, ci + 1);
		ils->coordIndex.set1Value(3 * ci + 2, -1);
	}
	lineSep->addChild(co);
	lineSep->addChild(ils);
	return lineSep;
}
	
SoSeparator* Painter::get1PointSep(Mesh* obj, int pnt, int drawWhat, float colorX, float colorY, float colorZ, bool showNeighbours)
{

	Mesh* mesh = obj;

	SoSeparator* pntSep = new SoSeparator;
	//material
	SoMaterial* mat = new SoMaterial;
	mat->diffuseColor.setValue(SbColor(colorX, colorY, colorZ)); //green
	pntSep->addChild(mat);
	SoDrawStyle* style = new SoDrawStyle;
	style->pointSize = 8.0f;
	pntSep->addChild(style);
	
	//shape
	SoVertexProperty* vp = new SoVertexProperty;
	vp->vertex.set1Value(0, 1 * mesh->verts[pnt]->coords[0], 1 * mesh->verts[pnt]->coords[1], 1 * mesh->verts[pnt]->coords[2]);
	SoPointSet* pSet = new SoPointSet;
	pSet->numPoints = 1;
	pSet->vertexProperty = vp;
	pntSep->addChild(pSet);
	if(showNeighbours)
		pntSep->addChild(paintNeighbours(mesh, pnt));

	return pntSep;
}

vector<int> Painter::getPath(int* prev,int source,int targ,int N){
	int u=targ;
	vector<int> path;
	if (prev[u] != NULL || u == source) {
		//printf("in if");
		while (u != -1) {
			path.push_back(u);
			u = prev[u];
		}
	}
	return path;
}

/// <summary>
/// Use this function for debugging, to check if we succesfully found the neighbours of the vertex.
/// </summary>
/// <param name="mesh"></param>
/// <param name="pnt"></param>
/// <returns></returns>
SoSeparator* Painter::paintNeighbours(Mesh* mesh, int pnt, bool showAllVerts) {
	SoSeparator* pntSep = new SoSeparator;
	//material
	SoMaterial* mat = new SoMaterial;
	mat->diffuseColor.setValue(SbColor(1.0f, 0.0f, 0.0f));
	pntSep->addChild(mat);
	SoDrawStyle* style = new SoDrawStyle;
	style->pointSize = 10.0f;
	pntSep->addChild(style);

	SoVertexProperty* vp = new SoVertexProperty;
	if(!showAllVerts)
		for (int i = 0; i < mesh->verts[pnt]->vertList.size(); i++) {
			//printf("Neighbour %d is %d", i, mesh->verts[pnt]->vertList[i]);
			vp->vertex.
				set1Value(i, mesh->verts[mesh->verts[pnt]->vertList[i]]->coords);
		}

	if(showAllVerts)
		for (int i = 0; i < mesh->verts.size(); i++) {
			//printf("Painting vertex %d\n", i);
			vp->vertex.
				set1Value(i, mesh->verts[i]->coords);
		}
	
	SoPointSet* pSet = new SoPointSet;
	if(!showAllVerts)
		pSet->numPoints = mesh->verts[pnt]->vertList.size();
	if(showAllVerts)
		pSet->numPoints = mesh->verts.size();
	pSet->vertexProperty = vp;
	pntSep->addChild(pSet);
	return pntSep;
}

SoSeparator* Painter::compareAndPaint(Mesh* baseMesh, Mesh* newMesh){
	SoSeparator* pntSep = new SoSeparator;
	float* distances = new float[baseMesh->verts.size()];
	float tempX, tempY, tempZ, tempDist, maxDist=0.0f;
	for (int i = 0; i < baseMesh->verts.size(); i++) {
		distances[i] = 1e5;
	}

	for (int i = 0; i < newMesh->verts.size(); i++) {

		for (int j = 0; j < baseMesh->tris.size(); j++) {
			tempX = (baseMesh->verts[baseMesh->tris[j]->v1i]->coords[0] +
				baseMesh->verts[baseMesh->tris[j]->v2i]->coords[0] +
				baseMesh->verts[baseMesh->tris[j]->v3i]->coords[0]) / 3.0f;

			tempY = (baseMesh->verts[baseMesh->tris[j]->v1i]->coords[1] +
				baseMesh->verts[baseMesh->tris[j]->v2i]->coords[1] +
				baseMesh->verts[baseMesh->tris[j]->v3i]->coords[1]) / 3.0f;

			tempZ = (baseMesh->verts[baseMesh->tris[j]->v1i]->coords[2] +
				baseMesh->verts[baseMesh->tris[j]->v2i]->coords[2] +
				baseMesh->verts[baseMesh->tris[j]->v3i]->coords[2]) / 3.0f;
			
			tempDist = sqrt( pow((newMesh->verts[i]->coords[0]-tempX),2) +
				pow((newMesh->verts[i]->coords[1] - tempY), 2) +
				pow((newMesh->verts[i]->coords[2] - tempZ), 2));

			if (tempDist < distances[i]) {
				distances[i] = tempDist;
			}
			if (tempDist > maxDist) {
				maxDist = tempDist;
			}
		}
	}

	for (int i = 0; i < newMesh->verts.size(); i++) {
		distances[i] /= maxDist;
		pntSep->addChild(get1PointSep(newMesh, newMesh->verts[i]->idx, 0, distances[i], 1 - distances[i], 0.0f, false));
	}
	return pntSep;


}