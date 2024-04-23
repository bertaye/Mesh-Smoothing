#include "Mesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void Mesh::loadOffTriMesh(const char* name)
{
	std::ifstream file(name);
	if (!file.is_open()) {
		std::cerr << "Failed to open file.\n";
		return;
	}

	std::string line;
	// Read and ignore the first line containing "OFF"
	std::getline(file, line);

	// Read the second line for vertex, triangle counts
	std::getline(file, line);
	std::istringstream iss(line);
	int nVerts, nTris, dummy;
	if (!(iss >> nVerts >> nTris >> dummy)) { // The last integer is not used
		std::cerr << "Failed to read header.\n";
		return;
	}

	float x, y, z;
	for (int i = 0; i < nVerts; ++i) {
		if (!getline(file, line)) {
			std::cerr << "Failed to read vertex data.\n";
			return;
		}
		std::istringstream vertStream(line);
		if (!(vertStream >> x >> y >> z)) {
			std::cerr << "Failed to parse vertex coordinates.\n";
			return;
		}
		addVertex(x, y, z);
	}

	int discard, v1, v2, v3;
	while (getline(file, line)) {
		std::istringstream triStream(line);
		if (!(triStream >> discard >> v1 >> v2 >> v3)) {
			std::cerr << "Failed to parse triangle indices.\n";
			return;
		}
		addTriangle(v1, v2, v3);
		setTriNormal(tris.size() - 1);
		makeVertsNeighbor(v1, v2);
		makeVertsNeighbor(v2, v3);
		makeVertsNeighbor(v1, v3);
	}
}

void Mesh::loadOffQuadMesh(const char* name) {
	FILE* fPtr;
	fopen_s(&fPtr,name, "r");
	char str[334];

	fscanf_s(fPtr, "%s", str,334);

	int nVerts, nQuads, n, i = 0;
	float x, y, z, w;

	fscanf_s(fPtr, "%d %d %d\n", &nVerts, &nQuads, &n);
	while (i++ < nVerts)
	{
		fscanf_s(fPtr, "%f %f %f", &x, &y, &z);
		addVertex(x, y, z);
	}

	while (fscanf_s(fPtr, "%d", &i) != EOF)
	{
		fscanf_s(fPtr, "%f %f %f %f", &x, &y, &z, &w);
		addQuad((int)x, (int)y, (int)z, (int)w);
		
		makeVertsNeighbor(x, y);
		makeVertsNeighbor(y, z);
		makeVertsNeighbor(z, w);
		makeVertsNeighbor(w, x);
	}

	fclose(fPtr);
}

void Mesh::CatmullClarkSubdiv() {
	float* tempVert = new float[3]; //[0]=x, [1]=y, [2]=z
	if (quads.size() < 1) {
		//printf("Cant perform Catmull-Clark subdivision on non-quad meshes");
		return;
	}

	//after creating new vertices we need to update coords of old verts.
	//we must know where should we stop.
	int oldVertsIdx = verts.size() - 1;

	////printf("Creating face points!\n");
	//create face points for each quad
	/*
		A______B
	    |  	   | 	
	    |  fp  |     Face Point: fp = (A+B+C+D)/4
	   C|______|D

	*/
	for (int i = 0; i < quads.size(); i++) {

		tempVert[0] = (verts[quads[i]->v1i]->coords[0] + verts[quads[i]->v2i]->coords[0]
			+ verts[quads[i]->v3i]->coords[0] + verts[quads[i]->v4i]->coords[0]) / 4.0f;

		tempVert[1] = (verts[quads[i]->v1i]->coords[1] + verts[quads[i]->v2i]->coords[1]
			+ verts[quads[i]->v3i]->coords[1] + verts[quads[i]->v4i]->coords[1]) / 4.0f;

		tempVert[2] = (verts[quads[i]->v1i]->coords[2] + verts[quads[i]->v2i]->coords[2]
			+ verts[quads[i]->v3i]->coords[2] + verts[quads[i]->v4i]->coords[2]) / 4.0f;

		addVertex(tempVert[0], tempVert[1], tempVert[2]);
		//////printf("Face point coordinats %.5g, %.5g, %.5g\n", tempVert[0], tempVert[1], tempVert[2]);

		//to decrease complexity we will store facepoint idx of each quad
		quads[i]->facePointIdx = verts.size() - 1;
	}

	//create edge points for each edge
	////printf("Creating edge points!\n");

	/*  _______E1______ 
		|      |      |
		|  F1  ep  F2 |  Edge Point: ep = (F1+F2+E1+E2)\4
		|______|______|
			   E2	
	*/
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i]->quadIdx.size() != 2) {
			//printf("This mesh NOT watertight and edges[%d] is at boundary.\n",i);
			edges[i]->midPointIdx = -1;
			continue;
		}

		//printf("edges[%d]->quadIdx.size() = %d\n",i, edges[i]->quadIdx.size());
		//printf("edges[%d]->quadIdx[0] = %d\n",i, edges[i]->quadIdx[0]);
		//printf("edges[%d]->quadIdx[1] = %d\n",i, edges[i]->quadIdx[1]);
		

		tempVert[0] = (verts[quads[edges[i]->quadIdx[0]]->facePointIdx]->coords[0] +
			verts[quads[edges[i]->quadIdx[1]]->facePointIdx]->coords[0] +
			verts[edges[i]->v1i]->coords[0] +
			verts[edges[i]->v2i]->coords[0]) / 4.0f;

		tempVert[1] = (verts[quads[edges[i]->quadIdx[0]]->facePointIdx]->coords[1] +
			verts[quads[edges[i]->quadIdx[1]]->facePointIdx]->coords[1] +
			verts[edges[i]->v1i]->coords[1] +
			verts[edges[i]->v2i]->coords[1]) / 4.0f;

		tempVert[2] = (verts[quads[edges[i]->quadIdx[0]]->facePointIdx]->coords[2] +
			verts[quads[edges[i]->quadIdx[1]]->facePointIdx]->coords[2] +
			verts[edges[i]->v1i]->coords[2] +
			verts[edges[i]->v2i]->coords[2]) / 4.0f;
		////printf("Edge point coordinats %.5g, %.5g, %.5g\n", tempVert[0], tempVert[1], tempVert[2]);
		addVertex(tempVert[0], tempVert[1], tempVert[2]);
		edges[i]->midPointIdx = verts.size() - 1;

	}

	//printf("*********Moving the old vertices now!\n");
	//move old verts to their new positions
	/*
	            |
		        |
		fp1	    ep1    fp2
				|
		        |
	  __ep2_____v_____ep4______        fp: face point, ep: edge point
		        |                       
				|                       v_new = (fp1+fp2+fp3+fp4)/4/n + 2*(ep1+ep2+ep3+ep4)/4/n + (n-3)v/n
        fp3 	ep3    fp4
				|
				|
	
	*/
	float tempX_fp=-1.0f,tempX_ep = -1.0f, tempY_fp = -1.0f, tempY_ep = -1.0f, tempZ_fp = -1.0f, tempZ_ep = -1.0f;
	int n=-1;
	for (int i = 0; i < oldVertsIdx+1; i++) {
		if (verts[i]->quadList.size() != verts[i]->edgeList.size()) {
			//printf("EdgeList.size() of %d th vertex is %d\n", i, n);
			//printf("QuadList.size() of %d th vertex is %d\n\n", i, verts[i]->quadList.size());
			continue;
		}

		//printf("Moving %d th vertex\n", i);
		n = verts[i]->edgeList.size();
		//printf("EdgeList.size() of %d th vertex is %d\n", i, n);
		//printf("QuadList.size() of %d th vertex is %d\n", i, verts[i]->quadList.size());
		tempX_fp = 0.0f;
		tempX_ep = 0.0f;

		tempY_fp = 0.0f;
		tempY_ep = 0.0f;

		tempZ_fp = 0.0f;
		tempZ_ep = 0.0f;
		for (int j = 0; j < n; j++) {
			//printf("quads.size() = %d\n", quads.size());
			//printf("verts[%d]->quadList[%d] = %d\n", i, j, verts[i]->quadList[j]);

			//printf("edges.size() = %d\n", edges.size());
			//printf("verts[%d]->edgeList[%d] = %d\n", i, j, verts[i]->edgeList[j]);

			tempX_fp += verts[quads[verts[i]->quadList[j]]->facePointIdx]->coords[0]/n;
			tempX_ep += verts[edges[verts[i]->edgeList[j]]->midPointIdx]->coords[0] / n;

			tempY_fp += verts[quads[verts[i]->quadList[j]]->facePointIdx]->coords[1] / n;
			tempY_ep += verts[edges[verts[i]->edgeList[j]]->midPointIdx]->coords[1] / n;

			tempZ_fp += verts[quads[verts[i]->quadList[j]]->facePointIdx]->coords[2] / n;
			tempZ_ep += verts[edges[verts[i]->edgeList[j]]->midPointIdx]->coords[2] / n;
		}
		tempVert[0] = ((float)tempX_fp / n) + ((float)2.0f * tempX_ep / n) + ((float)(n - 3) * verts[i]->coords[0] / (n));

		tempVert[1] = ((float)tempY_fp / n) + ((float)2.0f * tempY_ep / n) + ((float)(n - 3) * verts[i]->coords[1] / (n));

		tempVert[2] = ((float)tempZ_fp / n) + ((float)2.0f * tempZ_ep / n) + ((float)(n - 3) * verts[i]->coords[2] / (n));

		verts[i]->coords[0] = tempVert[0];
		verts[i]->coords[1] = tempVert[1];
		verts[i]->coords[2] = tempVert[2];

		
	}
	
	//printf("*********We will create new quads and edges now!\n");


	//now we need to create new quads and edges!
	vector<Edge*> catmullEdges;
	vector<Quad*> catmullQuads;

	//new quads will contain only 1 old vertex per quad!
	/*
		we will create quads by going through quadList of verties and then get the mutual edges of that quad & vertex
		because edges contain all information needed; -midpoint of edges and hence neighbouring midpoints to vertex! (ep1,ep2)
													  -quads of edges and hence neighbouring facepoints to vertex! (fp)
		and our new quad will be simply this-> v,ep1,fp,ep2
													  
	*/
	vector<int> verts_oldEdges;

	for (int i = 0; i < verts.size(); i++) {
		verts[i]->quadList.clear();
		verts[i]->vertList.clear();
	}

	for (int i = 0; i < oldVertsIdx+1; i++) {
		verts_oldEdges.push_back(verts[i]->edgeList.size());
	}

	int tempEdge1_2, tempEdge2_3, tempEdge3_4, tempEdge4_1;
	int tempQuadIdx=-1, tempEdge1=-1, tempEdge2=-1;
	int vert1_oldQuads, vert1_oldEdges, vert2_oldQuads, vert2_oldEdges,
		vert3_oldQuads, vert3_oldEdges, vert4_oldQuads, vert4_oldEdges;
	
	for (int i = 0; i < quads.size(); i++) {
		//printf("%d th quad\n", i);
		//printf("edges.size() = %d\n", edges.size());
		tempEdge1_2 = getEdgeFromVerts(quads[i]->v1i, quads[i]->v2i);
		//printf("tempEdge1_2 = %d\n", tempEdge1_2);
		tempEdge2_3 = getEdgeFromVerts(quads[i]->v2i, quads[i]->v3i);
		//printf("tempEdge2_3 = %d\n", tempEdge2_3);

		tempEdge3_4 = getEdgeFromVerts(quads[i]->v3i, quads[i]->v4i);
		//printf("tempEdge3_4 = %d\n", tempEdge3_4);

		tempEdge4_1 = getEdgeFromVerts(quads[i]->v4i, quads[i]->v1i);
		//printf("tempEdge4_1 = %d\n", tempEdge4_1);
		
		if (edges[tempEdge1_2]->midPointIdx == -1 || edges[tempEdge2_3]->midPointIdx == -1
			|| edges[tempEdge3_4]->midPointIdx == -1 || edges[tempEdge4_1]->midPointIdx == -1) {
			continue;
		}

		addCatmullClarkQuad(quads[i]->v1i, edges[tempEdge1_2]->midPointIdx, 
			quads[i]->facePointIdx, edges[tempEdge4_1]->midPointIdx,catmullQuads,catmullEdges, verts_oldEdges);

		addCatmullClarkQuad(quads[i]->v2i, edges[tempEdge2_3]->midPointIdx,
			quads[i]->facePointIdx, edges[tempEdge1_2]->midPointIdx, catmullQuads, catmullEdges, verts_oldEdges);

		addCatmullClarkQuad(quads[i]->v3i, edges[tempEdge3_4]->midPointIdx,
			quads[i]->facePointIdx, edges[tempEdge2_3]->midPointIdx, catmullQuads, catmullEdges, verts_oldEdges);

		addCatmullClarkQuad(quads[i]->v4i, edges[tempEdge4_1]->midPointIdx,
			quads[i]->facePointIdx, edges[tempEdge3_4]->midPointIdx, catmullQuads, catmullEdges, verts_oldEdges);

		
	}
	//printf("\n\n We will erase old verts \n\n");
	for (int i = 0; i < oldVertsIdx + 1; i++) {
		verts[i]->edgeList.erase(verts[i]->edgeList.begin(), verts[i]->edgeList.begin() +
			verts_oldEdges[i] );
	}
	
	copyEdges(catmullEdges);
	////printf("catmullEdges.size() = %d\n", catmullEdges.size());
	copyQuads(catmullQuads);
	catmullEdges.clear();
	catmullQuads.clear();
	delete[] tempVert;
	////printf("******Catmull Clark Subdivision Completed!\n");
}

void Mesh::Sqrt3SubDiv(bool checkNormals) {


	//TODO: 1 - First create face point for each triangles
	//TODO: 2 - Triangulate with this new point
	//TODO: 3 - Relax old vertices
	//TODO: 4 - Flip old edges
	int oldTrisIndex, oldEdgesIndex, oldVertsIdx;
	oldTrisIndex = tris.size();
	oldEdgesIndex = edges.size();
	oldVertsIdx = verts.size();
	float* tempVert;
	tempVert = new float[3];
	vector<Triangle*> sqrt3tris;
	vector<Edge*> sqrt3edges;



	float* an;
	an = new float[oldVertsIdx];
	float** pi = new float* [oldVertsIdx];
	for (int i = 0; i < oldVertsIdx; i++) {
		pi[i] = new float[3];
		pi[i][0] = 0.0f;
		pi[i][1] = 0.0f;
		pi[i][2] = 0.0f;
	}
	int n=-1;
	for (int i = 0; i < oldVertsIdx; i++) {
		n = verts[i]->vertList.size();
		an[i] = (float) (4 - 2 * cos(360.0f / n)) / 9.0f;
		//printf("an[%d] = %.6g\n", i, an[i]);
		for (int j = 0; j < n; j++) {
			pi[i][0] += verts[verts[i]->vertList[j]]->coords[0] / n;
			pi[i][1] += verts[verts[i]->vertList[j]]->coords[1] / n;
			pi[i][2] += verts[verts[i]->vertList[j]]->coords[2] / n;

		}
	}

	for (int i = 0; i < verts.size(); i++) {
		verts[i]->vertList.clear();
		verts[i]->edgeList.clear();
		verts[i]->triList.clear();
	}

	//Creating mid points
	//printf("***** Creating mid points\n\n");
	for (int i = 0; i < tris.size(); i++) {
		tempVert = calculateMidPoint(verts[tris[i]->v1i]->coords, verts[tris[i]->v2i]->coords);
		tempVert = calculateMidPoint(tempVert, verts[tris[i]->v3i]->coords, 2.0f, 1.0f);
		//printf("v1i : %.5g, %.5g, %.5g\nv2i : %.5g, %.5g, %.5g\nv3i : %.5g, %.5g, %.5g\n",
		//	verts[tris[i]->v1i]->coords[0], verts[tris[i]->v1i]->coords[1], verts[tris[i]->v1i]->coords[2],
		//	verts[tris[i]->v2i]->coords[0], verts[tris[i]->v2i]->coords[1], verts[tris[i]->v2i]->coords[2],
		//	verts[tris[i]->v3i]->coords[0], verts[tris[i]->v3i]->coords[1], verts[tris[i]->v3i]->coords[2]);
		//printf("Mid point of them : %.5g, %.5g, %.5g\n\n", tempVert[0], tempVert[1], tempVert[2]);
		addVertex(tempVert[0], tempVert[1], tempVert[2]);
		tris[i]->midPointIdx = verts.size() - 1;
	}
	/*         v1i 
			  / \
   		     /   \    mp = (v1i+v2i+v3i)/3
			/  mp \          and our new triangles will be v1i,v2i,mp (counter-clock wise)
		v2i/_______\v3i										v2i,v3i,mp
															v3i,v1i,mp

	*/

	//Adding new triangles
	//printf("***** Adding new triangles\n\n");
	for (int i = 0; i < oldTrisIndex; i++) {
		addSqrt3Triangle(tris[i]->v1i, tris[i]->v2i, tris[i]->midPointIdx, sqrt3tris, sqrt3edges);
		addSqrt3Triangle(tris[i]->v2i, tris[i]->v3i, tris[i]->midPointIdx, sqrt3tris, sqrt3edges);
		addSqrt3Triangle(tris[i]->v3i, tris[i]->v1i, tris[i]->midPointIdx, sqrt3tris, sqrt3edges);
	}


	/*               A				  A
					/|\				 / \
				   / | \			/   \
				X /t1|t2\ Y  ===> X/_____\Y			this is what we want.
				  \  |  /		   \     /
		           \ | /			\   /
				    \|/				 \ /
					 B				  B		
	
			First, we need to detect flippable edges. These are the edges who their vertices have more than one
			mutual triangles (in fact we need 2 mutual triangles exactly).

			But the major problem is finding the correct order (i.e. A, Y, X and  B, X, Y )
			
			We know that vertices of a triangle are already ordered. 
			(index_A+2)%3=(index_X+1)%3 = (index_B) for t1.

			So if we detect index of A, then we simply need to create a new triangle with
			(index_A,(index_A+1)%3,Y). This will ensure the ordering.

			Then for the last step we must clear the old data of this triangles and edges.

	*/
	int mutualTriCounter = 0;
	int* mutualTris = new int[2];
	int* tempTriVerts_1 = new int[3];
	int* tempTriVerts_2 = new int[3];
	int A=-1, B=-1, X=-1, Y=-1, AIdx_1=4, BIdx_1=4, AIdx_2=4,BIdx_2=4;
	int XIdx, YIdx;
	bool isFlippable = false;
	bool A_before_X = false;
	bool B_before_Y = false;

	//Flipping edges
	//printf("***** Flipping edges\n\n");
	for (int i = 0; i < sqrt3edges.size(); i++) {

		if (sqrt3edges[i]->v1i >= oldVertsIdx || sqrt3edges[i]->v2i >= oldVertsIdx) {
			continue;
		}

		mutualTriCounter = 0;
		//printf("Looking for mutual triangles of vertices %d and %d\n", sqrt3edges[i]->v1i, sqrt3edges[i]->v2i);
		for (int j = 0; j < verts[sqrt3edges[i]->v1i]->triList.size(); j++) {
			//printf("Vertex %d's triList[%d] = %d\n", sqrt3edges[i]->v1i,j,verts[sqrt3edges[i]->v1i]->triList[j]);
			for (int k = 0; k < verts[sqrt3edges[i]->v2i]->triList.size(); k++) {
				//printf("Vertex %d's triList[%d] = %d\n", sqrt3edges[i]->v2i, k, verts[sqrt3edges[i]->v2i]->triList[k]);
				if (verts[sqrt3edges[i]->v1i]->triList[j] == verts[sqrt3edges[i]->v2i]->triList[k]) {
					//printf("Mutual tri counter: %d\n", mutualTriCounter);
					mutualTris[mutualTriCounter] = verts[sqrt3edges[i]->v1i]->triList[j];
					mutualTriCounter++;
				}
			}

		}
		isFlippable = false;
		//printf("\nMutual triangle count for edge %d = %d\n", i, mutualTriCounter);
		if (mutualTriCounter == 2) {
			isFlippable = true;

			//printf("*******\n\nDot product of normals: %.11g\n\n",
				//sqrt3tris[mutualTris[0]]->normals[0] * sqrt3tris[mutualTris[1]]->normals[0] +
				//sqrt3tris[mutualTris[0]]->normals[1] * sqrt3tris[mutualTris[1]]->normals[1] +
				//sqrt3tris[mutualTris[0]]->normals[2] * sqrt3tris[mutualTris[1]]->normals[2]);

			if(checkNormals)
			for (int j = 0; j < 3; j++) {
				//printf("\n*MutualTri[0]->Normals[%d]=%.8g \n*MutualTri[1]->Normals[%d]=%.8g\n\n\n", j, 
					//sqrt3tris[mutualTris[0]]->normals[j],
					//j,
					//sqrt3tris[mutualTris[1]]->normals[j]);
				//printf("Differences of normals: %.9g\n", abs(sqrt3tris[mutualTris[0]]->normals[j] - sqrt3tris[mutualTris[1]]->normals[j]));
				if (abs(sqrt3tris[mutualTris[0]]->normals[j] - sqrt3tris[mutualTris[1]]->normals[j]) > 1) {
					isFlippable = false;
					//printf("\Mutual triangles have different normals!\n\n");
					break;
				}
			}
		}
			

		if (isFlippable) {

			//printf("\n\n** Edge %d is between vertices %d and %d and it is flippable!\n\n", i,sqrt3edges[i]->v1i,sqrt3edges[i]->v2i);
			//printf("** Edge %d's mutual triangles are: %d and %d\n\n",i, mutualTris[0], mutualTris[1]);


			A = sqrt3edges[i]->v1i;
			B = sqrt3edges[i]->v2i;
			
			tempTriVerts_1[0] = sqrt3tris[mutualTris[0]]->v1i;
			tempTriVerts_1[1] = sqrt3tris[mutualTris[0]]->v2i;
			tempTriVerts_1[2] = sqrt3tris[mutualTris[0]]->v3i;

			for (int j = 0; j < 3; j++) {
				if (tempTriVerts_1[j] != A && tempTriVerts_1[j] != B) {
					X = tempTriVerts_1[j];
					XIdx = j;
				}
				if (tempTriVerts_1[j] == A) {
					AIdx_1 = j;
				}
				if (tempTriVerts_1[j] == B) {
					BIdx_1 = j;
				}
			}


			if (XIdx > AIdx_1 || (XIdx == 0 && AIdx_1 == 2)) {
				A_before_X = true;
			}


			tempTriVerts_2[0] = sqrt3tris[mutualTris[1]]->v1i;
			tempTriVerts_2[1] = sqrt3tris[mutualTris[1]]->v2i;
			tempTriVerts_2[2] = sqrt3tris[mutualTris[1]]->v3i;

			for (int j = 0; j < 3; j++) {
				if (tempTriVerts_2[j] != A && tempTriVerts_2[j] != B) {
					Y = tempTriVerts_2[j];
					YIdx = j;
				}
				if (tempTriVerts_2[j] == A) {
					AIdx_2 = j;
				}
				if (tempTriVerts_2[j] == B) {
					BIdx_2 = j;
				}
			}

		

			//printf("\n\nTriangle %d will be replaced.....\n", mutualTris[0]);
			//printf("\nVertices before replacement: %d  %d  %d\n", sqrt3tris[mutualTris[0]]->v1i,
				//sqrt3tris[mutualTris[0]]->v2i,
				//sqrt3tris[mutualTris[0]]->v3i);

			//printf("\n\nTriangle %d will be replaced.....\n", mutualTris[1]);
			//printf("\nVertices before replacement: %d  %d  %d\n", sqrt3tris[mutualTris[1]]->v1i,
				//sqrt3tris[mutualTris[1]]->v2i,
				//sqrt3tris[mutualTris[1]]->v3i);


			if (!A_before_X) {
				replaceSqrt3Triangle(mutualTris[0], A, X, Y, sqrt3tris);
				replaceSqrt3Triangle(mutualTris[1], B, Y, X, sqrt3tris);
			}
			else {
				replaceSqrt3Triangle(mutualTris[0], A, Y, X, sqrt3tris);
				replaceSqrt3Triangle(mutualTris[1], B, X, Y, sqrt3tris);
			}

			//printf("Flipped...\n");
			replaceSqrt3Edge(i, X, Y, sqrt3edges);



			//printf("\n\nTriangle %d replaced.....\n", mutualTris[0]);
			//printf("\nVertices after replacement: %d  %d  %d\n", sqrt3tris[mutualTris[0]]->v1i,
				//sqrt3tris[mutualTris[0]]->v2i,
				//sqrt3tris[mutualTris[0]]->v3i);

			//printf("\n\nTriangle %d replaced.....\n", mutualTris[1]);
			//printf("\nVertices after replacement: %d  %d  %d\n", sqrt3tris[mutualTris[1]]->v1i,
				//sqrt3tris[mutualTris[1]]->v2i,
				//sqrt3tris[mutualTris[1]]->v3i);

		}

	}


	for (int i = 0; i < oldVertsIdx; i++) {
		verts[i]->coords[0] = (1 - an[i]) * verts[i]->coords[0] + an[i] * pi[i][0];
		verts[i]->coords[1] = (1 - an[i]) * verts[i]->coords[1] + an[i] * pi[i][1];
		verts[i]->coords[2] = (1 - an[i]) * verts[i]->coords[2] + an[i] * pi[i][2];
	}



	//edges.erase(edges.begin(), edges.begin() + oldEdgesIndex);
	//tris.erase(tris.begin(), tris.begin() + oldTrisIndex);
	//printf("\nWill copy pointer vectors...\n");
	copyTriangles(sqrt3tris);
	copyEdges(sqrt3edges);
	//printf("\nSuccesfully copied pointer vectors!\n");
	sqrt3edges.clear();
	sqrt3tris.clear();
	//printf("\nSuccesfully cleared pointer vectors!\n");
	for (int i = 0; i < oldVertsIdx; i++) {
		delete[] pi[i];
	}
	delete[] an;
	delete[] tempVert;
	//printf("\nDeleted TempVect\n");
	//delete[] mutualTris;
	//printf("\nDeleted Mutual Tris\n");
	delete[3] tempTriVerts_1;
	//printf("\nDeleted tempTriVerts_1\n");
	delete[3] tempTriVerts_2;
	//printf("\nDeleted tempTriVerts_2\n");

	//printf("\nSuccesfully deleted allocated memory!\n");
}

/// <summary>
/// Will create 3 quads for each triangle. Simple Catmull-Clark split.
/// </summary>
void Mesh::convertQuadbySplit() {

	if (quads.size() > 1)
		return;

	//We will clear the edges so we can fill our new edges.
	edges.clear();
	//we will fill our vertlist with new quad vertices
	for (int i = 0; i < verts.size(); i++) {
		verts[i]->vertList.clear();
		verts[i]->edgeList.clear();
		//we will use this for not creating the same vertices again and again, see below

		/*	v1
			/\ we will clear all neighbours so v1 and v2 are not neighbours anymore
		 v2/__\v3

			but after than we will split the v1-v2 edge and add new v4 vertex there. and we will repeat
			that for all edges.
			so the new structure will be like;
			  v1
			/   \
		   v4 v7 v6    and in here note that v1's new neighbours are v4 and v6
		v2/___v5___\v3 and v3's new neighbours are v6 and v5
					   thus when creating a new vertex in midpoints, wefirst need to check if v1 and v3
					   has any common neighbour. if so, that neighbour's idx will be our previously created
					   v6 index because we cleaned all neighbour list!
					   then we will not add a new vertex but just add v6 to our quad by obtaining that index.

		*/
	}

	int idxcenter, idx1_2, idx2_3, idx1_3;
	for (int i = 0; i < tris.size(); i++) {
		float* center, * mid1_2, * mid2_3, * mid1_3;
		center = new float[3];
		mid1_2 = new float[3];
		mid2_3 = new float[3];
		mid1_3 = new float[3];

		idx1_2 = isCommonNeighbourExist(tris[i]->v1i, tris[i]->v2i);

		//if v1i and v2i does not have any common vertex (means the are not splitted before)
		//we will begin the process. if they have common neighbour we will just add new quad
		if (idx1_2 == -1) {
			//middle point between first and second vertices of triangle (counter-clockwise)
			mid1_2 = calculateMidPoint(verts[tris[i]->v1i]->coords, verts[tris[i]->v2i]->coords);
			addVertex(mid1_2[0], mid1_2[1], mid1_2[2]);
			idx1_2 = verts.size() - 1;

		}



		idx2_3 = isCommonNeighbourExist(tris[i]->v2i, tris[i]->v3i);
		if (idx2_3 == -1) {

			//middle point between second and third vertices of triangle (counter-clockwise)
			mid2_3 = calculateMidPoint(verts[tris[i]->v2i]->coords, verts[tris[i]->v3i]->coords);
			addVertex(mid2_3[0], mid2_3[1], mid2_3[2]);
			idx2_3 = verts.size() - 1;

		}


		idx1_3 = isCommonNeighbourExist(tris[i]->v1i, tris[i]->v3i);

		if (idx1_3 == -1) {
			//middle point between first and second vertices of triangle (counter-clockwise)
			mid1_3 = calculateMidPoint(verts[tris[i]->v1i]->coords, verts[tris[i]->v3i]->coords);
			addVertex(mid1_3[0], mid1_3[1], mid1_3[2]);
			idx1_3 = verts.size() - 1;

		}

		//centers are unique for each triangle so we can continue safely

		center = calculateMidPoint(verts[tris[i]->v1i]->coords, verts[tris[i]->v2i]->coords);
		center = calculateMidPoint(center, verts[tris[i]->v3i]->coords, 2.0f, 1.0f);
		addVertex(center[0], center[1], center[2]);
		idxcenter = verts.size() - 1;


		//Now add quads. One triangle will be splitted to 3 quads.
		addQuad(tris[i]->v1i, idx1_2, idxcenter, idx1_3);
		addQuad(idx1_2, tris[i]->v2i, idx2_3, idxcenter); //order arrenged counter-clockwise
		addQuad(idx2_3, tris[i]->v3i, idx1_3, idxcenter); //order arrenged counter-clockwise

		//memory leakage prevention
		delete[3] center;
		delete[3] mid1_2;
		delete[3] mid2_3;
		delete[3] mid1_3;
	}

}

void Mesh::copyEdges(const vector<Edge*> edgev) {
	edges = edgev;
}
void Mesh::copyQuads(const vector<Quad*> quadv) {
	quads = quadv;
}
void Mesh::copyTriangles(const vector<Triangle*> triangleV) {
	tris = triangleV;
}


void Mesh::createCube(float sideLen)
{
	//coordinates
	float flbc[3] = {0, 0, 0}, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
			case 1:
				deltaX = sideLen;
				break;
			case 2:
				deltaZ = -sideLen;
				break;
			case 3:
				deltaX = 0;
				break;
			case 4:
				deltaZ = 0;
				deltaY = sideLen;
				break;
			case 5:
				deltaX = sideLen;
				break;
			case 6:
				deltaZ = -sideLen;
				break;
			default:
				deltaX = 0;;
				break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	addTriangle(0, 1, 5);
	addTriangle(0, 5, 4);
}

/// <summary>
/// creates cube with each surface as quad
/// </summary>
/// <param name="sideLength"></param>
void Mesh::createCubeQuad(float sideLength) {

	float flbc[3] = { 0, 0, 0 }, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
		case 1:
			deltaX = sideLength;
			break;
		case 2:
			deltaZ = -sideLength;
			break;
		case 3:
			deltaX = 0;
			break;
		case 4:
			deltaZ = 0;
			deltaY = sideLength;
			break;
		case 5:
			deltaX = sideLength;
			break;
		case 6:
			deltaZ = -sideLength;
			break;
		default:
			deltaX = 0;;
			break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addQuad(0, 1, 2, 3);
	addQuad(0, 1, 5, 4);
	addQuad(4, 5, 6, 7);
	addQuad(3, 2, 6, 7);
	addQuad(1, 2, 6, 5);
	addQuad(0, 3, 7, 4);


}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back( new Triangle(idx, v1, v2, v3) );

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	if (!makeVertsNeighbor(v1, v2)) {
		addEdge(v1, v2);
		////printf("Edges size: %d", edges.size());
		edges[edges.size()-1]->length = calculateEdgeLength(v1, v2);
	}
		

	if (!makeVertsNeighbor(v1, v3)) {
		addEdge(v1, v3);
		edges[edges.size()-1 ]->length = calculateEdgeLength(v1, v3);
	}

	if (!makeVertsNeighbor(v2, v3)) {
		addEdge(v2, v3);
		edges[edges.size()-1]->length = calculateEdgeLength(v2, v3);

	}

}

void Mesh::addSqrt3Triangle(int v1, int v2, int v3, vector<Triangle*>& sqrt3Triangles, vector<Edge*>& sqrt3Edge) {
	int idx = sqrt3Triangles.size();
	sqrt3Triangles.push_back(new Triangle(idx, v1, v2, v3));


	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);
	//printf("\nAded triangle %d to vertex %d \n", idx, v1);
	//printf("Aded triangle %d to vertex %d \n", idx, v2);
	//printf("Aded triangle %d to vertex %d \n\n", idx, v3);

	if (!makeVertsNeighbor(v1, v2)) {
		addSqrt3Edge(v1, v2, sqrt3Edge);
	}


	if (!makeVertsNeighbor(v1, v3)) {
		addSqrt3Edge(v1, v3, sqrt3Edge);
	}

	if (!makeVertsNeighbor(v2, v3)) {
		addSqrt3Edge(v2, v3, sqrt3Edge);
	}

	setSqrt3TriNormal(idx, sqrt3Triangles);
}

void Mesh::removeSqrt3Triangle(int triIdx, vector<Triangle*>& sqrt3Triangle) {

	printf("** Triangle requested for removal is %d and its vertices are %d, %d and %d\n", triIdx,
		sqrt3Triangle[triIdx]->v1i, sqrt3Triangle[triIdx]->v2i, sqrt3Triangle[triIdx]->v3i);

	printf("\nverts[tris[%d]->v1i]->triList.size() = %d\n\n",triIdx, verts[tris[triIdx]->v1i]->triList.size());

	for (int i = 0; i < verts[sqrt3Triangle[triIdx]->v1i]->triList.size(); i++) {
		printf("%d th loop for v1\n", i);
		if (verts[sqrt3Triangle[triIdx]->v1i]->triList[i] == triIdx) {
			printf("Removed the triangle from v1 \n");
			verts[sqrt3Triangle[triIdx]->v1i]->triList.erase(verts[sqrt3Triangle[triIdx]->v1i]->triList.begin() + i);
			break;
		}
	}

	for (int i = 0; i < verts[sqrt3Triangle[triIdx]->v2i]->triList.size(); i++) {
		printf("%d th loop for v2\n", i);

		if (verts[sqrt3Triangle[triIdx]->v2i]->triList[i] == triIdx) {
			printf("Will try to remove the triangle from %d \n",sqrt3Triangle[triIdx]->v2i);

			verts[sqrt3Triangle[triIdx]->v2i]->triList.erase(verts[sqrt3Triangle[triIdx]->v2i]->triList.begin() + i);

			printf("Removed the triangle from v2 \n");

			break;
		}

	}


	for (int i = 0; i < verts[sqrt3Triangle[triIdx]->v3i]->triList.size(); i++) {
		printf("%d th loop for v3\n", i);

		if (verts[sqrt3Triangle[triIdx]->v3i]->triList[i] == triIdx) {
			printf("Removed the triangle from v3 \n");

			verts[sqrt3Triangle[triIdx]->v3i]->triList.erase(verts[sqrt3Triangle[triIdx]->v3i]->triList.begin() + i);
			break;
		}

	}

	//sqrt3Triangle.erase(sqrt3Triangle.begin() + triIdx);
	return;	
}

void Mesh::replaceSqrt3Triangle(int triIdx, int v1, int v2, int v3, vector<Triangle*>& sqrt3Triangle) {

	//clear the old vertices triLists
	for (int i = 0; i < verts[sqrt3Triangle[triIdx]->v1i]->triList.size(); i++) {
		if (verts[sqrt3Triangle[triIdx]->v1i]->triList[i] == triIdx) {
			verts[sqrt3Triangle[triIdx]->v1i]->triList.erase(verts[sqrt3Triangle[triIdx]->v1i]->triList.begin() + i);
			break;
		}
	}

	for (int i = 0; i < verts[sqrt3Triangle[triIdx]->v2i]->triList.size(); i++) {
		if (verts[sqrt3Triangle[triIdx]->v2i]->triList[i] == triIdx) {
			verts[sqrt3Triangle[triIdx]->v2i]->triList.erase(verts[sqrt3Triangle[triIdx]->v2i]->triList.begin() + i);
			break;
		}

	}

	for (int i = 0; i < verts[sqrt3Triangle[triIdx]->v3i]->triList.size(); i++) {
		if (verts[sqrt3Triangle[triIdx]->v3i]->triList[i] == triIdx) {
			verts[sqrt3Triangle[triIdx]->v3i]->triList.erase(verts[sqrt3Triangle[triIdx]->v3i]->triList.begin() + i);
			break;
		}

	}

	//replace them
	sqrt3Triangle[triIdx]->v1i = v1;
	verts[v1]->triList.push_back(triIdx);

	sqrt3Triangle[triIdx]->v2i = v2;
	verts[v2]->triList.push_back(triIdx);

	sqrt3Triangle[triIdx]->v3i = v3;
	verts[v3]->triList.push_back(triIdx);

	setSqrt3TriNormal(triIdx, sqrt3Triangle);

}

bool Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->vertList.size(); i++)
		if (verts[v1i]->vertList[i] == v2i)
			return true;

	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addVertex(float x, float y, float z)
{

	////printf("New added vertex coordinats : %.5g, %.5g, %.5g \n", x, y, z);
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;
	////printf("Vertex num: %d coordinats : % f % f % f\n",idx, x, y, z);
	verts.push_back( new Vertex(idx, c) );
	
}

/// <summary>
/// During the Catmull Clark subdivision, we will store our new and temp edges with help of
/// this function.
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <param name="cmEdges"></param>
/// <param name="quadIdx"></param>
void Mesh::addCatmullClarkEdge(int v1, int v2 , vector<Edge*>& cmEdges, int quadIdx) {
	int idx = cmEdges.size();
	//printf("catmull edge index with %d created between: %d and %d with QUAD: %d\n",idx, v1, v2,quadIdx);
	cmEdges.push_back(new Edge(idx, v1, v2));
	cmEdges[idx]->quadIdx.push_back(quadIdx);

	
	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

void Mesh::addEdge(int v1, int v2, int quadIdx)
{
	int idx = edges.size();
	//printf("Edge created between: %d and %d with index of %d\n", v1, v2,idx);
	////printf("Edge %d has new quad %d \n", idx, quadIdx);
	edges.push_back( new Edge(idx, v1, v2) );
	if (quadIdx != -1) {
		edges[idx]->quadIdx.push_back(quadIdx);
	}
	verts[v1]->edgeList.push_back(idx);
	////printf("verts[%d]->edgelist added %d\n", v1, idx);
	verts[v2]->edgeList.push_back(idx);
	////printf("verts[%d]->edgelist added %d\n", v2, idx);

}

void Mesh::addSqrt3Edge(int v1, int v2, vector<Edge*>& sqrt3Edges) {
	int idx = sqrt3Edges.size();

	sqrt3Edges.push_back(new Edge(idx, v1, v2));

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

void Mesh::replaceSqrt3Edge(int edgeIdx, int v1, int v2, vector<Edge*> sqrt3Edges) {

	//clear old vertices edgeLists
	for (int i = 0; i < verts[sqrt3Edges[edgeIdx]->v1i]->edgeList.size(); i++) {
		if (verts[sqrt3Edges[edgeIdx]->v1i]->edgeList[i] == edgeIdx) {
			verts[sqrt3Edges[edgeIdx]->v1i]->edgeList.erase(verts[sqrt3Edges[edgeIdx]->v1i]->edgeList.begin() + i);
		}
	}

	for (int i = 0; i < verts[sqrt3Edges[edgeIdx]->v2i]->edgeList.size(); i++) {
		if (verts[sqrt3Edges[edgeIdx]->v2i]->edgeList[i] == edgeIdx) {
			verts[sqrt3Edges[edgeIdx]->v2i]->edgeList.erase(verts[sqrt3Edges[edgeIdx]->v2i]->edgeList.begin() + i);
		}
	}

	for (int i = 0; i < verts[sqrt3Edges[edgeIdx]->v1i]->vertList.size(); i++) {
		if (verts[sqrt3Edges[edgeIdx]->v1i]->vertList[i] == sqrt3Edges[edgeIdx]->v2i) {
			verts[sqrt3Edges[edgeIdx]->v1i]->vertList.erase(verts[sqrt3Edges[edgeIdx]->v1i]->vertList.begin() + i);
		}
	}

	for (int i = 0; i < verts[sqrt3Edges[edgeIdx]->v2i]->vertList.size(); i++) {
		if (verts[sqrt3Edges[edgeIdx]->v2i]->vertList[i] == sqrt3Edges[edgeIdx]->v1i) {
			verts[sqrt3Edges[edgeIdx]->v2i]->vertList.erase(verts[sqrt3Edges[edgeIdx]->v2i]->vertList.begin() + i);
		}
	}

	sqrt3Edges[edgeIdx]->v1i = v1;
	sqrt3Edges[edgeIdx]->v2i = v2;
	verts[v1]->edgeList.push_back(edgeIdx);
	verts[v2]->edgeList.push_back(edgeIdx);
	verts[v1]->vertList.push_back(v2);
	verts[v2]->vertList.push_back(v1);


}

void Mesh::removeSqrt3Edge(int edgeIdx, vector<Edge*>& sqrt3Edges) {

	for (int i = 0; i < verts[sqrt3Edges[edgeIdx]->v1i]->edgeList.size(); i++) {
		if (verts[sqrt3Edges[edgeIdx]->v1i]->edgeList[i] == edgeIdx) {
			verts[sqrt3Edges[edgeIdx]->v1i]->edgeList.erase(verts[sqrt3Edges[edgeIdx]->v1i]->edgeList.begin() + i);
		}
	}

	for (int i = 0; i < verts[sqrt3Edges[edgeIdx]->v2i]->edgeList.size(); i++) {
		if (verts[sqrt3Edges[edgeIdx]->v2i]->edgeList[i] == edgeIdx) {
			verts[sqrt3Edges[edgeIdx]->v2i]->edgeList.erase(verts[sqrt3Edges[edgeIdx]->v2i]->edgeList.begin() + i);
		}
	}

	//sqrt3Edges.erase(sqrt3Edges.begin() + edgeIdx);
	return;
}

float Mesh::calculateEdgeLength(int v1i, int v2i) {
	float* v1 = new float[3];
	float* v2 = new float[3];
	memcpy_s(v1, 3 * sizeof(float), verts[v1i]->coords, 3 * sizeof(float));
	memcpy_s(v2, 3 * sizeof(float), verts[v2i]->coords, 3 * sizeof(float));

	float distance = (
		pow(v1[0] - v2[0], 2) +
		pow(v1[1] - v2[1], 2) +
		pow(v1[2] - v2[2], 2));
	////printf("***\n\t v1_x: %f, v1_y: %f, v1_z: %f\n", v1[0], v1[1], v1[2]);
	////printf("***\n\t v2_x: %f, v2_y: %f, v2_z: %f\n", v2[0], v2[1], v2[2]);
	////printf("***\t\t distance: %f", distance);
	delete[3] v1;
	delete[3] v2;
	return sqrt(distance);

}

/// <summary>
/// Calculates mid point with respect to given coefficients.
/// One may use this function for calculating mid points ofany number of vertices with the help of the coefficients.
/// </summary>
/// <param name="c1"></param>
/// <param name="c2"></param>
/// <param name="c1_coeff"></param>
/// <param name="c2_coeff"></param>
/// <returns></returns>
float* Mesh::calculateMidPoint(float* c1, float* c2, float c1_coeff, float c2_coeff) {
	float* midP = new float[3];
	midP[0] = (c1_coeff*c1[0] + c2_coeff*c2[0]) /(c1_coeff+c2_coeff);
	midP[1] = (c1_coeff * c1[1] + c2_coeff * c2[1]) / (c1_coeff + c2_coeff);
	midP[2] = (c1_coeff * c1[2] + c2_coeff * c2[2]) / (c1_coeff + c2_coeff);
	return midP;
}

/// <summary>
/// Returns corresponding edge for given two vertices
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
int Mesh::getEdgeFromVerts(int v1, int v2) {
	for (int i = 0; i<verts[v1]->edgeList.size() ; i++) {
		

		if (edges[verts[v1]->edgeList[i]]->v1i == v1 && edges[verts[v1]->edgeList[i]]->v2i == v2 ||
			edges[verts[v1]->edgeList[i]]->v1i == v2 && edges[verts[v1]->edgeList[i]]->v2i == v1)

			return verts[v1]->edgeList[i];
	}
	return -1;
}

/// <summary>
/// During Catmull-Clark subdivision, we will have temporary edges list.
/// This function will return corresponding edge from that list, for given two vertices.
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <param name="cmEdge"></param>
/// <param name="verts_oldEdgeList"></param>
/// <returns></returns>
int Mesh::getCatmullClarkEdgeFromVerts(int v1, int v2, vector<Edge*>& cmEdge, vector<int>& verts_oldEdgeList) {
	//printf("We will get the edge for verts %d and %d \n",v1, v2);
	//printf("verts_oldEdgeList.size() = %d\n", verts_oldEdgeList.size());
	if (v1 < verts_oldEdgeList.size()) {
		for (int i = verts_oldEdgeList[v1]; i < verts[v1]->edgeList.size(); i++) {
			//printf("verts[%d].edgeList[%d] = %d\n", v1, i, verts[v1]->edgeList[i]);
			for (int j = 0; j < verts[v2]->edgeList.size(); j++) {
				//printf("verts[%d].edgeList[%d] = %d\n", v2, j, verts[v2]->edgeList[j]);


				if (verts[v1]->edgeList[i] == verts[v2]->edgeList[j]) {
					//printf("Edge found = %d\n", verts[v1]->edgeList[i]);
					return verts[v1]->edgeList[i];
				}


			}

		}
	}
	else if (v2 < verts_oldEdgeList.size()) {
		for (int i = 0; i < verts[v1]->edgeList.size(); i++) {
			//printf("verts[%d].edgeList[%d] = %d\n", v1, i, verts[v1]->edgeList[i]);
			for (int j = verts_oldEdgeList[v2]; j < verts[v2]->edgeList.size(); j++) {
				//printf("verts[%d].edgeList[%d] = %d\n", v2, j, verts[v2]->edgeList[j]);


				if (verts[v1]->edgeList[i] == verts[v2]->edgeList[j]) {
					//printf("Edge found = %d\n", verts[v1]->edgeList[i]);
					return verts[v1]->edgeList[i];
				}


			}

		}
	}
	else {
		for (int i = 0; i < verts[v1]->edgeList.size(); i++) {
			//printf("verts[%d].edgeList[%d] = %d\n", v1, i, verts[v1]->edgeList[i]);
			for (int j = 0; j < verts[v2]->edgeList.size(); j++) {
				//printf("verts[%d].edgeList[%d] = %d\n", v2, j, verts[v2]->edgeList[j]);


				if (verts[v1]->edgeList[i] == verts[v2]->edgeList[j]) {
					//printf("Edge found = %d\n", verts[v1]->edgeList[i]);
					return verts[v1]->edgeList[i];
				}


			}

		}
	}

	
	//printf("Edge not found!\n");
	return -1;
}

/// <summary>
/// Input parameters should be in counter-clockwise order.
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <param name="v3"></param>
/// <param name="v4"></param>
void Mesh::addQuad(int v1, int v2, int v3, int v4) {
	int idx = quads.size();
	quads.push_back(new Quad(idx, v1, v2, v3, v4));
	////printf("Add quad called for %d , %d , %d , %d\n", v1, v2, v3, v4);

	if (!makeVertsNeighbor(v1, v2))
		addEdge(v1, v2, idx);
	else {
		int temp = getEdgeFromVerts(v1, v2);
		//printf("These verts %d and %d already edge in edge: %d \n",v1,v2, temp);
		if (temp == -1) {
			//printf("temp == -1\n");
		}
		else {
			if (std::count(edges[temp]->quadIdx.begin(), edges[temp]->quadIdx.end(), idx) == 0)
			{
				edges[temp]->quadIdx.push_back(idx);
				////printf("Edge %d has new quad: %d \n", temp, idx);
			}
		}
	}

	if (!makeVertsNeighbor(v2, v3))
		addEdge(v2, v3, idx);
	else {
		int temp = getEdgeFromVerts(v2, v3);
		////printf("These verts %d and  %d are already edge in edge: %d \n",v2,v3, temp);
		if (temp == -1) {
			//printf("temp == -1\n");
		}
		else {
			if (std::count(edges[temp]->quadIdx.begin(), edges[temp]->quadIdx.end(), idx) == 0)
			{
				edges[temp]->quadIdx.push_back(idx);
				////printf("Edge %d has new quad: %d \n", temp, idx);
			}
		}
	}

	if (!makeVertsNeighbor(v3, v4))
		addEdge(v3, v4, idx);
	else {

		int temp = getEdgeFromVerts(v3, v4);
		////printf("These verts %d and %d are already edge in edge: %d \n",v3,v4, temp);

		if (temp == -1) {
			//printf("temp == -1\n");
		}
		else {
			if (std::count(edges[temp]->quadIdx.begin(), edges[temp]->quadIdx.end(), idx) == 0)
			{
				edges[temp]->quadIdx.push_back(idx);
				////printf("Edge %d has new quad: %d \n", temp, idx);
			}
		}
	}

	if (!makeVertsNeighbor(v1, v4))
		addEdge(v1, v4, idx);
	else {
		////printf(" %d , %d are already neighbour!\n", v1, v4);
		int temp = getEdgeFromVerts(v1, v4);
		////printf("These verts %d and %d already edge in edge: %d \n",v1,v4, temp);
		if (temp == -1) {
			//printf("temp == -1\n");
		}
		else {
			if (std::count(edges[temp]->quadIdx.begin(), edges[temp]->quadIdx.end(), idx) == 0)
			{
				edges[temp]->quadIdx.push_back(idx);
				////printf("Edge %d has new quad: %d \n", temp, idx);
			}
		}
	}

	verts[v1]->quadList.push_back(idx);
	verts[v2]->quadList.push_back(idx);
	verts[v3]->quadList.push_back(idx);
	verts[v4]->quadList.push_back(idx);
}

void Mesh::addCatmullClarkQuad(int v1, int v2, int v3, int v4, vector<Quad*>& cmQuads, vector<Edge*> &cmEdges, vector<int>& verts_oldEdgeList) {
	int idx = cmQuads.size();
	cmQuads.push_back(new Quad(idx, v1, v2, v3, v4));
	//printf("Quad will be created with verts %d , %d , %d , %d", v1, v2, v3, v4);
	if (!makeVertsNeighbor(v1, v2))
		addCatmullClarkEdge(v1, v2, cmEdges, idx);
	else {
		int temp = getCatmullClarkEdgeFromVerts(v1, v2,cmEdges, verts_oldEdgeList);
		if (temp != -1) {
			//printf("edge %d will be added new quad %d\n", temp, idx);
			cmEdges[temp]->quadIdx.push_back(idx);
		}
	}
	if (!makeVertsNeighbor(v2, v3))
		addCatmullClarkEdge(v2, v3, cmEdges, idx);
	else {
		int temp = getCatmullClarkEdgeFromVerts(v2, v3, cmEdges, verts_oldEdgeList);
		if (temp != -1) {
			//printf("edge %d will be added new quad %d\n", temp, idx);
			cmEdges[temp]->quadIdx.push_back(idx);
			//printf("edge %d added new quad %d\n", temp, idx);
		}
	}
	if (!makeVertsNeighbor(v3, v4))
		addCatmullClarkEdge(v3, v4, cmEdges, idx);
	else {
		int temp = getCatmullClarkEdgeFromVerts(v3, v4, cmEdges, verts_oldEdgeList);
		if (temp != -1) {
			//printf("edge %d will be added new quad %d\n", temp, idx);
			cmEdges[temp]->quadIdx.push_back(idx);
			//printf("edge %d added new quad %d\n", temp, idx);
		}
	}
	if (!makeVertsNeighbor(v1, v4))
		addCatmullClarkEdge(v1, v4, cmEdges, idx);
	else {
		int temp = getCatmullClarkEdgeFromVerts(v1, v4, cmEdges, verts_oldEdgeList);
		if (temp != -1) {
			//printf("edge %d will be added new quad %d\n", temp, idx);
			cmEdges[temp]->quadIdx.push_back(idx);
			//printf("edge %d added new quad %d\n", temp, idx);
		}
	}

	verts[v1]->quadList.push_back(idx);
	verts[v2]->quadList.push_back(idx);
	verts[v3]->quadList.push_back(idx);
	verts[v4]->quadList.push_back(idx);
}

int Mesh::isCommonNeighbourExist(int v1, int v2) {
	for (int i = 0; i < verts[v1]->vertList.size(); i++) {
		if (std::find(verts[v2]->vertList.begin(), verts[v2]->vertList.end(),verts[v1]->vertList[i]) != verts[v2]->vertList.end()) {
			return verts[v1]->vertList[i];
		}
	}
	return -1; //no neighbour detected
}

void Mesh::setTriNormal(int triIdx) {
	int v1, v2, v3;
	float Ux, Uy, Uz, Vx, Vy, Vz;
	v1 = tris[triIdx]->v1i;
	v2 = tris[triIdx]->v2i;
	v3 = tris[triIdx]->v3i;

	Ux = verts[v2]->coords[0] - verts[v1]->coords[0];
	Uy = verts[v2]->coords[1] - verts[v1]->coords[1];
	Uz = verts[v2]->coords[2] - verts[v1]->coords[2];

	Vx = verts[v3]->coords[0] - verts[v1]->coords[0];
	Vy = verts[v3]->coords[1] - verts[v1]->coords[1];
	Vz = verts[v3]->coords[2] - verts[v1]->coords[2];

	tris[triIdx]->normals[0] = Uy * Vz - Uz * Vy;
	tris[triIdx]->normals[1] = Uz * Vx - Ux * Vz;
	tris[triIdx]->normals[2] = Ux * Vy - Uy * Vx;
	return;
}

void Mesh::setSqrt3TriNormal(int triIdx, vector<Triangle*>& sqrt3Triangles) {
	int v1, v2, v3;
	float Ux, Uy, Uz, Vx, Vy, Vz;
	v1 = sqrt3Triangles[triIdx]->v1i;
	v2 = sqrt3Triangles[triIdx]->v2i;
	v3 = sqrt3Triangles[triIdx]->v3i;

	Ux = verts[v2]->coords[0] - verts[v1]->coords[0];
	Uy = verts[v2]->coords[1] - verts[v1]->coords[1];
	Uz = verts[v2]->coords[2] - verts[v1]->coords[2];

	Vx = verts[v3]->coords[0] - verts[v1]->coords[0];
	Vy = verts[v3]->coords[1] - verts[v1]->coords[1];
	Vz = verts[v3]->coords[2] - verts[v1]->coords[2];

	sqrt3Triangles[triIdx]->normals[0] = Uy * Vz - Uz * Vy;
	sqrt3Triangles[triIdx]->normals[1] = Uz * Vx - Ux * Vz;
	sqrt3Triangles[triIdx]->normals[2] = Ux * Vy - Uy * Vx;
	return;
}

void Mesh::removeTriFromVert(int v, int tri) {
	for (int i = 0; i < verts[v]->triList.size(); i++) {
		if (verts[v]->triList[i] == tri) {
			verts[v]->triList.erase(verts[v]->triList.begin() + i);
			return;
		}
	}

}