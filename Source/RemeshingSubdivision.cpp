#include "Mesh.h"
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <string>
#include <iostream>
#include "CmdConfigs.h"
#include <filesystem>
#include <fstream>

#define CHECK_FILE_EXTENSION(meshPath) \
do { \
	std::filesystem::path __MACRO_PATH__ = std::filesystem::path(meshPath);\
	bool __MACRO_EXTENSION_COMP__ = false; \
	for(auto& ext: SupportedFileFormats){ \
		if (ext == __MACRO_PATH__.extension()) { \
			__MACRO_EXTENSION_COMP__ = true; \
		}\
	}\
	if (!__MACRO_EXTENSION_COMP__) {\
		std::cout << "Unsupported file format: "+__MACRO_PATH__.extension().string() + ".\n" << "Supported formats are: \n";\
		for (auto& format : SupportedFileFormats) {\
			std::cout << "." + format + "\n";\
		}\
		return 0;\
	}\
} while (0);\

int main(int argc, char** argv)
{
	if (argc == 1 || argc == 2) {
		std::cout << "Options:\n";
		std::cout << catmull_clark_quad_mesh + "\t\tApplies Catmull-Clark algorithm to given quad mesh\n";
		std::cout << sqrt3_tri_mesh + "\t\tApplies Sqrt3 algorithm to given triangular mesh\n";
		return 0;
	}

	if (argc == 3) {
		uint64_t preSubdivVertexCount = 0;
		uint64_t postSubdivVertexCount = 0;
		std::string arg = argv[1];
		std::string filePath = argv[2];
		std::filesystem::path file(filePath);
		if (std::filesystem::exists(file) == false) {
			std::cout << "Could not find file:" << file << std::endl;
			return 0;
		}
		CHECK_FILE_EXTENSION(filePath);
		Mesh mesh;
		if (arg == catmull_clark_quad_mesh) {
			mesh.loadOffQuadMesh(file.string().c_str());
			preSubdivVertexCount = mesh.verts.size();
			mesh.CatmullClarkSubdiv();
			postSubdivVertexCount = mesh.verts.size();
			std::filesystem::path outputPath = file.parent_path() / (file.filename().string() + "_result.off");
			std::ofstream resultFile(outputPath);
			if (!resultFile.is_open()) {
				std::cerr << "Failed to open file for writing.\n";
				return 1;
			}

			resultFile << "OFF\n";
			resultFile << mesh.verts.size() << " " << mesh.quads.size() << " " << mesh.edges.size() << "\n";

			for (const auto& vert : mesh.verts) {
				resultFile << vert->coords[0] << " " << vert->coords[1] << " " << vert->coords[2] << "\n";
			}

			for (const auto& quad : mesh.quads) {
				resultFile << "4 " << quad->v1i << " " << quad->v2i << " " << quad->v3i << " " << quad->v4i << "\n";
			}
		}
		else if (arg == sqrt3_tri_mesh) {
			mesh.loadOffTriMesh(file.string().c_str());
			preSubdivVertexCount = mesh.verts.size();
			mesh.Sqrt3SubDiv();
			postSubdivVertexCount = mesh.verts.size();
			std::filesystem::path outputPath = file.parent_path() / (file.filename().string() + "_result.off");
			std::ofstream resultFile(outputPath);
			if (!resultFile.is_open()) {
				std::cerr << "Failed to open file for writing.\n";
				return 1;
			}

			resultFile << "OFF\n";
			resultFile << mesh.verts.size() << " " << mesh.tris.size() << " " << mesh.edges.size() << "\n";

			for (const auto& vert : mesh.verts) {
				resultFile << vert->coords[0] << " " << vert->coords[1] << " " << vert->coords[2] << "\n";
			}

			for (const auto& triangle : mesh.tris) {
				resultFile << "3 " << triangle->v1i << " " << triangle->v2i << " " << triangle->v3i << " " << "\n";
			}
		}
		std::cout<<"***************\nSubdivisions Completed with stats:\
			\nVertex Count before subdivision: " + std::to_string(preSubdivVertexCount) + "\n"<<\
			"Vertex Count after subdivision: " + std::to_string(postSubdivVertexCount) + "\n"<<\
			"***************\n";
	}
	return 0;
	
}