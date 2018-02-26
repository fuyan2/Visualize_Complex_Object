/***********************************************************
Starter code for Assignment 3

This code was originally written by Jack Wang for
CSC418, SPRING 2005

***********************************************************/







#include "raytracer.h"
#include "OBJ_Loader.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

using namespace std;

//set the Tree to global variable
Octree Tree(Point3D(0,0,0), Vector3D(2000, 2000, 2000));

int main(int argc, char* argv[])
{
	// Build your scene and setup your camera here, by calling
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the
	// assignment.
	Raytracer raytracer;
	int width = 320;
	int height = 240;

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold(Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
		Colour(0.628281, 0.555802, 0.366065),
		51.2,0.6);
	Material jade(Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
		Colour(0.316228, 0.316228, 0.316228),
		12.8,0.4);

	// Defines a point light source.
	raytracer.addLightSource(new PointLight(Point3D(0, 0, 5),
		Colour(0.9, 0.9, 0.9)));

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere1 = raytracer.addObject(new UnitSphere(), &gold);
	//SceneDagNode* sphere2 = raytracer.addObject(new UnitSphere(), &gold);
	SceneDagNode* plane = raytracer.addObject(new UnitSquare(), &jade);

	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere1, Vector3D(0, 0, -5));
	raytracer.rotate(sphere1, 'x', -45);
	raytracer.rotate(sphere1, 'z', 45);
	raytracer.scale(sphere1, Point3D(0, 0, 0), factor1);

	// raytracer.translate(sphere2, Vector3D(0, 0, 5));
	// raytracer.rotate(sphere2, 'x', 45);
	// raytracer.rotate(sphere2, 'z', -45);
	// raytracer.scale(sphere2, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -7));
	raytracer.rotate(plane, 'z', 45);
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	// Build the Octree
	objl::Loader ObjLoader;
	bool loadSuccess = ObjLoader.LoadFile("elephav623.obj");
	if(loadSuccess){
		cout << "load object success!\n";
		// Go through each loaded mesh and out its contents
		for (int i = 0; i < ObjLoader.LoadedMeshes.size(); i++)
		{
			// Copy one of the loaded meshes to be our current mesh
			objl::Mesh curMesh = ObjLoader.LoadedMeshes[i];
			// Build the Tree
			for (int j = 0; j < curMesh.Indices.size(); j += 3)
			{
				file << "T" << j / 3 << ": " << curMesh.Indices[j] << ", " << curMesh.Indices[j + 1] << ", " << curMesh.Indices[j + 2] << "\n";
				unsigned int In1 = curMesh.Indices[j];
				unsigned int In2 = curMesh.Indices[j+1];
				unsigned int In3 = curMesh.Indices[j+2];
				Vertex v1;
				Vertex v2;
				Vertex v3;
				v1.Position = Point3D(curMesh.Vertices[In1].Position.X, curMesh.Vertices[In1].Position.Y, curMesh.Vertices[In1].Position.Z);
				v2.Position = Point3D(curMesh.Vertices[In2].Position.X, curMesh.Vertices[In2].Position.Y, curMesh.Vertices[In2].Position.Z);
				v3.Position = Point3D(curMesh.Vertices[In3].Position.X, curMesh.Vertices[In3].Position.Y, curMesh.Vertices[In3].Position.Z);
				v1.mat = &gold;
				v2.mat = &gold;
				v3.mat = &gold;
				v1.normal = Vector3D(curMesh.Vertices[In1].Normal.X,curMesh.Vertices[In1].Normal.Y,curMesh.Vertices[In1].Normal.Z);
				v2.normal = Vector3D(curMesh.Vertices[In2].Normal.X,curMesh.Vertices[In2].Normal.Y,curMesh.Vertices[In2].Normal.Z);
				v3.normal = Vector3D(curMesh.Vertices[In3].Normal.X,curMesh.Vertices[In3].Normal.Y,curMesh.Vertices[In3].Normal.Z);
				Face f1(v1,v2,v3);
				Tree.insert(f1);
			}
		}
		file.close();
	}


	 cout << "finish build Octree\n";


	// Render the scene, feel free to make the image smaller for
	// testing purposes.
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");

	// Render it from a different point of view.
	Point3D eye2(4, 2, 5);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");

	return 0;
}

