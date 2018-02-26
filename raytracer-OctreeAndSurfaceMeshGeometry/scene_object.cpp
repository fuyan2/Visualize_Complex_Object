/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {

			// TODO: implement intersection code for UnitSquare, which is
			// defined on the xy-plane, with vertices (0.5, 0.5, 0),
			// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
			// (0, 0, 1).
			//
			// Your goal here is to fill ray.intersection with correct values
			// should an intersection occur.  This includes intersection.point,
			// intersection.normal, intersection.none, intersection.t_value.
			//
			// HINT: Remember to first transform the ray into object space
			// to simplify the intersection test.

			Ray3D rayModel;
			rayModel.origin = worldToModel * ray.origin;
			rayModel.dir = worldToModel * ray.dir;
			//Using the equation, point = ray.origin + t_value * ray.dir
			//and since the Square lies on the x-y Plane =>

			double t_value = (0 - rayModel.origin[2]) / rayModel.dir[2];

			if (t_value > 0){
				if (ray.intersection.none || t_value < ray.intersection.t_value){
				 Point3D point = rayModel.origin + t_value * rayModel.dir;
				 if (point[0]< 0.5 && point[0]>-0.5 && point[1]< 0.5 && point[1] > -0.5){
					 ray.intersection.point = modelToWorld * point;
					 ray.intersection.normal = transNorm(worldToModel, Vector3D(0.0, 0.0, 1.0));
					 ray.intersection.normal.normalize();
					 ray.intersection.t_value = t_value;
					 ray.intersection.none = false;
					 return true;
				 }
			 	}
			 }

	return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
			// Ray3D rayModel;
			// rayModel.origin = worldToModel * ray.origin;
			// rayModel.dir = worldToModel * ray.dir;
	// TODO: implement intersection code for UnitSphere, which is centred
	// on the origin.
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point,
	// intersection.normal, intersection.none, intersection.t_value.
	//
	// HINT: Remember to first transform the ray into object space
	// to simplify the intersection test.

	Ray3D rayModel;
	rayModel.origin = worldToModel * ray.origin;
	rayModel.dir = worldToModel * ray.dir;

	// equation of a sphere: (p-c) dot (p-c) = r^2
	//equation of ray: r = p0 + td
	//by substituting equation of sphere into ray
	//we get a quadratic equation
	//https://www.youtube.com/watch?v=bQKy3N4TshU&t=374s
	//This link derives the calculations for a,b and c
	Point3D center(0,0,0); // 0 being the center of the sphere

	double a = rayModel.dir.dot(rayModel.dir);
	double b = 2 * rayModel.dir.dot(((rayModel.origin - Point3D(0,0,0))));
	double c = ((rayModel.origin - Point3D(0,0,0)).dot((rayModel.origin - Point3D(0,0,0)))) - 1; // radius being 1

	double discriminant = b * b - 4 * a * c;

	double t_value;

	if (discriminant < 0)
		return false; // no intersection
	else if (discriminant == 0)
		t_value = -b/(2*a); // only one intersections
	else{
		t_value = ( -b - sqrt(discriminant) ) / (2 * a); // we don't care about the other root since it will be greater
	}

  if (ray.intersection.none || t_value < ray.intersection.t_value){
		if(t_value > 0){
			Point3D point = rayModel.origin + t_value * rayModel.dir;
			ray.intersection.point = modelToWorld * point;
			ray.intersection.normal = transNorm(worldToModel,point - center);
			ray.intersection.normal.normalize();
			ray.intersection.t_value = t_value;
			ray.intersection.none = false;
			return true;
		}
	}

	return false;
}

void UnitSquare::setPointArray(const Matrix4x4& modelToWorld, Material* mat){
	PointArray = new Vertex[10000];
	double phi = 0.1;
	double theta = 0.1;
	for(int i = 0; i < 10; i++){
		for(int j = 0; j < 10; j++){
			double x = i * phi;
			double y = j * theta;
			double z = 0.0;
		 	PointArray[i*10+j].Position = modelToWorld * Point3D(x, y, z);
		 	PointArray[i*10+j].normal = transNorm(modelToWorld, Vector3D(.0, .0, 1.0));
		 	PointArray[i*10+j].normal.normalize();
		 	PointArray[i*10+j].mat = mat;
		}
	}
}

void UnitSphere::setPointArray(const Matrix4x4& modelToWorld, Material* mat){
	PointArray = new Vertex[10000];
	double phi = 2 * M_PI / 100;
	double theta = 2 * M_PI / 100;
	for(int i = 0; i < 100; i++){
		for(int j = 0; j < 100; j++){
			double x = cos(phi * i) * sin(theta * j);
			double y = sin(phi * i) * sin(theta * j);
			double z = cos(theta * j);
		 	PointArray[i*100+j].Position = modelToWorld * Point3D(x, y, z);
		 	PointArray[i*100+j].normal = transNorm(modelToWorld, Vector3D(x, y, z));
		 	PointArray[i*100+j].normal.normalize();
		 	PointArray[i*100+j].mat = mat;
		}
	}
}
