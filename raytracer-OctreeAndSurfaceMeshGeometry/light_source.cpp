/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <algorithm>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray
	// is available.  So be sure that traverseScene() is called on the ray
	// before this function.
  Vector3D normal = ray.intersection.normal;


  // light vector
  Vector3D lightVec = _pos - ray.intersection.point;
  lightVec.normalize();

  //light Reflection vector
  //http://paulbourke.net/geometry/reflected/
  //calculation from above link
  Vector3D lightReflect =  2 * normal.dot(lightVec) * normal - lightVec;
  lightReflect.normalize();

  // view direction vector
  Vector3D viewVector = - ray.dir;
  viewVector.normalize();


  Colour Ia = (ray.intersection.mat->ambient) * _col_ambient;

  Colour Id = fmax(0, lightVec.dot(normal)) * (ray.intersection.mat->diffuse) * _col_diffuse;

  Colour Is = fmax(0, pow((viewVector.dot(lightReflect)), (ray.intersection.mat->specular_exp)))
                      * (ray.intersection.mat->specular) * _col_specular;

  Colour colour = Ia + Id + Is;

  ray.col = colour;

  ray.col.clamp();
}
