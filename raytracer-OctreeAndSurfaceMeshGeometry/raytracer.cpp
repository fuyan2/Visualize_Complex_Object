/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h,
		and the main function which specifies the
		scene to be rendered.

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent,
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;

	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}

	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;

	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation;
			angle = -angle;
		}
		else {
			node->invtrans = rotation*node->invtrans;
		}
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;

	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation;
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans;
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;

	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale;
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans;
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view,
		Vector3D up ) {
	Matrix4x4 mat;
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat;
}


void Raytracer::computeTransforms( SceneDagNode* node )
{
    SceneDagNode *childPtr;
    if (node->parent != NULL )
    {
        node->modelToWorld = node->parent->modelToWorld*node->trans;
        node->worldToModel = node->invtrans*node->parent->worldToModel;
    }
    else
    {
        node->modelToWorld = node->trans;
        node->worldToModel = node->invtrans;
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        computeTransforms(childPtr);
        childPtr = childPtr->next;
    }



}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
    SceneDagNode *childPtr;

    // Applies transformation of the current node to the global
    // transformation matrices.

    if (node->obj) {
        // Perform intersection.
        if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
            ray.intersection.mat = node->mat;
        }
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        traverseScene(childPtr, ray);
        childPtr = childPtr->next;
    }

}

void Raytracer::computeShading( Ray3D& ray ) {
    LightListNode* curLight = _lightSource;
    for (;;) {
        if (curLight == NULL) break;
        // Each lightSource provides its own shading function.

        // Implement shadows here if needed.
				Point3D light = curLight->light->get_position();
        Colour col_ambient = curLight->light->get_ambient();
        Vector3D lightdir = light - ray.intersection.point;
        double t = lightdir.length();
        lightdir.normalize();
				//0.00001 is needed to avoid black dots appearing on the plane due to
				//numerical errors which result in the shadowRay origin being below the
				//plane
        Ray3D shadow(ray.intersection.point + 0.00001*lightdir, lightdir);

        bool octree_enable = true;
	    if(octree_enable){
	    //	std::cout<<"traverseScene with OctTree\n";
	    	Tree.setIntersect(shadow);
	    }
	    else
        	traverseScene(_root, shadow);

        if (!shadow.intersection.none && t >= shadow.intersection.t_value) {
        	ray.col = ray.intersection.mat->ambient * col_ambient;
        	ray.col.clamp();
        }
				else{
					curLight->light->shade(ray);
				}
        curLight = curLight->next;
    }
}

void Raytracer::initPixelBuffer() {
    int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
    _rbuffer = new unsigned char[numbytes];
    std::fill_n(_rbuffer, numbytes,0);
    _gbuffer = new unsigned char[numbytes];
    std::fill_n(_gbuffer, numbytes,0);
    _bbuffer = new unsigned char[numbytes];
    std::fill_n(_bbuffer, numbytes,0);
}

void Raytracer::flushPixelBuffer( std::string file_name ) {
    bmp_write( file_name.c_str(), _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
    delete _rbuffer;
    delete _gbuffer;
    delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray , int depth) {
    Colour col(0.0, 0.0, 0.0);
    bool octree_enable = true;
    if(octree_enable){
    	//std::cout<<"traverseScene with OctTree\n";
    	Tree.setIntersect(ray);
    }
    else
		traverseScene(_root, ray);
    // Don't bother shading if the ray didn't hit
    // anything.
    if (!ray.intersection.none) {
			computeShading(ray);
			if (depth == 0)
				col = ray.col;
			if (ray.intersection.mat->reflection_coef && depth <=2 ){
				Vector3D reflectedRay =  2 * ray.intersection.normal.dot(-ray.dir) * ray.intersection.normal +ray.dir;
				reflectedRay.normalize();
				Ray3D Ray(ray.intersection.point,reflectedRay);
				col = col + (ray.intersection.mat->reflection_coef * shadeRay(Ray,depth+1));
			}
			else
				return ray.col;
    }

    // You'll want to call shadeRay recursively (with a different ray,
    // of course) here to implement reflection/refraction effects.
		col.clamp();
    return col;
}

void Raytracer::render( int width, int height, Point3D eye, Vector3D view,
        Vector3D up, double fov, std::string fileName ) {
    computeTransforms(_root);
    Matrix4x4 viewToWorld;
    _scrWidth = width;
    _scrHeight = height;
    double factor = (double(height)/2)/tan(fov*M_PI/360.0);

    initPixelBuffer();
    viewToWorld = initInvViewMatrix(eye, view, up);

    // Construct a ray for each pixel.
    for (int i = 0; i < _scrHeight; i++) {
        for (int j = 0; j < _scrWidth; j++) {
            // Sets up ray origin and direction in view space,
            // image plane is at z = -1.
            Point3D origin(0, 0, 0);
						Point3D imagePlane;
						imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
						imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
						imagePlane[2] = -1;

						// TODO: Convert ray to world space and call
						// shadeRay(ray) to generate pixel colour.
						Ray3D ray;

						ray.origin = viewToWorld * origin;
						ray.dir = viewToWorld * (imagePlane - origin);
						ray.dir.normalize();

						Colour col = shadeRay(ray,0);

						_rbuffer[i*width+j] = int(col[0]*255);
						_gbuffer[i*width+j] = int(col[1]*255);
						_bbuffer[i*width+j] = int(col[2]*255);
		}
	}

	flushPixelBuffer(fileName);
}
