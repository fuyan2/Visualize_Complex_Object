#ifndef Octree_H
#define Octree_H

#include <cstddef>
#include <vector>
#include "util.h"
#include <math.h>
using namespace std;
class Octree {
	// Physical position/size. This implicitly defines the bounding 
	// box of this node
	// The tree has up to eight children and can additionally store
	// a point, though in many applications only, the leaves will store data.
	Octree *children[8]; //! Pointers to child octants
	vector<Face> data;   //! Data point to be stored at a node
	int FaceNum;
	Point3D origin;         //! The physical center of this node
	Vector3D halfDimension;  //! Half the width/height/depth of this node
	public:
	Octree(const Octree& copy)
		: origin(copy.origin), halfDimension(copy.halfDimension), data(copy.data) {

		}
		
	Octree(const Point3D& origin, const Vector3D& halfDimension) 
		: origin(origin), halfDimension(halfDimension), data(NULL) {
			// Initially, there are no children
			FaceNum = 0;
			for(int i=0; i<8; ++i) 
				children[i] = NULL;
		}

	~Octree() {
		// Recursively destroy octants
		for(int i=0; i<8; ++i) 
			delete children[i];
	}

	// Determine which octant of the tree would contain 'point'
	int getOctantContainingPoint(const Point3D& point) const {
		int oct = 0;
		if(point[0] >= origin[0]) oct |= 4;
		if(point[1] >= origin[1]) oct |= 2;
		if(point[2] >= origin[2]) oct |= 1;
		return oct;
	}

	bool isLeafNode() const {
		// We are a leaf iff we have no children. Since we either have none, or 
		// all eight, it is sufficient to just check the first.
		return children[0] == NULL;
	}

	void insert(Face &face) {
		// If this node doesn't have a data point yet assigned 
		// and it is a leaf, then we're done!
		// cout <<point->Position<<" is going to be inserted\n";
		//cout<<"face is going to be insert\n";
		if(isLeafNode()) {
			if(FaceNum <= 10) {
				Face newFace = face;
				data.push_back(newFace);
				FaceNum++;
				//cout <<data->Position<<" is inserted\n";
				//cout<<"Face is inserted\n";
				return;
			} else {
				// Save this data point that was here for a later re-insert
				//cout <<data->Position<<" is removed for re-insert\n";
				vector<Face> oldFaces = data;
				//data.clear();
				//FaceNum = 0;
				// Split the current node and create new empty trees for each
				// child octant.
				for(int i=0; i<8; ++i) {
					// Compute new bounding box for this child
					Point3D newOrigin = origin;
					newOrigin[0] += halfDimension[0] * (i&4 ? 0.5 : -0.5);
					newOrigin[1] += halfDimension[1] * (i&2 ? 0.5 : -0.5);
					newOrigin[2] += halfDimension[2] * (i&1 ? 0.5 : -0.5);
					children[i] = new Octree(newOrigin, 0.5*halfDimension);
				}

				// Re-insert the old point, and insert this new point
				// (We wouldn't need to insert from the root, because we already
				// know it's guaranteed to be in this section of the tree)
				//cout<<"insert old Face\n";
				for(int i = 0; i < data.size(); i++){
					//cout << "insert old Face "<<i<<"\n";
					Vertex v1 = data[i][0];
					Vertex v2 = data[i][1];
					Vertex v3 = data[i][2];
					int node1 = getOctantContainingPoint(v1.Position);
					int node2 = getOctantContainingPoint(v2.Position);
					int node3 = getOctantContainingPoint(v3.Position);
					//children[node1]->insert(data[i]);
					//cout << "insert old Face "<<i<<" node \n";
					//Insert face to any Octant that contain at least one point of it
					if(node1 == node2 && node1 == node3){
						//cout << "insert old Face "<<i<<" all equal\n";
						children[node1]->insert(data[i]);
					}
					else if(node1 == node2){
						//cout << "insert old Face "<<i<<" 1 = 2 node1 is "<<node1 <<"node3 is "<<node3<<"\n";
						children[node1]->insert(data[i]);
						//children[node3]->insert(data[i]);
					}
					else if(node1 == node3){
						//cout << "insert old Face "<<i<<"1 = 3\n";
						children[node1]->insert(data[i]);
						children[node2]->insert(data[i]);
					}
					else if(node2 == node3){
						//cout << "insert old Face "<<i<<"2 = 3\n";
						children[node1]->insert(data[i]);
						children[node2]->insert(data[i]);
					}
					else{
						//cout << "insert old Face "<<i<<"all different\n";
						children[node1]->insert(data[i]);
						children[node2]->insert(data[i]);
						children[node3]->insert(data[i]);
					}
				}

				//cout<<"insert new Face\n";
				Vertex v1 = face[0];
				Vertex v2 = face[1];
				Vertex v3 = face[2];
				int node1 = getOctantContainingPoint(v1.Position);
				int node2 = getOctantContainingPoint(v2.Position);
				int node3 = getOctantContainingPoint(v3.Position);
				//Insert face to any Octant that contain at least one point of it
				if(node1 == node2 && node1 == node3){
					children[node1]->insert(face);
				}
				else if(node1 == node2){
					children[node1]->insert(face);
					children[node3]->insert(face);
				}
				else if(node1 == node3){
					children[node1]->insert(face);
					children[node2]->insert(face);
				}
				else if(node2 == node3){
					children[node1]->insert(face);
					children[node2]->insert(face);
				}
				else{
					children[node1]->insert(face);
					children[node2]->insert(face);
					children[node3]->insert(face);
				}
			}
		} else {
			// We are at an interior node. Insert recursively into the 
			// appropriate child octant
			//cout<<"insert new point\n";
			Vertex v1 = face[0];
			Vertex v2 = face[1];
			Vertex v3 = face[2];
			int node1 = getOctantContainingPoint(v1.Position);
			int node2 = getOctantContainingPoint(v2.Position);
			int node3 = getOctantContainingPoint(v3.Position);
			//Insert face to any Octant that contain at least one point of it
			if(node1 == node2 && node1 == node3){
				children[node1]->insert(face);
			}
			else if(node1 == node2){
				children[node1]->insert(face);
				children[node3]->insert(face);
			}
			else if(node1 == node3){
				children[node1]->insert(face);
				children[node2]->insert(face);
			}
			else if(node2 == node3){
				children[node1]->insert(face);
				children[node2]->insert(face);
			}
			else{
				children[node1]->insert(face);
				children[node2]->insert(face);
				children[node3]->insert(face);
			}
		}
	}

	bool rayIntersectNode(Ray3D &r) {
		double Tmin, Tmax, tymin, tymax, tzmin, tzmax; 
 		Vector3D invdir = Vector3D(1.0/r.dir[0], 1.0/r.dir[1], 1.0/r.dir[2]);
 		Point3D bounds[2];
 		bounds[0] = origin - halfDimension;
 		bounds[1] = origin + halfDimension;
 		int sign[3];
 		sign[0] = (invdir[0] < 0); 
        sign[1] = (invdir[1] < 0); 
        sign[2] = (invdir[2] < 0);
	    Tmin = (bounds[sign[0]][0] - r.origin[0]) * invdir[0]; 
	    Tmax = (bounds[1-sign[0]][0] - r.origin[0]) * invdir[0]; 
	    tymin = (bounds[sign[1]][1] - r.origin[1]) * invdir[1]; 
	    tymax = (bounds[1-sign[1]][1] - r.origin[1]) * invdir[1]; 	 	
	    if ((Tmin > tymax) || (tymin > Tmax)){
	        return false; 
	    }
	    if (tymin > Tmin) {
	        Tmin = tymin; 
	    }
	    if (tymax < Tmax) {
	        Tmax = tymax; 
	    }
	 
	    tzmin = (bounds[sign[2]][2] - r.origin[2]) * invdir[2]; 
	    tzmax = (bounds[1-sign[2]][2] - r.origin[2]) * invdir[2]; 	 
	    if ((Tmin > tzmax) || (tzmin > Tmax)) {
	        return false; 
	    }
	    if (tzmin > Tmin) {
	        Tmin = tzmin; 
	    }
	    if (tzmax < Tmax) {
	        Tmax = tzmax; 
	    }
	    return true; 
	}

bool rayTriangleIntersect(const Ray3D &ray, const Face &face, double &t, Point3D &IntersectPoint) 
{ 
    // compute plane's normal
    Point3D v0 = face[0].Position;
    Point3D v1 = face[1].Position;
    Point3D v2 = face[2].Position;
    Point3D orig = ray.origin;
    Vector3D dir = ray.dir;
    Vector3D v0v1 = v1 - v0; 
    Vector3D v0v2 = v2 - v0; 
    Vector3D N = v0v1.cross(v0v2); // N  
    // check if ray and plane are parallel ?
    double NdotRayDirection = N.dot(dir); 
    if (fabs(NdotRayDirection) < 0.01) // almost 0 
    	// they are parallel so they don't intersect ! 
        return false; 
 
    // compute d parameter
    Vector3D v0V = v0 - Point3D(0,0,0);
    double d = N.dot(v0V); 
 
    // compute t
    Vector3D origV = orig - Point3D(0,0,0);
    t = (N.dot(origV) + d) / NdotRayDirection; 
    // check if the triangle is in behind the ray
    if (t < 0) return false; // the triangle is behind 
 
    // compute the intersection point using equation 1
    Point3D P = orig + t * dir; 
    Vector3D C; // vector perpendicular to triangle's plane 
    // edge 0
    Vector3D edge0 = v1 - v0; 
    Vector3D vp0 = P - v0; 
    C = edge0.cross(vp0); 
    if (N.dot(C) < 0) {
    	return false; // P is on the right side 
    }
 
    // edge 1
    Vector3D edge1 = v2 - v1; 
    Vector3D vp1 = P - v1; 
    C = edge1.cross(vp1); 
    if (N.dot(C) < 0) {
    	return false; // P is on the right side 
    }
 
    // edge 2
    Vector3D edge2 = v0 - v2; 
    Vector3D vp2 = P - v2; 
    C = edge2.cross(vp2); 
    if (N.dot(C) < 0) {
    	return false; // P is on the right side; 
    }
 
 	IntersectPoint = P;
 	// this ray hits the triangle 
    return true; 
} 

void setIntersect(Ray3D &ray){
	if(isLeafNode()) {
		if(FaceNum != 0) {
			if(rayIntersectNode(ray)){
				double Tmin = 10000000000;
				double Ttemp = 0;
				Point3D IntersectPoint(0,0,0);
				//cout<<"find intersect octant\n";
				for(int i = 0; i < data.size(); i++){
					//cout <<"enter for loop\n";
					if(rayTriangleIntersect(ray, data[i], Ttemp, IntersectPoint)){
						//cout<<"find intersect triangle\n";
						if(Ttemp < Tmin){
							Tmin = Ttemp;
							if(ray.intersection.none || (Tmin < ray.intersection.t_value && Tmin > 0)){
					 			ray.intersection.point = IntersectPoint;
					 			Vector3D norm = data[i][0].normal + data[i][1].normal + data[i][2].normal;
					 			norm.normalize();
								ray.intersection.normal = norm;
								ray.intersection.t_value = Tmin;
								ray.intersection.mat = data[i][0].mat;
								ray.intersection.none = false;
							//	cout << "ray intersect "<<IntersectPoint<<"\n";
							}
						}
					}
				}
			}			
		}
		return;
	}else{
		if (rayIntersectNode(ray)) {
			for(int i=0; i<8; ++i) {
				children[i]->setIntersect(ray);
			}
		}
		return;
	}
	return;
}

};


#endif
