		//Below only apply when the view is along z direction
		Point3D center1 = origin;
		center1[2] = origin[2] - halfDimension[2];
		Point3D center2 = origin;
		center1[2] = origin[2] + halfDimension[2];
		double t_value1 = (center1[2] - ray.origin[2]) / ray.dir[2];	
		double t_value2 = (center2[2] - ray.origin[2]) / ray.dir[2];
		if (t_value1 > 0 || t_value2 > 0){
			Point3D point1 = ray.origin + t_value1 * ray.dir;
			Vector3D distance1 = point1 - center1;
			Point3D point2 = ray.origin + t_value2 * ray.dir;
			Vector3D distance2 = point2 - center2;
			if (distance1[0] <= halfDimension[0] && distance1[0]>-halfDimension[0] && distance1[1]< halfDimension[1] && distance1[1] > -halfDimension[1]){
				return true;
		 	}else if (distance2[0] <= halfDimension[0] && distance2[0]>-halfDimension[0] && distance2[1]< halfDimension[1] && distance2[1] > -halfDimension[1]){
				return true;
		 	}
		 	return false;
	 	}
		else
			return false;