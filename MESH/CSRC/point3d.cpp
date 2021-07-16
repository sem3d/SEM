
#include <cassert>
#include "point3d.h"


void PointsArray_3D::create_cube(double W, int N0)
{
    int NP = (N0+1)*(N0+1)*(N0+1);
    int count=0;
    pts.resize(NP);
    for(int k=0;k<=N0;++k) {
        double z=W*double(2*k-N0)/double(N0);
	for(int j=0;j<=N0;++j) {
	    double y=W*double(2*j-N0)/double(N0);
	    for(int i=0;i<=N0;++i) {
		double x=W*double(2*i-N0)/double(N0);
		Point3D_idx pti(x,y,z);
		Point3D& pt = pts[count];
		pts_map[pti] = count;
		pt.x = x;
		pt.y = y;
		pt.z = z;
		pt.level = 0;
		pt.app_data = 0;
		++count;
	    }
	}
    }
    assert(count==NP);
}

int PointsArray_3D::create_point(double x, double y, double z, int level)
{
    int id;
    Point3D_idx pti(x,y,z);
    pt_map_t::iterator it = pts_map.lower_bound(pti);

    if (it!=pts_map.end() && (pti==it->first)) {
	id = it->second;
	if (level>pts[id].level)
	    pts[id].level = level;
    } else {
	// New point
	id = pts.size();
	Point3D pt;
	pt.x = x;
	pt.y = y;
	pt.z = z;
	pt.level = level;
	pts.push_back(pt);
	pts_map.insert(it, pt_map_t::value_type(pti, id));
    }
    return id;
}

int PointsArray_3D::subpoint(double X[8], double Y[8], double Z[8], int level, double a, double b, double c)
{
    double x, y, z;

    x = (1-a)*(1-b)*(1-c)*X[0] + a*(1-b)*(1-c)*X[1] + a*b*(1-c)*X[2] + (1-a)*b*(1-c)*X[3] +
	(1-a)*(1-b)*c*X[4] + a*(1-b)*c*X[5] + a*b*c*X[6] + (1-a)*b*c*X[7];

    y = (1-a)*(1-b)*(1-c)*Y[0] + a*(1-b)*(1-c)*Y[1] + a*b*(1-c)*Y[2] + (1-a)*b*(1-c)*Y[3] +
	(1-a)*(1-b)*c*Y[4] + a*(1-b)*c*Y[5] + a*b*c*Y[6] + (1-a)*b*c*Y[7];

    z = (1-a)*(1-b)*(1-c)*Z[0] + a*(1-b)*(1-c)*Z[1] + a*b*(1-c)*Z[2] + (1-a)*b*(1-c)*Z[3] +
	(1-a)*(1-b)*c*Z[4] + a*(1-b)*c*Z[5] + a*b*c*Z[6] + (1-a)*b*c*Z[7];

    return create_point(x,y,z, level);
}

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
