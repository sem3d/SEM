/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// earthmesh.h : Generates a sphere meshed with hexaedrons and refined layers
#ifndef _POINT3D_H_
#define _POINT3D_H_

#include <map>
#include <vector>

const double PREC=1e-8;
const double DPREC=1./PREC;

struct Point3D {
    double x, y, z;
    int  level;
    int   app_data;
};

struct Point3D_idx
{
    Point3D_idx():x(0),y(0),z(0) {}
    Point3D_idx(const Point3D_idx& p):x(p.x), y(p.y), z(p.z) {}
    Point3D_idx(double _x, double _y, double _z):x(_x*DPREC), y(_y*DPREC), z(_z*DPREC) {}

    long x,y,z;

    bool operator==(const Point3D_idx& p) const {
	return x==p.x && y==p.y && z==p.z; }
    bool operator<(const Point3D_idx& p) const {
	return x<p.x || (x==p.x && y<p.y) || (x==p.x && y==p.y && z<p.z);;
    }
};


struct PointsArray_3D {
    typedef std::map<Point3D_idx,int> pt_map_t;
    pt_map_t pts_map; /// maps coordinates to index
    std::vector<Point3D> pts;           /// Point at index

    int create_point(double x, double y, double z, int level);
    int subpoint(double x[8], double y[8], double z[8], int level, double a, double b, double c);
    void create_cube(double W, int N0);

    Point3D* data() { return &pts[0]; }
    int level(int p) { return pts[p].level; }
    int npoints() const { return pts.size(); }
};

#endif
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
