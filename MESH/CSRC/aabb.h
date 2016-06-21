/* This file is part of SEM                                                */
/*                                                                         */
/* Copyright CEA, ECP, IPGP                                                */
/*                                                                         */
// mesh.h : Gestion maillage format SEM
#ifndef AABB_H
#define AABB_H

#include <limits>
/**  Axis Aligned Bounding Box
 */

struct Vec3
{
    Vec3() {}
    Vec3(double x, double y, double z) { v[0] = x; v[1] = y; v[2] = z; }
    double& operator[](int i) { return v[i]; }
    const double& operator[](int i) const { return v[i]; }
    double v[3];
};

class AABB {
public:
    AABB() {init_bounds();}
    AABB(const Vec3 &vmin,const Vec3 &vmax): min(vmin),max(vmax) {}

    AABB(AABB const& box)
        {
            min=box.min;max=box.max;
        }

    ~AABB() {};

    void init_bounds() {
        double maxval = std::numeric_limits<double>::max();
        double minval = -maxval;
        max = Vec3(minval, minval, minval);
        min = Vec3(maxval, maxval, maxval);
    }

    void update_bounds(const Vec3& pos) {
        if (pos[0] > max[0]) max[0] = pos[0];
        if (pos[1] > max[1]) max[1] = pos[1];
        if (pos[2] > max[2]) max[2] = pos[2];

        if (pos[0] < min[0]) min[0] = pos[0];
        if (pos[1] < min[1]) min[1] = pos[1];
        if (pos[2] < min[2]) min[2] = pos[2];
    }


    void set_bounds(const Vec3 &vmin,const Vec3 &vmax)
        {
            max[0]=vmax[0];max[1]=vmax[1];max[2]=vmax[2];
            min[0]=vmin[0];min[1]=vmin[1];min[2]=vmin[2];
        }

    double xmax() const { return max[0]; }
    double ymax() const { return max[1]; }
    double zmax() const { return max[2]; }
    double xmin() const { return min[0]; }
    double ymin() const { return min[1]; }
    double zmin() const { return min[2]; }
    double xl() const { return max[0] - min[0]; }
    double yl() const { return max[1] - min[1]; }
    double zl() const { return max[2] - min[2]; }


    bool contact(AABB &b);

    bool test_inside(const Vec3& p) const
        {
            for(int i=0;i<3;++i)
            {
                if (p[i]<min[i]) return false;
                if (p[i]>max[i]) return false;
            }
            return true;
        }
    Vec3 max;
    Vec3 min;
};

#endif

/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
