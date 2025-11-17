#ifndef GEOMETERY_H
#define GEOMETERY_H

struct point2d {
    double x;
    double y;
};

struct point3d {
    double x;
    double y;
    double z;
    unsigned char r = 255;
    unsigned char g = 255;
    unsigned char b = 255;
};

struct Sphere {
    double x0, y0, z0; // ÇòÐÄ×ø±ê
    double radius;      // °ë¾¶
};

#endif /* ifndef GEOMETERY_H */
