//
//  Vector.hpp
//  vector
//
//  Created by Alexey Karpov on 03.04.17.
//  Copyright © 2017 Alexey Karpov. All rights reserved.
//

#ifndef Vector_hpp
#define Vector_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

class Point;
class Segment;
class Line;
class Ray;
class Polygon;

const double EPS = 0.00000001;

class Vector {
private:
    double x;
    double y;
public:
    Vector(const double x = 0.0, const double y = 0.0);
    Vector(const Point &pnt);
    Vector(const double x1, const double y1, const double x2, const double y2);
    Vector (const Point &start, const Point &finish);
    Vector(const Vector &that);
    double getx() const;
    double gety() const;
    double abs() const;
    double operator * (const Vector &that) const;
    double vector_mult(const Vector &that) const;
    double triangle_area(const Vector &that) const;
    Vector operator * (double c) const;
    Vector operator / (double c) const;
    Vector operator + (const Vector &that) const;
    Vector operator - (const Vector &that) const;
    friend Vector operator * (const double c, const Vector &that);
    friend ostream& operator<<(ostream& res, const Vector &a);
    friend istream& operator>>(istream& in, Vector &a);
    Vector perpendikVector() const;
};



class Figure {
public:
    virtual void shift(const Vector &) = 0;
    virtual bool includesPoint(const Point &) const = 0;
    virtual bool crossSegment(const Segment &) const = 0;
};

class Point : public Figure
{
private:
    double x;
    double y;
public:
    Point(const double x, const double y);
    Point(const Point &that);
    Point();
    double dist(const Point& pnt) const;
    double getx() const;
    double gety() const;
    Point &operator=(const Point &that);
    void shift (const Vector &);
    bool includesPoint(const Point &) const;
    bool crossSegment(const Segment &) const;
    friend ostream &operator<<(ostream &out, const Point &point);
    friend istream &operator >> (istream &in, Point &point);
};

class Segment : public Figure {
    Point a;
    Point b;
public:
    Segment(const Point &a, const Point &b);
    Segment(const Segment &that);
    void shift(const Vector &);
    bool includesPoint(const Point &) const;
    bool crossSegment(const Segment &) const;
    Point geta() const;
    Point getb() const;
    Line getLine() const;
    double distanceToPoint(const Point &) const;
    double distanceToSegment(const Segment &) const;
};

class Line : public Figure {
    double a, b, c;
public:
    Line(double a, double b, double c);
    Line(const Point& pnt1, const Point & pnt2);
    double geta() const;
    double getb() const;
    double getc() const;
    Line(const Point& pnt1, const Vector & vec);
    
    Vector lineGetVector() const;
    bool ifCrossesLine(const Line &) const;
    double parallDist(const Line &) const;
    Point lineCrossLine(const Line &) const;
    Point lineGetPoint() const;
    Point proekciya(const Point &) const;
    double distanceToPoint(const Point &) const;
    
    void shift(const Vector&);
    bool includesPoint(const Point &) const;
    bool crossSegment(const Segment &) const;
};

class Ray : public Figure {
    Point a;
    Vector vec;
public:
    Ray(const Point &, const Vector &);
    Ray(const Point &, const Point &);
    
    Point rayGeta() const;
    Vector rayGetVector() const;
    void shift(const Vector &);
    bool includesPoint(const Point &) const;
    bool crossSegment(const Segment &) const;
    Line rayGetLine() const;
    double distanceToPoint(const Point &) const;
    
};

class Polygon : public Figure {
    int n;
    Point *pointMas;
public:
    Polygon(int n, Point *mas);
    void shift(const Vector &);
    bool includesPoint(const Point &) const;
    bool crossSegment(const Segment &) const;
    bool isProminent() const;
    double square() const;
    ~Polygon();
};

#endif /* Vector_hpp */
