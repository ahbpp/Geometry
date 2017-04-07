//
//  Vector.cpp
//  vector
//
//  Created by Alexey Karpov on 03.04.17.
//  Copyright © 2017 Alexey Karpov. All rights reserved.
//

#include "Vector.hpp"


Vector::Vector(const double x,const double y)
{
    this->x = x;
    this->y = y;
}
Vector::Vector(const Point &point)
{
    this->x = point.getx();
    this->y = point.gety();
}
Vector::Vector(const double x1, const double y1, const double x2, const double y2)
{
    this->x = x2 - x1;
    this->y = y2 - y1;
}
Vector::Vector (const Point &start, const Point &finish)
{
    this->x = finish.getx() - start.getx();
    this->y = finish.gety() - start.gety();
}
Vector::Vector(const Vector &that)
{
    this->x = that.x;
    this->y = that.y;
}
double Vector::getx() const
{
    return this->x;
}
double Vector::gety() const
{
    return this->y;
}
double Vector::abs() const
{
    return sqrt(this->x * this->x + this->y * this->y);
}
ostream& operator<<(ostream& res, const Vector &a)
{
    res << a.x << " " << a.y;
    return res;
}
istream& operator>>(istream& in, Vector &a)
{
    in >> a.x >> a.y;
    return in;
}
double Vector::operator * (const Vector &that) const
{
    double res = this->x * that.x + this->y * that.y;
    return res;
}
double Vector::vector_mult(const Vector &that) const
{
    double res = this->x * that.y - that.x * this->y;
    return res;
}
double Vector::triangle_area(const Vector &that) const
{
    return (vector_mult(that)) / 2;
}
Vector Vector::operator * (double c) const
{
    Vector res(this->x * c, this->y * c);
    return res;
}
Vector Vector::operator / (double c) const
{
    Vector res(this->x / c, this->y / c);
    return res;
}
Vector operator * (const double c, const Vector &that)
{
    return that * c;
}
Vector Vector::operator + (const Vector &that) const
{
    Vector res(this->x + that.x, this->y + that.y);
    return res;
}
Vector Vector::operator - (const Vector &that) const
{
    Vector res(this->x - that.x, this->y - that.y);
    return res;
}
Vector Vector::perpendikVector() const{
    return Vector(-y, x);
}

Point::Point()
{
    this->x = 0;
    this->y = 0;
}
Point::Point(const double x, const double y)
{
    this->x = x;
    this->y = y;
}
Point::Point(const Point &that)
{
    this->x = that.x;
    this->y = that.y;
}
double Point::getx() const
{
    return this->x;
}
double Point::gety() const
{
    return this->y;
}
Point &Point::operator=(const Point &that)
{
    this->x = that.x;
    this->y = that.y;
    return *this;
}
void Point::shift(const Vector & vec)
{
    this->x += vec.getx();
    this->y += vec.gety();
}
bool Point::includesPoint(const Point &that) const
{
    bool flag = false;
    if (this->x == that.x && this->y == that.y)
        flag = true;
    return flag;
}
bool Point::crossSegment(const Segment &that) const
{
    return that.includesPoint(*this);
}
double Point::dist(const Point& point) const {
    return sqrt((this->x - point.x)*(this->x - point.x) + (this->y - point.y)*(this->y - point.y));
}
ostream &operator<<(ostream &out, const Point &point)
{
    out << point.x << ' ' << point.y;
    return out;
}
istream &operator >> (istream &in, Point &point)
{
    in >> point.x >> point.y;
    return in;
}

Segment::Segment(const Point& a, const Point& b)
{
    this->a = a;
    this->b = b;
}
Segment::Segment(const Segment &that)
{
    this->a = that.a;
    this->b = that.b;
}
Point Segment::geta() const
{
    return this->a;
}
Point Segment::getb() const
{
    return this->b;
}
Line Segment::getLine() const {
    Line line(this->a, this->b);
    return line;
}
void Segment::shift(const Vector& vec) {
    this->a.shift(vec);
    this->b.shift(vec);
}

bool Segment::includesPoint(const Point& point) const {
    Ray ab(this->a, this->b);
    Ray ba(this->b, this->a);
    return ab.includesPoint(point) && ba.includesPoint(point);
}

bool Segment::crossSegment(const Segment& segment) const
{
    Line line1 = getLine();
    Line line2 = segment.getLine();
    double x1 = min(this->a.getx(), this->b.getx());
    double x2 = max(this->a.getx(), this->b.getx());
    double x3 = min(segment.geta().getx(), segment.getb().getx());
    double x4 = max(segment.geta().getx(), segment.getb().getx());
    
    double y1 = min(this->a.gety(), this->b.gety());
    double y2 = max(this->a.gety(), this->b.gety());
    double y3 = min(segment.geta().gety(), segment.getb().gety());
    double y4 = max(segment.geta().gety(), segment.getb().gety());
    
    bool projectionCross = ((x2 >= x3) && (x1 <= x4) && (y2 >= y3) && (y1 <= y4));
    
    if (line1.crossSegment(segment) && line2.crossSegment(*this) && projectionCross)
        return true;
    else
        return false;
}

double Segment::distanceToPoint(const Point& point) const {
    if (includesPoint(getLine().proekciya(point)))
        return getLine().distanceToPoint(point);
    else
        return min(point.dist(a), point.dist(b));
}

double Segment::distanceToSegment(const Segment& seg) const {
    return min( min(distanceToPoint(seg.a), distanceToPoint(seg.b)), min(seg.distanceToPoint(this->a), seg.distanceToPoint(this->b)) );
}

Line Ray::rayGetLine() const
{
    return Line(this->a, this->vec);
}
Line::Line(const Point& point1, const Point& point2) {
    
    if (point1.getx() == point2.getx()) {
        this->a = 1;
        this->b = 0;
        this->c = -point2.getx();
    }
    else {
        if (point1.gety() == point2.gety()) {
            this->a = 0;
            this->b = 1;
            this->c = -point2.gety();
        }
        else {
            this->a = 1 / (point1.getx() - point2.getx());
            this->b = 1 / (point2.gety() - point1.gety());
            this->c = point1.gety() / (point1.gety() - point2.gety()) - point1.getx() / (point1.getx() - point2.getx());
        }
    }
}
Line::Line(double a, double b, double c) : a(a), b(b), c(c) {}
Line::Line(const Point& pnt1, const Vector& vec) {
    Point pnt2 = pnt1;
    pnt2 = Point(pnt2.getx() + vec.getx(), pnt2.gety() + vec.gety());
    if (pnt1.getx() == pnt2.getx()) {
        this->a = 1;
        this->b = 0;
        this->c = -pnt2.getx();
    }
    else {
        if (pnt1.gety() == pnt2.gety()) {
            this->a = 0;
            this->b = 1;
            this->c = -pnt2.gety();
        }
        else {
            this->a = 1 / (pnt1.getx() - pnt2.getx());
            this->b = 1 / (pnt2.gety() - pnt1.gety());
            this->c = pnt1.gety() / (pnt1.gety() - pnt2.gety()) - pnt1.getx() / (pnt1.getx() - pnt2.getx());
        }
    }
}
double Line::geta() const {
    return this->a;
}
double Line::getb() const {
    return this->b;
}
double Line::getc() const {
    return this->c;
}

Point Line::lineGetPoint() const {
    if (a == 0)
        return Point(1, -this->c / this->b);
    
    if (b == 0)
        return Point(-this->c / this->a, 1);
    
    return Point(1, (-this->c - this->a) / this->b);
}

void Line::shift(const Vector& vec) {
    this->c = this->c - vec.getx() * this->a - vec.gety() * this->b;
    return;
}

bool Line::includesPoint(const Point& pnt) const
{
    double promez = this->a * pnt.getx() + this->b * pnt.gety();
    return ((promez <= -this->c + EPS) && (promez >= -this->c - EPS));
}

bool Line::crossSegment(const Segment& seg) const {
    Vector vec = lineGetVector();
    Point pnt = lineGetPoint();
    Vector vec1(pnt, seg.geta());
    Vector vec2(pnt, seg.getb());
    if (includesPoint(seg.geta()) || includesPoint(seg.getb()) || vec.vector_mult(vec1) * vec.vector_mult(vec2) <= 0)
        return true;
    else
        return false;
}


bool Line::ifCrossesLine(const Line& l1) const {
    if (a == 0) {
        if (l1.a == 0)
            return false;
        else
            return true;
    }
    if (b == 0) {
        if (l1.b == 0)
            return false;
        else
            return true;
    }
    if (a / l1.a == b / l1.b)
        return false;
    else
        return true;
}


Point Line::lineCrossLine(const Line& l1) const {
    double x, y;
    if (a == 0) {
        y = -c / b;
        x = (-l1.c - l1.b * y) / l1.a;
        Point pnt(x, y);
        return pnt;
    }
    if (l1.a == 0) {
        y = -l1.c / l1.b;
        x = (-c - b * y) / a;
        Point pnt(x, y);
        return pnt;
    }
    if (b == 0) {
        x = -c / a;
        y = (-l1.c - l1.a * x) / l1.b;
        Point pnt(x, y);
        return pnt;
    }
    if (l1.b == 0) {
        x = -l1.c / l1.a;
        y = (-c - a * x) / b;
        Point pnt(x, y);
        return pnt;
    }
    
    x = -(this->c * l1.b - l1.c * this->b) / (l1.b * this->a - this->b * l1.a);
    y = (-this->c - this->a * x) / this->b;
    Point pnt(x, y);
    return pnt;
}

double Line::parallDist(const Line& l1) const {
    
    if ( ((this->a * l1.b) == (this->b * l1.a)) && ((this->c * l1.b) == (this->b * l1.c)) && ((this->a * l1.c) == (this->c * l1.a)) )
        return 0;
    if (a == 0)
        return fabs(this->c / this->b - l1.c / l1.b);
    if (b == 0)
        return fabs(this->c / this->a - l1.c / l1.a);
    
    return fabs(this->c - l1.c * this->a / l1.a) / sqrt(this->a * this->a + this->b * this->b);
}

Vector Line::lineGetVector() const
{
    return Vector(-this->b, this->a);
}

Point Line::proekciya(const Point& pnt) const
{
    return lineCrossLine(Line(pnt, lineGetVector().perpendikVector()));
}

double Line::distanceToPoint(const Point& pnt) const {
    return pnt.dist(proekciya(pnt));
}

Ray::Ray(const Point& pnt, const Vector& v) : a(pnt), vec(v) {}

Ray::Ray(const Point& pnt1, const Point& pnt2) {
    this->a = pnt1;
    this->vec = Vector(pnt1, pnt2);
}
Point Ray::rayGeta() const
{
    return a;
}
Vector Ray::rayGetVector() const
{
    return vec;
}
void Ray::shift(const Vector& vec) {
    this->a.shift(vec);
}

bool Ray::includesPoint(const Point& pnt) const {
    Vector a_p (pnt, this->a);
    return this->vec.vector_mult(a_p) == 0 && this->vec * a_p <= 0;
}

bool Ray::crossSegment(const Segment& seg) const {
    if ((rayGetLine().crossSegment(seg) == true) && (includesPoint(rayGetLine().lineCrossLine(seg.getLine()))))
        return true;
    return false;
}

double Ray::distanceToPoint(const Point& pnt) const {
    if (includesPoint(rayGetLine().proekciya(pnt)))
        return rayGetLine().distanceToPoint(pnt);
    else
        return pnt.dist(this->a);
}


Polygon::Polygon(int n, Point *mas) {
    this->n = n;
    pointMas = new Point[this->n];
    for (int i = 0; i < this->n; i++)
    {
        pointMas[i] = mas[i];
    }
}
void Polygon::shift(const Vector& vec) {
    for (int i = 0; i < this->n; ++i) {
        pointMas[i].shift(vec);
    }
    return;
}

bool Polygon::includesPoint(const Point &pnt) const{
    
    if (Segment(pointMas[n - 1], pointMas[0]).includesPoint(pnt)) {
        return true;
    }
    for (int i = 0; i < n - 1; i++) {
        if (Segment(this->pointMas[i], this->pointMas[i + 1]).includesPoint(pnt)) {
            return true;
        }
    }
    double angleSum = 0;
    Vector one(pnt, this->pointMas[n - 1]);
    Vector two(pnt, this->pointMas[0]);
    int k = 1;
    if (one.vector_mult(two) < 0) {
        k = -1;
    }
    angleSum += (acos(one * two / one.abs() / two.abs())  * k);
    for (int i = 0; i < n - 1; i++) {
        one = two;
        two = Vector(this->pointMas[i + 1].getx() - pnt.getx(), pointMas[i + 1].gety() - pnt.gety());
        k = 1;
        if (one.vector_mult(two) < 0) {
            k = -1;
        }
        angleSum += (acos(one * two / one.abs() / two.abs())  * k);
    }
    return !(fabs(angleSum) < 1 + EPS);
    
}



bool Polygon::crossSegment(const Segment& seg) const {
    for (int i = 0; i < n - 1; ++i) {
        if (Segment(pointMas[i], pointMas[i + 1]).crossSegment(seg))
            return true;
    }
    if (Segment(this->pointMas[0], this->pointMas[n - 1]).crossSegment(seg))
        return true;
    if ( (includesPoint(seg.geta())) || (includesPoint(seg.getb())) )
        return true;
    
    return false;
}


bool Polygon::isProminent() const{
    Vector vec1(this->pointMas[n - 1], this->pointMas[0]);
    Vector vec2(this->pointMas[0], this->pointMas[1]);
    
    double mult = vec1.vector_mult(vec2);
    for (int i = 0; i < this->n - 1; i++) {
        if (i == this->n - 2) {
            vec1 = Vector(this->pointMas[n - 2], this->pointMas[n - 1]);
            vec2 = Vector(this->pointMas[n - 1], this->pointMas[0]);
        }
        else {
            vec1 = Vector(this->pointMas[i], this->pointMas[i + 1]);
            vec2 = Vector(this->pointMas[i + 1], this->pointMas[i + 2]);
        }
        if (mult == 0) {
            mult = vec1.vector_mult(vec2);
        }
        if (vec1.vector_mult(vec2) * mult < 0) {
            return false;
        }
    }
    return true;
}

double Polygon::square() const {
    double square = 0;
    square += Vector(this->pointMas[n - 1]).triangle_area(Vector(this->pointMas[n - 1], this->pointMas[0]));
    for (int i = 0; i < this->n - 1; ++i) {
        square += Vector(this->pointMas[i]).triangle_area(Vector(this->pointMas[i], this->pointMas[i + 1]));
    }
    return fabs(square);
}

Polygon::~Polygon()
{
    delete [] this->pointMas;
}

