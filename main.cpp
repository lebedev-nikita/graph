#include <iostream>
#include <cstdint>

#include <string>
#include <vector>
#include <unordered_map>

#include <cmath>

using namespace std;
#include "Gauss.h"
#include "Bitmap.h"

const uint32_t RED   = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE  = 0x00FF0000;

/*
    TODO: LATER:
      * Окончательно настроить сцены
      * Написать в README, какими источниками пользовался
*/
#define IMG_WIDTH 1024
#define IMG_HEIGHT 1024
#define FOV M_PI/3
#define BACKGROUND_COLOR Vec3d(1)
#define EPS 0.000000000001
#define REC_DEPTH 8


// Векторы:

class Vec3d {
public:
  double x, y, z;

  Vec3d(): x(0), y(0), z(0)
  {}
  Vec3d(double xx, double yy, double zz): x(xx), y(yy), z(zz)
  {}
  Vec3d(double xx): x(xx), y(xx), z(xx)
  {}
  inline Vec3d operator * (const double n) const
  { return Vec3d(this->x * n, this->y * n, this->z * n); }
  inline Vec3d operator * (const Vec3d &v) const
  { return Vec3d(this->x * v.x, this->y * v.y, this->z * v.z); }
  friend inline Vec3d operator * (const double &r, const Vec3d &v)
  { return Vec3d(v.x * r, v.y * r, v.z * r); }
  inline Vec3d operator - (const Vec3d &v) const
  { return Vec3d(this->x - v.x, this->y - v.y, this->z - v.z); }
  inline Vec3d operator - () const
  { return Vec3d(-x, -y, -z); }
  inline Vec3d operator + (const Vec3d &v) const
  { return Vec3d(this->x + v.x, this->y + v.y, this->z + v.z); }
  inline Vec3d operator += (const Vec3d &v)
  {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
  }
  friend inline double abs(const Vec3d &v)
  { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
};


// Вспомогательные функции:

inline double dotProduct(const Vec3d &a, const Vec3d &b)
{ return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }

  /* (В левом ортонормированном базисе) */
inline Vec3d vectProduct(const Vec3d &a, const Vec3d &b)
{ return Vec3d(a.z*b.y-a.y*b.z, a.x*b.z-a.z*b.x, a.y*b.x-a.x*b.y); }

inline Vec3d normalize(const Vec3d &vec)
{ return vec * (1 / sqrt(dotProduct(vec,vec))); }

inline double clamp(double num, double low = 0, double hi = 1)
{ return (num < low)?low:(num>hi)?hi:num; }

inline double cos(const Vec3d &a, const Vec3d &b)
{ return clamp(dotProduct(normalize(a), normalize(b)), -1, 1);}

inline Vec3d reflectRay(const Vec3d &I, const Vec3d &N)
{ return 2*N*dotProduct(I, N) - I; }

Vec3d refractRay(const Vec3d &I, const Vec3d &N, const double &ior)
{
    // считаем, что abs(I) == 1
    double cosi = dotProduct(I, N);
    double etai = 1, etat = ior;
    Vec3d n = N;
    if (cosi < 0) { cosi = -cosi; } else { swap(etai, etat); n= -N; }
    double eta = etai / etat;
    double k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrt(k)) * n;
}

bool solveQuadratic(double a, double b, double c, double &tMin, double &tMax)
{
  double discr = b*b - 4*a*c;
  if (discr < 0) {
    return false;
  }
  else if (discr == 0) {
    tMin = tMax = -b/(2*a);
    return true;
  }
  else {
    double t1, t2;
    double sqrtDiscr = sqrt(discr);
    t1 = (-b + sqrtDiscr) / (2*a);
    t2 = (-b - sqrtDiscr) / (2*a);

    if (t1 > t2)
      swap(t1, t2);

    tMin = t1;
    tMax = t2;
    return true;
  }
}

// Свет:

enum LightType { AMBIENT, POINT, DIRECTIONAL};

class Light {
public:
  Vec3d color;
  LightType type;

  Vec3d source = Vec3d(0,0,0);
  Vec3d direction = Vec3d(0,0,0);

  Light(const Vec3d &clr, LightType type) : color(clr), type(type) {}
};


// Объекты:

enum MaterialType { DIFFUSE, DIFFUSE_AND_GLOSSY, MIRROR, ALL_IN_ONE, GLASS };

class Object {
public:
  Vec3d color;
  MaterialType material;
  double kDiffuse;
  double kRefl;
  int specExp;
  double ior;

  Object(MaterialType m = DIFFUSE_AND_GLOSSY) : material(m), color(1)
  {
    if (m == DIFFUSE) {
      kDiffuse = 1;
      specExp = -1;
      kRefl = 0;
      ior = -1;
    }
    else if (m == DIFFUSE_AND_GLOSSY) {
      kDiffuse = 0.5;
      kRefl = 0.1;
      specExp = 25;
      ior = -1;
    }
    else if (m == MIRROR) {
      kDiffuse = 0;
      specExp = -1;
      kRefl = 1;
      ior = -1;
    }
    else if(m == GLASS) {
      kDiffuse = 0;
      specExp = -1;
      kRefl = 0;
      ior = 1.3;
    }
  }
  virtual ~Object() {}
  virtual bool intersect(const Vec3d &orig, const Vec3d &dir, double &tMin, double &tMax) = 0;
  virtual Vec3d getNormal(const Vec3d &point) const = 0;
};

class Parallelogram : public Object {
public:
  Vec3d corner;
  Vec3d edge1;
  Vec3d edge2;

  Parallelogram(const Vec3d &c, const Vec3d &e1, const Vec3d &e2, MaterialType m) :
    Object(m), corner(c), edge1(e1), edge2(e2) {}


  Vec3d getNormal(const Vec3d &point = Vec3d(0)) const
  {
    return -normalize(vectProduct(edge1, edge2));
  }

  // k * [dir] = ([corner] - [orig]) + a * [edge1] + b * [edge2]
  bool intersect(const Vec3d &orig, const Vec3d &dir, double &tMin, double &tMax)
  {
    Vec3d N = this->getNormal();

    Vec3d CO = corner - orig;
    Vec3d dirEdge1 = normalize(edge1);

    double system[3][3] = { dir.x, -edge1.x, -edge2.x,
                            dir.y, -edge1.y, -edge2.y,
                            dir.z, -edge1.z, -edge2.z};

    double f[3] = {CO.x, CO.y, CO.z};

    double k,a,b;
    solveSystem(system,f,k,a,b); // Gauss

    if(0 <= a && a <= 1 && 0 <= b && b <= 1)
    {
      tMin = tMax = k;
      return true;
    }
    else { return false; }
  }
};

class Sphere : public Object {
public:
  Vec3d center;
  double rad;
  double rad2;

  Sphere(const Vec3d &c, double r, MaterialType m) :
    center(c), rad(r), rad2(r*r), Object(m) {}

  bool intersect(const Vec3d &orig, const Vec3d &dir, double &tMin, double &tMax)
  {
    Vec3d OC = orig - center;
    double a = dotProduct(dir,dir);
    double b = 2 * dotProduct(OC, dir);
    double c = dotProduct(OC, OC) - rad2;
    bool res = solveQuadratic(a, b, c, tMin, tMax);
    return res;
  }

  Vec3d getNormal(const Vec3d &point) const
  {
    return normalize(point - center);
  }
};


// Основные функции:

/*
          (Light)
            |\    /|\    /|
              \    | N  /
               \   |   /
  dirToLight    \  |  /  reflDir
                 \ | /
                  \+/
                hitPoint
*/

bool trace (const Vec3d &orig,
            const Vec3d &dir,
            const vector<Object*> &objects,
            Object*& hitObject,
            Vec3d &hitPoint)
{
  hitObject = NULL;
  double closestDist = INFINITY;

  double tmin, tmax;

  for (int i = 0; i < objects.size(); i++)
  {
    if (objects[i]->intersect(orig, dir, tmin, tmax))
    {
      if (tmin > 0 && tmin < closestDist) {
        closestDist = tmin;
        hitObject = objects[i];
      }
      else if (tmax > 0 && tmax < closestDist) {
        closestDist = tmax;
        hitObject = objects[i];
      }
    }
  }

  if (hitObject != NULL)
  {
    hitPoint = orig + dir * closestDist;
    return true;
  }
  else
  {
    return false;
  }
}

bool hasIntersection(Vec3d orig,
                     Vec3d dir,
                     const vector<Object*> &objects,
                     LightType lightType)
{
  double tmin, tmax;

  if (lightType == AMBIENT)
    return false;

  for (int i = 0; i < objects.size(); i++)
  {
    if (objects[i]->intersect(orig, dir, tmin, tmax))
    {
      if (lightType == DIRECTIONAL && (tmin > 0 || tmax > 0))
        return true;
      else if (lightType == POINT && ((tmin > 0 && tmin < 1) || (tmax > 0 && tmax < 1)))
        return true;
    }
  }
  return false;
}

void computeShadows(const Vec3d &hitPoint,
                    const Vec3d &N,
                    bool inShadow[],
                    const vector<Object*> &objects,
                    const vector<Light*> &lights)
{
  Vec3d dirToLight, toLight;

  for (int i = 0; i < lights.size(); i++)
  {
    switch (lights[i]->type) {
      case AMBIENT:
        inShadow[i] = false;
        continue;
      case POINT:
        toLight = lights[i]->source - hitPoint;
        break;
      case DIRECTIONAL:
        toLight = -lights[i]->direction;
        break;
    }
    inShadow[i] = hasIntersection(hitPoint, toLight, objects, lights[i]->type);
  }
}

Vec3d computeSpecular (const Vec3d &hitPoint,
                       const Vec3d &minusDir,
                       Vec3d N,
                       int specExp,
                       const vector<Light*> &lights,
                       bool inShadow[])
{
  Vec3d retColor(0);
  Vec3d dirToLight;

  for (int i = 0; i < lights.size(); i++)
  {
    if (!inShadow[i])
    {
      switch (lights[i]->type) {
        case AMBIENT:
          continue;
        case POINT:
          dirToLight = normalize(lights[i]->source - hitPoint);
          break;
        case DIRECTIONAL:
          dirToLight = normalize(-lights[i]->direction);
          break;
      }

      if (specExp != -1) {
        Vec3d reflDir = normalize(reflectRay(dirToLight, N));
        double cosA = cos(reflDir, minusDir);
        if (cosA > 0) {
          retColor += lights[i]->color * pow(cosA, specExp);
        }
      }
    }
  }

  return retColor;
}

Vec3d computeDiffuse (const Vec3d &hitPoint,
                      const Vec3d &N,
                      vector<Light*> &lights,
                      bool inShadow[])
{
  Vec3d retColor(0);
  Vec3d dirToLight;

  for (int i = 0; i < lights.size(); i++)
  {
    if (!inShadow[i])
    {
      switch (lights[i]->type) {
        case AMBIENT:
          retColor += lights[i]->color;
          continue;
        case POINT:
          dirToLight = normalize(lights[i]->source - hitPoint);
          break;
        case DIRECTIONAL:
          dirToLight = normalize(-lights[i]->direction);
          break;
      }

      retColor += lights[i]->color * clamp(abs(dotProduct(dirToLight, N)));
    }
  }
  return retColor;
}

Vec3d castRay(Vec3d orig,
              Vec3d dir,
              vector<Object*> objects,
              vector<Light*> lights,
              int recDepth)
{
  Vec3d retColor = BACKGROUND_COLOR;

  Object* hitObject;
  Vec3d hitPoint, reflectedDir, minusDir, N;
  Vec3d objColor;
  double kD, kR, ior;
  int specExp;
  bool inShadow[lights.size()];
  Vec3d reflectionOrig;
  Vec3d refractionOrig;
  double k;

  if (recDepth == 0) { return retColor; }

  if (trace(orig, dir, objects, hitObject, hitPoint))
  {
    objColor = hitObject->color;
    specExp  = hitObject->specExp;
    kD = hitObject->kDiffuse;
    kR = hitObject->kRefl;
    ior = hitObject->ior;
    N  = hitObject->getNormal(hitPoint);
    reflectedDir = reflectRay(-normalize(dir), N);
    reflectionOrig = (dotProduct(dir, N) < 0)? (hitPoint + EPS * N): (hitPoint - EPS * N);

    computeShadows(reflectionOrig, N, inShadow, objects, lights);

    retColor = Vec3d(0);

    if (kD) // DIFFUSE_AND_GLOSSY
    {
      retColor += kD * computeDiffuse(reflectionOrig, N, lights, inShadow);

      if (specExp != -1 && (1-kD))
        retColor += (1-kD)  * computeSpecular(reflectionOrig, -normalize(dir),
                                              N, specExp, lights, inShadow);
    }
    else if (ior == -1) // MIRROR
    {
        retColor += castRay(reflectionOrig, reflectedDir, objects, lights, recDepth-1);
    }
    else // GLASS
    {
      if (kR)
      {
        retColor += kR * castRay(reflectionOrig, reflectedDir,
                                 objects, lights, recDepth-1);
      }

      if (1-kR)
      {
        Vec3d refractionDirection = normalize(refractRay(dir, N, hitObject->ior));
        refractionOrig = (dotProduct(refractionDirection, N) < 0)?
                         (hitPoint - N * EPS): (hitPoint + N * EPS);

        retColor += (1-kR) * castRay(refractionOrig, refractionDirection,
                                     objects, lights, recDepth-1);
      }
    }

    retColor = retColor * objColor;
  }

  return retColor;
}

void createScene(vector<Object*> &objects, vector<Light*> &lights)
{

  { // spheres
    Sphere* mirrorSph1 = new Sphere(Vec3d(3,-3,13), 2, MIRROR);
    mirrorSph1->color = Vec3d(0.9);
    objects.push_back(mirrorSph1);

    Sphere* sph2 = new Sphere(Vec3d(3,0,13), 1.25, DIFFUSE);
    sph2->color = Vec3d(1);
    objects.push_back(sph2);

    Sphere* sph3 = new Sphere(Vec3d(3,1.85,13), 0.75, DIFFUSE);
    sph3->color = Vec3d(1);
    objects.push_back(sph3);

    Sphere* glassSph = new Sphere(Vec3d(1.5,-1,10), 1, GLASS);
    glassSph->color = Vec3d(1);
    objects.push_back(glassSph);

    Sphere* normalSph = new Sphere(Vec3d(-3.5,-3.5,10), 1.5, DIFFUSE_AND_GLOSSY);
    normalSph->color = Vec3d(1);
    objects.push_back(normalSph);
  }

  { // back wall
    Parallelogram* pgBack = new Parallelogram(Vec3d(6,-5.5,-5), Vec3d(-12,0,0),
                                              Vec3d(0,11,0), DIFFUSE_AND_GLOSSY);
    pgBack->color = Vec3d(0,0,1);
    objects.push_back(pgBack);
  }

  { // front wall
    Parallelogram* pgFront = new Parallelogram(Vec3d(-6,-6,15), Vec3d(12,0,0),
                                               Vec3d(0,12,0), DIFFUSE_AND_GLOSSY);
    pgFront->color = Vec3d(1,1,0.3);
    objects.push_back(pgFront);
  }

  { // right wall, window 3x3
    // far
    Parallelogram* pgRight1 = new Parallelogram(Vec3d(5,6,16), Vec3d(0,-12,0),
                                               Vec3d(0,0,-5), DIFFUSE_AND_GLOSSY);
    pgRight1->color = Vec3d(0,1,1);
    objects.push_back(pgRight1);
    // top
    Parallelogram* pgRight2 = new Parallelogram(Vec3d(5,6,11), Vec3d(0,-4,0),
                                                Vec3d(0,0,-2), DIFFUSE_AND_GLOSSY);
    pgRight2->color = Vec3d(0,1,1);
    objects.push_back(pgRight2);
    // bottom
    Parallelogram* pgRight3 = new Parallelogram(Vec3d(5,-1,11), Vec3d(0,-5,0),
                                                Vec3d(0,0,-2), DIFFUSE_AND_GLOSSY);
    pgRight3->color = Vec3d(0,1,1);
    objects.push_back(pgRight3);
    // near
    Parallelogram* pgRight4 = new Parallelogram(Vec3d(5,6,9), Vec3d(0,-12,0),
                                               Vec3d(0,0,-15), DIFFUSE_AND_GLOSSY);
    pgRight4->color = Vec3d(0,1,1);
    objects.push_back(pgRight4);
  }

  { // left wall, window 3x3
    // far
    Parallelogram* pgLeft1 = new Parallelogram(Vec3d(-5,6,16), Vec3d(0,-12,0),
                                               Vec3d(0,0,-5), DIFFUSE_AND_GLOSSY);
    pgLeft1->color = Vec3d(1,0,1);
    objects.push_back(pgLeft1);
    // top
    Parallelogram* pgLeft2 = new Parallelogram(Vec3d(-5,6,11), Vec3d(0,-4,0),
                                                Vec3d(0,0,-2), DIFFUSE_AND_GLOSSY);
    pgLeft2->color = Vec3d(1,0,1);
    objects.push_back(pgLeft2);
    // bottom
    Parallelogram* pgLeft3 = new Parallelogram(Vec3d(-5,-1,11), Vec3d(0,-5,0),
                                                Vec3d(0,0,-2), DIFFUSE_AND_GLOSSY);
    pgLeft3->color = Vec3d(1,0,1);
    objects.push_back(pgLeft3);
    // near
    Parallelogram* pgLeft4 = new Parallelogram(Vec3d(-5,6,9), Vec3d(0,-12,0),
                                               Vec3d(0,0,-15), DIFFUSE_AND_GLOSSY);
    pgLeft4->color = Vec3d(1,0,1);
    objects.push_back(pgLeft4);
  }

  { // bottom wall
    Parallelogram* pgBot = new Parallelogram(Vec3d(-6,-5,16), Vec3d(12,0,0),
                                             Vec3d(0,0,-22), DIFFUSE_AND_GLOSSY);
    pgBot->color = Vec3d(1,0,0);
    objects.push_back(pgBot);
  }

  { // top wall, window 3x3
    Parallelogram* pgTop1 = new Parallelogram(Vec3d(6,5,16), Vec3d(-12,0,0),
                                              Vec3d(0,0,-14.5), DIFFUSE_AND_GLOSSY);
    pgTop1->color = Vec3d(0,1,0);
    objects.push_back(pgTop1);

    Parallelogram* pgTop2 = new Parallelogram(Vec3d(6,5,-1.5), Vec3d(-12,0,0),
                                              Vec3d(0,0,-4.5), DIFFUSE_AND_GLOSSY);
    pgTop2->color = Vec3d(0,1,0);
    objects.push_back(pgTop2);

    Parallelogram* pgTop3 = new Parallelogram(Vec3d(6,5,1.5), Vec3d(-4.5,0,0),
                                              Vec3d(0,0,-3), DIFFUSE_AND_GLOSSY);
    pgTop3->color = Vec3d(0,1,0);
    objects.push_back(pgTop3);

    Parallelogram* pgTop4 = new Parallelogram(Vec3d(-1.5,5,1.5), Vec3d(-4.5,0,0),
                                              Vec3d(0,0,-3), DIFFUSE_AND_GLOSSY);
    pgTop4->color = Vec3d(0,1,0);
    objects.push_back(pgTop4);
  }

  Parallelogram *pg1 = new Parallelogram(Vec3d(-3.53,-1.8,11), Vec3d(0,0,-2),
  Vec3d(7,-3.5,0), DIFFUSE_AND_GLOSSY);
  pg1->color = Vec3d(0,1,1);
  objects.push_back(pg1);


  { // lights
    Light* ambientLight = new Light(0.2, AMBIENT);
    lights.push_back(ambientLight);

    Light* rightLight = new Light(0.4, POINT);
    rightLight->source = Vec3d(8,1,9.5);
    lights.push_back(rightLight);

    Light* leftLight = new Light(0.3, POINT);
    leftLight->source = Vec3d(-8,1,9.5);
    lights.push_back(leftLight);

    Light* l3 = new Light(0.5, POINT);
    l3->source = Vec3d(0,4.5,0);
    lights.push_back(l3);
  }
}


void doEverything(uint32_t image[IMG_HEIGHT][IMG_WIDTH])
{
  double xScale = tan(FOV / 2.0);
  double yScale = (xScale * IMG_HEIGHT) / IMG_WIDTH;

  Vec3d camera(0,0,0);

  // Initialize objects:

  vector<Object*> objects;
  vector<Light*> lights;
  createScene(objects, lights);

  // Fill image:

  double zz, xx, yy;
  Vec3d dir, color;
  uint32_t red, green, blue;

  for (int y = -IMG_HEIGHT/2; y < IMG_HEIGHT/2; y++) {
    for (int x = -IMG_WIDTH/2; x < IMG_WIDTH/2; x++) {
      zz = 1;
      xx = (x + 0.5) / IMG_WIDTH  * 2*zz * xScale;
      yy = (y + 0.5) / IMG_HEIGHT * 2*zz * yScale;
      dir = normalize(Vec3d(xx, yy, zz));

      color = castRay(camera, dir, objects, lights, REC_DEPTH);

      red = uint32_t(clamp(color.x) * 255);
      green = uint32_t(clamp(color.y) * 255) << 8;
      blue = uint32_t(clamp(color.z) * 255) << 16;

      image[y+IMG_HEIGHT/2][x+IMG_HEIGHT/2] = red + green + blue;
    }
  }

  // Clean pointers:
  for (Light* l : lights)
    delete l;
  for (Object* o : objects)
    delete o;
}

int main(int argc, const char** argv)
{
  unordered_map<string, string> cmdLineParams;

  for(int i=0; i<argc; i++)
  {
    string key(argv[i]);

    if(key.size() > 0 && key[0]=='-')
    {
      if(i != argc-1) // not last argument
      {
        cmdLineParams[key] = argv[i+1];
        i++;
      }
      else
        cmdLineParams[key] = "";
    }
  }

  string outFilePath = "zout.bmp";
  if(cmdLineParams.find("-out") != cmdLineParams.end())
    outFilePath = cmdLineParams["-out"];

  int sceneId = 0;
  if(cmdLineParams.find("-scene") != cmdLineParams.end())
    sceneId = atoi(cmdLineParams["-scene"].c_str());

  uint32_t color = 9999;
  if(sceneId == 1)
    color = RED;
  else if(sceneId == 2)
    color = RED | GREEN;
  else if(sceneId == 3)
    color = BLUE;
  uint32_t image[IMG_HEIGHT][IMG_WIDTH];
  for (int i = 0; i < IMG_HEIGHT; i++)
    for (int j = 0; j < IMG_WIDTH; j++)
      image[i][j] = color;

  doEverything(image); // MY CODE

  SaveBMP(outFilePath.c_str(), (uint32_t*)image, IMG_WIDTH, IMG_HEIGHT);

  cout << "end." << endl;
  return 0;
}
