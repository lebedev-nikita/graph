#include <iostream>
#include <cstdint>

#include <string>
#include <vector>
#include <unordered_map>

#include <cmath>

#include "Bitmap.h"

using namespace std;

const uint32_t RED   = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE  = 0x00FF0000;

#define IMG_WIDTH 512
#define IMG_HEIGHT 512
#define FOV M_PI/3
#define BACKGROUND_COLOR Vec3d(0,0.7,1)

// Векторы

class Vec3d {
public:
  double x, y, z;

  Vec3d(double xx, double yy, double zz): x(xx), y(yy), z(zz)
  {}
  Vec3d operator * (const double n) const
  { return Vec3d(this->x * n, this->y * n, this->z * n); }
  Vec3d operator - (const Vec3d v) const
  { return Vec3d(this->x - v.x, this->y -v.y, this->z - v.z); }
};

// Вспомогательные функции

double dotProduct(const Vec3d &v, const Vec3d &u)
{ return (v.x * u.x) + (v.y * u.y) + (v.z * u.z); }

Vec3d normalize(const Vec3d &vec)
{
  double len = sqrt(dotProduct(vec, vec));
  return vec * (1/len);
}

double clamp(double num) {
  if (num < 0) return 0;
  if (num > 1) return 1;
  else         return num;
}

bool solveQuadratic(double a, double b, double c, double tMin, double tMax)
{
  double discr = b*b - 4*a*c;
  if (discr < 0) {
    return false;
  }
  else if (discr == 0) {
    tMin = -b/(2*a);
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

// Объекты

class Object {
public:
  Vec3d color;

  Object() : color(1,0,0) {};
  virtual ~Object() {}
  virtual bool intersect(Vec3d orig, Vec3d dir, double &tMin, double &tMax) = 0;
};

class Sphere : public Object {
public:
  Vec3d center;
  double rad;
  double rad2;

  Sphere(const Vec3d &c, double r) : center(c), rad(r), rad2(r*r)
  {}

  bool intersect(Vec3d orig, Vec3d dir, double &tMin, double &tMax)
  {
    Vec3d OC = center - orig;
    double a = dotProduct(dir,dir);
    double b = 2 * dotProduct(OC, dir);
    double c = dotProduct(OC,OC) - rad2;
    return solveQuadratic(a, b, c, tMin, tMax);
  }
};

// Основные функции

bool trace (Vec3d orig, Vec3d dir, vector<Object*> objects, Object*& hitObject)
{
  hitObject = NULL;
  double closestDist = INFINITY;

  Object* obj;
  double t1, t2;

  for (int i = 0; i < objects.size(); i++)
  {
    if (objects[i]->intersect(orig, dir, t1, t2))
    {
      if (t1 > 0 && t1 < closestDist) {
        closestDist = t1;
        hitObject = objects[i];
      }
      else if (t2 > 0 && t2 < closestDist) {
        closestDist = t2;
        hitObject = objects[i];
      }
    }
  }

  return hitObject != NULL;
}

Vec3d castRay(Vec3d orig, Vec3d dir, vector<Object*> objects)
{
  Vec3d color = BACKGROUND_COLOR;

  Object* hitObject;
  if (trace(orig, dir, objects, hitObject)) {
    color = hitObject->color;
  }

  return color;
}


void doEverything(uint32_t image[IMG_HEIGHT][IMG_WIDTH])
{
  double xScale = tan(FOV / 2.0);
  double yScale = (xScale * IMG_HEIGHT) / IMG_WIDTH;

  Vec3d camera(0,0,0);
  uint32_t color;

  vector<Object*> objects;

  Sphere* sph1 = new Sphere(Vec3d(-3,0,16), 1);
  sph1->color = Vec3d(1,0,0);
  Sphere* sph2 = new Sphere(Vec3d(0,0,16), 1);
  sph2->color = Vec3d(0,1,0);
  Sphere* sph3 = new Sphere(Vec3d(3,0,16), 1);
  sph3->color = Vec3d(0,0,1);

  objects.push_back(sph1);
  objects.push_back(sph2);
  objects.push_back(sph3);



  for (int y = -IMG_HEIGHT/2; y < IMG_HEIGHT/2; y++) {
    for (int x = -IMG_WIDTH/2; x < IMG_WIDTH/2; x++) {
      double xx = (x + 0.5) / IMG_WIDTH * xScale;
      double yy = (y + 0.5) / IMG_HEIGHT * yScale;
      double zz = 1;
      Vec3d dir = normalize(Vec3d(xx, yy, zz));

      Vec3d color = castRay(camera, dir, objects);

      uint32_t red = uint32_t(clamp(color.x) * 255);
      uint32_t green = uint32_t(clamp(color.y) * 255) << 8;
      uint32_t blue = uint32_t(clamp(color.z) * 255) << 16;

      image[y+IMG_HEIGHT/2][x+IMG_HEIGHT/2] = red + green + blue;
    }
  }
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
