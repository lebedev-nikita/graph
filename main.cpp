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

/*
    TODO:
      * настроить так, чтобы все работало при разных ширине и высоте
        (попробовать image = malloc(...) ?)
      * POINT: сделать, чтобы свет отражался и от внутренней поверхности сферы,
        если она окружает камеру
*/
#define IMG_WIDTH 1024
#define IMG_HEIGHT 1024
#define FOV M_PI/3
#define BACKGROUND_COLOR Vec3d(0,0.7,1)

// Векторы

class Vec3d {
public:
  double x, y, z;

  Vec3d(): x(0), y(0), z(0)
  {}
  Vec3d(double xx, double yy, double zz): x(xx), y(yy), z(zz)
  {}
  Vec3d operator * (const double n) const
  { return Vec3d(this->x * n, this->y * n, this->z * n); }
  Vec3d operator * (const Vec3d &v) const
  { return Vec3d(this->x * v.x, this->y * v.y, this->z * v.z); }
  Vec3d operator - (const Vec3d &v) const
  { return Vec3d(this->x - v.x, this->y - v.y, this->z - v.z); }
  Vec3d operator - () const
  { return Vec3d(-x, -y, -z); }
  Vec3d operator + (const Vec3d &v) const
  { return Vec3d(this->x + v.x, this->y + v.y, this->z + v.z); }
  Vec3d operator += (const Vec3d &v)
  {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
  }
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

// Свет

enum LightType { AMBIENT, POINT, DIRECTIONAL};

class Light {
public:
  Vec3d color;
  LightType lightType;

  Vec3d source;
  Vec3d direction;

  Light(const Vec3d &clr, LightType type) :
    color(clr), lightType(type), source(0,0,0), direction(0,0,0) {}
};

// Объекты

class Object {
public:
  Vec3d color;

  Object() : color(1,0,0) {};
  virtual ~Object() {}
  virtual bool intersect(const Vec3d &orig, const Vec3d &dir, double &tMin, double &tMax) = 0;
  virtual Vec3d getNormal(const Vec3d &point) = 0;
};

class Sphere : public Object {
public:
  Vec3d center;
  double rad;
  double rad2;

  Sphere(const Vec3d &c, double r) : center(c), rad(r), rad2(r*r)
  {}

  bool intersect(const Vec3d &orig, const Vec3d &dir, double &tMin, double &tMax)
  {
    Vec3d OC = orig - center;
    double a = dotProduct(dir,dir);
    double b = 2 * dotProduct(OC, dir);
    double c = dotProduct(OC, OC) - rad2;
    bool res = solveQuadratic(a, b, c, tMin, tMax);
    return res;
  }

  Vec3d getNormal(const Vec3d &point)
  {
    return normalize(point - center);
  }
};

// Основные функции

Vec3d computeIllumination (const Vec3d &hitPoint, const Vec3d &N, vector<Light*> &lights)
{
  Vec3d illumination(0,0,0);
  Vec3d lightDir;

  for (int i = 0; i < lights.size(); i++)
  {
    switch (lights[i]->lightType)
    {
      case AMBIENT:
        illumination += lights[i]->color;
        break;
      case POINT:
        lightDir = normalize(hitPoint - lights[i]->source);
        illumination += lights[i]->color * clamp(dotProduct(-lightDir, N));
        break;
      case DIRECTIONAL:
        // TODO
        break;
    }
  }
  return illumination;
}

bool trace (Vec3d orig,
            Vec3d dir,
            vector<Object*> &objects,
            Object*& hitObject,
            Vec3d &hitPoint)
{
  hitObject = NULL;
  double closestDist = INFINITY;

  Object* obj;
  double tMin, tMax;

  for (int i = 0; i < objects.size(); i++)
  {
    if (objects[i]->intersect(orig, dir, tMin, tMax))
    {
      if (tMin > 0 && tMin < closestDist) {
        closestDist = tMin;
        hitObject = objects[i];
      }
      else if (tMax > 0 && tMax < closestDist) {
        closestDist = tMax;
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

Vec3d castRay(Vec3d orig, Vec3d dir, vector<Object*> objects, vector<Light*> lights)
{
  Vec3d color = BACKGROUND_COLOR;

  Object* hitObject;
  Vec3d hP(0,0,0);

  if (trace(orig, dir, objects, hitObject, hP))
  {
    color = hitObject->color * computeIllumination(hP, hitObject->getNormal(hP), lights);
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

  Sphere* sph1 = new Sphere(Vec3d(-3,0,16), 2);
  sph1->color = Vec3d(1,0,0);
  Sphere* sph2 = new Sphere(Vec3d(0,0,16), 2);
  sph2->color = Vec3d(0,1,0);
  Sphere* sph3 = new Sphere(Vec3d(3,0,16), 2);
  sph3->color = Vec3d(0,0,1);
  // Sphere* sph4 = new Sphere(Vec3d(0,0,0), 16);
  // sph4->color = Vec3d(0,0.7,0.9);

  objects.push_back(sph1);
  objects.push_back(sph2);
  objects.push_back(sph3);
  // objects.push_back(sph4);

  vector<Light*> lights;

  // TODO: проработать, чтобы сумма освещения не была больше 1
  Light* light1 = new Light(Vec3d(0.2, 0.2, 0.2), AMBIENT);
  Light* light2 = new Light(Vec3d(0.3, 0.3, 0.3), POINT);
  light2->source = Vec3d(-4,3,10);
  Light* light4 = new Light(Vec3d(0.3, 0.3, 0.3), POINT);
  light4->source = Vec3d(4,3,10);
  // Light* light3 = new Light(Vec3d(0.2, 0.2, 0.2), DIRECTIONAL);
  // light3->direction = Vec3d(1, 4, 4);

  lights.push_back(light1);
  lights.push_back(light2);
  lights.push_back(light4);
  // lights.push_back(light3);


  /*
          xx
        ______    xx = zz * tg(a)
        |   /     xx = IMG_WIDTH/2 * ?
     zz |⏜/
        |a/       ? = zz * tg(a) * 2 / IMG_WIDTH
        |/        xx = (IMG_WIDTH/2) * zz * tg(a) * 2 / IMG_WIDTH
  */


  for (int y = -IMG_HEIGHT/2; y < IMG_HEIGHT/2; y++) {
    for (int x = -IMG_WIDTH/2; x < IMG_WIDTH/2; x++) {
      double zz = 1;
      double xx = (x + 0.5) / IMG_WIDTH  * 2*zz * xScale;
      double yy = (y + 0.5) / IMG_HEIGHT * 2*zz * yScale;
      Vec3d dir = normalize(Vec3d(xx, yy, zz));

      Vec3d color = castRay(camera, dir, objects, lights);

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
