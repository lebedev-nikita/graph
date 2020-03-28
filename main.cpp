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
    TODO: NEXT:
      * добавить преломление
      * добавить меши
      * распараллелить
    TODO: LATER:
      * сгладить пиксели
      * настроить так, чтобы все работало при разных ширине и высоте
        (попробовать image = malloc(...) ?)
      * POINT: сделать, чтобы свет отражался и от внутренней поверхности сферы,
        если она окружает камеру
*/
#define IMG_WIDTH 1024
#define IMG_HEIGHT 1024
#define FOV M_PI/3
#define BACKGROUND_COLOR Vec3d(1)
#define EPS 0.0000001
#define REC_DEPTH 3


/* Векторы: */

class Vec3d {
public:
  double x, y, z;

  Vec3d(): x(0), y(0), z(0)
  {}
  Vec3d(double xx, double yy, double zz): x(xx), y(yy), z(zz)
  {}
  Vec3d(double xx): x(xx), y(xx), z(xx)
  {}
  Vec3d operator * (const double n) const
  { return Vec3d(this->x * n, this->y * n, this->z * n); }
  Vec3d operator * (const Vec3d &v) const
  { return Vec3d(this->x * v.x, this->y * v.y, this->z * v.z); }
  friend Vec3d operator * (const float &r, const Vec3d &v)
  { return Vec3d(v.x * r, v.y * r, v.z * r); }
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


/* Вспомогательные функции: */

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

Vec3d reflectRay(const Vec3d &R, const Vec3d &N) {
    return 2*N*dotProduct(R, N) - R;
}


/* Свет: */

enum LightType { AMBIENT, POINT, DIRECTIONAL};

class Light {
public:
  Vec3d color;
  LightType type;

  Vec3d source = Vec3d(0,0,0);
  Vec3d direction = Vec3d(0,0,0);

  Light(const Vec3d &clr, LightType type) : color(clr), type(type) {}
};


/* Объекты: */

enum MaterialType { DIFFUSE, DIFFUSE_AND_GLOSSY, MIRROR, ALL_IN_ONE };

class Object {
public:
  Vec3d color = Vec3d(1,1,1);
  MaterialType material = DIFFUSE_AND_GLOSSY;
  double kDiffuse;
  double kSpecular;
  double kReflection;
  int specExp;

  Object(MaterialType m = DIFFUSE_AND_GLOSSY)
  {
    material = m;
    if (m == DIFFUSE) {
      kDiffuse = 1;
      kSpecular = 0;
      kReflection = 0;
      specExp = -1;
    }
    else if (m == DIFFUSE_AND_GLOSSY) {
      kDiffuse = 0.5;
      kSpecular = 1;
      kReflection = 0;
      specExp = 25;
    }
    else if (m == MIRROR) {
      kDiffuse = 0;
      kSpecular = 0;
      kReflection = 1;
      specExp = -1;
    }
    else if (m == ALL_IN_ONE) {
      kDiffuse = 0.5;
      kSpecular = 100;
      kReflection = 1;
      specExp = 25;
    }
  }
  virtual ~Object() {}
  virtual bool intersect(const Vec3d &orig, const Vec3d &dir, double &tMin, double &tMax) = 0;
  virtual Vec3d getNormal(const Vec3d &point) const = 0;
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


/* Основные функции: */

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
            Vec3d &hitPoint,
            double &tMin,
            double &tMax)
{
  hitObject = NULL;
  double closestDist = INFINITY;

  Object* obj;
  double tmin, tmax;

  for (int i = 0; i < objects.size(); i++)
  {
    if (objects[i]->intersect(orig, dir, tmin, tmax))
    {
      if (tmin > 0 && tmin < closestDist) {
        closestDist = tmin;
        tMin = tmin;
        tMax = tmax;
        hitObject = objects[i];
      }
      else if (tmax > 0 && tmax < closestDist) {
        closestDist = tmax;
        tMin = tmin;
        tMax = tmax;
        hitObject = objects[i];
      }
    }
  }

  if (hitObject != NULL)
  {
    hitPoint = orig + dir * (closestDist - EPS);
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
                       const Vec3d &dir,
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
        double cosA = clamp(dotProduct(reflDir, -dir));
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

      retColor += lights[i]->color * clamp(dotProduct(dirToLight, N));
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
  Vec3d hitPoint(0,0,0);
  double tMin, tMax;
  Vec3d objColor;
  int specExp;
  double kD;
  double kS;
  double kR;
  Vec3d N;
  Vec3d reflectedDir;
  bool inShadow[lights.size()];

  if (recDepth == 0) return BACKGROUND_COLOR;

  if (trace(orig, dir, objects, hitObject, hitPoint, tMin, tMax))
  {
    objColor = hitObject->color;
    specExp  = hitObject->specExp;
    kD = hitObject->kDiffuse;
    kS = hitObject->kSpecular;
    kR = hitObject->kReflection;
    N  = hitObject->getNormal(hitPoint);
    reflectedDir = reflectRay(-normalize(dir), N);

    computeShadows(hitPoint, N, inShadow, objects, lights);

    retColor = Vec3d(0);
    double k = 1;

    if (kD)
      retColor += objColor * k * kD * computeDiffuse(hitPoint, N, lights, inShadow);
    k *= clamp(1-kD);

    if (specExp != -1 && kS)
      retColor += objColor * kS  * computeSpecular(hitPoint, dir, N, specExp, lights, inShadow);
    k *= clamp(1-kS);

    if (kR)
      retColor += objColor * k * kR * castRay(hitPoint, reflectedDir, objects, lights, recDepth-1);
    k *= clamp(1-kR);
  }

  return retColor;
}


void doEverything(uint32_t image[IMG_HEIGHT][IMG_WIDTH])
{
  double xScale = tan(FOV / 2.0);
  double yScale = (xScale * IMG_HEIGHT) / IMG_WIDTH;

  Vec3d camera(0,0,0);
  uint32_t color;

  /* Initialize objects: */
  vector<Object*> objects;

  Sphere* sph1 = new Sphere(Vec3d(-3,0,16), 1.4, DIFFUSE_AND_GLOSSY);
  sph1->color = Vec3d(1,0,1);
  objects.push_back(sph1);

  Sphere* sph2 = new Sphere(Vec3d(0,0,16), 1.4, DIFFUSE_AND_GLOSSY);
  sph2->color = Vec3d(1,1,0);
  objects.push_back(sph2);

  Sphere* sph3 = new Sphere(Vec3d(3,0,16), 1.4, DIFFUSE_AND_GLOSSY);
  sph3->color = Vec3d(0,1,1);
  objects.push_back(sph3);

  Sphere* bigSph = new Sphere(Vec3d(0,-160,16), 156, DIFFUSE_AND_GLOSSY);
  bigSph->color = Vec3d(1);
  objects.push_back(bigSph);

  /* Initialize lights: */
  vector<Light*> lights;

  Light* ambientLight = new Light(0.2, AMBIENT);
  lights.push_back(ambientLight);

  Light* topLight = new Light(0.5, POINT);
  topLight->source = Vec3d(0, 10, 16);
  lights.push_back(topLight);

  Light* directionalLight = new Light(Vec3d(1), DIRECTIONAL);
  directionalLight->direction = normalize(Vec3d(-1,-0.5,0));
  lights.push_back(directionalLight);

  Light* rightTopLight = new Light(0.4, POINT);
  rightTopLight->source = Vec3d(10,10,16);
  lights.push_back(rightTopLight);

  // Light* light4 = new Light(0.5, POINT);
  // light4->source = Vec3d(10, 10, 16);
  // lights.push_back(light4);


  /* Fill image: */

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

      Vec3d color = castRay(camera, dir, objects, lights, REC_DEPTH);

      uint32_t red = uint32_t(clamp(color.x) * 255);
      uint32_t green = uint32_t(clamp(color.y) * 255) << 8;
      uint32_t blue = uint32_t(clamp(color.z) * 255) << 16;

      image[y+IMG_HEIGHT/2][x+IMG_HEIGHT/2] = red + green + blue;
    }
  }

  /* Clean pointers: */
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
