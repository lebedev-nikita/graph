#include <iostream>
#include <cstdint>

#include <string>
#include <vector>
#include <unordered_map>

#include <cmath>
// #include <omp.h>
// #include "/usr/local/opt/libomp/include/omp.h"

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
#define BACKGROUND_COLOR Vec3d(0)
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
  friend Vec3d operator * (const double &r, const Vec3d &v)
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
  friend inline double abs(const Vec3d &v)
  { return sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
};


/* Вспомогательные функции: */

inline double dotProduct(const Vec3d &v, const Vec3d &u)
{ return (v.x * u.x) + (v.y * u.y) + (v.z * u.z); }

inline Vec3d normalize(const Vec3d &vec)
{ return vec * (1 / sqrt(dotProduct(vec,vec))); }

inline double clamp(double num, double low = 0, double hi = 1)
{ return (num < low)?low:(num>hi)?hi:num; }

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

double kFresnel(const Vec3d &I, const Vec3d &N, const double &ior)
{
    double kr;
    double cosi = clamp(dotProduct(I, N), -1, 1);
    double etai = 1, etat = ior;
    if (cosi > 0) {  std::swap(etai, etat); }
    // Compute sini using Snell's law
    double sint = etai / etat * ((cosi * cosi < 1)?sqrt(1 - cosi * cosi):0);
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        double cost = (sint * sint < 1)?sqrtf(1 - sint * sint):0;
        cosi = abs(cosi);
        double Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        double Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // kt = 1 - kr;
    return kr;
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

enum MaterialType { DIFFUSE, DIFFUSE_AND_GLOSSY, MIRROR, ALL_IN_ONE, GLASS };

class Object {
public:
  Vec3d color = Vec3d(1,1,1);
  MaterialType material = DIFFUSE_AND_GLOSSY;
  double kDiffuse;
  double kSpecular;
  double kReflection;
  int specExp;
  double ior;

  Object(MaterialType m = DIFFUSE_AND_GLOSSY)
  {
    material = m;
    if (m == DIFFUSE) {
      kDiffuse = 1;
      kSpecular = 0;
      kReflection = 0;
      specExp = -1;
      ior = -1;
    }
    else if (m == DIFFUSE_AND_GLOSSY) {
      kDiffuse = 0.5;
      kSpecular = 1;
      kReflection = 0;
      specExp = 25;
      ior = -1;
    }
    else if (m == MIRROR) {
      kDiffuse = 0;
      kSpecular = 0;
      kReflection = 1;
      specExp = -1;
      ior = -1;
    }
    else if (m == ALL_IN_ONE) {
      kDiffuse = 0.5;
      kSpecular = 0.5;
      kReflection = 1;
      specExp = 25;
      ior = -1;
    }
    else if(m == GLASS) {
      kDiffuse = 0;
      kSpecular = 0;
      kReflection = 0.2;
      specExp = -1;
      ior = 1.3;

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
    hitPoint = orig + dir * closestDist;
    // hitPoint = orig + dir * (closestDist - EPS);
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
        double cosA = clamp(dotProduct(reflDir, minusDir), -1, 1);
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
  Vec3d hitPoint, reflectedDir, minusDir, N;
  Vec3d objColor;
  double kD, kS, kR, ior;
  double tMin, tMax;
  int specExp;
  bool inShadow[lights.size()];
  Vec3d reflectionOrig;
  Vec3d refractionOrig;
  double k;

  if (recDepth == 0) { return retColor; }

  if (trace(orig, dir, objects, hitObject, hitPoint, tMin, tMax))
  {
    objColor = hitObject->color;
    specExp  = hitObject->specExp;
    kD = hitObject->kDiffuse;
    kS = hitObject->kSpecular;
    kR = hitObject->kReflection;
    ior = hitObject->ior;
    N  = hitObject->getNormal(hitPoint);
    reflectedDir = reflectRay(-normalize(dir), N);
    if (dotProduct(dir, N) < 0) {
      reflectionOrig = hitPoint + EPS * N;
      refractionOrig = hitPoint - EPS * N;
    }
    else {
      reflectionOrig = hitPoint - EPS * N;
      refractionOrig = hitPoint + EPS * N;
    }


    computeShadows(reflectionOrig, N, inShadow, objects, lights);

    retColor = Vec3d(0);
    k = 1;

    if (kD)
      retColor += k * kD * computeDiffuse(reflectionOrig, N, lights, inShadow);
    k *= clamp(1-kD);

    if (specExp != -1 && kS)
      retColor += kS  * computeSpecular(reflectionOrig, -normalize(dir) , N, specExp, lights, inShadow);
    k *= clamp(1-kS);

    if (false);
    else {
      double kF = kFresnel(dir, N, hitObject->ior);
      if (kR)
      retColor += k * kR * castRay(reflectionOrig, reflectedDir, objects, lights, recDepth-1);
      k *= clamp(1-kR);

      if (ior != -1) {
        Vec3d refractionDirection = normalize(refractRay(dir, N, hitObject->ior));
        // Vec3d refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1, 1);
        // double kr;
        // fresnel(dir, N, hitObject->ior, kr);
        // hitColor = refractionColor * (1 - kr);

        retColor += k * castRay(refractionOrig, refractionDirection, objects, lights, recDepth-1);

      }
    }

    retColor = retColor * objColor;
  }

  return retColor;
}


void doEverything(uint32_t image[IMG_HEIGHT][IMG_WIDTH])
{
  double xScale = tan(FOV / 2.0);
  double yScale = (xScale * IMG_HEIGHT) / IMG_WIDTH;

  Vec3d camera(0,0,0);

  /* Initialize objects: */
  vector<Object*> objects;

  Sphere* sph1 = new Sphere(Vec3d(-3,0,16), 1.4, GLASS);
  sph1->color = Vec3d(1,0,1);
  objects.push_back(sph1);

  Sphere* sph2 = new Sphere(Vec3d(0,0,16), 1.4, GLASS);
  sph2->color = Vec3d(1,1,0);
  objects.push_back(sph2);

  Sphere* sph3 = new Sphere(Vec3d(3,0,16), 1.4, GLASS);
  sph3->color = Vec3d(0,1,1);
  objects.push_back(sph3);

  Sphere* backSph = new Sphere(Vec3d(0, 0, 30), 5, ALL_IN_ONE);
  backSph->color = Vec3d(1);
  objects.push_back(backSph);

  Sphere* bigSph = new Sphere(Vec3d(0,-160,16), 156, GLASS);
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


  /* Fill image: */

  /*
          xx
        ______    xx = zz * tg(a)
        |   /     xx = IMG_WIDTH/2 * ?
     zz |⏜/
        |a/       ? = zz * tg(a) * 2 / IMG_WIDTH
        |/        xx = (IMG_WIDTH/2) * zz * tg(a) * 2 / IMG_WIDTH
  */
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
