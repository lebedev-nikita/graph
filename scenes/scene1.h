/* Initialize objects: */
vector<Object*> objects;

Sphere* sph1 = new Sphere(Vec3d(-3,0,16), 1.4, MIRROR);
sph1->color = Vec3d(1,0,1);
objects.push_back(sph1);

Sphere* sph2 = new Sphere(Vec3d(0,0,16), 1.4, GLASS);
sph2->color = Vec3d(1,1,0);
objects.push_back(sph2);

Sphere* sph3 = new Sphere(Vec3d(3,0,16), 1.4, GLASS);
sph3->color = Vec3d(0,1,1);
objects.push_back(sph3);

Sphere* backSph = new Sphere(Vec3d(0, 0, 24), 4, DIFFUSE_AND_GLOSSY);
backSph->color = Vec3d(1,0,1);
objects.push_back(backSph);

Sphere* bigSph = new Sphere(Vec3d(0,0,16), 17, DIFFUSE_AND_GLOSSY);
bigSph->color = Vec3d(0,1,0);
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
