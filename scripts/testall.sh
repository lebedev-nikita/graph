./clean.sh
mkdir tmp
cd tmp
mv zout.bmp old.bmp
cmake ..
make
cd ..

clang++ -Xpreprocessor -fopenmp -lomp -Wc++11-extensions main.cpp Bitmap.cpp  -o tests/rt

time tmp/rt
time tests/rt
