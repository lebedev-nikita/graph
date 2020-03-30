# ./clean.sh
mkdir tmp
cd tmp
mv zout.bmp old.bmp
cmake ..
make
time ./rt
cd ..
./op.sh
