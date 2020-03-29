# scripts/clean.sh
mkdir tmp
cd tmp
mv zout.bmp old.bmp
cmake ..
make
./rt
cd ..
./op.sh
