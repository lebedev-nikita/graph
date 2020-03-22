# ./clean.sh
mkdir tmp
cd tmp
mv zout.bmp old.bmp
cmake ..
make
./rt
open old.bmp
sleep 0.2
open zout.bmp
cd ..
