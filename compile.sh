g++ -c -o GETDecoder.o GETDecoder.cc `root-config --cflags --glibs` -Wall -O2
g++ -c -o GETPad.o GETPad.cc `root-config --cflags --glibs` -Wall -O2
g++ -o analyzeTPC GETDecoder.o GETPad.o analyzeTPC.cc `root-config --cflags --glibs` -fopenmp -Wall -O2
