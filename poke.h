// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cinttypes>
#include <deque>
#include <queue>
#include <limits>
#include <math.h> 
#include <unordered_map>

using namespace std;

class Poke {
    private:
    
    int totalPokemons;
    double totalDistance = 0;

    // for cin in readInputFile
    int xCoordinate;
    int yCoordinate;


    public:

    struct Modes{
        bool mst = false;
        bool fasttsp = false;
        bool opttsp = false;
    };

    // struct just for verticies 
    struct Vertices{
        int x;
        int y;
        bool landPokemon;
        bool coastPokemon;
        bool seaPokemon;
        // enum or char instead of bool
    };
    // store it into a vector of the verticies struct parent is the index of the vector
    vector<Vertices> allVertices;

    struct Pokemons{
        double distance = numeric_limits<double>::infinity();
        int parent = -1;
        bool visited = false;
    };
    vector<Pokemons> allPokemons;

    void readInputFile(bool mst1, bool fasttsp1, bool opttsp1);

    // PartA
    //--------------------------------------------------------------------------------------------
    
    void terrainType (Vertices &vertex);

    void getMode(int argc, char * argv[], Modes &mode1);

    int minDistanceVertex (vector<Pokemons> const &allPokemons1);

    void linearPrimsAlgorithm (vector<Vertices> const &allVertices1);

    double euclideanDistance (Vertices vert1, Vertices vert2);

    void outPut ();
    //--------------------------------------------------------------------------------------------

    // PartB

    vector<int> tour;

    void fasttsp(vector<Vertices> const &allVertices2);

    double bEuclideanDistance (Vertices vert1, Vertices vert2);

    void tspOutPut();

    //--------------------------------------------------------------------------------------------

    // PartC

    vector<int> path;

    // full permutation, start as part B's path
    vector<int> currBestVec;

    //  start at 0
    double currBestDistance = 0;

    // upper bound (TOTAL) weight start as part B's
    double upperBound;

    double grayDistance;

    void genPerms(int permLength);

    bool promising(vector<int> const &path1, int permLength1);

    void partCPrims (vector<Vertices> allVertices1, int permLength1);

    void partB(vector<Vertices> allVertices1);

    void optOutput();
    

};