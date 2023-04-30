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
#include <unordered_map>
#include <limits>
#include <math.h> 
#include "poke.h"

using namespace std;


void Poke :: readInputFile(bool mst1, bool fasttsp1, bool opttsp1){
    Vertices verts;

    // read in total number of Pokemons
    cin >> totalPokemons;
    allVertices.reserve(totalPokemons);

    for (int i = 0; i < totalPokemons; ++i){
        //  read in x and y coordinates of the pokemons
        cin >> xCoordinate >> yCoordinate;
        
        // reset coordinates
        verts.x = xCoordinate;
        verts.y = yCoordinate;

        // reset bool values;
        verts.landPokemon = false;
        verts.seaPokemon = false;
        verts.coastPokemon = false;
        // set the terrain type of the vertices
        terrainType(verts);
        // store each vertex(x,y) into the vector of all vertices
        allVertices.push_back(verts);
        //path.push_back((uint32_t) i);
        
    }

    if (mst1 == true){
        linearPrimsAlgorithm(allVertices);
        outPut();
    }
    else if (fasttsp1 == true){
        fasttsp(allVertices);
        tspOutPut();
    }
    else if (opttsp1 == true){
        partB(allVertices);
        genPerms(1);
        optOutput();
    }
}

// PART A
//------------------------------------------------------------------------------------------------------------------------------------

// make a function that returns classification(Land, coast, sea) of a vertex
void Poke :: terrainType (Vertices &vertex){
    // land: x > 0 || y > 0
    if ( (vertex.x > 0 || vertex.y > 0) ){
        vertex.landPokemon = true;
    }
    // sea: x < 0 && y < 0
    else if ( (vertex.x < 0 && vertex.y < 0) ){
        vertex.seaPokemon = true;
    }
    // else: coast
    else{
        vertex.coastPokemon = true;
    }

    // ***** error checking ******
    // ***** if coast count is still == 0 when theres only sea and land-> cerr ******

}
// *Took out the &*
void Poke :: linearPrimsAlgorithm (vector<Vertices> const &allVertices1){
    Pokemons pokemon;
    // store pokemon objects into vector of allPokemons
    allPokemons.resize(allVertices1.size());
    
    allPokemons[0].distance = 0;

    for (int i = 0; i < totalPokemons; ++i){
        // find and store the min distance vertex's index: A first time around the loop
        int tempIndex = minDistanceVertex(allPokemons);
        
        allPokemons[tempIndex].visited = true; // set A (first Vertex) as visited

        totalDistance += sqrt(allPokemons[tempIndex].distance);

        // MST
        for (int j = 0; j < totalPokemons; ++j){ // loop through all pokemons
            if (allPokemons[j].visited == true){
                continue;
            }

            // check the distance between A and all other unvisited pokemons and calculate the distance in between
            double newDistance = euclideanDistance(allVertices1[tempIndex], allVertices1[j]);
            
            // if the new distance is < that unvisited pokemon's distance, set that pokemons distance to the new distance
            // and set that unvisited pokemon's parent to A
            if (newDistance < allPokemons[j].distance){
                allPokemons[j].distance = newDistance;
                allPokemons[j].parent = tempIndex;
            }
        }
    }

}

// returns index of vertex with smallest d value among vertices whose visted = false
int Poke :: minDistanceVertex (vector<Pokemons> const &allPokemons1){

    int minDistancePokemon = -1; // since distance cant be < 0
    double minDist = numeric_limits<double>::infinity(); // Temp infinity variable
    
    for (int i = 0; i < (int)allPokemons1.size(); ++i){
        if (allPokemons1[i].visited == false && allPokemons1[i].distance < minDist){
            minDist = allPokemons1[i].distance;
            minDistancePokemon = i;
        }
    }

    return minDistancePokemon;
}
 
double Poke :: euclideanDistance (Vertices vert1, Vertices vert2){

    double euclideanDistance = 0;

    if ( (vert1.landPokemon == true && vert2.seaPokemon == true) || 
          (vert1.seaPokemon == true && vert2.landPokemon == true)  ) {
        euclideanDistance = numeric_limits<double>::infinity();
        return euclideanDistance;
    }
    else {
        double xdifference = (double)(vert1.x - vert2.x);
        double ydifference = (double)(vert1.y - vert2.y);

        euclideanDistance = ( (xdifference) * (xdifference) ) + ( (ydifference) * (ydifference) );
        return euclideanDistance;
    }

    //     euclideanDistance = sqrt( ((xdifference) * (xdifference)) + ((ydifference) * (ydifference)) );

}

void Poke :: outPut(){
    cout << totalDistance << "\n";
    
    for (int j = 1; j < (int)allPokemons.size(); ++j){
        if (j < allPokemons[j].parent){
            cout << j << " " << allPokemons[j].parent << "\n";
        }
        else if (j > allPokemons[j].parent){
            cout << allPokemons[j].parent << " " << j << "\n";
        }
    }
}

//---------------------------------------------------------------------------------------------------------------//

// PartB
void Poke :: fasttsp(vector<Vertices> const &allVertices2){

    // store pokemon objects into vector of allPokemons
    // allPokemons.resize(allVertices2.size());
    // initialize vertex(Pokemon) i (initial/A)
    // Identify vertex j and set its parent as i, visted = true, distance from j and i 
    // insert index of pokemons from allPokemons in the tour path

    tour.reserve(totalPokemons);

    int i = 0;
    int j = 1;
    int k = 2;

    tour.push_back(i);
    tour.push_back(j);
    tour.push_back(k);

    totalDistance += 2 * sqrt((bEuclideanDistance(allVertices2[i], allVertices2[j])));

    // triangle A,B,C's total weight
    totalDistance += sqrt(bEuclideanDistance(allVertices2[i], allVertices2[k]));
    totalDistance += sqrt(bEuclideanDistance(allVertices2[k], allVertices2[j]));
    totalDistance -= sqrt(bEuclideanDistance(allVertices2[i], allVertices2[j]));

    double calcDistance = 0.0;
    double minDistance = 0.0;
    int minIndex = 0;

    // step 3. arbitrarily select k
    for (int p = 3; p < totalPokemons; ++p){
        minDistance = numeric_limits<double>::infinity();
        k = p;
        for (int t = 0; t < (int)tour.size(); ++t){
            calcDistance = 0;
            if (t != (int)tour.size()-1){
                i = tour[t];
                j = tour[t+1];
            }
            else {
                i = tour[t];
                j = 0;
            }

            calcDistance += sqrt(bEuclideanDistance(allVertices2[i], allVertices2[k]));
            calcDistance += sqrt(bEuclideanDistance(allVertices2[k], allVertices2[j]));
            calcDistance -= sqrt(bEuclideanDistance(allVertices2[i], allVertices2[j]));

            if (calcDistance < minDistance){
                minDistance = calcDistance;
                minIndex = t;
            }
        }
        auto it = tour.begin() + minIndex + 1;

        tour.insert(it, k);

        totalDistance += minDistance;
    }

    upperBound = totalDistance;
}

double Poke :: bEuclideanDistance(Vertices vert1, Vertices vert2){
    double bEuclideanDistance = 0;

    double xdifference = (double)(vert1.x - vert2.x);
    double ydifference = (double)(vert1.y - vert2.y);

    bEuclideanDistance = ( (xdifference) * (xdifference) ) + ( (ydifference) * (ydifference) );
    
    return bEuclideanDistance;
}

void Poke :: tspOutPut(){
    cout << totalDistance << "\n";
    for (int i = 0; i < (int)tour.size(); ++i){
        cout << tour[i] << " ";
    }
    cout << "\n";
}

//--------------------------------------------------------------------------------------------

// Part C

// set upperbound and currBestVec before calling genperms using part B
void Poke :: partB(vector<Vertices> allVertices1){
    // set upperBound
    fasttsp(allVertices1);
    
    path = tour;
    upperBound = totalDistance;
    currBestVec = tour;
    totalDistance = 0;
    
}

void Poke :: genPerms(int permLength){
    if (permLength == (int)path.size()){
    // Do something with the path: every pokemon you need to visit
        // add the currBestPath + loop from end to finish
        double startToEndDist = bEuclideanDistance(allVertices[path[0]], allVertices[path[totalPokemons - 1]]);
        currBestDistance += sqrt(startToEndDist);
        // check if that is < upperBound, and if it is < update the upperBound
        if (currBestDistance < upperBound){
            upperBound = currBestDistance;
            currBestVec = path;
        }
        // subtract the loop from end to finish
        currBestDistance -= sqrt(startToEndDist);
        return; // update() 
    }

    // promising: mst estimate/prune part
    if (!promising(path, permLength)){
        return;
    }

    for (int i = permLength; i < (int)path.size(); ++i){
        swap(path[permLength], path[i]);
        // add path length
        currBestDistance += sqrt(bEuclideanDistance(allVertices[path[permLength -1]], allVertices[path[permLength]]));
        genPerms(permLength + 1);
        // remove path length
        currBestDistance -= sqrt(bEuclideanDistance(allVertices[path[permLength -1]], allVertices[path[permLength]]));
        swap(path[permLength], path[i]);
    }
}

// black is the currLength, blue calculate the distance between start and ending vertex of curr subtour to what ever vertex not in current tour 0, permlength - 1
// gray is unvisited vertices calculate w mst.
// initialize currpath, currlength with part B called once
bool Poke :: promising(vector<int> const &path1, int permLength1){
    double lowerBoundDistance = 0;

    // edge cases:
    if ((int)path1.size() - permLength1 < 5){
        return true;
    }

    if (currBestDistance > upperBound){
        return false;
    }

    // call Prims function
    partCPrims(allVertices, permLength1);

    allPokemons.clear();
        
    // compare start pokemon and the permlength-1 pokemon at allVertices[curren
    // double tempDistance1 = 0;
    // double tempDistance2 = 0;
    double minDistance1 = numeric_limits<double>::infinity();
    double minDistance2 = numeric_limits<double>::infinity();
    int start0 = path1[0];
    int permLengthminus1 = path1[permLength1 - 1];

    // loop through path1 and compare the distances to find the minDistances.
    for (int i = permLength1; i < (int)path1.size(); ++i){
        int pathi = path1[i];
        double tempDistance1 = bEuclideanDistance(allVertices[start0], allVertices[pathi]);
        double tempDistance2 = bEuclideanDistance(allVertices[permLengthminus1], allVertices[pathi]);
        if (tempDistance1 < minDistance1){
            minDistance1 = tempDistance1;
        }
        if(tempDistance2 < minDistance2){
            minDistance2 = tempDistance2;
        }
    }

    // calculate the lowerbound distance.
    lowerBoundDistance = grayDistance + currBestDistance + sqrt(minDistance1) + sqrt(minDistance2);

    grayDistance = 0;

    if (lowerBoundDistance < upperBound){
        return true;
    }
    else {
        return false;
    }
}

void Poke :: partCPrims(vector<Vertices> allVertices1, int permLength1){
    //Pokemons pokemon;
    // store pokemon objects into vector of allPokemons
    grayDistance = 0;
    allPokemons.resize(allVertices1.size());
    int pathPermlength = path[permLength1];
    
    allPokemons[pathPermlength].distance = 0;

    for (int i = permLength1; i < (int)path.size(); ++i){
        // find and store the min distance vertex's index: A first time around the loop
        // new mindistanceVertex
        int tempIndex;
        int minDistancePokemon = -1;
        double minDist = numeric_limits<double>::infinity();
    
        for (int k = permLength1; k < (int)path.size(); ++k){
            int pathk = path[k];
            if (allPokemons[pathk].visited == false && allPokemons[pathk].distance < minDist){
                minDist = allPokemons[pathk].distance;
                minDistancePokemon = pathk;
            }
        }
        
        tempIndex = minDistancePokemon;
        
        allPokemons[tempIndex].visited = true; // set A as visited

        grayDistance += sqrt(allPokemons[tempIndex].distance);

        // MST
        for (int j = permLength1; j < (int)path.size(); ++j){ // loop through all pokemons
            int pathj = path[j];
            // auto& pt = allVertices1[tempIndex];
            // auto& pj = allVertices1[pathj];
            if (allPokemons[pathj].visited == true){
                continue;
            }
            // check the distance between A and all other unvisited pokemons and calculate the distance in between
            double newDistance = bEuclideanDistance(allVertices1[tempIndex], allVertices1[pathj]);
            
            // if the new distance is < that unvisited pokemon's distance, set that pokemons distance to the new distance
            // and set that unvisited pokemon's parent to A
            if (newDistance < allPokemons[pathj].distance){
                allPokemons[pathj].distance = newDistance;
                allPokemons[pathj].parent = tempIndex;
            }
        }
    }
}

void Poke :: optOutput (){
    cout << upperBound << "\n";
    for (int i = 0; i < (int)currBestVec.size(); ++i){
        cout << currBestVec[i] << " ";
    }
}

