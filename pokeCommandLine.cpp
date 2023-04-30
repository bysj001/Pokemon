// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "poke.h"

using namespace std;

void Poke :: getMode(int argc, char * argv[], Poke::Modes &mainMode) {
    bool modeSpecified = false;
    std::string modeArgument;
    
    // These are used with getopt_long()
    opterr = false; // Let us handle all error output for command line options
    int choice;
    int option_index = 0;

    option long_options[] = {
        { "help", no_argument,      nullptr, 'h' },
        { "mode", required_argument, nullptr, 'm' },
        { nullptr, 0,                 nullptr, '\0' }
    };

    // TODO: Fill in the double quotes, to match the mode and help options.
    while ((choice = getopt_long(argc, argv, "hm:", long_options, &option_index)) != -1) {
        switch (choice) {
        
          case 'h':{
            cout << "Helpful message " << "\n";
            exit(0);
            break;
          }
          case 'm':{
            modeArgument = optarg;
            

            if (modeArgument != "MST" && modeArgument != "FASTTSP" && modeArgument != "OPTTSP"){
                cerr << "Invalid mode " << "\n";
                exit(1);
            }
            else if (modeArgument == "MST"){
                mainMode.mst = true;
                modeSpecified = true;
            }
            else if (modeArgument == "FASTTSP"){
                mainMode.fasttsp = true;
                modeSpecified = true;
            }
            else if (modeArgument == "OPTTSP"){
                mainMode.opttsp = true;
                modeSpecified = true;
            }

            break;
          }
          default:{
            cerr << "Invalid command line option" << endl;
            exit(1);
          }
        } // switch
    } // while

    if (!modeSpecified) {
        cerr << "No mode specified " << endl;
        exit(1);
    } // if

} // getMode()



int main(int argc, char *argv[]) {
    // This should be in all of your projects, speeds up I/O
    ios_base::sync_with_stdio(false);

    cout << std::setprecision(2);
    cout << std::fixed;

    Poke pokeObject;
    Poke :: Modes mainMode;

    pokeObject.getMode(argc, argv, mainMode);

    pokeObject.readInputFile(mainMode.mst, mainMode.fasttsp, mainMode.opttsp);
    
    return 0;
}