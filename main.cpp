/// ------------------------------------------
/// @file main.cpp
///
/// @brief Start point for nodal analysis program
/// ------------------------------------------

#include <iostream>

#include "inc/Complex.h"
#include "inc/Matrix.h"
#include "inc/Nodal_Analysis.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Arguments: filepath" << endl;
        return EXIT_FAILURE;
    }

    std::string inpFile(argv[1]);

    Nodal_Analysis_DC_t analysis = readDCAnalysisFile(inpFile);

    cout << "Addmitance mat: " << endl << analysis.conductance_mat << endl;
    cout << "Net currents: " << endl << analysis.net_currents << endl;

    auto results = DCNodalAnalysis(analysis);

    cout << "Voltages:" << endl;
    for (auto res : results)
    {
        cout << res.first << ": " << res.second << endl;
    }

    return EXIT_SUCCESS;
}