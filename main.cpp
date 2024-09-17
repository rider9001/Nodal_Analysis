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
    if (argc != 3)
    {
        cout << "Arguments: [type A/D] [filepath]" << endl;
        return EXIT_FAILURE;
    }

    std::string inpFile(argv[2]);

    std::string anaylsis_type(argv[1]);

    if (anaylsis_type == "A")
    {
        Nodal_Analysis_AC_t analysis = readACAnalysisFile(inpFile);

        cout << "Addmitance mat: " << endl << analysis.admittance_mat << endl;
        cout << "Net currents: " << endl << analysis.net_currents << endl;

        auto results = ACNodalAnalysis(analysis);

        cout << "Voltages:" << endl;
        for (auto res : results)
        {
            cout << res.first << ": " << res.second << endl;
        }
    }
    else if (anaylsis_type == "D")
    {
        Nodal_Analysis_DC_t analysis = readDCAnalysisFile(inpFile);

        cout << "Addmitance mat: " << endl << analysis.conductance_mat << endl;
        cout << "Net currents: " << endl << analysis.net_currents << endl;

        auto results = DCNodalAnalysis(analysis);

        cout << "Voltages:" << endl;
        for (auto res : results)
        {
            cout << res.first << ": " << res.second << endl;
        }
    }
    else
    {
        cout << "Unknown analysis type: " + anaylsis_type << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}