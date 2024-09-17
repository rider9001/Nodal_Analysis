/// ------------------------------------------
/// @file Nodal_Analysis.h
///
/// @brief Header for nodal analysis functions
/// ------------------------------------------
#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include "Matrix.h"
#include "Complex.h"

/// @brief All whitespace chars for comparing
const std::string whitespace(" \r\n\t\v\f");

/// @brief GND is a reserved node name
const std::string ground_node_name("GND");

/// @brief Key:
/// I: current source
/// V: voltage source
/// R: resistor
/// C: capacitor
/// L: inductor
const std::vector<char> valid_component_symbols({'I','V','R','L','C'});

/// @brief Stores the needed matricies and net names required for a DC analysis
struct Nodal_Analysis_DC_t
{
    /// @brief Names of the node names used in analysis
    /// Order of the names corrisponds to both row on the conductance matrix
    /// and net current on the net_currents list.
    std::vector<std::string> node_names;

    /// @brief (n,n) Matrix of conductances between nodes
    Matrix<double> conductance_mat;

    /// @brief (n, 1) Matrix of net currents on each node
    Matrix<double> net_currents;
};

/// @brief Stores the needed matricies and net names required for a AC analysis
struct Nodal_Analysis_AC_t
{
    /// @brief Names of the node names used in analysis
    /// Order of the names corrisponds to both row on the conductance matrix
    /// and net current on the net_currents list.
    std::vector<std::string> node_names;

    /// @brief (n,n) Matrix of admittances between nodes
    Matrix<Complex_P_t> admittance_mat;

    /// @brief (n, 1) Matrix of net current phasors on each node
    Matrix<Complex_P_t> net_currents;
};


///--------------------------------------------------------
/// @brief Uses conductance matrix and net currents to calculate
/// the voltage at all nodes
///
/// @param node_info conductance and current matricies and net names
///
/// @return list of pairs of net names and calculated voltages
std::vector<std::pair<std::string, double>> DCNodalAnalysis(const Nodal_Analysis_DC_t& node_info);

///--------------------------------------------------------
/// @brief Uses the admittance matrix and net currents to calculate voltages for all nodes
///
/// @param node_info admittance and current matricies and net names
///
/// @return List of pairs of node names and voltage phasors
std::vector<std::pair<std::string, Complex_P_t>> ACNodalAnalysis(const Nodal_Analysis_AC_t& node_info);

///--------------------------------------------------------
/// @brief Splits string into vector using single char delimiter
///
/// @param str string to split
/// @param delim delimiter character
///
/// @return vector of split strings
std::vector<std::string> split(const std::string& str, const char& delim);

///--------------------------------------------------------
/// @brief Reads a DC analysis file and compiles components into a conductance/net current matrix
///
/// @param filename local path of file to read
///
/// @return Compiled DC nodal analysis data
Nodal_Analysis_DC_t readDCAnalysisFile(const std::string& filename);