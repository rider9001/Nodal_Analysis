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
/// @brief Parses a text file into a vector, blanking any comment lines
///
/// @param filename name of file to parse
///
/// @return parsed file contents
std::vector<std::string> parseTextContent(const std::string& filename);

///--------------------------------------------------------
/// @brief Reads a DC analysis file and compiles components into a conductance/net current matrix
///
/// @param filename local path of file to read
///
/// @return Compiled DC nodal analysis data
Nodal_Analysis_DC_t readDCAnalysisFile(const std::string& filename);

///--------------------------------------------------------
/// @brief Reads a DC analysis file and compiles components into a conductance/net current matrix
///
/// @param filename local path of file to read
///
/// @return Compiled AC nodal analysis data
Nodal_Analysis_AC_t readACAnalysisFile(const std::string& filename);

///--------------------------------------------------------
/// @brief Decodes a phasor from a string in the form [mag],[phase]
///
/// @param phasorStr string containing phasor in form
///
/// @return complex phasor
Complex_P_t decodePhasor(const std::string& phasorStr);

///--------------------------------------------------------
/// @brief Adds a given admittance to the admittance matrix given in mat
///
/// @tparam T type of admittance (pure real, complex)
///
/// @param mat matrix to add admittance to
/// @param admittance admittance to add
/// @param node1 node 1 of the connected component, -1 indicates ground
/// @param node2 node 2 of the connected component
template <typename T>
void addAdmittance(const Matrix<T>& mat, const T& admittance, const int& node1, const int& node2);

///--------------------------------------------------------
/// @brief Converts a component value string into a double value
/// e.g: 20k -> 20,000, 10m -> 0.001
///
/// @param comp string to convert to value
///
/// @return resoved value of the component string
double convertCompToValue(const std::string& comp);