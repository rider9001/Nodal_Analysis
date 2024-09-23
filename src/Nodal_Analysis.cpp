/// ------------------------------------------
/// @file Nodal_Analysis.cpp
///
/// @brief Source for nodal analysis functions
/// ------------------------------------------

#include "../inc/Nodal_Analysis.h"

///--------------------------------------------------------
std::vector<std::pair<std::string, double>> DCNodalAnalysis(const Nodal_Analysis_DC_t& node_info)
{
    Matrix<double> voltRes = node_info.conductance_mat.inverse() % node_info.net_currents;

    std::vector<std::pair<std::string, double>> nodeResults;
    for (size_t i = 0; i < voltRes.getRowCount(); i++)
    {
        nodeResults.push_back({node_info.node_names.at(i), voltRes.get(i,0)});
    }

    return nodeResults;
}

///--------------------------------------------------------
std::vector<std::pair<std::string, Complex_P_t>> ACNodalAnalysis(const Nodal_Analysis_AC_t& node_info)
{
    Matrix<Complex_P_t> voltRes = node_info.admittance_mat.inverse() % node_info.net_currents;

    std::vector<std::pair<std::string, Complex_P_t>> nodeResults;
    for (size_t i = 0; i < voltRes.getRowCount(); i++)
    {
        nodeResults.push_back({node_info.node_names.at(i), voltRes.get(i,0)});
    }

    return nodeResults;
}

///--------------------------------------------------------
std::vector<std::string> split(const std::string& str, const char& delim)
{
    std::stringstream ss(str);
    std::string token;
    std::vector<std::string> tokens;
    while (getline(ss, token, delim))
    {
        tokens.push_back(token);
    }
    return tokens;
}

///--------------------------------------------------------
std::vector<std::string> parseTextContent(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<std::string> fileLines;

    std::string line;
    while (std::getline(file, line))
    {
        line.erase(0, line.find_first_not_of(whitespace));

        // Skip lines commented out or blank lines
        if (line.substr(0, 2) != "//" and line.size() != 0)
        {
            fileLines.push_back(line);
        }
        else
        {
            // Pushing blank strings onto list so line numbers are correct
            fileLines.push_back("");
        }
    }

    return fileLines;
}

///--------------------------------------------------------
template<typename T>
void addAdmittance(Matrix<T>& mat, const T& admittance, const int& node1, const int& node2)
{
    if (node1 != -1)
    {
        T newVal = mat.get(node1, node1) + admittance;
        mat.set(node1, node1, newVal);

        // apply negative to col of opposite node if not ground
        if (node2 != -1)
        {
            newVal = mat.get(node1, node2) - admittance;
            mat.set(node1, node2, newVal);
        }
    }

    if (node2 != -1)
    {
        T newVal = mat.get(node2, node2) + admittance;
        mat.set(node2, node2, newVal);

        // apply negative to col of opposite node if not ground
        if (node1 != -1)
        {
            newVal = mat.get(node2, node1) - admittance;
            mat.set(node2, node1, newVal);
        }
    }
}

///--------------------------------------------------------
Nodal_Analysis_DC_t readDCAnalysisFile(const std::string& filename)
{
    std::vector<std::string> fileLines = parseTextContent(filename);

    auto it = std::find_if_not(fileLines.begin(), fileLines.end(),
                [](const std::string& x) { return x.empty();});
    if (it == fileLines.end())
    {
        throw std::invalid_argument("File has no content");
    }

    // First non-empty line should be a space-seperated list of the names of all nodes
    std::vector<std::string> node_names = split(fileLines.at(std::distance(fileLines.begin(), it)), ' ');
    if (std::find(node_names.begin(), node_names.end(), ground_node_name) != node_names.end())
    {
        throw std::invalid_argument("GND is a reserved node name and cannot be in the node list");
    }

    Nodal_Analysis_DC_t analysis{
        node_names,
        Matrix<double>(node_names.size(), node_names.size()),
        Matrix<double>(node_names.size(), 1)
        };

    // Start at the first line after the net names
    for (size_t i = std::distance(fileLines.begin(), it) + 1; i < fileLines.size(); i++)
    {
        // skip empty lines
        if (fileLines.at(i).empty())
        {
            continue;
        }

        // Each line should be in the following form:
        // [Symbol char] [component value] [Node1] [Node2]
        auto lineSplit = split(fileLines.at(i), ' ');

        if (lineSplit.size() != 4)
        {
            throw std::invalid_argument("Bad component command (line " + std::to_string(i+1) + ")");
        }

        if (lineSplit.at(0).size() != 1 or
            std::find(valid_component_symbols.begin(), valid_component_symbols.end(), lineSplit.at(0)[0]) == valid_component_symbols.end())
        {
            throw std::invalid_argument("Symbol: " + lineSplit.at(0) + " is not a valid symbol {I,V,R,L,C} (line " + std::to_string(i+1)+ ")");
        }

        char symbol = lineSplit.at(0)[0];
        // Write conversion for 10k -> 10,000 etc
        double magnitude = convertCompToValue(lineSplit.at(1));
        std::pair<std::string, std::string> nodes_connected{lineSplit.at(2), lineSplit.at(3)};

        if (symbol == 'L' or symbol == 'C')
        {
            throw std::invalid_argument("Symbol: " + std::to_string(symbol) + " is not allowed in DC analysis {I,V,R} (line " + std::to_string(i+1) + ")");
        }

        // If first node is groud on a direction agnostic component, swap nodes to make sure calculation in correct magnitude
        if (symbol == 'R' and nodes_connected.first == ground_node_name)
        {
            swap(nodes_connected.first, nodes_connected.second);
        }

        int node_idx_1, node_idx_2;
        if (nodes_connected.first == ground_node_name)
        {
            node_idx_1 = -1;
        }
        else
        {
            auto it = std::find(node_names.begin(), node_names.end(), nodes_connected.first);
            if (it == node_names.end())
            {
                throw std::invalid_argument("Node name: " + nodes_connected.first +
                " is not found in the initial node name delcaration (line " + std::to_string(i+1) + ")");
            }

            node_idx_1 = std::distance(node_names.begin(), it);
        }

        if (nodes_connected.second == ground_node_name)
        {
            node_idx_2 = -1;
        }
        else
        {
            auto it = std::find(node_names.begin(), node_names.end(), nodes_connected.second);
            if (it == node_names.end())
            {
                throw std::invalid_argument("Node name: " + nodes_connected.first +
                " is not found in the initial node name delcaration (line " + std::to_string(i+1) + ")");
            }

            node_idx_2 = std::distance(node_names.begin(), it);
        }

        switch(symbol)
        {
            case 'I':
                // set the net current values for both node columns in the net currents matrix
                // only if the node is not ground
                if (node_idx_1 != -1)
                {
                    double newVal = analysis.net_currents.get(node_idx_1, 0) + magnitude;
                    analysis.net_currents.set(node_idx_1, 0, newVal);
                }

                if (node_idx_2 != -1)
                {
                    double newVal = analysis.net_currents.get(node_idx_2, 0) - magnitude;
                    analysis.net_currents.set(node_idx_2, 0, newVal);
                }
                break;

            case 'V':
                throw std::invalid_argument("V is not implemented yet");
                break;

            case 'R':
                // 1 / magnitude is conductance
                addAdmittance<double>(analysis.conductance_mat, (1/magnitude), node_idx_1, node_idx_2);
                break;
        }
    }

    return analysis;
}

///--------------------------------------------------------
Nodal_Analysis_AC_t readACAnalysisFile(const std::string& filename)
{
    std::vector<std::string> fileLines = parseTextContent(filename);

    auto it = std::find_if_not(fileLines.begin(), fileLines.end(),
                [](const std::string& x) { return x.empty();});
    if (it == fileLines.end())
    {
        throw std::invalid_argument("File has no content");
    }

    // First non-empty line should be a space-seperated list of the names of all nodes
    std::vector<std::string> node_names = split(fileLines.at(std::distance(fileLines.begin(), it)), ' ');
    if (std::find(node_names.begin(), node_names.end(), ground_node_name) != node_names.end())
    {
        throw std::invalid_argument("GND is a reserved node name and cannot be in the node list");
    }

    double freq = -1;
    // second non empty line has the frequency
    it = std::find_if_not(++it, fileLines.end(),
                [](const std::string& x) { return x.empty();});
    if (it == fileLines.end())
    {
        throw std::invalid_argument("Frequnecy should be stated on line after netnames");
    }
    freq = convertCompToValue(fileLines.at(std::distance(fileLines.begin(), it)));

    if (freq <= 0)
    {
        throw std::invalid_argument("Freq must be greater than 0");
    }

    Nodal_Analysis_AC_t analysis{
        node_names,
        Matrix<Complex_P_t>(node_names.size(), node_names.size()),
        Matrix<Complex_P_t>(node_names.size(), 1)
        };

    // Start at the first line after the net names and freq (+2)
    for (size_t i = std::distance(fileLines.begin(), it) + 2; i < fileLines.size(); i++)
    {
        // skip empty lines
        if (fileLines.at(i).empty())
        {
            continue;
        }

        // Each line should be in the following form, phase only used on voltage/current sources:
        // [Symbol char] [component magnitude,phase] [Node1] [Node2]
        auto lineSplit = split(fileLines.at(i), ' ');

        if (lineSplit.size() != 4)
        {
            throw std::invalid_argument("Bad component command (line " + std::to_string(i+1) + ")");
        }

        if (lineSplit.at(0).size() != 1 or
            std::find(valid_component_symbols.begin(), valid_component_symbols.end(), lineSplit.at(0)[0]) == valid_component_symbols.end())
        {
            throw std::invalid_argument("Symbol: " + lineSplit.at(0) + " is not a valid symbol {I,V,R,L,C} (line " + std::to_string(i+1)+ ")");
        }

        char symbol = lineSplit.at(0)[0];
        std::pair<std::string, std::string> nodes_connected{lineSplit.at(2), lineSplit.at(3)};

        // If first node is groud on a direction agnostic component, swap nodes to make sure calculation in correct magnitude
        if (symbol == 'R' and nodes_connected.first == ground_node_name)
        {
            swap(nodes_connected.first, nodes_connected.second);
        }

        int node_idx_1, node_idx_2;
        if (nodes_connected.first == ground_node_name)
        {
            node_idx_1 = -1;
        }
        else
        {
            auto it = std::find(node_names.begin(), node_names.end(), nodes_connected.first);
            if (it == node_names.end())
            {
                throw std::invalid_argument("Node name: " + nodes_connected.first +
                " is not found in the initial node name delcaration (line " + std::to_string(i+1) + ")");
            }

            node_idx_1 = std::distance(node_names.begin(), it);
        }

        if (nodes_connected.second == ground_node_name)
        {
            node_idx_2 = -1;
        }
        else
        {
            auto it = std::find(node_names.begin(), node_names.end(), nodes_connected.second);
            if (it == node_names.end())
            {
                throw std::invalid_argument("Node name: " + nodes_connected.first +
                " is not found in the initial node name delcaration (line " + std::to_string(i+1) + ")");
            }

            node_idx_2 = std::distance(node_names.begin(), it);
        }

        if (symbol == 'I')
        {
            // set the net current values for both node columns in the net currents matrix
            // only if the node is not ground
            Complex_P_t phasor = decodePhasor(lineSplit.at(1));

            if (node_idx_1 != -1)
            {
                Complex_P_t newVal = analysis.net_currents.get(node_idx_1, 0) + phasor;
                analysis.net_currents.set(node_idx_1, 0, newVal);
            }

            if (node_idx_2 != -1)
            {
                Complex_P_t newVal = analysis.net_currents.get(node_idx_2, 0) - phasor;
                analysis.net_currents.set(node_idx_2, 0, newVal);
            }
        }
        else if (symbol == 'V')
        {
            throw std::invalid_argument("V is not implemented yet");
        }
        else if (symbol == 'R')
        {
            // 1 / magnitude is addmittance
            Complex_P_t res_admittance{1 / convertCompToValue(lineSplit.at(1)), 0};
            addAdmittance<Complex_P_t>(analysis.admittance_mat, res_admittance, node_idx_1, node_idx_2);
        }
        else if (symbol == 'C')
        {
            Complex_C_t cap_admittance{0, 2 * M_PI * freq * convertCompToValue(lineSplit.at(1))};
            addAdmittance<Complex_P_t>(analysis.admittance_mat, cartToPolar(cap_admittance), node_idx_1, node_idx_2);
        }
        else if (symbol == 'L')
        {
            Complex_C_t ind_admittance{0, 1 / (2 * M_PI * freq * convertCompToValue(lineSplit.at(1)))};
            addAdmittance<Complex_P_t>(analysis.admittance_mat, cartToPolar(ind_admittance), node_idx_1, node_idx_2);
        }
        else
        {
            // should be caught by prevoius check, but keeping this here for completeness
            throw std::invalid_argument("Unkwon symbol: " + std::to_string(symbol));
        }
    }

    return analysis;
}

///--------------------------------------------------------
Complex_P_t decodePhasor(const std::string& phasorStr)
{
    std::vector<std::string> phasorSplit = split(phasorStr, ',');

    if (phasorSplit.size() == 1)
    {
        return Complex_P_t{convertCompToValue(phasorSplit.at(0))};
    }
    else if (phasorSplit.size() == 2)
    {
        return Complex_P_t{convertCompToValue(phasorSplit.at(0)), convertCompToValue(phasorSplit.at(1))};
    }
    else
    {
        throw std::invalid_argument(phasorStr + " is not a valid phasor");
    }
}

///--------------------------------------------------------
double convertCompToValue(const std::string& comp)
{
    std::string workingComp = comp;
    // remove all whitespace
    size_t first_non_space = workingComp.find_first_not_of(whitespace);
    workingComp.erase(0, first_non_space);

    size_t last_non_space = workingComp.find_last_not_of(whitespace);
    workingComp.erase(last_non_space + 1);

    // extract multipler from string
    char last = workingComp.at(workingComp.size() - 1);

    if (isdigit(last))
    {
        return stod(workingComp);
    }

    // extract number from string
    workingComp = workingComp.substr(0, workingComp.size() - 1);

    // Check that rest of component value is completely numeric
    for (char c : workingComp)
    {
        // Decimal place dot is only allowed char
        if (!isdigit(c) and c != '.')
        {
            throw std::invalid_argument("Only one modifier may be used on a component value: " + comp);
        }
    }

    double num = stod(workingComp);

    switch(last)
    {
        case 'p':
            // Pico
            return num * pow(10, -12);

        case 'n':
            // Nano
            return num * pow(10, -9);

        case 'u':
            // Micro
            return num * pow(10, -6);

        case 'm':
            // Milli
            return num * pow(10, -3);

        case 'k':
            // Kilo
            return num * pow(10, 3);

        case 'M':
            // Mega
            return num * pow(10, 6);

        case 'G':
            // Giga
            return num * pow(10, 9);

        default:
            // unrecognized multiplier
            throw std::invalid_argument(comp + " cannot be evaluated, " + last + " is not a recognized multiplier.");
    }
}