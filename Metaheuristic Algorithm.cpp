/**
 * @file Metaheuristic ALgorithm.cpp
 * @author Paria Fakhrzad (fakhrzap@mcmaster.ca)
 * @brief 
 * @version 2
 * @date 2021-12-09
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <exception>
#include <sstream>
#include <cmath>
#include <random>
#include <set>
#include <array>
#include <algorithm>
#include <chrono>
using namespace std;

/**
 * @brief  This struct stores all possible solutions for each try
 */
struct solutions
{
    uint64_t ID;                                 // The number of particle in one population
    vector<vector<uint64_t>> solution;           // Matrix of possible allocation
    double cost_amount;                          // Amount of minimize objective function
    double fitness_amount;                       // Amount of maximize objective function
    vector<double> fulfillment_percentage;       // To keep percentage of shipment for each node
    double rank;                                 // Rank based on nondominated algorithm
    vector<vector<uint64_t>> best_self_position; // The last try with best rank
    vector<double> best_fulfillment_percentage;  // To keep percentage of best shipment in memory
    vector<vector<double>> velocity;             // Velocity of particle movement
};

// List of particles based on population for solution with order with first initialization
vector<solutions> particle_list = {{1, {{1}}, 0, 0, {1}, 1, {{1}}, {1}, {{0}}}};
// List of particles based on population for solution without order with first initialization
vector<solutions> particle_list2 = {{1, {{1}}, 0, 0, {1}, 1, {{1}}, {1}, {{0}}}};

uint64_t population;     // Population/Swarm size
uint16_t Algorithm1 = 1; // For giving as input argument to function to know algorithm with ordering should be run
uint16_t Algorithm2 = 2; // For giving as input argument to function to know algorithm without ordering should be run

/**
 * @brief This class  will read the path scheduling model parameteres
 * @param model File name
 * @param node_x A control variable
 * @param node_y A control variable
 * @param stock Stock in each node
 * @param partial The boolian variable for allowing partial shipment
 * @param node_count Number of nodes in graph
 * @param truck_count Number of trucks
 * @param truck_capacity  Capacity per truck
 * @param rate  Power consumption rate
 * @param max_iteration  Maximum allowed iteration, the default would be 50
 * @return Double -1 if the reading is not successful, 0 if was successful
 */
class input
{
public:
    input() {}
    double read_dataset(const string &model, vector<double> &node_x, vector<double> &node_y, vector<double> &stock, vector<string> &partial, uint64_t &node_count,
                        uint64_t &truck_count,
                        double &truck_capacity,
                        double &rate,
                        double &max_iteration)
    {
        ifstream dataset(model, ios::in);
        if (dataset.is_open())
        {
            while (!dataset.eof() && exit == 0) // To check end of file
            {
                while (getline(dataset, line))
                {
                    if ((line_number > 0) && exit == 0)
                    {
                        stringstream extractLine(line);
                        string element;
                        uint32_t element_number = 0;

                        while (getline(extractLine, element, ','))
                        {
                            if (exit == 0)
                            {
                                try
                                {
                                    // Number of nodes
                                    if (line_number == 1 && element_number == 1)
                                    {
                                        (node_count = stol(element));
                                    }
                                    // Number of trucks
                                    else if (line_number == 2 && element_number == 1)
                                    {
                                        truck_count = stol(element);
                                    }
                                    // Truck capacity
                                    else if (line_number == 3 && element_number == 1)
                                    {
                                        truck_capacity = stod(element);
                                    }
                                    // Truck power consumption per 1km
                                    else if (line_number == 4 && element_number == 1)
                                    {
                                        rate = stod(element);
                                    }
                                    // Maximum iteration
                                    else if (line_number == 5 && element_number == 1)
                                    {
                                        max_iteration = stod(element);
                                        if (max_iteration == 0)
                                        {
                                            max_iteration = 50; // If user doesn't enter the max iteration it will be 50
                                        }
                                    }

                                    // Reading (x,y) and stock and partial allowed
                                    // x
                                    else if (element_number == 6)
                                    {
                                        if (node_x.size() - 1 < node_count)
                                        {
                                            node_x.push_back(stod(element));
                                        }

                                        else if (node_x.size() - 1 >= node_count)
                                        {
                                            if (line_number >= 5)
                                            {
                                                exit = 1;
                                            }
                                        }
                                    }
                                    // y
                                    else if (element_number == 7)
                                    {
                                        if (node_y.size() - 1 < node_count)
                                        {
                                            node_y.push_back(stod(element));
                                        }

                                        else if (node_y.size() - 1 >= node_count)
                                        {
                                            if (line_number >= 5)
                                            {
                                                exit = 1;
                                            }
                                        }
                                    }

                                    // Stock
                                    else if (element_number == 8)
                                    {
                                        if (stock.size() < node_count)
                                        {
                                            stock.push_back(stod(element));
                                        }

                                        else if (stock.size() >= node_count)
                                        {
                                            if (line_number >= 5)
                                            {
                                                exit = 1;
                                            }
                                        }
                                    }

                                    // Partial pick up
                                    else if (element_number == 9)
                                    {
                                        if (partial.size() < node_count)
                                        {
                                            if (element == "No" || element == "Yes")
                                            {
                                                partial.push_back(element);
                                            }
                                            else
                                            {
                                                throw invalid_argument("You entered wrong Yes/No");
                                            }
                                        }

                                        else if (partial.size() >= node_count)
                                        {
                                            if (line_number >= 5)
                                            {
                                                exit = 1;
                                            }
                                        }
                                    }
                                    if (node_count == 0 || truck_count == 0 || truck_capacity == 0 || rate == 0)
                                    {
                                        throw invalid_argument("The algorithm cannot work with 0 amount for parameters");
                                    }
                                }
                                catch (const invalid_argument &e)
                                {
                                    if (element_number == 6 && ((line_number - node_count) == 0))
                                    {
                                        cout << "Error 3: Number of (x,y) positions should be the same as number of nodes " << '\n';
                                        throw invalid_argument("Update dataset!");
                                    }
                                    else if (element == "")
                                    {
                                        cout << "Error 2: The cell in line " << line_number + 1 << " and column " << element_number + 1 << " is mandatory" << '\n';
                                        throw invalid_argument("Update dataset!");
                                    }
                                    else
                                    {
                                        cout << e.what() << '\n';
                                        cout << "Error 4: in line " << line_number + 1 << " and column " << element_number + 1 << " the type is not correct" << '\n';
                                        throw invalid_argument("Update dataset!");
                                    }
                                }
                                element_number++;
                            }
                        }
                    }
                    line_number++;
                }
            }
        }
        else
        {
            cout << "The input file opened with error" << '\n';
            dataset.close();
            return -1;
        }
        dataset.close();
        return 0;
    }

private:
    string line;              // Reading each line of file
    uint64_t line_number = 0; // Counting the line number
    uint32_t exit = 0;        // Condition of exiting the file
};

//--------------
// PSO Algorithm
//--------------

/**
 * @brief The PSO is main class to run the algorithm of finding the best solution for path scheduling optimization model
 * 
 */

class PSO
{
public:
    // Constructor to get main parameters for PSO algorithm
    PSO() {}

    // Member function of particle initialization
    void particle_initialization(const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> stock, const vector<string> partial, const vector<vector<double>> distance);

    // Member function of particle ranking
    void particle_ranking(uint16_t Algorithm_number);

    // Member function of best_global
    void best_global_particle(uint16_t Algorithm_number, double &Try, vector<vector<uint64_t>> &BestGlobal, double &Best_cost, double &Best_fitness);

    // Member function of particle movement
    void Particle_movement(uint16_t Algorithm_number, vector<vector<uint64_t>> &BestGlobal, const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> &stock, const vector<string> &partial);

    // Member function of particle objective calculation
    void particle_objective(uint16_t Algorithm_number, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock);

    // Member function of particle best memory position
    void best_memory_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance);

private:
};

/**
 * @brief This class will compare the position of particles in one population to rank them based on two objectives then based on the best rank decide to move the particles
 * 
 */
class NonDominated : public PSO
{
public:
    void dominated_list(vector<double> &rank, vector<double> &cost_amount, vector<double> &fitness_amount, uint64_t &population)
    {
        // Counting the dominant number for each particle
        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t w = 1; w < population; w++)
            {
                if (cost_amount[i - 1] > cost_amount[w - 1] || fitness_amount[i - 1] < fitness_amount[w - 1])
                {
                    dominate_number++;
                }
            }
            dominate_num.push_back(dominate_number);
            dominate_number = 0;
        }

        // Ranking each particle
        for (uint64_t i = 1; i < population; i++)
        {
            if (i == 1)
            {
                best_rank = 1;
            }
            for (uint64_t y = 1; y < population; y++)

            {
                if (dominate_num[i - 1] > dominate_num[y - 1])
                {
                    best_rank++;
                }
            }
            rank.push_back(best_rank);
            best_rank = 1;
        }
    }

private:
    double best_rank = 0;
    double dominate_number = 0;
    vector<double> dominate_num;
};

/**
 * @brief This class will calculate all attributes of a particle in a population 
 * 
 */
class particles : public PSO
{
public:
    particles() {}
    // Member function to calculate the cost objective
    double cost_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<vector<double>> &, const uint64_t &, const uint64_t &);

    // Member function to calculate the fitness objective
    double fitness_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<double> &stock, const uint64_t &node_count);

    // Member function to set best-self position for particle
    void best_self_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance);

private:
    vector<vector<uint64_t>> position;
    vector<vector<uint64_t>> velocity;
    vector<vector<uint64_t>> best_position;
    double best_cost;
    double best_fitness;
};

/**
 * @brief  This FUnction will calculate the cost objective F(x) based on the solution 
 * 
 * @param Algorithm_number If it is 1 so function will calculate based on algorithm1, and if 2 the formula is related to algorithm2
 * @param memory If the value is 1 it mean that objective needs to be calculated based on memory of particle
 * @param ID  Index of particle in a population
 * @param distance  Distance matrix for all nodes
 * @param node_count  Number of nodes in model
 * @param truck_count Number of truck in model
 * @return double  The amount of cost objective
 */
double particles::cost_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count)
{
    double cost = 0;
    // Cost function based on algorithm1
    if (Algorithm_number == 1)
    {

        // Calculate cost for current particle in algorithm1
        uint64_t last_allocation = 0;
        if (memory == 0)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    if (particle_list[ID].solution[j][k] == 1)
                    {
                        cost = cost + distance[last_allocation][k + 1];
                        last_allocation = k + 1;
                    }
                }
                last_allocation = 0;
            }
            return cost;
        }

        // Calculate cost for best memory of particle
        else
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    if (particle_list[ID].best_self_position[j][k] == 1)
                    {
                        cost = cost + distance[last_allocation][k + 1];
                        last_allocation = k + 1;
                    }
                }
                last_allocation = 0;
            }
            return cost;
        }
    }

    // Cost function in algorithm2
    else
    {
        // Calculate cost for current particle in algorithm1
        if (memory == 0)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    // Check for first allocated node
                    if (particle_list2[ID].solution[j][k] == 1)
                    {
                        cost = cost + distance[0][k + 1];
                    }
                    else if (particle_list2[ID].solution[j][k] > 1)
                    {
                        for (uint64_t l = 0; l < node_count; l++)
                        {
                            // Finding the previouse code
                            if (particle_list2[ID].solution[j][k] - particle_list2[ID].solution[j][l] == 1)
                            {
                                cost = cost + distance[k + 1][l + 1];
                                break;
                            }
                        }
                    }
                }
            }
            return cost;
        }

        // Calculate cost for best memory of particle
        else
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    // Check for first allocated node
                    if (particle_list2[ID].best_self_position[j][k] == 1)
                    {
                        cost = cost + distance[0][k + 1];
                    }
                    else if (particle_list2[ID].best_self_position[j][k] > 1)
                    {
                        for (uint64_t l = 0; l < node_count; l++)
                        {
                            // Finding the previouse code
                            if (particle_list2[ID].best_self_position[j][k] - particle_list2[ID].best_self_position[j][l] == 1)
                            {
                                cost = cost + distance[k + 1][l + 1];
                                break;
                            }
                        }
                    }
                }
            }
            return cost;
        }
    }
}
/**
 * @brief This function will calculate fitness objective G(x) for particles
 * 
 * @param Algorithm_number If it is 1 so function will calculate based on algorithm1, and if 2 the formula is related to algorithm2
 * @param memory If the value is 1 it mean that objective needs to be calculated based on memory of particle
 * @param ID Index of particle in a population
 * @param stock Available stock in each node
 * @param node_count Number of nodes in model
 * @return double 
 */
double particles::fitness_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<double> &stock, const uint64_t &node_count)
{
    try
    {
        double fulfillment_den = 0; // Denominator
        double fulfillment_num = 0; // Numerator
        double fulfillment = 0;     // The calculated value of objective

        // Calculating fitness objective for algorithm 1 that order of nodes is important
        if (Algorithm_number == 1)
        {

            // Calculate fitness for current particle
            if (memory == 0)
            {
                for (uint64_t i = 0; i < node_count; i++)
                {
                    fulfillment_num += stock[i] * particle_list[ID].fulfillment_percentage[i];
                    fulfillment_den += stock[i];
                }
                fulfillment = fulfillment_num / fulfillment_den;
                return fulfillment;
            }

            // Calculating fitness for best memory of particle
            else
            {
                for (uint64_t i = 0; i < node_count; i++)
                {
                    fulfillment_num += stock[i] * particle_list[ID].best_fulfillment_percentage[i];
                    fulfillment_den += stock[i];
                }
                fulfillment = fulfillment_num / fulfillment_den;
                return fulfillment;
            }
        }
        // Calculating fitness objective based on algorithm2 that is without ordering
        else if(Algorithm_number==2)
        {
            // Calculate fitness for current particle in algorithm2
            if (memory == 0)
            {
                for (uint64_t i = 0; i < node_count; i++)
                {
                    fulfillment_num += stock[i] * particle_list2[ID].fulfillment_percentage[i];
                    fulfillment_den += stock[i];
                }
                fulfillment = fulfillment_num / fulfillment_den;
                return fulfillment;
            }

            // Calculating fitness for best memory of particle in algorithm2
            else
            {
                for (uint64_t i = 0; i < node_count; i++)
                {
                    fulfillment_num += stock[i] * particle_list2[ID].best_fulfillment_percentage[i];
                    fulfillment_den += stock[i];
                }
                fulfillment = fulfillment_num / fulfillment_den;
                return fulfillment;
            }
        }
    }
    catch (invalid_argument &e)
    {
        cout << "Exitflag-2: NAN value";
    }
}

/**
 * @brief This function will update the best memory of one particle in each iteration
 * 
 * @param Algorithm_number If it is 1 so function will calculate based on algorithm1, and if 2 the formula is related to algorithm2
 * @param node_count  Number of nodes in model
 * @param truck_count  Number of trucks in model
 * @param stock   Stock in each node for pick up
 * @param distance Distance matrix for nodes in a graph
 */
void particles::best_self_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance)
{
    particles best_self_cost;
    particles best_self_fitness;

    // Updating best self memory for each particle in algorithm1
    if (Algorithm_number == 1)
    {
        for (uint64_t i = 1; i < population; i++)
        {

            double self_cost = best_self_cost.cost_function(Algorithm_number, 1, i, distance, node_count, truck_count); // Second argument shows it is memory
            double self_fitness = best_self_fitness.fitness_function(Algorithm_number, 1, i, stock, node_count);

            if (particle_list[i].cost_amount < self_cost && particle_list[i].fitness_amount > self_fitness)
            {
                particle_list[i].best_self_position = particle_list[i].solution;
                particle_list[i].best_fulfillment_percentage = particle_list[i].fulfillment_percentage;
            }
        }
    }

    // Updating best self memory for each particle in algorithm2
    else
    {
        for (uint64_t i = 1; i < population; i++)
        {

            double self_cost = best_self_cost.cost_function(Algorithm_number, 1, i, distance, node_count, truck_count); // Second argument shows it is a memory
            double self_fitness = best_self_fitness.fitness_function(Algorithm_number, 1, i, stock, node_count);        // Second argument shows it is a memory

            if (particle_list2[i].cost_amount < self_cost && particle_list2[i].fitness_amount > self_fitness)
            {
                particle_list2[i].best_self_position = particle_list2[i].solution;
                particle_list2[i].best_fulfillment_percentage = particle_list2[i].fulfillment_percentage;
            }
        }
    }
}

/**
 * @brief  This function will randomly initialize the particles in solution space and calculate all section of a solution
 * 
 * @param node_count  This is number of nodes that has been entered in input file
 * @param truck_count  This is number of available truck daily that has been entered in input file
 * @param truck_capacity Capacity of each truck
 * @param stock  Available stock in each node for picking
 * @param partial This is a boolian variable that shows partial shipment is allowed or not
 * @param distance This is a matrix that keeps the distance between two nodes
 */
void PSO::particle_initialization(const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> stock, const vector<string> partial, const vector<vector<double>> distance)
{
    // Control parameters for initializing
    double filled_capacity = 0;
    set<uint64_t> allocated_set;
    uint64_t NumberOfNodes_truck = 0;

    // Random allocation between 0 and 1
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int16_t> uid(0, 1);

    // Main steps
    vector<vector<uint64_t>> first_pop(truck_count, vector<uint64_t>(node_count)); // Just for initialization
    vector<vector<double>> first_pop2(truck_count, vector<double>(node_count));    // Just for initialization
    vector<double> fulfillment_percentage_truck(node_count);

    // For filling the struct vector
    for (uint64_t i = 1; i < population; i++)
    {
        particle_list.push_back({i, first_pop, 0, 1, fulfillment_percentage_truck, 0, first_pop, fulfillment_percentage_truck, first_pop2});
        particle_list2.push_back({i, first_pop, 0, 1, fulfillment_percentage_truck, 0, first_pop, fulfillment_percentage_truck, first_pop2});
    }

    for (uint64_t i = 1; i < population; i++)
    {
        particle_list[i].ID = i;
        particle_list2[i].ID = i;
        for (uint64_t j = 0; j < truck_count; j++)
        {
            NumberOfNodes_truck = 0;
            for (uint64_t k = 0; k < node_count; k++)
            {
                uint64_t cell = uid(mt);
                if (cell == 0)
                {
                    particle_list[i].solution[j][k] = cell;
                    particle_list[i].best_self_position[j][k] = cell;
                    particle_list2[i].solution[j][k] = cell;
                    particle_list2[i].best_self_position[j][k] = cell;
                }
                else
                {
                    if (allocated_set.count(k) == 0)
                    {
                        if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "Yes"))
                        {
                            NumberOfNodes_truck++;
                            particle_list[i].solution[j][k] = cell;
                            particle_list[i].best_self_position[j][k] = cell;
                            particle_list2[i].solution[j][k] = NumberOfNodes_truck;
                            particle_list2[i].best_self_position[j][k] = NumberOfNodes_truck;
                            particle_list[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                            particle_list[i].best_fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                            filled_capacity += (double)cell * (truck_capacity - filled_capacity);
                            allocated_set.insert(k);
                        }
                        else if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "No"))
                        {
                            particle_list[i].solution[j][k] = 0;
                            particle_list[i].best_self_position[j][k] = 0;
                            particle_list2[i].solution[j][k] = 0;
                            particle_list2[i].best_self_position[j][k] = 0;
                        }
                        else if ((truck_capacity - filled_capacity) == 0)
                        {
                            particle_list[i].solution[j][k] = 0;
                            particle_list[i].best_self_position[j][k] = 0;
                            particle_list2[i].solution[j][k] = 0;
                            particle_list2[i].best_self_position[j][k] = 0;
                        }
                        else
                        {
                            NumberOfNodes_truck++;
                            particle_list[i].solution[j][k] = cell;
                            particle_list[i].best_self_position[j][k] = cell;
                            particle_list2[i].solution[j][k] = NumberOfNodes_truck;
                            particle_list2[i].best_self_position[j][k] = NumberOfNodes_truck;
                            filled_capacity += (double)particle_list[i].solution[j][k] * stock[k];
                            particle_list[i].fulfillment_percentage[k] = 1;
                            particle_list[i].best_fulfillment_percentage[k] = 1;
                            allocated_set.insert(k);
                        }
                    }
                    else
                    {
                        particle_list[i].solution[j][k] = 0;
                        particle_list[i].best_self_position[j][k] = 0;
                        particle_list2[i].solution[j][k] = 0;
                        particle_list2[i].best_self_position[j][k] = 0;
                    }
                }
            }
            filled_capacity = 0;

            // based on shortest distance assign nodes for each truck when the order of nodes is not important
            vector<double> min_distance;        // A vector to keep the distances
            set<uint64_t> allocated_node_truck; // Allocated nodes number to each truck
            uint64_t NumberNode = 0; // For just preventing duplicated node number in a path
            allocated_node_truck.clear();
            for (uint64_t k = 0; k < node_count; k++)
            {
                if (particle_list2[i].solution[j][k] > 0)
                {
                  //  cout << distance[0][k + 1] << '\n';
                    vector<double>::const_iterator i = find(min_distance.begin(), min_distance.end(), distance[0][k + 1]);
                    if (i == min_distance.end())
                    {
                        min_distance.push_back(distance[0][k + 1]);
                        NumberNode++;
                    }
                }
            }
            sort(min_distance.begin(), min_distance.end());
            NumberOfNodes_truck=NumberNode;
            NumberNode=0;
            for (uint64_t k2 = 0; k2 < NumberOfNodes_truck; k2++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    if (particle_list[i].solution[j][k] > 0)
                    {
                     //   cout << particle_list2[i].solution[j][k] << '\t';
                        if (min_distance[k2] == distance[0][k + 1])
                        {
                            if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "Yes"))
                            {
                                NumberNode++;
                                particle_list2[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                particle_list2[i].best_fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                filled_capacity += truck_capacity - filled_capacity;
                                particle_list2[i].solution[j][k] = NumberNode;
                                particle_list2[i].best_self_position[j][k] = NumberNode;
                                allocated_node_truck.insert(NumberNode);
                             //   cout << particle_list2[i].solution[j][k] << '\t';
                            }
                            else if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "No"))
                            {
                                particle_list2[i].solution[j][k] = 0;
                                particle_list2[i].best_self_position[j][k] = 0;
                            }
                            else if ((truck_capacity - filled_capacity) == 0)
                            {
                                particle_list2[i].solution[j][k] = 0;
                                particle_list2[i].best_self_position[j][k] = 0;
                            }
                            else
                            {
                                NumberNode++;
                                particle_list2[i].solution[j][k] = NumberNode;
                                particle_list2[i].best_self_position[j][k] = NumberNode;
                                filled_capacity += stock[k];
                                particle_list2[i].fulfillment_percentage[k] = 1;
                                particle_list2[i].best_fulfillment_percentage[k] = 1;
                                allocated_node_truck.insert(NumberNode);
                              //  cout << particle_list2[i].solution[j][k] << '\t';
                            }
                        }
                    }
                }
            }
            // This part will check if all path included the correct node counts
            for (uint64_t node = node_count; node > 0; node--)
            {
                if (allocated_node_truck.count(node) > 0)
                {
                    for (uint64_t node_inverse = 1; node_inverse < node; node_inverse++)
                    {
                        // Check the in this path the previouse node has been exists
                        if ((node - node_inverse != 0) && (allocated_node_truck.count(node - node_inverse) == 0))
                        {
                            uint64_t update_node = node_inverse;
                            while (update_node >= 1)
                            {
                                for (uint64_t k = 0; k < node_count; k++)
                                {
                                    // cout << particle_list2[i].solution[j][k] << '\t';
                                    if (particle_list2[i].solution[j][k] == node - update_node + 1)
                                    {
                                        particle_list2[i].solution[j][k] = node - update_node;
                                        particle_list2[i].best_self_position[j][k] = node - update_node;
                                        allocated_node_truck.erase(node - update_node + 1);
                                        allocated_node_truck.insert(node - update_node);
                                        //   cout << particle_list2[i].solution[j][k] << '\t';
                                        break;
                                    }
                                }
                                update_node--;
                            }
                        }
                    }
                }
            }
            filled_capacity = 0;
        }
        allocated_set.clear();
    }

    // Here we fill the cost and fitness objectives of initial particles
    for (uint64_t i = 1; i < population; i++)
    {
        particles particle;
        particle_list[i].cost_amount = particle.cost_function(Algorithm1, 0, i, distance, node_count, truck_count);
        particle_list[i].fitness_amount = particle.fitness_function(Algorithm1, 0, i, stock, node_count);
        particle_list2[i].cost_amount = particle.cost_function(Algorithm2, 0, i, distance, node_count, truck_count);
        particle_list2[i].fitness_amount = particle.fitness_function(Algorithm2, 0, i, stock, node_count);
    }
}

/**
 * @brief This function will specify the rank for each particle based on population
 * 
 * @param Algorithm_number If it is 1 so function will calculate based on algorithm1, and if 2 the formula is related to algorithm2
 */
void PSO::particle_ranking(uint16_t Algorithm_number)
{
    NonDominated first_rank;
    vector<double> rank_list;
    vector<double> cost_amount;
    vector<double> fitness_amount;

    // Rank particles based on algorithm1
    if (Algorithm_number == 1)
    {
        for (uint64_t i = 1; i < population; i++)
        {
            cost_amount.push_back(particle_list[i].cost_amount);
            fitness_amount.push_back(particle_list[i].fitness_amount);
        }

        first_rank.dominated_list(rank_list, cost_amount, fitness_amount, population);
        for (uint64_t i = 1; i < population; i++)
        {
            particle_list[i].rank = rank_list[i - 1];
        }
    }

    // Rank particles based on algorithm2
    else if (Algorithm_number == 2)
    {
        for (uint64_t i = 1; i < population; i++)
        {
            cost_amount.push_back(particle_list2[i].cost_amount);
            fitness_amount.push_back(particle_list2[i].fitness_amount);
        }

        first_rank.dominated_list(rank_list, cost_amount, fitness_amount, population);
        for (uint64_t i = 1; i < population; i++)
        {
            particle_list2[i].rank = rank_list[i - 1];
        }
    }
    else
    {
        throw invalid_argument("Error: Yhe algorithm number should be specified");
    }
}
/**
 * @brief This function will specify the best global particle in a population and update it in each iteration
 * 
 * @param Algorithm_number If it is 1 so function will calculate based on algorithm1, and if 2 the formula is related to algorithm2
 * @param Try  Number of iteration
 * @param BestGlobal A martix that keeps the particle with best ranking in population
 * @param Best_cost The amount of cost objective 
 * @param Best_fitness  The amount of fitness objective
 */

void PSO::best_global_particle(uint16_t Algorithm_number, double &Try, vector<vector<uint64_t>> &BestGlobal, double &Best_cost, double &Best_fitness)
{
    try
    {
        // Calculating the best global particle for algorithm1
        if (Algorithm_number == 1)
        {
            // If this is the first time for filling best global and no need to comparison with memory
            if (Try == 0)
            {
                uint64_t rank_count = 0;
                for (uint64_t i = 1; i < population; i++)
                {
                    if (particle_list[i].rank == 1)
                    {
                        rank_count++;

                        // First filling the targets for next compare
                        if (rank_count == 1)
                        {
                            Best_cost = particle_list[i].cost_amount;
                            Best_fitness = particle_list[i].fitness_amount;
                            BestGlobal = particle_list[i].solution;
                        }
                        else
                        {
                            if (particle_list[i].cost_amount < Best_cost && particle_list[i].fitness_amount > Best_fitness)
                            {
                                Best_cost = particle_list[i].cost_amount;
                                Best_fitness = particle_list[i].fitness_amount;
                                BestGlobal = particle_list[i].solution;
                            }
                        }
                    }
                }
            }
            // IF try is more than 0 and need to compare with memory of best global
            else
            {
                for (uint64_t i = 1; i < population; i++)
                {
                    if (particle_list[i].rank == 1)
                    {
                        // Comparing to old best_global
                        if (particle_list[i].cost_amount < Best_cost && particle_list[i].fitness_amount > Best_fitness)
                        {
                            Best_cost = particle_list[i].cost_amount;
                            Best_fitness = particle_list[i].fitness_amount;
                            BestGlobal = particle_list[i].solution;
                        }
                        else if ((((Best_cost - particle_list[i].cost_amount) / Best_cost) > .4))
                        {
                            if ((Best_fitness - particle_list[i].fitness_amount) / Best_fitness < .3)
                            {
                                Best_cost = particle_list[i].cost_amount;
                                Best_fitness = particle_list[i].fitness_amount;
                                BestGlobal = particle_list[i].solution;
                            }
                        }
                        else
                        {
                            if (((particle_list[i].fitness_amount - Best_fitness) / Best_fitness) > .4)
                            {
                                if (((particle_list2[i].cost_amount - Best_cost) / Best_cost) < .3)
                                {
                                    Best_cost = particle_list[i].cost_amount;
                                    Best_fitness = particle_list[i].fitness_amount;
                                    BestGlobal = particle_list[i].solution;
                                }
                            }
                        }
                    }
                }
            }
        }
        // Calculating the best global particle for algorithm 2 that the object is particle_list2
        else if (Algorithm_number == 2)
        {
            // Calculating the best global for first time and without needing to compare with the memory
            if (Try == 0)
            {
                uint64_t rank_count = 0;
                for (uint64_t i = 1; i < population; i++)
                {
                    if (particle_list2[i].rank == 1)
                    {
                        rank_count++;

                        // First filling the targets for next compare
                        if (rank_count == 1)
                        {
                            Best_cost = particle_list2[i].cost_amount;
                            Best_fitness = particle_list2[i].fitness_amount;
                            BestGlobal = particle_list2[i].solution;
                        }
                        else
                        {
                            if (particle_list2[i].cost_amount < Best_cost && particle_list2[i].fitness_amount > Best_fitness)
                            {
                                Best_cost = particle_list2[i].cost_amount;
                                Best_fitness = particle_list2[i].fitness_amount;
                                BestGlobal = particle_list2[i].solution;
                            }
                        }
                    }
                }
            }

            // If the try is bigger than 0 and best global needs to be compared with memory
            else
            {
                for (uint64_t i = 1; i < population; i++)
                {
                    if (particle_list2[i].rank == 1)
                    {
                        // Comparing to memory of best global
                        if (particle_list2[i].cost_amount < Best_cost && particle_list2[i].fitness_amount > Best_fitness)
                        {
                            Best_cost = particle_list2[i].cost_amount;
                            Best_fitness = particle_list2[i].fitness_amount;
                            BestGlobal = particle_list2[i].solution;
                        }
                        else if (((Best_cost - particle_list2[i].cost_amount) / Best_cost) > .4)
                        {
                            if ((Best_fitness - particle_list2[i].fitness_amount) / Best_fitness < .3)
                            {
                                Best_cost = particle_list2[i].cost_amount;
                                Best_fitness = particle_list2[i].fitness_amount;
                                BestGlobal = particle_list2[i].solution;
                            }
                        }
                        else
                        {
                            if (((particle_list2[i].fitness_amount - Best_fitness) / Best_fitness) > .4)
                            {
                                if (((particle_list2[i].cost_amount - Best_cost) / Best_cost) < .3)
                                {
                                    Best_cost = particle_list2[i].cost_amount;
                                    Best_fitness = particle_list2[i].fitness_amount;
                                    BestGlobal = particle_list2[i].solution;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch (invalid_argument &e)
    {
        cout << "Error: calculation error with inf number";
        throw invalid_argument("Update dataset!");
    }
};

/**
 * @brief This function will calculate the velocity for each particle based on itself memory and public memory then move particle to the optimal solution
 * 
 * @param Algorithm_number If it is 1 so function will calculate based on algorithm1, and if 2 the formula is related to algorithm2
 * @param BestGlobal A matrix that keeps the best global particle
 * @param node_count Number of nodes in model
 * @param truck_count Number of trucks in model
 * @param truck_capacity Capacity of each truck 
 * @param stock Available stock  for each node to pick up
 * @param partial A boolean variable that check if stock of the node can be picked partially or not
 */
void PSO::Particle_movement(uint16_t Algorithm_number, vector<vector<uint64_t>> &BestGlobal, const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> &stock, const vector<string> &partial)
{
    uint64_t new_element;          // New calculated element in solution matrix
    set<uint64_t> allocated_set;   // A set of assigned nodes in each solution, for checking that each node be in just one truck
    double filled_capacity = 0;    // Filled capacity of each truck per steps
    double max_element;            // Maximum amount that a solution element can have
    double min_element;            // Minimum amount that a solution element can have
    double max_velocity;           // Maximum movement that each particle can have
    double min_velocity;           // Minimum movement that each particle can have
    double velocity_range;         // Just for checking that is in the range or not
    static double w1 = 1;          // The weight of old velocity of algorithm1
    static double w2 = 1;          // The weight of old velocity of algorithm2
    double wdamp = 0.99;           // The coefficient of changing velocity for each iteration
    double c1 = 2;                 // Personal learning coefficient
    double c2 = 2;                 // Global learning coefficient
    double r1, r2;                 // A random parameter that is [0, 1]
    vector<double> velocity_truck; // This is just for creating each row of matrix
    random_device random;          // For getting random value of {0,1}
    mt19937 mt(random());
    uniform_real_distribution<double> uid(0.01, 0.99);

    // Movign particles in algorithm1 when the order of nodes is important and should be based on input dataset
    if (Algorithm_number == 1)
    {
        max_element = 1;
        min_element = 0;
        max_velocity = max_element - min_element;
        min_velocity = -max_velocity;

        // Calculating each particle new velocity
        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    r1 = uid(mt);
                    r2 = uid(mt);
                    velocity_range = ((w1 * particle_list[i].velocity[j][k]) + c1 * r1 * ((double)particle_list[i].best_self_position[j][k] - (double)particle_list[i].solution[j][k]) + c2 * r2 * ((double)BestGlobal[j][k] - (double)particle_list[i].solution[j][k]));
                    if (velocity_range > max_velocity)
                    {
                        velocity_range = max_element;
                    }
                    if (velocity_range < min_velocity)
                    {
                        velocity_range = -max_element;
                    }
                    particle_list[i].velocity[j][k] = velocity_range;
                }
            }
        }
        w1 = wdamp * w1;

        // Calculating new particle positions based on new velocity
        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                filled_capacity = 0;
                for (uint64_t k = 0; k < node_count; k++)
                {
                    if ((double)particle_list[i].solution[j][k] + particle_list[i].velocity[j][k] < 0.2)
                    {
                        new_element = 0;
                    }
                    else
                    {
                        new_element = 1;
                    }
                    if (new_element == 0)
                    {
                        particle_list[i].solution[j][k] = new_element;
                    }
                    else
                    {
                        if (allocated_set.count(k) == 0)
                        {
                            if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "Yes"))
                            {
                                particle_list[i].solution[j][k] = new_element;
                                filled_capacity += (truck_capacity - filled_capacity);
                                particle_list[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                allocated_set.insert(k);
                            }
                            else if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "No"))
                            {
                                particle_list[i].solution[j][k] = 0;
                            }
                            else if ((truck_capacity - filled_capacity) == 0)
                            {
                                particle_list[i].solution[j][k] = 0;
                            }

                            else
                            {
                                particle_list[i].solution[j][k] = new_element;
                                filled_capacity += stock[k];
                                particle_list[i].fulfillment_percentage[k] = 1;
                                allocated_set.insert(k);
                            }
                        }
                        else
                        {
                            particle_list[i].solution[j][k] = 0;
                        }
                    }
                }
            }
            allocated_set.clear();
        }
    }
    // Moving particles based on algorithm2 when the order of nodes is not important
    else if (Algorithm_number == 2)
    {
        max_element = (double)node_count;
        min_element = 0;
        max_velocity = max_element - min_element;
        min_velocity = -max_velocity;

        // Calculating each particle new velocity based on algorithm2
        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    r1 = uid(mt);
                    r2 = uid(mt);
                    velocity_range = ((w2 * particle_list2[i].velocity[j][k]) + c1 * r1 * ((double)particle_list2[i].best_self_position[j][k] - (double)particle_list2[i].solution[j][k]) + c2 * r2 * ((double)BestGlobal[j][k] - (double)particle_list2[i].solution[j][k]));
                    if (velocity_range > max_velocity)
                    {
                        velocity_range = max_element;
                    }
                    if (velocity_range < min_velocity)
                    {
                        velocity_range = -min_element;
                    }
                    particle_list2[i].velocity[j][k] = velocity_range;
                }
            }
        }
        w2 = wdamp * w2;

        // Calculating new particle positions based on new velocity
        uint64_t NumberOfNodes_truck;       // Number of nodes in each truck, update in each iteration
        set<uint64_t> allocated_node_truck; // Allocated nodes number to each truck

        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                filled_capacity = 0;
                NumberOfNodes_truck = 0;
                for (uint64_t k = 0; k < node_count; k++)
                {
                    // For calculating new element of each solution matrix
                    if ((double)particle_list2[i].solution[j][k] + particle_list2[i].velocity[j][k] > max_element)
                    {
                        new_element = (uint64_t)max_element;
                    }
                    else if ((double)particle_list2[i].solution[j][k] + particle_list2[i].velocity[j][k] < min_element)
                    {
                        new_element = 0;
                    }
                    else
                    {
                        // To round new element
                        new_element = particle_list2[i].solution[j][k] + (uint64_t)particle_list2[i].velocity[j][k];
                        for (double round = min_element; round <= max_element; round++)
                        {
                            double check = (double)particle_list2[i].solution[j][k] + particle_list2[i].velocity[j][k];
                            if (((0 < check - round) && (check - round < 0.5)) || ((-0.5 < check - round) && (check - round < 0)))
                            {
                                new_element = (uint64_t)round;
                                break;
                            }
                        }
                    }
                    if (new_element == 0)
                    {
                        particle_list2[i].solution[j][k] = new_element;
                    }
                    else
                    {
                        if (allocated_set.count(k) == 0) // Check if the node has already been assigned in this particle or not
                        {
                            for (uint64_t node = 0; node < node_count; node++)
                            {
                                if ((new_element - node != 0) && (allocated_node_truck.count(new_element - node) == 0)) // Check if this path number has already been assigned to truck
                                {
                                    new_element = new_element - node;
                                    if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "Yes"))
                                    {
                                        NumberOfNodes_truck++;

                                        particle_list2[i].solution[j][k] = new_element;
                                        filled_capacity += (truck_capacity - filled_capacity);
                                        particle_list2[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                        allocated_set.insert(k);
                                        allocated_node_truck.insert(new_element);
                                    }
                                    else if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "No"))
                                    {
                                        particle_list2[i].solution[j][k] = 0;
                                    }

                                    else if ((truck_capacity - filled_capacity) == 0)
                                    {
                                        particle_list2[i].solution[j][k] = 0;
                                    }

                                    else
                                    {
                                        NumberOfNodes_truck++;
                                        particle_list2[i].solution[j][k] = new_element;
                                        filled_capacity += stock[k];
                                        particle_list2[i].fulfillment_percentage[k] = 1;
                                        allocated_set.insert(k);
                                        allocated_node_truck.insert(new_element);
                                    }
                                    break;
                                }
                                else if (allocated_node_truck.count(new_element + node + 1) == 0) // Check if this path number has already been assigned to truck
                                {
                                    new_element = new_element + node + 1;
                                    if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "Yes"))
                                    {
                                        NumberOfNodes_truck++;

                                        particle_list2[i].solution[j][k] = new_element;
                                        filled_capacity += (truck_capacity - filled_capacity);
                                        particle_list2[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                        allocated_set.insert(k);
                                        allocated_node_truck.insert(new_element);
                                    }
                                    else if (((truck_capacity - filled_capacity) > 0) && ((truck_capacity - filled_capacity) < stock[k]) && (partial[k] == "No"))
                                    {
                                        particle_list2[i].solution[j][k] = 0;
                                    }

                                    else if ((truck_capacity - filled_capacity) == 0)
                                    {
                                        particle_list2[i].solution[j][k] = 0;
                                    }

                                    else
                                    {
                                        NumberOfNodes_truck++;
                                        particle_list2[i].solution[j][k] = new_element;
                                        filled_capacity += stock[k];
                                        particle_list2[i].fulfillment_percentage[k] = 1;
                                        allocated_set.insert(k);
                                        allocated_node_truck.insert(new_element);
                                    }
                                    break;
                                }
                            }
                        }
                        else
                        {
                            particle_list2[i].solution[j][k] = 0;
                        }
                    }
                }
                // Check the capacity of truck
                if ((truck_capacity - filled_capacity) > 0)
                {
                    for (uint64_t node_index = 0; node_index < node_count; node_index++)
                    {
                        if ((allocated_set.count(node_index) == 0) && ((stock[node_index]) < (truck_capacity - filled_capacity)))
                        {
                            particle_list2[i].solution[j][node_index] = NumberOfNodes_truck + 1;
                            filled_capacity += stock[node_index];
                            particle_list2[i].fulfillment_percentage[node_index] = 1;
                            allocated_set.insert(node_index);
                            allocated_node_truck.insert(NumberOfNodes_truck + 1);
                            NumberOfNodes_truck++;
                        }
                    }
                }
                // This part will check if all path included the correct node counts
                for (uint64_t node = node_count; node > 0; node--)
                {
                    if (allocated_node_truck.count(node) > 0)
                    {
                        for (uint64_t node_inverse = 1; node_inverse < node; node_inverse++)
                        {
                            // Check the in this path the previouse node has been exists
                            if ((node - node_inverse != 0) && (allocated_node_truck.count(node - node_inverse) == 0))
                            {
                                uint64_t update_node = node_inverse;
                                while (update_node >= 1)
                                {
                                    for (uint64_t k = 0; k < node_count; k++)
                                    {
                                        if (particle_list2[i].solution[j][k] == node - update_node + 1)
                                        {
                                            particle_list2[i].solution[j][k] = node - update_node;
                                            allocated_node_truck.erase(node - update_node + 1);
                                            allocated_node_truck.insert(node - update_node);
                                            break;
                                        }
                                    }
                                    update_node--;
                                }
                            }
                        }
                    }
                }
                allocated_node_truck.clear();
            }
            allocated_set.clear();
        }
    }
};

/**
 * @brief This is just a member function in PSO that call the main functions in another class
 * 
 */
void PSO::particle_objective(uint16_t Algorithm_number, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock)
{
    particles objective;

    // To update objective of algorithm1 new particles
    if (Algorithm_number == 1)
    {
        for (uint64_t i = 1; i < population; i++)
        {
            particle_list[i].cost_amount = objective.cost_function(Algorithm_number, 1, i, distance, node_count, truck_count);
            particle_list[i].fitness_amount = objective.fitness_function(Algorithm_number, 1, i, stock, node_count);
        }
    }

    // To update objective of algorithm2 new particles
    else
    {
        for (uint64_t i = 1; i < population; i++)
        {
            particle_list2[i].cost_amount = objective.cost_function(Algorithm_number, 1, i, distance, node_count, truck_count);
            particle_list2[i].fitness_amount = objective.fitness_function(Algorithm_number, 1, i, stock, node_count);
        }
    }
}

void PSO::best_memory_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance)
{
    particles best_memory;
    best_memory.best_self_particle(Algorithm_number, node_count, truck_count, stock, distance);
}
/**
 * @brief This class will print the result of algorithms and save the solution in a file
 * 
 */
class output
{
public:
    output(){};
    void print(const double &Try_o, const string Algorithm_o, vector<vector<uint64_t>> &BestGlobal, double Best_cost, double Best_fitness, uint64_t node_count, uint64_t truck_count, double rate)
    {
        Algorithm_name = Algorithm_o;

        cout << "Iteration : " << Try_o << '\n';
        cout << "Algorithm : " << Algorithm_o << '\n';
        cout << "Min F(x)  : " << Best_cost * rate << '\n';
        cout << "Max G(x)  : " << Best_fitness << '\n';
        cout << '\n';
        cout << "Best vehicle allocations :" << '\n';

        for (uint64_t j = 0; j < truck_count; j++)
        {
            cout << "Truck " << j + 1 << " :" << '\n';
            cout << '\t';
            for (uint64_t k = 0; k < node_count; k++)
            {
                cout << BestGlobal[j][k] << '\t';
            }
            cout << '\n';
        }
        cout << '\n';
    }
    uint64_t file(ostream &OutputFile, const double &Try_o, vector<vector<uint64_t>> &BestGlobal, uint64_t node_count, uint64_t truck_count)
    {
        OutputFile << Algorithm_name << " Try:" << Try_o << ",";
        for (uint64_t k = 0; k < node_count; k++)
        {
            OutputFile << "node" << k + 1 << ",";
        }
        OutputFile << endl;
        for (uint64_t j = 0; j < truck_count; j++)
        {
            OutputFile << "Truck: " << j + 1 << ",";
            for (uint64_t k = 0; k < node_count; k++)
            {
                OutputFile << BestGlobal[j][k] << ",";
            }
            OutputFile << endl;
        }
        return 0;
    }

private:
    string Algorithm_name;
};

/**
 * @brief This class is used the same as Lecture note CSE701- link in document
 * 
 */
class timer
{
public:
    void start()
    {
        start_time = chrono::steady_clock ::now();
    }

    void end()
    {
        end_time = chrono::steady_clock ::now();
        duration = end_time - start_time;
        cout << duration.count() << '\n';
        cout << "-------------------------------------------------------" << '\n';
    }

private:
    // Parameters for checking the time taken by functions
    chrono::time_point<chrono::steady_clock> start_time;
    chrono::time_point<chrono::steady_clock> end_time;
    chrono::duration<double> duration;
};

/**
 * @brief This is the main block of program that will run all steps of algorithm
 * 
 * @return The return 0 shows that program run successfully and return -1 will close program with printing the error description
 */
int main()
{
    // Input parameters from dataset
    uint64_t node_count;    // Number of nodes in network
    uint64_t truck_count;   // Number of available trucks
    double truck_capacity;  // Capacity of a truck
    double rate;            // Power consumption rate for trucks
    double max_iteration;   // Maximum allowed iteration
    vector<double> node_x;  // x_axis of position of nodes
    node_x.push_back(0);    // The origin x
    vector<double> node_y;  // y_axis of position of nodes
    node_y.push_back(0);    // The origin y
    vector<double> stock;   // Stocks available in each node
    vector<string> partial; // Partial shipment available in each node
    timer duration;         // TO calculate the take taken by each function

    // Reading the model parameters from dataset
    string model = "Model.csv";
    try
    {
        cout << "-------------------------------------------------------" << '\n';
        cout << " Reading File....    elapsed_time: ";
        duration.start();
        input input_reading; // An object of class reading dataset
        input_reading.read_dataset(model, node_x, node_y, stock, partial, node_count, truck_count, truck_capacity, rate, max_iteration);
        duration.end();
    }
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
        return -1;
    }

    // Calculating the population for starting it is a multiple of number of nodes

    //population = (node_count * node_count);
    population = 2;
    if (population < 2)
    {
        cout << "Error: number of population is not" << '\n';
        return -1;
    }

    // Print the main parameters of the model
    cout << "node count     : " << node_count << '\n';
    cout << "truck count    : " << truck_count << '\n';
    cout << "truck capacity : " << truck_capacity << '\n';
    cout << "rate           : " << rate << '\n';
    cout << "max_iteration  : " << max_iteration << '\n';
    cout << "population     : " << population << '\n';

    // The main contsiners that will be used for the algorithm
    vector<vector<double>> nodes;    // Matrix of position of nodes
    vector<vector<double>> distance; // Matrix of distance between nodes
    vector<double> dis_node;         // This shows distance vector for each node

    // Calculating distance matrix between nodes by Euclidean method
    for (uint64_t i = 0; i < node_count + 1; i++)
    {
        for (uint64_t j = 0; j < node_count + 1; j++)
        {
            dis_node.push_back(sqrt(((node_x[j] - node_x[i]) * (node_x[j] - node_x[i]) + (node_y[j] - node_y[i]) * (node_y[j] - node_y[i]))));
        }
        distance.push_back(dis_node);
        dis_node.clear();
    }

    //------------------
    // Run the Algorithm1
    //------------------
    /**
     * @brief Running PSO Algorithm for Algorithm number1 that order of nodes is important and should be the same as input file
     * 
     */
    try
    {
        PSO run; // Define an object of PSO algorithm
        static double Try = 0;

        // Step 1- initialization the particles by random
        cout << "-------------------------------------------------------" << '\n';
        cout << "Population initialization...    elapsed_time: ";
        duration.start();
        run.particle_initialization(node_count, truck_count, truck_capacity, stock, partial, distance);

        // Velocity initialization-first we assume that velocity is zero for all particles
        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    particle_list[i].velocity[j][k] = 0;
                    particle_list2[i].velocity[j][k] = 0;
                }
            }
        }

        // Step 2- put rank for initial particle based on non-domination sorting algorithm
        run.particle_ranking(Algorithm1);

        // Step 3- initial  best global particle
        vector<vector<uint64_t>> BestGlobal_Algorithm1;
        double Best_cost_Algorithm1;
        double Best_fitness_Algorithm1;
        PSO initial_BestGlobal;

        initial_BestGlobal.best_global_particle(Algorithm1, Try, BestGlobal_Algorithm1, Best_cost_Algorithm1, Best_fitness_Algorithm1);
        if (Best_cost_Algorithm1 == 0)
        {
            cout << "----------------------------------" << '\n';
            cout << "Error : There is issue in calculating the distance, check the nodes (x,y)";
            return (-1);
        }
        duration.end();

        // Check for exceeding to the maximum number of iteration
        cout << '\n'
             << '\n'
             << "Running PSO algorithm1...        elapsed_time: ";
        duration.start();
        while (Try < max_iteration)
        {
            // Updating iteration number
            Try = Try + 1;

            // Step4- particles based on best velocity of self-memory and global memory will move
            run.Particle_movement(Algorithm1, BestGlobal_Algorithm1, node_count, truck_count, truck_capacity, stock, partial);

            // Step5- objective calculation
            run.particle_objective(Algorithm1, distance, node_count, truck_count, stock);

            // Step6- particle best-self position
            run.best_memory_particle(Algorithm1, node_count, truck_count, stock, distance);

            // Step7- particle ranking
            run.particle_ranking(Algorithm1);

            // Step8- run.best_global();
            run.best_global_particle(Algorithm1, Try, BestGlobal_Algorithm1, Best_cost_Algorithm1, Best_fitness_Algorithm1);
            if (Best_cost_Algorithm1 == 0)
            {
                cout << "----------------------------------" << '\n';
                cout << " Error : There is issue in calculating the distance, check the nodes (x,y)";
                return (-1);
            }
        }
        duration.end();

        // Publishing the result for  algorithm1
        string OutputFile = "NSPSO_Output.csv";
        ofstream data_output;
        string Algorithm1_name = "NS-PSO-with order";
        output result;
        data_output.open(OutputFile, ios::out);
        if (!data_output.is_open())
        {
            cout << "Error: opening output file is not successfull!";
            data_output.close();
            return -1;
        }
        result.print(Try, Algorithm1_name, BestGlobal_Algorithm1, Best_cost_Algorithm1, Best_fitness_Algorithm1, node_count, truck_count, rate);
        result.file(data_output, Try, BestGlobal_Algorithm1, node_count, truck_count);

        //------------------
        // Run the Algorithm2
        //------------------

        /**
         * @brief Running PSO Algorithm for second solution that order of nodes is not important
         * 
         */
        PSO run2; // Define an object of PSO algorithm
        static double Try2 = 0;
        // Step1 had been already done in algorithm1 because we want to have same initial particle to compare the algorithms so no need to do it again

        // Step 2- Put rank for initial particle based on non-domination sorting algorithm
        run2.particle_ranking(Algorithm2);

        // Step3- initial global best particle position
        vector<vector<uint64_t>> BestGlobal_Algorithm2;
        double Best_cost_Algorithm2;
        double Best_fitness_Algorithm2;
        PSO initial_BestGlobal2;
        initial_BestGlobal2.best_global_particle(Algorithm2, Try2, BestGlobal_Algorithm2, Best_cost_Algorithm2, Best_fitness_Algorithm2);

        // Check for excedding to the maximum number of iteration
        cout << "-------------------------------------------------------" << '\n';
        cout << "Running PSO algorithm2...        elapsed_time: ";
        duration.start();
        while (Try2 < max_iteration)
        {
            // Updating iteration number
            Try2 = Try2 + 1;

            // Step4- particles based on best velocity of self-memory and global memory will move
            run2.Particle_movement(Algorithm2, BestGlobal_Algorithm2, node_count, truck_count, truck_capacity, stock, partial);

            // Step5- objective calculation
            run2.particle_objective(Algorithm2, distance, node_count, truck_count, stock);

            // Step6- particle best-self position
            run2.best_memory_particle(Algorithm2, node_count, truck_count, stock, distance);

            // Step7- particle ranking
            run2.particle_ranking(Algorithm2);

            // Step8- run.best_global();
            run2.best_global_particle(Algorithm2, Try2, BestGlobal_Algorithm2, Best_cost_Algorithm2, Best_fitness_Algorithm2);
        }
        duration.end();
        // Publishing the result for  algorithm2
        string Algorithm2_name = "NS-PSO-without order";
        if (!data_output.is_open())
        {
            cout << "Error: opening output file is not successfull!";
            data_output.close();
            return -1;
        }
        result.print(Try2, Algorithm2_name, BestGlobal_Algorithm2, Best_cost_Algorithm2, Best_fitness_Algorithm2, node_count, truck_count, rate);
        result.file(data_output, Try2, BestGlobal_Algorithm2, node_count, truck_count);
        data_output.close();
    }
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
        return -1;
    }

    return 0;
}
