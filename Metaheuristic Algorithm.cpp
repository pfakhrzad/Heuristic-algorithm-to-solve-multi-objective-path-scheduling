/**
 * @file Project3.cpp
 * @author "Paria" "Fakhrzad" (fakhrzap@mcmaster.ca)
 * @brief 
 * @version 0.1
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
using namespace std;

//functions prototype
double read_dataset(const string &, vector<double> &, vector<double> &, vector<double> &, vector<string> &, uint64_t &, uint64_t &, double &, double &, double &);

/**
 * @brief  This struct stores all possible solutions for each try
 */
struct solutions
{
    uint64_t ID;                                 //the number of particle in one population
    vector<vector<uint64_t>> solution;           //matrix of possible allocation
    double cost_amount;                          //amount of minimize objective function
    double fitness_amount;                       //amount of maximize objective function
    vector<double> fulfillment_percentage;       // to keep percentage of shipment for each node
    double rank;                                 //rank based on nondominated algorithm
    vector<vector<uint64_t>> best_self_position; // the last try with best rank
    vector<double> best_fulfillment_percentage;  //to keep percentage of best shipment in memory
    vector<vector<double>> velocity;             //velocity of particle movement
};

//list of particles based on population for solution with order with first initialization
vector<solutions> particle_list = {{1, {{1}}, 0, 0, {1}, 1, {{1}}, {1}, {{0}}}};
uint16_t Algorithm1 = 1; // For giving as input argument to function to know algorithm with ordering should be run
uint16_t Algorithm2 = 2; // For giving as input argument to function to know algorithm without ordering should be run

//list of particles based on population for solution without order with first initialization
vector<solutions> particle_list2 = {{1, {{1}}, 0, 0, {1}, 1, {{1}}, {1}, {{0}}}};

uint64_t population; // population/Swarm size

/**
 * @brief The PSO is main class to run the algorithm of finding the best solution for optimization model
 * 
 */

class PSO
{
public:
    //constructor to get main parameters for PSO algorithm
    PSO() {}

    //member function of particle initialization
    void particle_initialization(const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> stock, const vector<string> partial, const vector<vector<double>> distance);

    //member function of particle ranking
    void particle_ranking(uint16_t Algorithm_number);

    //member function of best_global
    void best_global_particle(uint16_t Algorithm_number, double &Try, vector<vector<uint64_t>> &BestGlobal, double &Best_cost, double &Best_fitness);

    //member function of particle movement
    void Particle_movement(uint16_t Algorithm_number, vector<vector<uint64_t>> &BestGlobal, const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> &stock, const vector<string> &partial);

    //member function of particle objective calculation
    void particle_objective(uint16_t Algorithm_number, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock);

    // member function of particle best memory position
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
    //member function to calculate the cost objective
    double cost_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<vector<double>> &, const uint64_t &, const uint64_t &);

    //member function to calculate the fitness objective
    double fitness_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<double> &stock, const uint64_t &node_count);

    //member function to set best-self position for particle
    void best_self_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance);

private:
    vector<vector<uint64_t>> position;
    vector<vector<uint64_t>> velocity;
    vector<vector<uint64_t>> best_position;
    double best_cost;
    double best_fitness;
};

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

        //calculate cost for best memory of particle
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

    // cost function in algorithm2
    else
    {
        // Calculate cost for current particle in algorithm1
        if (memory == 0)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    //check for first allocated node
                    if (particle_list2[ID].solution[j][k] == 1)
                    {
                        cost = cost + distance[0][k + 1];
                    }
                    else if (particle_list2[ID].solution[j][k] > 1)
                    {
                        for (uint64_t l = 0; l < node_count; l++)
                        {
                            //finding the previouse code
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

        //calculate cost for best memory of particle
        else
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    //check for first allocated node
                    if (particle_list2[ID].best_self_position[j][k] == 1)
                    {
                        cost = cost + distance[0][k + 1];
                    }
                    else if (particle_list2[ID].best_self_position[j][k] > 1)
                    {
                        for (uint64_t l = 0; l < node_count; l++)
                        {
                            //finding the previouse code
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

double particles::fitness_function(uint16_t Algorithm_number, uint64_t memory, uint64_t &ID, const vector<double> &stock, const uint64_t &node_count)
{
    try
    {
        double fulfillment_den = 0; //denominator
        double fulfillment_num = 0; //numerator
        double fulfillment = 0;

        //calculating fitness objective for algorithm 1 that order of nodes is important
        if (Algorithm_number == 1)
        {

            //calculate fitness for current particle
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

            //calculating fitness for best memory of particle
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
        else
        {
            //calculate fitness for current particle in algorithm2
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

            //calculating fitness for best memory of particle in algorithm2
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

void particles::best_self_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance)
{
    for (uint64_t i = 1; i < population; i++)
    {
        particles best_self_cost;
        particles best_self_fitness;

        double self_cost = best_self_cost.cost_function(Algorithm1, 1, i, distance, node_count, truck_count); // first argument shows it is memory
        double self_fitness = best_self_fitness.fitness_function(Algorithm1, 1, i, stock, node_count);

        if (particle_list[i].cost_amount < self_cost && particle_list[i].fitness_amount > self_fitness)
        {
            particle_list[i].best_self_position = particle_list[i].solution;
            particle_list[i].best_fulfillment_percentage = particle_list[i].fulfillment_percentage;
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
    //control parameters for initializing
    double filled_capacity = 0;
    set<uint64_t> allocated_set;
    uint64_t NumberOfNodes_truck = 0;

    //random allocation between 0 and 1
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int16_t> uid(0, 1);

    //main steps
    vector<vector<uint64_t>> first_pop(truck_count, vector<uint64_t>(node_count)); //just for initialization
    vector<vector<double>> first_pop2(truck_count, vector<double>(node_count));    //just for initialization
    vector<double> fulfillment_percentage_truck(node_count);

    //for filling the struct vector
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
                        if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "Yes")
                        {
                            NumberOfNodes_truck++;
                            particle_list[i].solution[j][k] = cell;
                            particle_list[i].best_self_position[j][k] = cell;
                            particle_list2[i].solution[j][k] = NumberOfNodes_truck;
                            particle_list2[i].best_self_position[j][k] = NumberOfNodes_truck;
                            if (particle_list[i].solution[j][k] == 1)
                            {
                                particle_list[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                particle_list[i].best_fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                particle_list2[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                particle_list2[i].best_fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                            }
                            filled_capacity += (double)cell * (truck_capacity - filled_capacity);

                            allocated_set.insert(k);
                        }
                        else if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "No")
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
                            if (particle_list[i].solution[j][k] == 1)
                            {
                                particle_list[i].fulfillment_percentage[k] = 1;
                                particle_list[i].best_fulfillment_percentage[k] = 1;
                                particle_list2[i].fulfillment_percentage[k] = 1;
                                particle_list2[i].best_fulfillment_percentage[k] = 1;
                                allocated_set.insert(k);
                            }
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
        }

        fulfillment_percentage_truck.clear();
        allocated_set.clear();
    }

    //Here we fill the cost and fitness objectives of initial particles
    for (uint64_t i = 1; i < population; i++)
    {
        particles particle;
        particle_list[i].cost_amount = particle.cost_function(Algorithm1, 0, i, distance, node_count, truck_count);
        particle_list[i].fitness_amount = particle.fitness_function(Algorithm1, 0, i, stock, node_count);
        particle_list2[i].cost_amount = particle.cost_function(Algorithm2, 0, i, distance, node_count, truck_count);
        particle_list2[i].fitness_amount = particle.fitness_function(Algorithm2, 0, i, stock, node_count);
    }
    /*

    //will be removed
    cout << "-------------------------------------------------------" << '\n';
    for (uint64_t i = 1; i < population; i++)
    {
        cout << "particle " << i << " solution1" << '\n';
        for (uint64_t j = 0; j < truck_count; j++)
        {
            for (uint64_t k = 0; k < node_count; k++)
            {

                cout << particle_list[i].solution[j][k] << '\t';
            }
            cout << '\n';
        }
        cout << "particle " << i << " cost and fitness" << '\n';
        cout << particle_list[i].cost_amount << '\t' << particle_list[i].fitness_amount;
        cout << '\n';
    }

    //will be removed
    cout << "-------------------------------------------------------" << '\n';
    for (uint64_t i = 1; i < population; i++)
    {
        cout << "particle " << i << " solution2" << '\n';
        for (uint64_t j = 0; j < truck_count; j++)
        {
            for (uint64_t k = 0; k < node_count; k++)
            {

                cout << particle_list2[i].solution[j][k] << '\t';
            }
            cout << '\n';
        }
        cout << "particle " << i << " cost and fitness" << '\n';
        cout << particle_list2[i].cost_amount << '\t' << particle_list2[i].fitness_amount;
        cout << '\n';
    }*/
}

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
    else
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
}

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
                        else if ((((Best_cost - particle_list[i].cost_amount) / Best_cost) > .5))
                        {
                            if ((Best_fitness - particle_list[i].fitness_amount) / Best_fitness < .2)
                            {
                                Best_cost = particle_list[i].cost_amount;
                                Best_fitness = particle_list[i].fitness_amount;
                                BestGlobal = particle_list[i].solution;
                            }
                        }
                        else
                        {
                            if (((particle_list[i].fitness_amount - Best_fitness) / Best_fitness) > .5)
                            {
                                if (((particle_list2[i].cost_amount - Best_cost) / Best_cost) < .2)
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
        else
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
                        else if (((Best_cost - particle_list2[i].cost_amount) / Best_cost) > .5)
                        {
                            if ((Best_fitness - particle_list2[i].fitness_amount) / Best_fitness < .2)
                            {
                                Best_cost = particle_list2[i].cost_amount;
                                Best_fitness = particle_list2[i].fitness_amount;
                                BestGlobal = particle_list2[i].solution;
                            }
                        }
                        else
                        {
                            if (((particle_list2[i].fitness_amount - Best_fitness) / Best_fitness) > .5)
                            {
                                if (((particle_list2[i].cost_amount - Best_cost) / Best_cost) < .2)
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

void PSO::Particle_movement(uint16_t Algorithm_number, vector<vector<uint64_t>> &BestGlobal, const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> &stock, const vector<string> &partial)
{
    // Movign particles in algorithm1 when the order of nodes is important and should be based on input dataset
    if (Algorithm_number == 1)
    {
        double max_variable = 1;
        double min_variable = 0;

        // Calculating each particle new velocity
        double max_velocity = 1 * (max_variable - min_variable);
        double min_velocity = -max_velocity;
        double velo; //just for checking that is in the range or not
        static double w = 1;
        double wdamp = 0.99;
        double c1 = 2;
        double c2 = 2; //personal and global learning coefficient
        double r1, r2; //a random parameter that is [0, 1]
        vector<vector<double>> velocity_new;
        vector<double> velocity_truck; // this is just for creating each row of matrix
        random_device random;
        mt19937 mt(random());
        uniform_real_distribution<double> uid(0, 1);

        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    r1 = uid(mt);
                    r2 = uid(mt);
                    velo = ((w * particle_list[i].velocity[j][k]) + c1 * r1 * ((double)particle_list[i].best_self_position[j][k] - (double)particle_list[i].solution[j][k]) + c2 * r2 * ((double)BestGlobal[j][k] - (double)particle_list[i].solution[j][k]));

                    if (velo > max_velocity || velo < min_velocity)
                    {
                        velo = 0;
                    }
                    velocity_truck.push_back(velo);
                }
                particle_list[i].velocity.push_back(velocity_truck);
                velocity_truck.clear();
            }
        }

        w = wdamp * w;
        
        //calculating new particle positions based on new velocity
        uint64_t new_variable;       // just for checking that new variable is between range
        set<uint64_t> allocated_set; //for checking that each node be in just one truck
        double filled_capacity = 0;  // for capacity

        for (uint64_t i = 1; i < population; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                for (uint64_t k = 0; k < node_count; k++)
                {
                    if ((double)particle_list[i].solution[j][k] + particle_list[i].velocity[j][k] < 0.5)
                    {
                        new_variable = 0;
                    }
                    else
                    {
                        new_variable = 1;
                    }
                    if (new_variable == 0)
                        particle_list[i].solution[j][k] = new_variable;
                    else
                    {
                        if (allocated_set.count(k) == 0)
                        {
                            if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "Yes")
                            {
                                particle_list[i].solution[j][k] = new_variable;
                                filled_capacity += (truck_capacity - filled_capacity);
                                particle_list[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                allocated_set.insert(k);
                            }
                            else if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "No")
                            {
                                particle_list[i].solution[j][k] = 0;
                            }
                            else
                            {
                                particle_list[i].solution[j][k] = new_variable;
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
    else
    {
    }
};

void PSO::particle_objective(uint16_t Algorithm_number, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock)
{
    particles objective;
    for (uint64_t i = 1; i < population; i++)
    {
        particle_list[i].cost_amount = objective.cost_function(Algorithm1, 1, i, distance, node_count, truck_count);
        particle_list[i].fitness_amount = objective.fitness_function(Algorithm1, 1, i, stock, node_count);
    }
}

void PSO::best_memory_particle(uint16_t Algorithm_number, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance)
{
    particles best_memory;
    best_memory.best_self_particle(Algorithm_number, node_count, truck_count, stock, distance);
}

class output
{
public:
    output(){};
    void print(const double &Try_o, const string Algorithm_o, vector<vector<uint64_t>> &BestGlobal, double Best_cost, double Best_fitness, uint64_t node_count, uint64_t truck_count)
    {
        cout << "----------------------------------------------------" << '\n';
        cout << "----------------------------------------------------" << '\n';
        cout << "Iteration : " << Try_o << '\n';
        cout << "Algorithm : " << Algorithm_o << '\n';
        cout << "Min F(x)  : " << Best_cost << '\n';
        cout << "Max G(x)  : " << Best_fitness << '\n';
        cout << "----------------------------------------------------" << '\n';
        cout << "----------------------------------------------------" << '\n';
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
    }
    uint64_t file(ostream &OutputFile, const double &Try_o, vector<vector<uint64_t>> &BestGlobal, uint64_t node_count, uint64_t truck_count)
    {
        OutputFile << "Result" << Try_o << ",";
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
};

/**
 * @brief This is the main block of program that will run all steps of algorithm
 * 
 * @return int the return 0 shows that program run successfully and return -1 will close program with printing the error description
 */
int main()
{
    // input parameters from dataset
    uint64_t node_count = 1;   //number of nodes in network
    uint64_t truck_count = 0;  //number of available trucks
    double truck_capacity;     //Capacity of a truck
    double rate;               //power consumption rate for trucks
    double max_iteration = 50; //maximum allowed iteration
    vector<double> node_x;     //x of position of nodes
    node_x.push_back(0);       //the origin x
    vector<double> node_y;     //y of position of nodes
    node_y.push_back(0);       //the origin y
    vector<double> stock;      //stocks available in each node
    vector<string> partial;    //partial shipment available in each node

    //Reading the model parameters from dataset
    string model = "Model.csv";
    try
    {
        read_dataset(model, node_x, node_y, stock, partial, node_count, truck_count, truck_capacity, rate, max_iteration);
    }
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
        return -1;
    }

    // Calculating the population for starting it is a multiple of number of nodes
    population = (node_count * node_count);

    // Print the main parameters of the model
    cout << "-------------------------------------------------------" << '\n';
    cout << "node count     :" << node_count << '\n';
    cout << "truck count    :" << truck_count << '\n';
    cout << "truck capacity :" << truck_capacity << '\n';
    cout << "rate           :" << rate << '\n';
    cout << "max_iteration  :" << max_iteration << '\n';
    cout << "population     :" << population << '\n';

    // The main contsiners that will be used for the algorithm
    vector<vector<double>> nodes;    //matrix of position of nodes
    vector<vector<double>> distance; //matrix of distance between nodes
    vector<double> dis_node;         //This shows distance vector for each node

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

    /**
     * @brief Running PSO Algorithm for Algorithm number1 that order of nodes is important and should be the same as input file
     * 
     */
    PSO run; // Define an object of PSO algorithm
    static double Try = 0;

    // Step 1- initialization the particles by random
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

    // Check for exceeding to the maximum number of iteration
    while (Try < max_iteration - 1)
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

    // Publishing the result for  algorithm1
    string OutputFile = "NSPSO_Output.csv";
    ofstream data_output;
    string Algorithm = "NS-PSO-with order";
    output result;
    data_output.open(OutputFile, ios::out);
    if (!data_output.is_open())
    {
        cout << "Error: opening output file is not successfull!";
        data_output.close();
        return -1;
    }
    result.print(Try, Algorithm, BestGlobal_Algorithm1, Best_cost_Algorithm1, Best_fitness_Algorithm1, node_count, truck_count);
    result.file(data_output, Try, BestGlobal_Algorithm1, node_count, truck_count);

    /**
 * @brief Running PSO Algorithm for second solution that order of nodes is not important
 * 
 */
    PSO run2; // Define an object of PSO algorithm
    static double Try2 = 0;
    // Step1 had been already done in algorithm1 because we want to have same initial particle to compare the algorithms so no need to do it again

    // Step 2- Put rank for initial particle based on non-domination sorting algorithm
    run.particle_ranking(Algorithm2);

    // Step3- initial global best particle position
    vector<vector<uint64_t>> BestGlobal_Algorithm2;
    double Best_cost_Algorithm2;
    double Best_fitness_Algorithm2;
    PSO initial_BestGlobal2;
    initial_BestGlobal2.best_global_particle(Algorithm2, Try2, BestGlobal_Algorithm2, Best_cost_Algorithm2, Best_fitness_Algorithm2);

    // Check for excedding to the maximum number of iteration
    while (Try2 < max_iteration)
    {
        // updating iteration number
        Try2 = Try2 + 1;

        // Step4- particles based on best velocity of self-memory and global memory will move
        run.Particle_movement(Algorithm2, BestGlobal_Algorithm2, node_count, truck_count, truck_capacity, stock, partial);

        // Step5- objective calculation
        run.particle_objective(Algorithm2, distance, node_count, truck_count, stock);

        // Step6- particle best-self position
        run.best_memory_particle(Algorithm2, node_count, truck_count, stock, distance);

        // Step7- particle ranking
        run.particle_ranking(Algorithm2);

        // Step8- run.best_global();
        run.best_global_particle(Algorithm2, Try, BestGlobal_Algorithm2, Best_cost_Algorithm2, Best_fitness_Algorithm2);
    }
    data_output.close();
    return 0;
}

/**
 * @brief This function will read the path scheduling model parameteres
 * @param model file name
 * @param node_x a control variable
 * @param node_y a control variable
 * @param stock stock in each node
 * @param partial the boolian variable for allowing partial shipment
 * @param node_count number of nodes in graph
 * @param truck_count number of trucks
 * @param truck_capacity  capacity per truck
 * @param rate  power consumption rate
 * @param max_iteration  maximum allowed iteration, the default would be 50
 * @return double 
 */
double read_dataset(const string &model, vector<double> &node_x, vector<double> &node_y, vector<double> &stock, vector<string> &partial, uint64_t &node_count,
                    uint64_t &truck_count,
                    double &truck_capacity,
                    double &rate,
                    double &max_iteration)
{
    ifstream dataset(model, ios::in);
    if (dataset.is_open())
    {
        string line;              //reading each line of file
        uint64_t line_number = 0; //counting the line number
        while (!dataset.eof())    //to check end of file
        {
            while (getline(dataset, line))
            {
                //cout<< line<<'\n'; will be removed
                if ((line_number > 0) && (node_count > line_number - 1))
                {
                    stringstream extractLine(line);
                    string element;
                    uint32_t element_number = 0;

                    while (getline(extractLine, element, ','))
                    {
                        //  cout<< element<<'\n'; will be removed
                        try
                        {
                            //reading main parameters
                            if (line_number == 1 && element_number == 1)
                            {
                                //istringstream nu(element);
                                (node_count = stol(element));
                            }
                            else if (line_number == 2 && element_number == 1)
                            {
                                truck_count = stol(element);
                            }
                            else if (line_number == 3 && element_number == 1)
                            {
                                truck_capacity = stod(element);
                            }
                            else if (line_number == 4 && element_number == 1)
                            {
                                rate = stod(element);
                            }
                            else if (line_number == 5 && element_number == 1)
                            {
                                max_iteration = stol(element);
                            }

                            //reading x,y, stock and partial allowed
                            else if (element_number == 6)
                            {
                                node_x.push_back(stod(element));
                            }
                            else if (element_number == 7)
                            {
                                node_y.push_back(stod(element));
                            }
                            else if (element_number == 8)
                            {
                                stock.push_back(stod(element));
                            }
                            else if (element_number == 9)
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