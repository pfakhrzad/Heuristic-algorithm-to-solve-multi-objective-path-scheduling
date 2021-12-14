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
 * @brief  This struct stores all possible solutions for each iterator and will be updated 
 */
struct solutions
{
    uint64_t ID;                                 //the number of particle in one population
    vector<vector<uint16_t>> solution;           //matrix of possible allocation
    double cost_amount;                          //amount of minimize objective function
    double fitness_amount;                       //amount of maximize objective function
    vector<double> fulfillment_percentage;       // to keep percentage of shipment for each node
    double rank;                                 //rank based on nondominate algorithm
    vector<vector<uint16_t>> best_self_position; // the last try with best rank
    double velocity;                             //velocity of particle movement
};
uint64_t population; // population/Swarm size

//list of particles
vector<solutions> particle_list = {{1, {{1}}, 0, 0, {1}, 1, {{1}},1}};

//This is for checking segmentation fault and will be removed
array<solutions, 20> particle_list2;

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
    void particle_ranking();

    //member function pf best_global
    void best_global_initialization(vector<vector<uint16_t>> &BestGlobal, double &Best_cost, double &Best_fitness, const uint64_t &node_count, const uint64_t &truck_count);

    //member function of particle movement
    void Particle_movement(const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> stock, const vector<string> partial, const vector<vector<double>> distance)
    {
    }

private:
};

/**
 * @brief This class will compare the position of particles in one population to rank them based on two objectives then based on the best rank decide to move the particles
 * 
 */
class NonDominate : public PSO
{
public:
    void dominated_list(vector<double> &rank, vector<double> &cost_amount, vector<double> &fitness_amount, uint64_t &population)
    {
        //counting the dominate number for each particle
        for (uint64_t u = 1; u < population; u++)
        {
            for (uint64_t w = 1; w < population; w++)
            {
                if (cost_amount[u] > cost_amount[w] || fitness_amount[u] < fitness_amount[w])
                {
                    dominate_number++;
                }
            }
            dominate_num.push_back(dominate_number);
            dominate_number = 0;
        }

        //ranking each particle
        for (uint64_t x = 1; x < population; x++)
        {
            if (x == 1)
            {
                best_rank = 1;
            }
            for (uint64_t y = 1; y < population; y++)
            {
                if (dominate_num[x] > dominate_num[y])
                {
                    best_rank++;
                }
            }
            rank.push_back(best_rank);
        }
    }

private:
    //  double best_cost = 0;
    //  double best_fitness = 1;
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
    double cost_function(const vector<vector<uint16_t>> &solution, const vector<vector<double>> &, const uint64_t &, const uint64_t &);

    //member function to calculate the fitness objective
    double fitness_function(const vector<double> &stock, const vector<double> &fulfillment_percentage, const uint64_t &node_count);

private:
    vector<vector<uint64_t>> position;
    vector<vector<uint64_t>> velocity;
    vector<vector<uint64_t>> best_position;
    double best_cost;
    double best_fitness;
};

double particles::cost_function(const vector<vector<uint16_t>> &solution, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count)
{
    double cost = 0;
    uint16_t last_allocation = 0;
    for (uint16_t i = 0; i < truck_count; i++)
    {
        for (uint16_t j = 0; j < node_count; j++)
        {
            if (solution[i][j] == 1)
            {
                cost = cost + distance[last_allocation][j + 1];
                last_allocation = j + 1;
            }
        }
        last_allocation = 0;
    }
    return cost;
}

double particles::fitness_function(const vector<double> &stock, const vector<double> &fulfillment_percentage, const uint64_t &node_count)
{
    double fulfillment_den = 0;
    double fulfillment_num = 0;
    double fulfillment = 0;
    for (uint16_t i = 0; i < node_count; i++)
    {
        fulfillment_num += stock[i] * fulfillment_percentage[i];
        fulfillment_den += stock[i];
    }
    fulfillment = fulfillment_num / fulfillment_den;
    return fulfillment;
}

void PSO::particle_initialization(const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> stock, const vector<string> partial, const vector<vector<double>> distance)
{
    //control parameters for initializing
    double filled_capacity = 0;
    set<uint64_t> allocated_set;

    //calculating the population for starting
    population = ((node_count * node_count) / 2);

    //random allocation between 0 and 1
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int16_t> uid(0, 1);

    //main steps
    vector<vector<uint16_t>> first_pop(truck_count, vector<uint16_t>(node_count));
    vector<double> fulfillment_percentage_truck(node_count);
    vector<uint16_t> solution_truck(node_count);

    //for filling the struct vector
    for (uint64_t i = 1; i < population; i++)
    {
        //  particle_list[i].ID=i;
        //  particle_list[i].solution=first_pop;
        //  particle_list[i].cost_amount=0;
        //  particle_list[i].fitness_amount=1;
        //  particle_list[i].fulfillment_percentage=fulfillment_percentage_truck;
        //  particle_list[i].rank=1;
        particle_list.push_back({i, first_pop, 0, 1, fulfillment_percentage_truck, 0, first_pop});
    }

    for (uint64_t i = 1; i < population; i++)
    {
        for (uint64_t node = 0; node < node_count; node++)
        {
            fulfillment_percentage_truck.push_back(0);
        }
        particle_list[i].ID = i;

        for (uint64_t j = 0; j < truck_count; j++)
        {
            for (uint64_t k = 0; k < node_count; k++)
            {
                uint16_t cell = uid(mt);
                if (cell == 0)
                    solution_truck[k] = cell;
                else
                {
                    if (allocated_set.count(k) == 0)
                    {
                        if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "Yes")
                        {
                            solution_truck[k] = cell;
                            filled_capacity += solution_truck[k] * (truck_capacity - filled_capacity);
                            if (solution_truck[k] == 1)
                                fulfillment_percentage_truck[k] = ((truck_capacity - filled_capacity) / stock[k]);
                            allocated_set.insert(k);
                        }
                        else if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "No")
                        {
                            solution_truck[k] = 0;
                        }
                        else
                        {
                            solution_truck[k] = cell;
                            filled_capacity += solution_truck[k] * stock[k];
                            if (solution_truck[k] == 1)
                            {
                                fulfillment_percentage_truck[k] = 1;
                                allocated_set.insert(k);
                            }
                        }
                    }
                    else
                    {
                        solution_truck[k] = 0;
                    }
                }
            }
            particle_list[i].solution.push_back(solution_truck);
        }
        particle_list[i].fulfillment_percentage = (fulfillment_percentage_truck); // this is wrong
        //particle_list[i].fulfillment_percentage.push_back(fulfillment_percentage_truck)
        fulfillment_percentage_truck.clear();
        solution_truck.clear();
        allocated_set.clear();
    }

    //Here we fill the cost and fitness objectives of initial particles
    for (uint64_t k = 1; k < population; k++)
    {
        vector<vector<uint16_t>> first_particles;
        vector<uint16_t> first_particles_rows;
        vector<double> first_fulfillment_percentage;
        for (uint16_t s = 0; s < truck_count; s++)
        {
            for (uint16_t l = 0; l < node_count; l++)
            {
                first_particles_rows.push_back(particle_list[k].solution[s][l]);
                first_fulfillment_percentage.push_back(particle_list[k].fulfillment_percentage[l]);
            }
            first_particles.push_back(first_particles_rows);
        }
        particles particle;
        particle_list[k].cost_amount = particle.cost_function(first_particles, distance, node_count, truck_count);
        particle_list[k].fitness_amount = particle.fitness_function(stock, first_fulfillment_percentage, node_count);
        first_particles_rows.clear();
        first_particles.clear();
        first_fulfillment_percentage.clear();
    }
}

void PSO::particle_ranking()
{
    NonDominate first_rank;
    vector<double> rank_list;
    vector<double> cost_amount;
    vector<double> fitness_amount;
    for (uint64_t k = 1; k < population; k++)
    {
        cost_amount.push_back(particle_list[k].cost_amount);
        fitness_amount.push_back(particle_list[k].fitness_amount);
    }
    first_rank.dominated_list(rank_list, cost_amount, fitness_amount, population);
    for (uint64_t k = 1; k < population; k++)
    {
        particle_list[k].rank = rank_list[k - 1];
    }
}

void PSO::best_global_initialization(vector<vector<uint16_t>> &BestGlobal, double &Best_cost, double &Best_fitness, const uint64_t &node_count, const uint64_t &truck_count)
{
    int rank_count = 0;
    for (uint64_t i = 1; i < population; i++)
    {

        if (particle_list[i].rank == 1)
        {
            rank_count++;

            //first filling the targets for next compare
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
template <typename T>

class output
{
public:
    output(const uint64_t &Try_o, const string Algorithm_o, const T cost_function_o, const T fitness_function_o, const vector<vector<T>> &best_solution_o, int16_t &ExitFlag_o, const int64_t node_count_o, const int64_t truck_count_o)
    {
        Try = Try_o;
        Algorithm = Algorithm_o;
        cost_function = cost_function_o;
        fitness_function = fitness_function_o;
        best_solution = best_solution_o;
        ExitFlag = ExitFlag_o;
        node_count = node_count_o;
        truck_count = truck_count_o;
    }

    void print()
    {
        cout << "iteration :" << Try << '\n';
        cout << "Algorithm:" << Algorithm << '\n';
        cout << "min F(x) :" << cost_function << '\n';
        cout << "max G(x) :" << fitness_function << '\n';
        cout << "Allocations :" << '\n';

        for (uint64_t i = 0; i < truck_count; i++)
        {
            for (uint64_t j = 0; j < truck_count; j++)
            {
                cout << best_solution[i][j] << '\t';
            }
            cout << '\n';
        }
    }

private:
    uint64_t Try;
    string Algorithm;
    T cost_function;
    T fitness_function;
    vector<vector<T>> best_solution;
    int16_t ExitFlag;
    int64_t node_count;
    int64_t truck_count;
};

class AntColony
{
public:
private:
};

//main block of program
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
    read_dataset(model, node_x, node_y, stock, partial, node_count, truck_count, truck_capacity, rate, max_iteration);

    //The main contsiners that will be used for the algorithm
    vector<vector<double>> nodes;    //matrix of position of nodes
    vector<vector<double>> distance; //matrix of distance between nodes
    vector<double> dis_node;         //This shows distance vector for each node
    uint64_t c1, c2 = 2;             //personal and global learning coefficient

    //define the velocity for particle movement

    //Calculating distance matrix between nodes by Euclidean method
    for (uint16_t i = 0; i < node_count + 1; i++)
    {
        for (uint16_t j = 0; j < node_count + 1; j++)
        {
            dis_node.push_back(sqrt(((node_x[j] - node_x[i]) * (node_x[j] - node_x[i]) + (node_y[j] - node_y[i]) * (node_y[j] - node_y[i]))));
        }
        distance.push_back(dis_node);
        dis_node.clear();
    }

    //initialization the particles in a population
    PSO initial;
    initial.particle_initialization(node_count, truck_count, truck_capacity, stock, partial, distance);

    //put rank for initial particle based on non-domination algorithm
    PSO initial_rank;
    initial_rank.particle_ranking();

    //initial global best particle position
    static double Try = 0;
    vector<vector<uint16_t>> BestGlobal;
    double Best_cost;
    double Best_fitness;
    PSO initial_BestGlobal;
    initial_BestGlobal.best_global_initialization(BestGlobal, Best_cost, Best_fitness, node_count, truck_count);

    // Running PSO Algorithm
    while (Try < max_iteration)
    {
        Try = Try + 1;
        PSO run;
        run.Particle_movement(node_count, truck_count, truck_capacity, stock, partial, distance);

        //   run.objective;
        //   run.particle_ranking();
        //   run.best_global();
    }

    // Running Ant Colony Algorithm

    //Publishing the result
    string Algorithm = "non-dominated PSO";
    //  output result(Try, Algorithm, cost_function, itness_function, best_solution, ExitFlag, node_count, truck_count);

    return 0;
}

/**
 * @brief This function will read the optimization problem parameteres
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
                if ((line_number > 0) && (node_count > line_number - 1))
                {
                    stringstream extractLine(line);
                    string element;
                    uint32_t element_number = 0;

                    while (getline(extractLine, element, ','))
                    {
                        try
                        {
                            //reading main parameters
                            if (line_number == 1 && element_number == 1)
                            {
                                //istringstream nu(element);
                                node_count = stoi(element);
                            }
                            else if (line_number == 2 && element_number == 1)
                            {
                                truck_count = stoi(element);
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
                                max_iteration = stod(element);
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
                                    throw invalid_argument("you entered wrong Yes/No");
                                }
                            }
                        }
                        catch (const invalid_argument &e)
                        {
                            cout << "error 4: in line " << line_number << " and column " << element_number << " the type is not correct";
                            return -1;
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