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
    vector<vector<uint64_t>> solution;           //matrix of possible allocation
    double cost_amount;                          //amount of minimize objective function
    double fitness_amount;                       //amount of maximize objective function
    vector<double> fulfillment_percentage;       // to keep percentage of shipment for each node
    double rank;                                 //rank based on nondominated algorithm
    vector<vector<uint64_t>> best_self_position; // the last try with best rank
    vector<double> best_fulfillment_percentage;  //to keep percentage of best shipment in memory
    vector<vector<double>> velocity;             //velocity of particle movement
};
uint64_t population; // population/Swarm size

//list of particles
vector<solutions> particle_list = {{1, {{1}}, 0, 0, {1}, 1, {{1}}, {1}, {{0}}}};

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

    //member function of best_global
    void best_global_particle(double &Try, vector<vector<uint64_t>> &BestGlobal, double &Best_cost, double &Best_fitness);

    //member function of particle movement
    void Particle_movement(vector<vector<uint64_t>> &BestGlobal, const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> &stock, const vector<string> &partial);

    //member function of particle objective calculation
    void particle_objective(const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock);

    // member function of particle best memory position
    void best_memory_particle(const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance);

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
        //counting the dominant number for each particle
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

        //ranking each particle
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
    double cost_function(uint64_t memory, uint64_t &ID, const vector<vector<double>> &, const uint64_t &, const uint64_t &);

    //member function to calculate the fitness objective
    double fitness_function(uint64_t memory, uint64_t &ID, const vector<double> &stock, const uint64_t &node_count);

    //member function to set best-self position for particle
    void best_self_particle(const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance);

private:
    vector<vector<uint64_t>> position;
    vector<vector<uint64_t>> velocity;
    vector<vector<uint64_t>> best_position;
    double best_cost;
    double best_fitness;
};

double particles::cost_function(uint64_t memory, uint64_t &ID, const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count)
{
    double cost = 0;
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

double particles::fitness_function(uint64_t memory, uint64_t &ID, const vector<double> &stock, const uint64_t &node_count)
{
    try
    {
        double fulfillment_den = 0; //denominator
        double fulfillment_num = 0; //numerator
        double fulfillment = 0;
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
        else //needs to be modified for best self position
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
    catch (invalid_argument &e)
    {
        cout << "Exitflag-2: NAN value";
    }
}

void particles::best_self_particle(const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance)
{
    for (uint64_t i = 1; i < population; i++)
    {
        particles best_self_cost;
        particles best_self_fitness;

        double self_cost = best_self_cost.cost_function(1, i, distance, node_count, truck_count); // first argument shows it is memory
        double self_fitness = best_self_fitness.fitness_function(1, i, stock, node_count);

        if (particle_list[i].cost_amount < self_cost && particle_list[i].fitness_amount > self_fitness)
        {
            particle_list[i].best_self_position = particle_list[i].solution;
            particle_list[i].best_fulfillment_percentage = particle_list[i].fulfillment_percentage;
        }
    }
}

void PSO::particle_initialization(const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> stock, const vector<string> partial, const vector<vector<double>> distance)
{
    //control parameters for initializing
    double filled_capacity = 0;
    set<uint64_t> allocated_set;

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
    }

    for (uint64_t i = 1; i < population; i++)
    {
        particle_list[i].ID = i;
        for (uint64_t j = 0; j < truck_count; j++)
        {
            for (uint64_t k = 0; k < node_count; k++)
            {
                uint64_t cell = uid(mt);
                if (cell == 0)
                {
                    particle_list[i].solution[j][k] = cell;
                    particle_list[i].best_self_position[j][k] = cell;
                }
                else
                {
                    if (allocated_set.count(k) == 0)
                    {
                        if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "Yes")
                        {
                            particle_list[i].solution[j][k] = cell;
                            particle_list[i].best_self_position[j][k] = cell;
                            if (particle_list[i].solution[j][k] == 1)
                            {
                                particle_list[i].fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                                particle_list[i].best_fulfillment_percentage[k] = ((truck_capacity - filled_capacity) / stock[k]);
                            }
                            filled_capacity += (double)cell * (truck_capacity - filled_capacity);

                            allocated_set.insert(k);
                        }
                        else if ((truck_capacity - filled_capacity) < stock[k] && partial[k] == "No")
                        {
                            particle_list[i].solution[j][k] = 0;
                            particle_list[i].best_self_position[j][k] = 0;
                        }
                        else
                        {
                            particle_list[i].solution[j][k] = cell;
                            particle_list[i].best_self_position[j][k] = cell;
                            filled_capacity += (double)particle_list[i].solution[j][k] * stock[k];
                            if (particle_list[i].solution[j][k] == 1)
                            {
                                particle_list[i].fulfillment_percentage[k] = 1;
                                particle_list[i].best_fulfillment_percentage[k] = 1;
                                allocated_set.insert(k);
                            }
                        }
                    }
                    else
                    {
                        particle_list[i].solution[j][k] = 0;
                        particle_list[i].best_self_position[j][k] = 0;
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
        particle_list[i].cost_amount = particle.cost_function(0, i, distance, node_count, truck_count);
        particle_list[i].fitness_amount = particle.fitness_function(0, i, stock, node_count);
    }
}

void PSO::particle_ranking()
{
    NonDominated first_rank;
    vector<double> rank_list;
    vector<double> cost_amount;
    vector<double> fitness_amount;
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

void PSO::best_global_particle(double &Try, vector<vector<uint64_t>> &BestGlobal, double &Best_cost, double &Best_fitness)
{
    if (Try == 0)
    {
        uint64_t rank_count = 0;
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
    else
    {
        for (uint64_t i = 1; i < population; i++)
        {
            if (particle_list[i].rank == 1)
            {
                //comparing to old best_global
                if (particle_list[i].cost_amount < Best_cost && particle_list[i].fitness_amount > Best_fitness)
                {
                    Best_cost = particle_list[i].cost_amount;
                    Best_fitness = particle_list[i].fitness_amount;
                    BestGlobal = particle_list[i].solution;
                }
                else if (particle_list[i].cost_amount < Best_cost && (((Best_cost - particle_list[i].cost_amount) / Best_cost) > .5))
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
                    if (((particle_list[i].fitness_amount - Best_fitness) / Best_fitness) > .7)
                    {
                        Best_cost = particle_list[i].cost_amount;
                        Best_fitness = particle_list[i].fitness_amount;
                        BestGlobal = particle_list[i].solution;
                    }
                }
            }
        }
    }
};

void PSO::Particle_movement(vector<vector<uint64_t>> &BestGlobal, const uint64_t &node_count, const uint64_t &truck_count, const double &truck_capacity, const vector<double> &stock, const vector<string> &partial)
{
    double max_variable = 1;
    double min_variable = 0;

    //calculating each particle new velocity
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
};

void PSO::particle_objective(const vector<vector<double>> &distance, const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock)
{
    particles objective;
    for (uint64_t i = 1; i < population; i++)
    {
        particle_list[i].cost_amount = objective.cost_function(1, i, distance, node_count, truck_count);
        particle_list[i].fitness_amount = objective.fitness_function(1, i, stock, node_count);
    }
}

void PSO::best_memory_particle(const uint64_t &node_count, const uint64_t &truck_count, const vector<double> &stock, const vector<vector<double>> &distance)
{
    particles best_memory;
    best_memory.best_self_particle(node_count, truck_count, stock, distance);
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

//will be added
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
    try
    {
        read_dataset(model, node_x, node_y, stock, partial, node_count, truck_count, truck_capacity, rate, max_iteration);
    }
    catch (const invalid_argument &e)
    {
        cout << "Error: " << e.what() << '\n';
        return -1;
    }
    //calculating the population for starting
    population = ((node_count * node_count) / 2) + 1;

    // print the main parameters of the model
    cout << "-------------------------------------------------------" << '\n';
    cout << "node count     :" << node_count << '\n';
    cout << "truck count    :" << truck_count << '\n';
    cout << "truck capacity :" << truck_capacity << '\n';
    cout << "rate           :" << rate << '\n';
    cout << "max_iteration  :" << max_iteration << '\n';
    cout << "population     :" << population << '\n';

    //The main contsiners that will be used for the algorithm
    vector<vector<double>> nodes;    //matrix of position of nodes
    vector<vector<double>> distance; //matrix of distance between nodes
    vector<double> dis_node;         //This shows distance vector for each node

    //Calculating distance matrix between nodes by Euclidean method
    for (uint64_t i = 0; i < node_count + 1; i++)
    {
        for (uint64_t j = 0; j < node_count + 1; j++)
        {
            dis_node.push_back(sqrt(((node_x[j] - node_x[i]) * (node_x[j] - node_x[i]) + (node_y[j] - node_y[i]) * (node_y[j] - node_y[i]))));
        }
        distance.push_back(dis_node);
        dis_node.clear();
    }

    // Running PSO Algorithm
    PSO run; //define an object of PSO algorithm
    static double Try = 0;

    // Step 1- initialization the particles in a population
    run.particle_initialization(node_count, truck_count, truck_capacity, stock, partial, distance);

    //velocity initialization
    for (uint64_t i = 1; i < population; i++)
    {
        for (uint64_t j = 0; j < truck_count; j++)
        {
            for (uint64_t k = 0; k < node_count; k++)
            {
                particle_list[i].velocity[j][k] = 0;
            }
        }
    }

    //put rank for initial particle based on non-domination algorithm
    run.particle_ranking();

    //initial global best particle position
    vector<vector<uint64_t>> BestGlobal;
    double Best_cost;
    double Best_fitness;
    PSO initial_BestGlobal;
    initial_BestGlobal.best_global_particle(Try, BestGlobal, Best_cost, Best_fitness);

    while (Try < max_iteration)
    {
        //check the iteration
        Try = Try + 1;

        //Step1- particles based on best velocity of self-memory and global memory will move
        run.Particle_movement(BestGlobal, node_count, truck_count, truck_capacity, stock, partial);

        // Step2- objective calculation
        run.particle_objective(distance, node_count, truck_count, stock);

        // Step3- particle best-self position
        run.best_memory_particle(node_count, truck_count, stock, distance);

        // Step4- particle ranking
        run.particle_ranking();

        // Step5- run.best_global();
        run.best_global_particle(Try, BestGlobal, Best_cost, Best_fitness);
        //will be removed
    }

    // Running Ant Colony Algorithm
    //It will be added

    //Publishing the result
    string OutputFile = "NS-PSO.csv";
    ofstream data_output;
    string Algorithm = "NS-PSO";
    output result;

    data_output.open(OutputFile, ios::out);
    if (!data_output.is_open())
    {
        cout << "Error: opening output file is not successfull!";
        data_output.close();
        return -1;
    }
    result.print(Try, Algorithm, BestGlobal, Best_cost, Best_fitness, node_count, truck_count);
    result.file(data_output, Try, BestGlobal, node_count, truck_count);
    data_output.close();
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