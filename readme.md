A non-dominated sorting PSO Algorithm in C++ to solve multiple objective Optimization Problems
==============

***The vehicle route and allocation optimization model is used for this programm***

**Author:** *Paria Fakhrzad*

# 1- Introduction 
In real world whenever we encounter shipment planning in supply chain management, we have problems as how to schedule/allocate vehicles to target customers/vendors locations. This is an logistic distribution route optimization problem. The goal is finding the best allocation of vehicles with the minimum power consumption and maximum customers fulfillment.

In this program we will consider three major parts:
* **First:** designing  the intelligent optimization model whereas be able to accept different cases. The model can be multiple objective, can have cost functions and fitness functions, with constraints or without constraints.

* **Second:** solving the problem by and non-dominated sorting meta-heuristic algorithm that here we choose PSO and GA.

* **Third:** compare the results and give output

# 2- Optimization model
In this program we have a complex model that cannot be solved by exact algorithms and we need to apply evolutionary algorithms. There are many of them in application of researches during last decade. The algorithm that we consider is particle swarm optimization(PSO) that the main concept comes from the movement of birds and fishes for searching a common target. tests

![Algorithm](Algorithm.png)

Here we are going to apply these two methods together for a multiple objective optimization model. 

## 2-1- Parameter setting 
### 2-1-1 Input parameters
In input file these parameters should be entered by users:

* node_count: number of nodes in distribution network
* truck_count: number of trucks available daily
* $N_{i}=(x_{i},y_{i})$  $1<i<node\_count$ : a matrix of having position of nodes
* o : Starting point of trucks
* truck_capacity : capacity of trucks
* $s_{i}$ : volume of available stock in node i
* f: power consumption ratio of unit distance
* $ps_{i}$ : a binary parameter to show if partial shipment is allowed in node i or not

### 2-1-2 Model variable
The solution in This model will result in having below variables:
* $x_{ij}$ : a binary variable that show truck j has stop in node i or not $1<i<node\_count$  and $1<j<truck\_count$
* $F(x)$ : cost function of the model
* $G(x)$ : fitness function of the model

### 2-1-3 model set or parameter definitions
These sets or parameters will be defined during the algorithm steps
* Distance_matrix :  a i-dimensional matrix of distances between nodes
* $d_{ik}$ distance between node i and k , $d_{ik}=d_{ki}$
* $P_{i}$ the fulfilled percentage of each node

## 2-2- Objective function
Cost function here in this model is minimizing the sum of power consumption for all vehicles. the formula is:

$$F(x) =min \sum_{j=1}^{truk\_count}\sum_{i=1}^{node\_count} fd_{ik}$$

Fitness function in this model is maximize the total percentage of customer fulfillment

$$G(x) =max \sum_{i=1}^{node\_count} P_{i}$$

## 2-3- Hypothesis of model 
Since this model is complex so we assume below assumption to make it simple:
- The condition of trucks can meet all stocks requirements
- The capacity of all trucks are same
- The power consumption of trucks are same
- There is no limitation in weight of stocks
- The to tal volume of stocks should not be more than trucks capacity otherwise will see the error
- We assume that all stocks and all nodes have same important index and weight so there is no priority that which node should be services first
- We assume that the path between all nodes are possible 
- One truck cannot pass one point twice
- If the position of two or more nodes are the same and the partial shipment is allowed they stock be collected in ones
# 3- Methodology

## 3-1- Input 
There is one file that input parameters should be entered and the program reads from this file. Need to notice to below points for having clear and trusted input:
   1. Maximum iteration is optional if you don't write it the program will consider 50 as default
   2. The orange cells in below picture are mandatory and if you don't enter any number in them you will receive the error2
   3. The number of node positions should be the same as number of nodes otherwise you will receive error3
   4. If you enter other variable type you will see error4
   5. Changing other cells doesn't have any effect in the result

![input](input.png)

## 3-2- 

# 4- Result