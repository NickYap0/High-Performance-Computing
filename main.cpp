/**
 * @file main.cpp
 * @author Nicholas Yap (nicholas.yap19@imperial.ac.uk)
 * @brief   High Performance Computing Parallel Program Assignment;
 * Code utilises the shallow-water equations to simulate a Tsunami using 4 different test cases
 * @date 2023-08-07
 * 
 */
#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <cblas.h>
#include <omp.h>
#include "ShallowWater.h"

using namespace std;
namespace po = boost::program_options;

/**
 * @brief main function, called by MakeFile
 */
int main(int argc, char* argv[]){
    namespace po = boost::program_options;

    // Reading parameters from command line 
    po::options_description desc("Allowed Options");

    // declare arguments
    desc.add_options()
        ("dt",  po::value<double>()->default_value(0.1),  "Time-step to use.")
        ("T",   po::value<double>()->default_value(80.0),     "Total integration time.")
        ("Nx",  po::value<int>()->default_value(100),     "Number of grid points in x")
        ("Ny",  po::value<int>()->default_value(100),     "Number of grid points in y")
        ("ic",  po::value<int>()->default_value(1),  "Test case number")
        ("np",  po::value<int>()->default_value(9), "Number of Threads used")
        ("calc",po::value<int>()->default_value(1), "Loop-BLAS Selector")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    cout<<"Is this running?"<<endl;
    
    //double dx = 1.0, dy = 1.0;

    // store value in variables
    double dt = vm["dt"].as<double>();
    double T = vm["T"].as<double>();
    int Nx = vm["Nx"].as<int>();
    int Ny = vm["Ny"].as<int>();
    int ic = vm["ic"].as<int>();
    int np = vm["np"].as<int>();
    int calc = vm["calc"].as<int>();

    omp_set_num_threads(np);

    const int chunk = Nx;

    //running code
    ShallowWater tsunami(dt,T,ic,Nx,Ny,calc,chunk);
    
    //setting initial conditions
    tsunami.SetInitialConditions();

    //Time Integration
    tsunami.TimeIntegrate();

    //Saving u,v,h to txt file
    tsunami.WriteTXTFile();

    return 0;
}