#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cblas.h>
#include <cstring>
#include <chrono>

#include "ShallowWater.h"
#include <omp.h>

using namespace std;

#define SW ShallowWater

/**
 * @brief Construct a new SW::ShallowWater object
 * 
 * @param dt    Time Step
 * @param T     Total time
 * @param ic    Test case initial condition indicator (Test Cases specified in function SW::SetInitialConditions)
 * @param Nx    Number of grid points in the x-direction
 * @param Ny    Number of grid points in the y-direction
 * @param calc  Choice of Loops(1) or BLAS(2) based calculations
 * @param chunk Defining chunk size to split workload for loop based parallelisim
 * @param u     x-component of velocity
 * @param v     y-component of velocity
 * @param h     Surface height
 * @param M     Stencil Matrix storage for BLAS application on x-derivatives
 * @param M_tps Stencil Matrix storage for BLAS application on y-derivatives
 */
SW::ShallowWater(double dt, double T, int ic, int Nx, int Ny, int calc, int chunk)
{ 
    this->dt = dt;
    this->T = T;
    this->ic = ic;
    this->Nx = Nx;
    this->Ny = Ny;
    this->calc = calc;
    this->chunk = chunk;

    // allocate memory for matrix
    u = new double[Nx*Ny];
    v = new double[Nx*Ny];
    h = new double[Nx*Ny];
    M = new double[Nx*Ny];
    M_tps = new double[Nx*Ny];
}


SW::~ShallowWater() {
    if(u) delete[] u;
    if(v) delete[] v;
    if(h) delete[] h;
    if(M) delete[] M;
    if(M_tps) delete[] M_tps;
}

/**
 * @brief Sets the initial conditions for the test cases, and initialises the spatial derivative terms for BLAS implementation
 *
 * Test Case 1: Waves Propagating in x starting at x = 50.
 * Test Case 2: Waves Propagating in y starting at y = 50.
 * Test Case 3: Single droplet at x = 50, y = 50.
 * Test Case 4: Two droplets at x = 25, y = 25 and x = 75, y = 75.
 * 
 */
void SW::SetInitialConditions(){
    // Setting u and v initial conditions 
    memset(u, 0.0, (Nx*Ny)*sizeof(double));
    memset(v, 0.0, (Nx*Ny)*sizeof(double));
    
    switch (ic) {
        case 1:
            cout<<"Test Case 1: Waves Propagating in x starting at x = 50."<<endl;
            for (int i = 0; i<Nx ; i++){ //i = x
                for (int j = 0; j<Ny ; j++){ // j = y
                    h[ i + Ny*j ] = 10.0 + (double)(exp(-(i*dx-50.0)*(i*dx-50.0)/25.0));
                }
            }
            break;
        case 2:
            cout<<"Test Case 2: Waves Propagating in y starting at y = 50."<<endl;
            for (int i = 0; i<Ny ; i++){ // i = y
                for (int j = 0; j<Nx ; j++){ // j = x
                    h[i * Ny +  j ] =  10.0 + (double)(exp(-(i*dy-50.0)*(i*dy-50.0)/25.0));
                }
            }
            break;
        case 3:
            cout<<"Test Case 3: Single droplet at x = 50, y = 50."<<endl;
            for (int i = 0; i<Nx ; i++){ //i=x
                for (int j = 0; j<Ny ; j++){ //j=y
                    h[ i + Ny * j ] = 10.0 + (double)(exp(-((i*dx-50.0)*(i*dx-50.0)+(j*dy-50.0)*(j*dy-50.0))/25.0));
                }
            }
            break;
        case 4:
            cout<<"Test Case 4: Two droplets at x = 25, y = 25 and x = 75, y = 75."<<endl;
            for (int i = 0; i<Nx ; i++){//i=x
                for (int j = 0; j<Ny ; j++){//j=y
                    h[ i + Ny * j ] = 10.0 + (double)(exp(-((i*dx-25.0)*(i*dx-25.0)+(j*dy-25.0)*(j*dy-25.0))/25.0)) + (double)(exp(-((i*dx-75.0)*(i*dx-75.0)+(j*dy-75.0)*(j*dy-75.0))/25.0));
                }
            }
            break;
    
    }

    switch (calc) {
    case 1:
        break;
    case 2:
        StencilMatrix(M,M_tps);
        break;
    }

}

/**
 * @brief Generates a stencil matrix for BLAS based evaluation of x and y derivatives. Matrix generated is as follows
 * 0        0.75        -0.15       0.016667        0       0       0.........       -0.01667        0.15        -0.75
 * -0.75    0           0.75        -0.15           0.01667 0       0.........       0              -0.01667      0.15
 * 0.15     -0.75       0           0.75            -0.15   0.01667 0.........       0               0        -0.01667
 * -0.01667 0.15        -0.75       0               0.75    -0.15   0.01667...       0               0               0
 * .                                                                .........
 * .
 * .
 * .
 * .
 * 0.01667  0           0           0           0           0..... -0.01667   0.15      -0.75     0      0.75     -0.15
 * -0.15    0.01667     0           0           0           0.....  0        -0.01667   0.15      -0.75     0      0.75
 * 0.75     -0.15       0.01667     0           0           0.....  0         0        -0.01667   0.15      -0.75     0
 * @param Mat       Stencil Matrix storage for BLAS application on x-derivatives
 * @param Mat_tps   Stencil Matrix storage for BLAS application on y-derivatives
 */
void SW::StencilMatrix(double* Mat, double* Mat_tps) {
    //stencil vars
    const double m1 = 3.0 / 4.0;
    const double m2 = -3.0 / 20.0;
    const double m3 = 1.0/60.0;
    
    //filling diagonals
    //Upper Diagonal
    for (int  i = 0; i<Nx-1 ; i++){
        Mat[i * Nx + i + 1] = m1;
        Mat[i * Nx + i + 2] = m2;
        Mat[i * Nx + i + 3] = m3;
    }
    //Lower diagonal
    for (int i = 1; i<Nx ; i++){
        Mat[i * Nx + i - 1] = -m1;
        Mat[i * Nx + i - 2] = -m2;
        Mat[i * Nx + i - 3] = -m3;
    }

    //Filling Upper and lower triangles to enforce periodic BC
    //Upper
    Mat[Nx - 3] = -m3;
    Mat[Nx - 2] = -m2;
    Mat[Nx - 1] = -m1;
    Mat[2 * Nx - 1] = -m2;
    Mat[2 * Nx - 2] = -m3;
    Mat[3 * Nx - 1] = -m3;
    // Lower
    Mat[Nx * Ny - Ny]         = m1;//Last Row
    Mat[Nx * Ny - Ny + 1]     = m2;
    Mat[Nx * Ny - Ny + 2]     = m3;
    Mat[Nx * Ny - 2 * Ny]     = m2;//Second Last Row
    Mat[Nx * Ny - 2 * Ny + 1] = m3;
    Mat[Nx * Ny - 3 * Ny]     = m3;//Third Last Row

    //Creating a transpose of the Stencil Matrix to evaluate y-direction differentials
    for(int j = 0; j<Ny; j++){
        for(int i = 0; i<Nx; i++){
            //M_tps[i * Nx + j] = M[j * Ny + i];
            Mat_tps[i * Nx + j] = -Mat[i * Nx + j];
        }
    }
}

/**
 * @brief Loop based evaluation of the x-differentials
 * 
 * @param ddx Output derivative
 * @param pos Input vector/height
 */
void SW::xEval(double* in, double* out) {
    const double m1 = 0.75;
    const double m2 = -0.15;
    const double m3 = 1.0/60.0;

    #pragma omp parallel for collapse(2)//schedule(static,chunk)
    for (int i = 0; i < Ny; i++) {//col
        for (int j = 0; j < Nx; j++) {//row
            out[i+Ny*j] = (((j >= 3) ? in[(j-3)*Ny+i] : 0.0) * ((-1.0) * m3) + 
                              ((j >= 2) ? in[i+Ny*(j-2)] : 0.0) * ((-1.0) * m2) + 
                              ((j >= 1) ? in[i+Ny*(j-1)] : 0.0) * ((-1.0) * m1) + 
                              ((j < Nx - 1) ? in[i+Ny*(j+1)] : 0.0) * (m1) + 
                              ((j < Nx - 2) ? in[i+Ny*(j+2)] : 0.0) * (m2) +
                              ((j < Nx - 3) ? in[i+Ny*(j+3)] : 0.0) * (m3) +
                              //setting periodic BCs
                              ((j == 0) ? (in[i+Ny*(Nx-3)] * ((-1.0) * m3) + in[(Nx-2)*Ny+i]*((-1.0) * m2) + in[(Nx-1)*Ny+i]*((-1.0) * m1)): 0.0) +
                              ((j == 1) ? (in[(Nx-2)*Ny+i]*((-1.0) * m3) + in[(Nx-1)*Ny+i]*((-1.0) * m2)) : 0.0) +
                              ((j == 2) ? (in[(Nx-1)*Ny+i]*((-1.0) * m3)) : 0.0) +
                              ((j == Nx-1) ? (in[i]*(m1) + in[Ny+i]*(m2) + in[2*Ny+i]*(m3)) : 0.0) +
                              ((j == Nx-2) ? (in[i]*(m2) + in[Ny+i]*(m3)) : 0.0) +
                              ((j == Nx-3) ? (in[i]*(m3)) : 0.0)
                              )/ dx;
        }
    }
}

/**
 * @brief Loop based evaluation of the y-differentials
 * 
 * @param ddy Output derivative
 * @param pos Input vector/height
 */
void SW::yEval(double* in, double* out) {
    const double m1 = 3.0/4.0;
    const double m2 = -3.0/20.0;
    const double m3 = 1.0/60.0;

    #pragma omp parallel for collapse(2)//schedule(static,chunk)
    for (int i = 0; i < Nx; i++) {//row
        for (int j = 0; j < Ny; j++) {//col
            out[i*Nx+j] = (((j >= 3) ? in[i*Ny+(j-3)] : 0.0) * ((-1.0) * m3) + 
                              ((j >= 2) ? in[i*Nx+(j-2)] : 0.0) * ((-1.0) * m2) + 
                              ((j >= 1) ? in[i*Nx+(j-1)] : 0.0) * ((-1.0) * m1) + 
                              ((j < Nx - 1) ? in[i*Nx+(j+1)] : 0.0) * (m1) + 
                              ((j < Nx - 2) ? in[i*Nx+(j+2)] : 0.0) * (m2) +
                              ((j < Nx - 3) ? in[i*Nx+(j+3)] : 0.0) * (m3) +
                              //Setting periodic BCs
                              ((j == 0) ? (in[(i+1)*Ny-3]*((-1.0) * m3) + in[(i+1)*Ny-2]*((-1.0) * m2) + in[(i+1)*Ny-1]*((-1.0) * m1)): 0.0) +
                              ((j == 1) ? (in[(i+1)*Ny-2]*((-1.0) * m3) + in[(i+1)*Ny-1]*((-1.0) * m2)) : 0.0) +
                              ((j == 2) ? (in[(i+1)*Ny-1]*((-1.0) * m3)) : 0.0) +
                              ((j == Nx-1) ? (in[i*Ny] * (m1) + in[i*Ny+1] * (m2) + in[i*Ny+2] * (m3)) : 0.0) +
                              ((j == Nx-2) ? (in[i*Ny] * (m2) + in[i*Ny+1] * (m3)) : 0.0) +
                              ((j == Nx-3) ? (in[i*Ny]*(m3)) : 0.0)
                              ) / dy;
        }
    }
}

/**
 * @brief Function that evaluates the "Right Hand Side" of the 2-D ShallowWater eqations, summing them up to form time-derivatives
 * 
 * @param u1        Input u-velocity component
 * @param v1        Input v-velocity component
 * @param h1        Input surface height component
 * @param u_calc    Output k-value for u-velocity component
 * @param v_calc    Output k-value for v-velocity component
 * @param h_calc    Output k-value for surface height component
 */
void SW::RHSCalc(double* u1, double* v1, double* h1, double* u_calc, double* v_calc, double* h_calc) {

    const double g = 9.81;

    double* DUDX = new double[Nx*Ny];
    double* DUDY = new double[Nx*Ny];
    double* DVDX = new double[Nx*Ny];
    double* DVDY = new double[Nx*Ny];
    double* DHDX = new double[Nx*Ny];
    double* DHDY = new double[Nx*Ny];

    double* ududx = new double[Nx*Ny];
    double* udvdx = new double[Nx*Ny];
    double* udhdx = new double[Nx*Ny];
    double* vdudy = new double[Nx*Ny];
    double* vdvdy = new double[Nx*Ny];
    double* vdhdy = new double[Nx*Ny];
    double* hdiff  = new double[Nx*Ny];
    
    switch (calc) {
    case 1:
        #pragma omp parallel
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                xEval(u1, DUDX);

                #pragma omp section
                xEval(v1, DVDX);

                #pragma omp section
                xEval(h1, DHDX);

                #pragma omp section
                yEval(u1, DUDY);

                #pragma omp section
                yEval(v1, DVDY);

                #pragma omp section
                yEval(h1, DHDY);
            }
        }

        #pragma omp parallel for schedule(static, chunk)
        for (int i = 0; i < Nx*Ny; i++) {
            u_calc[i] = -(u1[i] * DUDX[i]) -(v1[i] * DUDY[i]) -(g * DHDX[i]);
            v_calc[i] = -(u1[i] * DVDX[i]) -(v1[i] * DVDY[i]) -(g * DHDY[i]);
            h_calc[i] = -(u1[i] * DHDX[i]) -(v1[i] * DHDY[i]) -(h1[i] * (DUDX[i] + DVDY[i]));
        }
        break;
    case 2:
        #pragma omp parallel
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0/dx, M, Nx, u1, Ny, 0.0, DUDX, Nx);

                #pragma omp section
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0/dx, M, Nx, v1, Ny, 0.0, DVDX, Nx);

                #pragma omp section
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0/dx, M, Nx, h1, Ny, 0.0, DHDX, Nx);

                #pragma omp section
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0/dy, M_tps, Nx, u1, Ny, 0.0, DUDY, Nx);

                #pragma omp section
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0/dy, M_tps, Nx, v1, Ny, 0.0, DVDY, Nx);

                #pragma omp section
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nx, Ny, Nx, 1.0/dy, M_tps, Nx, h1, Ny, 0.0, DHDY, Nx);
            }
        }

        //Conducting elementwise multiplication
        #pragma omp parallel
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, u1, 1, DUDX, 1, 0.0, ududx, 1);

                #pragma omp section
                cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, u1, 1, DVDX, 1, 0.0, udvdx, 1);

                #pragma omp section
                cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, u1, 1, DHDX, 1, 0.0, udhdx, 1);

                #pragma omp section
                cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, v1, 1, DUDY, 1, 0.0, vdudy, 1);

                #pragma omp section
                cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, v1, 1, DVDY, 1, 0.0, vdvdy, 1);

                #pragma omp section
                cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, v1, 1, DHDY, 1, 0.0, vdhdy, 1);

                #pragma omp section
                {
                    cblas_daxpy(Nx*Ny, 1.0, DUDX, 1, DVDY, 1);
                    cblas_dgbmv(CblasColMajor, CblasNoTrans, Nx*Ny, Nx*Ny, 0, 0, 1.0, h1, 1, DVDY, 1, 0.0, hdiff, 1);
                }
            }
        }

        // sum terms up for dudt, dvdt, dhdt
        
        cblas_daxpy(Nx*Ny, 1.0, ududx, 1, vdudy, 1);
        cblas_daxpy(Nx*Ny, g, DHDX, 1, vdudy, 1);
        cblas_daxpy(Nx*Ny, -1.0, vdudy, 1, u_calc, 1);

        cblas_daxpy(Nx*Ny, 1.0, udvdx, 1, vdvdy, 1);
        cblas_daxpy(Nx*Ny, g, DHDY, 1, vdvdy, 1);
        cblas_daxpy(Nx*Ny, -1.0, vdvdy, 1, v_calc, 1);

        cblas_daxpy(Nx*Ny, 1.0, udhdx, 1, vdhdy, 1);
        cblas_daxpy(Nx*Ny, 1.0, hdiff, 1, vdhdy, 1);
        cblas_daxpy(Nx*Ny, -1.0, vdhdy, 1, h_calc, 1);

        break;
    }

    delete[] DUDX;
    delete[] DVDX;
    delete[] DHDX;
    delete[] DUDY;
    delete[] DVDY;
    delete[] DHDY;

    delete[] ududx;
    delete[] udvdx;
    delete[] udhdx;
    delete[] vdudy;
    delete[] vdvdy;
    delete[] vdhdy;
    delete[] hdiff;
}

/**
 * @brief Integrating the RHS to find spatial (u, v, h) terms using 4th order Runge-Kutta
 * 
 */
void SW::TimeIntegrate(){
    double* k1u = new double [Nx*Ny];
    double* k2u = new double [Nx*Ny];
    double* k3u = new double [Nx*Ny];
    double* k4u = new double [Nx*Ny];

    double* k1v = new double [Nx*Ny];
    double* k2v = new double [Nx*Ny];
    double* k3v = new double [Nx*Ny];
    double* k4v = new double [Nx*Ny];

    double* k1h = new double [Nx*Ny];
    double* k2h = new double [Nx*Ny];
    double* k3h = new double [Nx*Ny];
    double* k4h = new double [Nx*Ny];

    double* u1temp = new double [Nx*Ny];
    double* v1temp = new double [Nx*Ny];
    double* h1temp = new double [Nx*Ny];

    double* u2temp = new double [Nx*Ny];
    double* v2temp = new double [Nx*Ny];
    double* h2temp = new double [Nx*Ny];

    double* u3temp = new double [Nx*Ny];
    double* v3temp = new double [Nx*Ny];
    double* h3temp = new double [Nx*Ny];

    const double a = dt /2.0;
    const double b = dt / 6.0;
    const double c = b * 2.0;
    //cout<<"Not broken yet"<<endl;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    for(double i = 0.0; i < T; i += dt){
        //resetting k terms to 0
        memset(k1u, 0.0, (Nx*Ny)*sizeof(double));
        memset(k2u, 0.0, (Nx*Ny)*sizeof(double));
        memset(k3u, 0.0, (Nx*Ny)*sizeof(double));
        memset(k4u, 0.0, (Nx*Ny)*sizeof(double));

        memset(k1v, 0.0, (Nx*Ny)*sizeof(double));
        memset(k2v, 0.0, (Nx*Ny)*sizeof(double));
        memset(k3v, 0.0, (Nx*Ny)*sizeof(double));
        memset(k4v, 0.0, (Nx*Ny)*sizeof(double));

        memset(k1h, 0.0, (Nx*Ny)*sizeof(double));
        memset(k2h, 0.0, (Nx*Ny)*sizeof(double));
        memset(k3h, 0.0, (Nx*Ny)*sizeof(double));
        memset(k4h, 0.0, (Nx*Ny)*sizeof(double));

        //evaluating k1
        RHSCalc(u,v,h,k1u,k1v,k1h);
        #pragma omp parallel for schedule(static,chunk)
        for(int i = 0; i< Nx*Ny ; i++){
            u1temp[i] = u[i] + (a * k1u[i]);
            v1temp[i] = v[i] + (a * k1v[i]);
            h1temp[i] = h[i] + (a * k1h[i]);
        }

        //evaluating k2
        RHSCalc(u1temp,v1temp,h1temp,k2u,k2v,k2h);
        #pragma omp parallel for schedule(static,chunk)
        for(int i = 0; i< Nx*Ny ; i++){
            u2temp[i] = u[i] + (a * k2u[i]);
            v2temp[i] = v[i] + (a * k2v[i]);
            h2temp[i] = h[i] + (a * k2h[i]);
        }

        //evaluating k3
        RHSCalc(u2temp,v2temp,h2temp,k3u,k3v,k3h);
        #pragma omp parallel for schedule(static,chunk)
        for(int i = 0; i< Nx*Ny ; i++){
            u3temp[i] = u[i] + (dt * k3u[i]);
            v3temp[i] = v[i] + (dt * k3v[i]);
            h3temp[i] = h[i] + (dt * k3h[i]);
        }

        //evaluating k4
        RHSCalc(u3temp,v3temp,h3temp,k4u,k4v,k4h);
        #pragma omp parallel for schedule(static,chunk)
        for(int i = 0; i < Nx*Ny ; i++){
            u[i] += (b * (k1u[i] + k4u[i])) + (c * (k2u[i] + k3u[i]));
            v[i] += (b * (k1v[i] + k4v[i])) + (c * (k2v[i] + k3v[i]));
            h[i] += (b * (k1h[i] + k4h[i])) + (c * (k2h[i] + k3h[i]));
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    
    cout << "Time Spent (sec) = " <<  (chrono::duration_cast<chrono::microseconds>(end - begin).count()) /1000000.0  <<endl;

    cout<<"Evaluation for Test Case "<< ic << " complete using "<<((calc == 1) ? "loops based calculations.": "BLAS based calculations")<<endl;

    delete[] k1u;
    delete[] k1v;
    delete[] k1h;
    delete[] k2u;
    delete[] k2v;
    delete[] k2h;
    delete[] k3u;
    delete[] k3v;
    delete[] k3h;
    delete[] k4u;
    delete[] k4v;
    delete[] k4h;
    delete[] u1temp;
    delete[] v1temp;
    delete[] h1temp;
    delete[] u2temp;
    delete[] v2temp;
    delete[] h2temp;
    delete[] u3temp;
    delete[] v3temp;
    delete[] h3temp;
}

/**
 * @brief Writes u, v, h to a txt file for gnuplot
 * 
 */
void SW::WriteTXTFile(){
    // Open the output file stream
    ofstream outfile;
    outfile.open("output.txt");
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < Nx; j++){
                outfile << j*dx << " " << i*dy << " " << u[i*Ny+j] << " " << v[i*Nx+j] << " " << h[i*Nx+j]  << endl;
            }
            outfile << endl;
        }
    outfile.close();
}