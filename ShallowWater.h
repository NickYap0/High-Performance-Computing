#ifndef SHALLOWWATER_H
#define SHALLOWWATER_H

class ShallowWater
{
    private:
            //Calculating vars
            int     ic;
            double  dt;
            double  T;
            int     Nx;
            int     Ny;
            int     calc;
            int     chunk;
            const double  dx = 1.0;
            const double  dy = 1.0;

            //Storage vars
            double* M;
            double* u;
            double* v;
            double* h;
            double* M_tps;
            

    public:
            ShallowWater(double dt, double T, int ic, int Nx, int Ny, int calc, int chunk);
            ~ShallowWater();
            //void printMatrix(double* matrix, int n);
            void StencilMatrix(double* Mat, double* Mat_tps);
            void SetInitialConditions();
            void xEval(double* ddx, double* pos);
            void yEval(double* ddy, double* pos);
            void RHSCalc(double* u1, double* v1, double* h1,  double* u_calc, double* v_calc, double* h_calc);
            void TimeIntegrate();
            void WriteTXTFile();
        
};
#endif