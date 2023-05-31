#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;

int Nx = 201;
const int Ny = 101;

const double Lx = 0.1, Ly = 0.05;
const double rho = 1000, nu = 1e-6;
const double P_max = 0.5;
const double t_end = 50.0;
const double dt_min = 1.e-3;
const double courant = 0.01;
const double dt_out = 0.5;

vector< vector<double> > P, P_old, u, u_old, v, v_old, PPrhs;
double dx, dy, dt, t;


// current process id
int id;

// number of rows computed by current process
int numRow;

// dimensions of the processor grid
int px, py;

// current process coordinates in the processor grid
int ix, iy;

// number of rows and columns computed by the current process
int numRow, numCol;

// tag for u, v
int uTag = 1;
int vTag = 2;

// total ranks
int total_ranks = 0;

// defined MPI type
MPI_Datatype MPI_type;

MPI_Request requests[8];

// number of requests
int cnt;


// Modifications to handle the parallel output of grid variables. 
// File names modified to include the process ID (id) to avoid 
// conflicts when multiple processes write to the same file.
void grids_to_file(int out)
{

    //Write the output for a single time step to file
    stringstream fname;
    fstream f1;

    // file name also contains process id information
    fname << "./out/P" << "_" << out << "-R_" << id << ".dat";
    f1.open(fname.str().c_str(), ios_base::out);

    int start = 1;

    int end = Nx - 1;

    // Ensures each process writes its portion of the grid 
    // Excluding boundary rows that are computed by neighboring processes.
    if (id == 0) {

        start = 0;
    }

    if (id == total_ranks - 1) {
        end = Nx;
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Ny; j++)
            f1 << P[i][j] << "\t";
        f1 << endl;
    }
    f1.close();
    fname.str("");
    fname << "./out/u" << "_" << out << "-R_" << id << ".dat";
    f1.open(fname.str().c_str(), ios_base::out);
    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Ny; j++)
            f1 << u[i][j] << "\t";
        f1 << endl;
    }
    f1.close();
    fname.str("");
    fname << "./out/v" << "_" << out << "-R_" << id << ".dat";
    f1.open(fname.str().c_str(), ios_base::out);
    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Ny; j++)
            f1 << v[i][j] << "\t";
        f1 << endl;
    }
    f1.close();
}

//update Nx and resize grid-related vectors
void setup(void)
{

    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);

    // find suitable dimensions for the processor grid
    px = int(sqrt(total_ranks));
    py = total_ranks / px;

    // determine the position of the current process in the processor grid
    ix = id % px;
    iy = id / px;

    // determine the number of rows and columns for the current process
    numRow = Nx / px;
    numCol = Ny / py;

    // adjust for possibly uneven division
    if (ix < Nx % px) ++numRow;
    if (iy < Ny % py) ++numCol;

    // add two to account for up and down, left and right boundary rows and columns
    numRow += 2;
    numCol += 2;

    P.resize(numRow, vector<double>(numCol, 0.0));
    P_old.resize(numRow, vector<double>(numCol, 0.0));
    u.resize(numRow, vector<double>(numCol, 0.0));
    u_old.resize(numRow, vector<double>(numCol, 0.0));
    v.resize(numRow, vector<double>(numCol, 0.0));
    v_old.resize(numRow, vector<double>(numCol, 0.0));
    PPrhs.resize(numRow, vector<double>(numCol, 0.0));

    if (id == 0) {
        for (int j = 0; j < Ny; j++)
            P[0][j] = P_max;
    }

    P_old = P;

    t = 0.0;

    // don't send the two corner elements
    int block_length = Ny - 2;
    MPI_Aint displacement = 0;
    MPI_Datatype typelist = MPI_DOUBLE;

    // create "MPI_Type_create_struct"
    MPI_Type_create_struct(1, &block_length, &displacement, &typelist, &MPI_type);
    MPI_Type_commit(&MPI_type);
}


void calculate_ppm_RHS_central(void)
{
    // first compute non-boundary values
    for (int i = 2; i < Nx - 2; i++) {
        for (int j = 1; j < Ny - 1; j++)
        {
            PPrhs[i][j] = rho / dt * ((u[i + 1][j] - u[i - 1][j]) / (2. * dx) + (v[i][j + 1] - v[i][j - 1]) / (2. * dy));
        }
    }

    // wait for receving boundary values
    MPI_Waitall(cnt, requests, MPI_STATUSES_IGNORE);

    // compute up boundary values
    int i = 1;
    for (int j = 1; j < Ny - 1; j++)
    {
        PPrhs[i][j] = rho / dt * ((u[i + 1][j] - u[i - 1][j]) / (2. * dx) + (v[i][j + 1] - v[i][j - 1]) / (2. * dy));
    }

    // compute down boundary values
    i = Nx - 2;
    for (int j = 1; j < Ny - 1; j++)
    {
        PPrhs[i][j] = rho / dt * ((u[i + 1][j] - u[i - 1][j]) / (2. * dx) + (v[i][j + 1] - v[i][j - 1]) / (2. * dy));
    }

}


void set_pressure_BCs(void)
{
    for (int i = 0; i < Nx; i++)
    {
        P[i][0] = P[i][1];
        P[i][Ny - 1] = P[i][Ny - 2];
    }

    // only computed by last process
    if (id == total_ranks - 1) {
        for (int j = Ny / 2; j < Ny; j++)
            P[Nx - 1][j] = P[Nx - 2][j];
    }
}

int pressure_poisson_jacobi(double rtol = 1.e-5)
{
    double tol = 10. * rtol;
    int it = 0;

    //store sum_val and tol locally computed within the function.
    double resid_vals[2];

    while (tol > rtol)
    {

        swap(P, P_old);

        double sum_val = 0.0;
        tol = 0.0;
        it++;

        //Jacobi iteration
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++)
            {
                P[i][j] = 1.0 / (2.0 + 2.0 * (dx * dx) / (dy * dy)) * (P_old[i + 1][j] + P_old[i - 1][j] +
                    (P_old[i][j + 1] + P_old[i][j - 1]) * (dx * dx) / (dy * dy)
                    - (dx * dx) * PPrhs[i][j]);

                sum_val += fabs(P[i][j]);
                tol += fabs(P[i][j] - P_old[i][j]);
            }
        }

        resid_vals[0] = sum_val;
        resid_vals[1] = tol;

        //Compute the residual sum from all processes. 
        MPI_Request request;
        MPI_Iallreduce(MPI_IN_PLACE, resid_vals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request);

        cnt = 0;

        if (id > 0) {

            // send data to up neighbor
            MPI_Isend(&P[1][1], 1, MPI_type, id - 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);

            // receive data from neighbor
            MPI_Irecv(&P[0][1], 1, MPI_type, id - 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);
        }

        if (id < total_ranks - 1) {

            // send data to down neighbor
            MPI_Isend(&P[numRow][1], 1, MPI_type, id + 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);

            // receive data from down neighbor
            MPI_Irecv(&P[numRow + 1][1], 1, MPI_type, id + 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);
        }

        set_pressure_BCs();

        // wait for sum residual finish
        MPI_Wait(&request, MPI_STATUS_IGNORE);

        sum_val = resid_vals[0];
        tol = resid_vals[1];
        tol = tol / max(1.e-10, sum_val);

        // wait sending and receiving finish
        MPI_Waitall(cnt, requests, MPI_STATUSES_IGNORE);

    }

    return it;
}

void calculate_intermediate_velocity(void)
{
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++)
        {
            //viscous diffusion
            u[i][j] = u_old[i][j] + dt * nu * ((u_old[i + 1][j] + u_old[i - 1][j] - 2.0 * u_old[i][j]) / (dx * dx) + (u_old[i][j + 1] + u_old[i][j - 1] - 2.0 * u_old[i][j]) / (dy * dy));
            v[i][j] = v_old[i][j] + dt * nu * ((v_old[i + 1][j] + v_old[i - 1][j] - 2.0 * v_old[i][j]) / (dx * dx) + (v_old[i][j + 1] + v_old[i][j - 1] - 2.0 * v_old[i][j]) / (dy * dy));
            //advection - upwinding
            if (u_old[i][j] > 0.0)
            {
                u[i][j] -= dt * u_old[i][j] * (u_old[i][j] - u_old[i - 1][j]) / dx;
                v[i][j] -= dt * u_old[i][j] * (v_old[i][j] - v_old[i - 1][j]) / dx;
            }
            else
            {
                u[i][j] -= dt * u_old[i][j] * (u_old[i + 1][j] - u_old[i][j]) / dx;
                v[i][j] -= dt * u_old[i][j] * (v_old[i + 1][j] - v_old[i][j]) / dx;
            }
            if (v_old[i][j] > 0.0)
            {
                u[i][j] -= dt * v_old[i][j] * (u_old[i][j] - u_old[i][j - 1]) / dy;
                v[i][j] -= dt * v_old[i][j] * (v_old[i][j] - v_old[i][j - 1]) / dy;
            }
            else
            {
                u[i][j] -= dt * v_old[i][j] * (u_old[i][j + 1] - u_old[i][j]) / dy;
                v[i][j] -= dt * v_old[i][j] * (v_old[i][j + 1] - v_old[i][j]) / dy;
            }
        }
    }
}

void set_velocity_BCs(void)
{
    // only computed by first process
    if (id == 0) {
        for (int j = 0; j < Ny; j++)
            u[0][j] = u[1][j];
    }

    // only computed by last process
    if (id == total_ranks - 1) {
        for (int j = 0; j < Ny / 2; j++)
            u[Nx - 1][j] = u[Nx - 2][j];
    }

}

// send and receive u, v boundary values to and from neighbours
void send_receive_u_v() {

    cnt = 0;

    if (id > 0) {

        // send u and v to up neighbour
        MPI_Isend(&u[1][1], 1, MPI_type, id - 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);
        MPI_Isend(&v[1][1], 1, MPI_type, id - 1, vTag, MPI_COMM_WORLD, &requests[cnt++]);

        // receive u and v from up neighbour
        MPI_Irecv(&u[0][1], 1, MPI_type, id - 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);
        MPI_Irecv(&v[0][1], 1, MPI_type, id - 1, vTag, MPI_COMM_WORLD, &requests[cnt++]);
    }

    if (id < total_ranks - 1) {

        // send u and v to down neighbour
        MPI_Isend(&u[numRow][1], 1, MPI_type, id + 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);
        MPI_Isend(&v[numRow][1], 1, MPI_type, id + 1, vTag, MPI_COMM_WORLD, &requests[cnt++]);

        // receive u and v from down neighbour
        MPI_Irecv(&u[numRow + 1][1], 1, MPI_type, id + 1, uTag, MPI_COMM_WORLD, &requests[cnt++]);
        MPI_Irecv(&v[numRow + 1][1], 1, MPI_type, id + 1, vTag, MPI_COMM_WORLD, &requests[cnt++]);
    }
}

double project_velocity(void)
{
    double vmax = 0.0;
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++)
        {
            u[i][j] = u[i][j] - dt * (1. / rho) * (P[i + 1][j] - P[i - 1][j]) / (2. * dx);
            v[i][j] = v[i][j] - dt * (1. / rho) * (P[i][j + 1] - P[i][j - 1]) / (2. * dy);

            double vel = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);

            vmax = max(vmax, vel);
        }
    }

    // send and receive boundary values
    send_receive_u_v();

    set_velocity_BCs();

    return vmax;
}


void solve_NS(void)
{
    double vel_max = 0.0;
    int time_it = 0;
    int its;
    int out_it = 0;
    double t_out = dt_out;

    grids_to_file(out_it);

    while (t < t_end)
    {
        if (vel_max > 0.0)
        {
            dt = min(courant * min(dx, dy) / vel_max, dt_min);
        }
        else dt = dt_min;

        MPI_Request request;

        // start to compute the minimum time step
        MPI_Iallreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &request);

        time_it++;
        swap(u, u_old);
        swap(v, v_old);

        MPI_Wait(&request, MPI_STATUS_IGNORE);

        t += dt;

        // Calculate intermediate velocity fields.
        calculate_intermediate_velocity();

        // Exchange boundary velocity values with neighboring MPI processes.
        send_receive_u_v();

        calculate_ppm_RHS_central();

        its = pressure_poisson_jacobi(1.e-5);
        vel_max = project_velocity();

        // start to compute max velocity
        MPI_Iallreduce(MPI_IN_PLACE, &vel_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &request);

        MPI_Wait(&request, MPI_STATUS_IGNORE);

        MPI_Waitall(cnt, requests, MPI_STATUSES_IGNORE);

        if (t >= t_out)
        {
            out_it++;
            t_out += dt_out;
            if (id == 0) {
                cout << time_it << ": " << t << " Jacobi iterations: " << its << " vel_max: " << vel_max << endl;
            }

            grids_to_file(out_it);
        }
    }
}

int main(int argc, char* argv[])
{


    MPI_Init(&argc, &argv);

    // total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);

    // current process id
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // number of rows for each process
    vector<int> counts(total_ranks);

    int total_rows = Nx - 2;
    int avg = total_rows / total_ranks;
    int remain = total_rows % total_ranks;

    for (int r = 0; r < total_ranks; r++) {
        counts[r] = avg;
    }

    for (int r = 0; r < remain; r++) {
        counts[r]++;
    }

    // number of rows will be computed by this process
    numRow = counts[id];

    setup();
    solve_NS();

    // free memory of mpi type
    MPI_Type_free(&MPI_type);
    MPI_Finalize();

    return 0;
}