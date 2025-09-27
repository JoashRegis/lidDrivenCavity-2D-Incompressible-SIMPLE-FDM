#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <fstream>

using namespace std;

// number of grids
const int x = 21;
const int y = 21;

// dimension of the cavity
const double lx = 1;
const double ly = 1;

// simulation properties
const double rho = 1;
const double mu = 0.01;

const double U_lid = 1;

// solver properties
const double dt = 0.005;
const double T = 5;
const int max_iter = 100;
const double residual_tolerance = 1e-6;
const double alpha_u = 0.7;
const double alpha_v = 0.7;
const double alpha_p = 0.3;

const double dx = lx / (x - 1);
const double dy = ly / (y - 1);

void setBoundaryConditions(vector<vector<double>>& u, vector<vector<double>>& v);

void solveMomentumPredictor(vector<vector<double>>& u_star, vector<vector<double>>& v_star, const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& p, vector<vector<double>>& u_old, vector<vector<double>>& v_old);

void solvePressureCorrection(vector<vector<double>>& p_prime, vector<vector<double>>& u_star, vector<vector<double>>& v_star);

int main() {

    vector<vector<double>> u(x, vector<double>(y + 1, 0));
    vector<vector<double>> u_old(x, vector<double>(y + 1, 0));
    vector<vector<double>> v(x + 1, vector<double>(y, 0));
    vector<vector<double>> v_old(x + 1, vector<double>(y, 0));
    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> v_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));

    setBoundaryConditions(u, v);

    for (double t = 0; t <= T; t += dt) {
        u_old = u;
        v_old = v;

        for (int iter = 0; iter < max_iter; iter++) {
            solveMomentumPredictor(u_star, v_star, u, v, p, u_old, v_old);

            solvePressureCorrection(p_prime, u_star, v_star);
        }
    }

    return 0;
}

void setBoundaryCondition(vector<vector<double>>& u, vector<vector<double>>& v) {
    // setting normal velocity on top and bottom
    for (int j = 0; j <= y; j++) {
        v[0][j] = 0;
        v[x-1][j] = 0;
    }

    // setting normal velocity for left and right walls
    for (int i = 0; i <= x; i++) {
        u[i][0] = 0;
        u[i][y-1] = 0;
    }

    // setting tangential velocity for top and bottom wall
    for (int i = 1; i < x-1; i++) {
        u[i][0] = -u[i][1];
        u[i][y] = 2 * U_lid - u[i][y-1];
    }

    // setting tangential velocity for left and right wall
    for (int j = 1; j < y-1; j++) {
        v[0][j] = -v[1][j];
        v[x][j] = -v[x-1][j];
    }
}

void solveMomentumPredictor(vector<vector<double>>& u_star, vector<vector<double>>& v_star, const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& p, vector<vector<double>>& u_old, vector<vector<double>>& v_old) {

    for (int i = 1; i < x-1; i++) {
        for (int j = 1; j < y; j++) {
            double u_avg_e = 0.5 * (u[i][j] + u[i+1][j]);
            double u_avg_w = 0.5 * (u[i][j] + u[i-1][j]);
            double v_avg_n = 0.5 * (v[i][j] + v[i][j+1]);
            double v_avg_s = 0.5 * (v[i][j] + v[i][j-1]);
            double u_avg_n = 0.5 * (u[i][j] + u[i][j+1]);
            double u_avg_s = 0.5 * (u[i][j] * u[i][j-1]);

            double u_conv = rho * (u_avg_e*u_avg_e - u_avg_w*u_avg_w) / dx + rho * (v_avg_n*u_avg_n - v_avg_s*u_avg_s) / dy;

            double u_pres_grad = (p[i+1][j] - p[i][j]) / dx;

            double u_diff = mu * ((u[i+1][j] - 2 * u[i][j] + u[i-1][j]) / (2 * dx) + (u[i][j+1] - 2 * u[i][j] + u[i][j-1]) / (2 * dy));

            double time_term_u = rho * u_old[i][j] * dt;

            double u_source_term = u_diff - u_conv - u_pres_grad + time_term_u;

            double a_P_u = rho / dt + 2 * mu / (dx * dx) + 2 * mu / (dy * dy) + rho * (abs(u_avg_e) + abs(u_avg_w)) / dx + rho * (abs(v_avg_n) + abs(v_avg_s)) / dy;

            u_star[i][j] = u_source_term /a_P_u;
        }
    }

    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y-1; j++) {
            double u_avg_e = 0.5 * (u[i][j] + u[i+1][j]);
            double u_avg_w = 0.5 * (u[i][j] + u[i-1][j]);
            double v_avg_n = 0.5 * (v[i][j] + v[i][j+1]);
            double v_avg_s = 0.5 * (v[i][j] + v[i][j-1]);
            double v_avg_e = 0.5 * (v[i][j] + v[i+1][j]);
            double v_avg_w = 0.5 * (v[i][j] * v[i-1][j]);

            double v_conv = rho * (v_avg_n*v_avg_n - v_avg_s*v_avg_s) / dy + rho * (v_avg_e*u_avg_e - v_avg_w*u_avg_e) / dx;

            double v_pres_grad = (p[i][j+1] - p[i][j]) / dy;

            double v_diff = mu * ((v[i+1][j] - 2 * v[i][j] + v[i-1][j]) / (2 * dx) + (v[i][j+1] - 2 * v[i][j] + v[i][j-1]) / (2 * dy));

            double time_term_v = rho * v_old[i][j] * dt;

            double v_source_term = v_diff - v_conv - v_pres_grad + time_term_v;

            double a_P_v = rho / dt + 2 * mu / (dx * dx) + 2 * mu / (dy * dy) + rho * (abs(u_avg_e) + abs(u_avg_w)) / dx + rho * (abs(v_avg_n) + abs(v_avg_s)) / dy;

            u_star[i][j] = v_source_term / a_P_v;
        }
    }
}

void solvePressureCorrection(vector<vector<double>>& p_prime, vector<vector<double>>& u_star, vector<vector<double>>& v_star) {

    for (auto& row : p_prime) {
        fill(row.begin(), row.end(), 0);
    }

    const double p_tol = 1e-5;
    const int p_iter = 100;

    for (int iter = 0; iter < p_iter; iter++) {
        double max_error = 0;

        for (int i = 1; i < x; i++) {
            for (int j = 1; j < y; j++) {

                double source = rho * ((u_star[i][j] - u_star[i-1][j]) / dx + (v_star[i][j] - v_star[i][j-1]) / dy);

                double a = 2 * dt * ((1 / (dx * dx)) + (1 / (dy * dy)));
                double b = dt / (dx * dx);
                double c = dt / (dy * dy);

                double p_prime_old = p_prime[i][j];
                
                p_prime[i][j] = ((b * p_prime[i+1][j]) - (b * p_prime[i-1][j]) - (c * p_prime[i][j+1]) - (c * p_prime[i][j-1]) - source) / a;

                max_error = max(max_error, abs(p_prime_old - p_prime[i][j]));
            }
        }

        for (int i = 1; i < x; i++) {
            p_prime[i][0] = p_prime[i][1];
            p_prime[i][x] = p_prime[i][x-1];
        }

        for (int j = 1; j < y; j++) {
            p_prime[0][j] = p_prime[1][j];
            p_prime[x][j] = p_prime[x-1][j];
        }

        if (max_error < p_tol) break;
    }
}

void write_results(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p) {
    std::ofstream outfile("results.csv");
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results.csv for writing." << std::endl;
        return;
    }

    // Write CSV header
    outfile << "X,Y,U,V,P\n";

    // Output data at cell centers
    for (int j = 1; j < y; ++j) {
        for (int i = 1; i < x; ++i) {
            double x_pos = (i - 0.5) * dx;
            double y_pos = (j - 0.5) * dy;

            // Interpolate velocities from faces to cell center
            double u_center = 0.5 * (u[i-1][j] + u[i][j]);
            double v_center = 0.5 * (v[i][j-1] + v[i][j]);
            double p_center = 0.25 * (p[i][j] + p[i+1][j] + p[i][j+1] + p[i+1][j+1]);

            outfile << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "\n";
        }
    }

    outfile.close();
    std::cout << "\nResults written to results.csv" << std::endl;
    std::cout << "You can now open the CSV file in any spreadsheet software or use the provided python script to view the data." << std::endl;
}