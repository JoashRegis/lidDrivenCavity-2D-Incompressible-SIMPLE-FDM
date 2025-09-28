#include <vector>

using namespace std;

// number of grids
const int x = 128;
const int y = 128;

// dimension of cavity
const double lx = 1;
const double ly = 1;

// simulation parameters
const double mu = 0.01;
const double rho = 1;
const double u_lid = 1;

// solver parameters
const double dx = lx / (x - 1);
const double dy = lx / (y - 1);
const double T = 5;
const double dt = 0.0078;
const double alpha_p = 0.7;

void setBoundaryConditions(vector<vector<double>>& u, vector<vector<double>>& v) {
    for (int i = 1; i < x; i++) {
        u[i][0] = 0;
        u[i][y-1] = 0;
    }

    for (int j = 1; j < y; j++) {
        v[0][j] = 0;
        v[x-1][j] = 0;
    }

    for (int j = 0; j < y; j++) {
        u[0][j] = u_lid;
        u[1][j] = u_lid;
        u[x-1][j] = 0;
        u[x][j] = 0;
    }

    for (int i = 0; i < x; i++) {
        v[i][0] = 0;
        v[i][1] = 0;
        v[i][y-1] = 0;
        v[i][y] = 0;
    }
}

void SolveMomentum(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& p, vector<vector<double>>& u_star, vector<vector<double>>& v_star) {
    //x-momentum equation
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y-1; j++){
            double u_conv = (((u[i][j+1] * u[i][j+1]) - (u[i][j-1] * u[i][j-1]) / (2 * dx))) + ((u[i][j] + u[i-1][j]) * (v[i-1][j+1] + v[i-1][j]) - (u[i][j] + u[i+1][j]) * (v[i][j+1] + v[i][j])) / (4 * dy);

            double u_pres_grad = (1 / rho) * (p[i][j+1] - p[i][j]) / dx;

            double u_diff = (mu / rho) * ((u[i][j+1] - 2 * u[i][j] + u[i][j-1]) / (dx * dx) + (u[i-1][j] - 2 * u[i][j] + u[i+1][j]) / (dy * dy));

            double u_rhs = u_diff - u_pres_grad - u_conv;

            u_star[i][j] = (u_rhs * dt) + u[i][j];
        }
    }

}

int main() {

    vector<vector<double>> u(x + 1, vector<double>(y, 0));
    vector<vector<double>> v(x, vector<double>(y + 1, 0));
    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> v_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));

    setBoundaryConditions(u, v);

    for (int time = 0; time < T; time += dt) {

    }

    return 0;
}
