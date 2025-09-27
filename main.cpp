//#include <iostream>
#include <vector>

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

const int max_iter = 20000;
const double residual_tolerance = 1e-6;
const double alpha_u = 0.7;
const double alpha_v = 0.7;
const double alpha_p = 0.3;

const double dx = lx / (x - 1);
const double dy = ly / (y - 1);

void setBoundaryConditions(vector<vector<double>>& u, vector<vector<double>>& v);

int main() {

    vector<vector<double>> u(x, vector<double>(y + 1, 0));
    vector<vector<double>> v(x + 1, vector<double>(y, 0));
    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> v_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));

    setBoundaryConditions(u, v);

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
