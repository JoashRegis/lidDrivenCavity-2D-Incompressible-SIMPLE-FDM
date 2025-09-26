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

    return 0;
}

void setBoundaryCondition(vector<vector<double>>& u, vector<vector<double>>& v) {
    // setting normal velocity for left and right wall
    for (int i = 0; i < x; i++) {
        u[0][i] = 0;
        u[x - 1][i] = 0;
    }

    // setting normal velocity for bottom ball and top surface
    for (int i = 0; i < y; i++) {
        v[i][0] = 0;
        v[i][y - 1] = 0;
    }
}
