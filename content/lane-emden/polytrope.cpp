#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

struct State {
    double theta;
    double dtheta_dxi;
};

class Polytrope {

public:

    std::vector<double> xi;
    std::vector<State> sol;

    int n;

    Polytrope(const int _n=1.5) {
        assert(_n > 0);
        n = _n;
    }

    State rhs(const double x, const State& s) {

        State f{0.0};

        f.theta = s.dtheta_dxi;

        if (x == 0.0) {
            f.dtheta_dxi = (2.0/3.0) - std::pow(s.theta, n);
        } else {
            f.dtheta_dxi = -2.0 * s.dtheta_dxi / x - std::pow(s.theta, n);
        }

        return f;
    }

    int npts() {
        return xi.size();
    }

    void integrate(const double h0=1.e-2, const double tol=1.e-12) {

        // initial conditions

        State y{1.0, 0.0};

        // storage for the stages

        State tmp{0.0};

        auto h = h0;

        double _xi = 0.0;

        while (h > tol) {

            // 4th order RK integration

            auto k1 = rhs(_xi, y);

            tmp.theta = y.theta + 0.5 * h * k1.theta;
            tmp.dtheta_dxi = y.dtheta_dxi + 0.5 * h * k1.dtheta_dxi;

            auto k2 = rhs(_xi + 0.5*h, tmp);

            tmp.theta = y.theta + 0.5 * h * k2.theta;
            tmp.dtheta_dxi = y.dtheta_dxi + 0.5 * h * k2.dtheta_dxi;

            auto k3 = rhs(_xi + 0.5*h, tmp);

            tmp.theta = y.theta + h * k3.theta;
            tmp.dtheta_dxi = y.dtheta_dxi + h * k3.dtheta_dxi;

            auto k4 = rhs(_xi + h, tmp);

            y.theta += (1./6.) * h * (k1.theta + 2.0*k2.theta + 2.0*k3.theta + k4.theta);
            y.dtheta_dxi += (1./6.) * h * (k1.dtheta_dxi + 2.0*k2.dtheta_dxi + 2.0*k3.dtheta_dxi + k4.dtheta_dxi);
            

            _xi += h;

            // set the new stepsize--our systems is always convex
            // (theta'' < 0), so the intersection of theta' with the
            // x-axis will always be a conservative estimate of the
            // radius of the star.  Make sure that the stepsize does
            // not take us past that.

            double R_est = _xi - y.theta/y.dtheta_dxi;

            if (_xi + h > R_est) {
                h = -y.theta/y.dtheta_dxi;
            }

            // store the solution:

            xi.push_back(_xi);
            sol.push_back(y);

        }
    }

};


int main() {

    Polytrope p(1.5);

    p.integrate();

    for (int n = 0; n < p.npts(); ++n) {
        std::cout << std::setw(20) << p.xi[n] << std::setw(20) << p.sol[n].theta << std::endl;
    }
}