#include "diamond.h"

#define _USE_MATH_DEFINES
#define REL_TOL 1e-8
#define N_XTRA_ARGS 4

#include <fsolve.h>
#include <math.h>
#include <stdio.h>

typedef void (*Fcn)(int, double*, double*, int, double*);

void arc_constraints(int n, double* x, double* fvec, int nXtraArgs, double* xtraArgs)
{
    double alpha = xtraArgs[0] * M_PI / 180.0;
    double beta = xtraArgs[1] * M_PI / 180.0;
    double l = xtraArgs[2];
    double r = xtraArgs[3];

    fvec[0] = (x[0] - x[2]) * (x[0] - x[2]) + (tan(alpha) * x[0] - x[3]) * (tan(alpha) * x[0] - x[3]) - r * r;

    fvec[1] = (x[1] - x[2]) * (x[1] - x[2]) + (tan(beta) * (l - x[1]) - x[3]) * (tan(beta) * (l - x[1]) - x[3]) - r * r;

    fvec[2] = x[0] * (x[2] - x[0]) + tan(alpha) * x[0] * (x[3] - tan(alpha) * x[0]);

    fvec[3] = (l - x[1]) * (x[1] - x[2]) - tan(beta) * (l - x[1]) * (tan(beta) * (l - x[1]) - x[3]);
}

int calculate_arc_parameters(double* x_guess, Diamond* diamond)
{
    int status;
    Fcn fcn = arc_constraints;

    double fvec[N_ARC_PARAMS];
    double xtraArgs[N_XTRA_ARGS];

    xtraArgs[0] = diamond->alpha;
    xtraArgs[1] = diamond->beta;
    xtraArgs[2] = diamond->l;
    xtraArgs[3] = diamond->r;

    status = fsolve(fcn, N_ARC_PARAMS, x_guess, fvec, REL_TOL, N_XTRA_ARGS, xtraArgs);

    diamond->x1 = x_guess[0];
    diamond->x2 = x_guess[1];
    diamond->cx = x_guess[2];
    diamond->cy = x_guess[3];

    return 0;
}
