#include <stdio.h>
#include <math.h>

#define VARS		(250)	/* max # of variables	     */
#define RHO_BEGIN	(0.5)	/* stepsize geometric shrink */
#define EPSMIN		(1E-3)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */

/* global variables */
int	funevals = 0;
int V[] = { 700, 200, 500, 150, 800 };
int K[] = { 5, 5, 20, 3, 4 };
int S[] = { 15, 4, 10, 2, 20 };

double
f(x, n)
double	x[VARS];
int		n;
{
	double x1, x2, x3, x4, x5, fun;
	funevals++;
	x1 = x[0];
	x2 = x[1];
	x3 = x[2];
	x4 = x[3];
	x5 = x[4];
	fun = (V[0] * K[0] / x1) + (V[1] * K[1] / x2) + (V[2] * K[2] / x3) + (V[3] * K[3] / x4) + (V[4] * K[4] / x5) +
		0.5 * (S[0] * x1 + S[1] * x2 + S[2] * x3 + S[3] * x4 + S[4] * x5);
	return fun;
}

/* given a point, look for a better one nearby, one coord at a time */
double
best_nearby(delta, point, prevbest, nvars)
double	delta[VARS], point[VARS];
double	prevbest;
int		nvars;
{
	double	z[VARS];
	double	minf, ftmp;
	int		i;
	minf = prevbest;
	for (i = 0; i < nvars; i++)
		z[i] = point[i];
	for (i = 0; i < nvars; i++) {
		z[i] = point[i] + delta[i];
		ftmp = f(z, nvars);
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i] = point[i] + delta[i];
			ftmp = f(z, nvars);
			if (ftmp < minf)
				minf = ftmp;
			else
				z[i] = point[i];
		}
	}
	for (i = 0; i < nvars; i++)
		point[i] = z[i];
	return (minf);
}

int
hooke(nvars, startpt, endpt, rho, epsilon, itermax)
double	startpt[VARS], endpt[VARS];
int		nvars, itermax;
double	rho, epsilon;
{
	double	delta[VARS];
	double	newf, fbefore, steplength, tmp;
	double	xbefore[VARS], newx[VARS];
	int		i, j, keep;
	int		iters, iadj;
	for (i = 0; i < nvars; i++) {
		newx[i] = xbefore[i] = startpt[i];
		delta[i] = fabs(startpt[i] * rho);
		if (delta[i] == 0.0)
			delta[i] = rho;
	}
	iadj = 0;
	steplength = rho;
	iters = 0;
	fbefore = f(newx, nvars);
	newf = fbefore;
	while ((iters < itermax) && (steplength > epsilon)) {
		iters++;
		iadj++;
		printf("\nAfter %d funevals, f(x) =  %.4f at\n", funevals, fbefore);
		for (j = 0; j < nvars; j++)
			printf("\tx[%d] = %2.4f\n", j, xbefore[j]);
		/* find best new point, one coord at a time */
		for (i = 0; i < nvars; i++) {
			newx[i] = xbefore[i];
		}
		newf = best_nearby(delta, newx, fbefore, nvars);
		/* if we made some improvements, pursue that direction */
		keep = 1;
		while ((newf < fbefore) && (keep == 1)) {
			iadj = 0;
			for (i = 0; i < nvars; i++) {
				/* firstly, arrange the sign of delta[] */
				if (newx[i] <= xbefore[i])
					delta[i] = 0.0 - fabs(delta[i]);
				else
					delta[i] = fabs(delta[i]);
				/* now, move further in this direction */
				tmp = xbefore[i];
				xbefore[i] = newx[i];
				newx[i] = newx[i] + newx[i] - tmp;
			}
			fbefore = newf;
			newf = best_nearby(delta, newx, fbefore, nvars);
			/* if the further (optimistic) move was bad*/
			if (newf >= fbefore)
				break;
			/* make sure that the differences between the new */
			/* and the old points are due to actual */
			/* displacements; beware of roundoff errors that */
			/* might cause newf < fbefore */
			keep = 0;
			for (i = 0; i < nvars; i++) {
				keep = 1;
				if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i])))
					break;
				else
					keep = 0;
			}
		}
		if ((steplength >= epsilon) && (newf >= fbefore)) {
			steplength = steplength * rho;
			for (i = 0; i < nvars; i++) {
				delta[i] *= rho;
			}
		}
	}
	for (i = 0; i < nvars; i++)
		endpt[i] = xbefore[i];
	return (iters);
}

main()
{
	double	startpt[VARS], endpt[VARS];
	int		nvars, itermax;
	double	rho, epsilon;
	int		i, j;

	/* starting guess for rosenbrock test function */
	nvars = 5;
	startpt[0] = 1.0;	endpt[0] = 1.0;
	startpt[1] = 1.0;	endpt[1] = 1.0;
	startpt[2] = 1.0;	endpt[2] = 1.0;
	startpt[3] = 1.0;	endpt[3] = 1.0;
	startpt[4] = 1.0;	endpt[4] = 1.0;

	itermax = IMAX;
	rho = RHO_BEGIN;
	epsilon = EPSMIN;
	j = hooke(nvars, startpt, endpt, rho, epsilon, itermax);
	printf("\n\n\nNUMBER OF ITERATIONS %d\n", j);
	for (i = 0; i < nvars; i++)
		printf("\tx[%d] = %2.4f \n", i, endpt[i]);
}