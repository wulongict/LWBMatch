/*
	Copyright (c) <year> <copyright holders>

	Permission is hereby granted, free of charge, to any person
	obtaining a copy of this software and associated documentation
	files (the "Software"), to deal in the Software without
	restriction, including without limitation the rights to use,
	copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following
	conditions:

	The above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
	OTHER DEALINGS IN THE SOFTWARE.

	*/

#include "lowess.h"
#include <math.h>
#include <iostream>
 double fmax2(double x, double y)
{
	if (x > y) return x;
	else return y;
}

 int imax2(int x, int y)
{
	if (x > y) return x;
	else return y;
}
 int imin2(int x, int y)
{
	if (x > y) return y;
	else return x;
}

#define psort_body						\
    bool nalast=1;					\
    int L, R, i, j;						\
								\
    for (L = lo, R = hi; L < R; ) {				\
	v = x[k];						\
	for(i = L, j = R; i <= j;) {				\
		    while (TYPE_CMP(x[i], v, nalast) < 0) i++;			\
				    while (TYPE_CMP(v, x[j], nalast) < 0) j--;			\
	    if (i <= j) { w = x[i]; x[i++] = x[j]; x[j--] = w; }\
		}							\
	if (j < k) L = i;					\
	if (k < i) R = j;					\
	    }


bool ISNAN(double x){
	if (fabs(x) < 1e-8) return 0;
	else return 1;
}
int rcmp(double x, double y, bool nalast)
{
	int nax = ISNAN(x), nay = ISNAN(y);
	if (nax && nay)	return 0;
	if (nax)		return nalast ? 1 : -1;
	if (nay)		return nalast ? -1 : 1;
	if (x < y)		return -1;
	if (x > y)		return 1;
	return 0;
}

void rPsort2(double *x, int lo, int hi, int k)
{
	double v, w;
#define TYPE_CMP rcmp
	psort_body
#undef TYPE_CMP
}

void rPsort(double *x, int n, int k)
{
	rPsort2(x, 0, n - 1, k);
}
double fsquare(double x)
{
	return x * x;
}

double fcube(double x)
{
	return x * x * x;
}

void lowest(double *x, double *y, int n, double *xs, double *ys,
	int nleft, int nright, double *w,
	bool userw, double *rw, bool *ok)
{
    //cout << "in lowest" << endl;
	int nrt, j;
	double a, b, c, h, h1, h9, r, range;

	x--;
	y--;
	w--;
	rw--;

	range = x[n] - x[1];
	h = fmax2(*xs - x[nleft], x[nright] - *xs);
	h9 = 0.999*h;
	h1 = 0.001*h;

	/* sum of weights */

	a = 0.;
	j = nleft;
	while (j <= n) {

		/* compute weights */
		/* (pick up all ties on right) */

		w[j] = 0.;
		r = fabs(x[j] - *xs);
		if (r <= h9) {
			if (r <= h1)
				w[j] = 1.;
			else
				w[j] = fcube(1. - fcube(r / h));
			if (userw)
				w[j] *= rw[j];
			a += w[j];
		}
		else if (x[j] > *xs)
			break;
		j = j + 1;
	}

	/* rightmost pt (may be greater */
	/* than nright because of ties) */

	nrt = j - 1;
	if (a <= 0.)
		*ok = 0;
	else {
		*ok = 1;

		/* weighted least squares */
		/* make sum of w[j] == 1 */

		for (j = nleft; j <= nrt; j++)
			w[j] /= a;
		if (h > 0.) {
			a = 0.;

			/*  use linear fit */
			/* weighted center of x values */

			for (j = nleft; j <= nrt; j++)
				a += w[j] * x[j];
			b = *xs - a;
			c = 0.;
			for (j = nleft; j <= nrt; j++)
				c += w[j] * fsquare(x[j] - a);
			if (sqrt(c) > 0.001*range) {
				b /= c;

				/* points are spread out */
				/* enough to compute slope */

				for (j = nleft; j <= nrt; j++)
					w[j] *= (b*(x[j] - a) + 1.);
			}
		}
		*ys = 0.;
		for (j = nleft; j <= nrt; j++)
			*ys += w[j] * y[j];
	}
    //cout << "out of lowest" << endl;
}

// commentby lwu
// Read this part carefully
void clowess( int n,	double f, int nsteps, double delta, LOWESSData &ld	)
{
	int i, iter, j, last, m1, m2, nleft, nright, ns;
	bool ok;
	double alpha, c1, c9, cmad, cut, d1, d2, denom, r, sc;
    //----------------------------------------------------------------------------
    //vector <double> ys, rw,res;
    cout << "[Info] n=" << n<< endl;

    vector<double> lowess_x = ld.getM_lowess_x(), lowess_y = ld.getM_lowess_y();
    cout << "[Info] size = " << lowess_x.size() << endl;
    double * y = new double[n+2];
    double * x = new double[n+2];
    double * ys = new double[n+2];
    double * rw = new double[n+2];
    double * res = new double[n+2];

    double *xtmp = x, *ytmp = y, *ystmp=ys, *rwtmp = rw, *restmp = res;

    for (int idx = 0; idx < n; ++idx) {
        x[idx] = lowess_x[idx];
        y[idx] = lowess_y[idx];
        ys[idx] = 0;
        rw[idx] = 0;
        res[idx] = 0;
    }

    //---------------------------------------------------------------------------


	if (n < 2) {
		ys[0] = y[0]; return;
	}

	/* nleft, nright, last, etc. must all be shifted to get rid of these: */
	x--;
	y--;
	ys--;


	/* at least two, at most n points */
	ns = imax2(2, imin2(n, (int)(f*n + 1e-7)));
    cout << "ns = " << ns << endl;


	/* robustness iterations */

	iter = 1;
	while (iter <= nsteps + 1) {
        cout << "[Info] iter ----------------------------------" << iter << endl;
		nleft = 1;
		nright = ns;
		last = 0;	/* index of prev estimated point */
		i = 1;		/* index of current point */

		for (;;) {
			if (nright < n) {

				/* move nleft,  nright to right */
				/* if radius decreases */

				d1 = x[i] - x[nleft];
				d2 = x[nright + 1] - x[i];

				/* if d1 <= d2 with */
				/* x[nright+1] == x[nright], */
				/* lowest fixes */

				if (d1 > d2) {

					/* radius will not */
					/* decrease by */
					/* move right */

					nleft++;
					nright++;
					continue;
				}
			}

			/* fitted value at x[i] */
            //cout <<"[Info] x, y= " << x[1] << " " << y[1] << " center is " << x[i]<< endl;
			lowest(&x[1], &y[1], n, &x[i], &ys[i],
				nleft, nright, res, iter > 1, rw, &ok);
			if (!ok) ys[i] = y[i];

			/* all weights zero */
			/* copy over value (all rw==0) */
             //cout << "in block a" << endl;
			if (last < i - 1) {
				denom = x[i] - x[last];

				/* skipped points -- interpolate */
				/* non-zero - proof? */

				for (j = last + 1; j < i; j++) {
					alpha = (x[j] - x[last]) / denom;
					ys[j] = alpha*ys[i] + (1. - alpha)*ys[last];
				}
			}
            //cout << "out block a" << endl;

			/* last point actually estimated */
			last = i;

			/* x coord of close points */
            //cout << "in block b" << endl;
			cut = x[last] + delta;
			for (i = last + 1; i <= n; i++) {
				if (x[i] > cut)
					break;
				if (x[i] == x[last]) {
					ys[i] = ys[last];
					last = i;
				}
			}
            //cout << "out block b" << endl;
			i = imax2(last + 1, i - 1);
			if (last >= n)
				break;
		}
		/* residuals */
		for (i = 0; i < n; i++)
			res[i] = y[i + 1] - ys[i + 1];

		/* overall scale estimate */
		sc = 0.;
		for (i = 0; i < n; i++) sc += fabs(res[i]);
		sc /= n;

		/* compute robustness weights */
		/* except last time */

		if (iter > nsteps)
			break;
		/* Note: The following code, biweight_{6 MAD|Ri|}
		   is also used in stl(), loess and several other places.
		   --> should provide API here (MM) */
		for (i = 0; i < n; i++)
			rw[i] = fabs(res[i]);

		/* Compute   cmad := 6 * median(rw[], n)  ---- */
		/* FIXME: We need C API in R for Median ! */
		m1 = n / 2;
		/* partial sort, for m1 & m2 */
		rPsort(rw, n, m1);
		if (n % 2 == 0) {
			m2 = n - m1 - 1;
			rPsort(rw, n, m2);
			cmad = 3.*(rw[m1] + rw[m2]);
		}
		else { /* n odd */
			cmad = 6.*rw[m1];
		}

		if (cmad < 1e-7 * sc) /* effectively zero */
			break;
		c9 = 0.999*cmad;
		c1 = 0.001*cmad;
		for (i = 0; i < n; i++) {
			r = fabs(res[i]);
			if (r <= c1)
				rw[i] = 1.;
			else if (r <= c9)
				rw[i] = fsquare(1. - fsquare(r / cmad));
			else
				rw[i] = 0.;
		}
		iter++;
	}
    //------------------------------------------------
    //Release the resouces
    //cout << "Final block" << endl;
    vector<double> lowess_rw, lowess_res, lowess_ys;
    lowess_res.assign(restmp,restmp+n);
    lowess_rw.assign(rwtmp,rwtmp+n);
    lowess_ys.assign(ystmp,ystmp+n);

    ld.setM_lowess_res(lowess_res);
    ld.setM_lowess_rw(lowess_rw);
    ld.setM_lowess_ys(lowess_ys);

    // delete [] pointers
    cout << "delete pointers" << endl;
    cout << x << " " <<  xtmp << endl;
    delete [] xtmp;
    delete [] ytmp;
    delete [] ystmp;
    //cout << "First half" << endl;
    delete [] rwtmp;
    delete [] restmp;
    cout << "Done" << endl;

    //-----------------------------------------------
}

void LOWESS::lowess(int n,double f, int nsteps, double delta,LOWESSData & ld)
{
	cout << "[Info] Running LOWESS method... " << endl;
	clowess( n, f, nsteps, delta, ld);
}
