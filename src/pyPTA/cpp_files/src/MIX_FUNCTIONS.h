#pragma once

#include "PENG_ROBINSON.h"
#include "ADSORPTION_POTENTIAL_FUNCTIONS.h"
#include "GLOBAL_FUNCTIONS.h"


#include <stdio.h>
#include <iostream>
using namespace std;
#include <vector>
#include <numeric>
#include <cmath>
#include <cfloat>

const double tol = 1e-10;

struct eq { double Pz;	vector<double> x; };



struct eq calcEqProps(double Pz_old, vector<double> xold, vector<double> fb,
                      vector<double> f_eps, double T, vector<double> Pc,
                      vector<double> Tc, vector<double> Zc, vector<double> c_1,
                      vector<double> w);


double MPTAdfxMIX(vector<double> fb,
                  vector<double> f_eps,
                  vector<double> phiz,
                  double Pz);

vector<double> CALCnewX(vector<double> fb,
                        vector<double> f_eps,
                        vector<double> phiz,
                        double Pz);

double MPTAfxMIX(vector<double> fb,
                 vector<double> f_eps,
                 vector<double> phiz,
                 double Pz);

double MPTAMIX_CALCPZ(vector<double> x, double Pz, vector<double> fb,
                      vector<double> f_eps, double T, vector<double> Pc,
                      vector<double> Tc, vector<double> Zc, vector<double> c_1,
                      vector<double> w);
double MPTAd2fxMIX(vector<double> fb, vector<double> f_eps, vector<double> phiz, double Pz);


vector<double> M_LEE_NCALC(vector<double> eps_k, double L, double A,
                            size_t Ncamadas,
                            double Pbulk,
							double T,
							vector<double> ybulk,
                            vector<double> Pc,
							vector<double> Tc,
							vector<double> Zc,
                            vector<double> c_1,
							vector<double> w,
							vector<double> sigma_ff,
                   			double sigma_ss,
                   			double rho_atomos){
	size_t i;
	size_t j;
	size_t ncomp = eps_k.size();
	double eps;
	double PSI_z, PSI_lz;
	vector<double> x = ybulk;



	struct PRMIX results = calcmix(ybulk, Pbulk, T, Pc, Tc, Zc, c_1, w);
	vector<double> fb = results.fug;
	double densbulk = results.dens;

	double Pz = Pbulk;
	double minint = 0.7;
	double minSigma = minVal(sigma_ff);
	double minimal_space = minint * minSigma;
	vector<double> z = linespace(L / 2.0, minimal_space, Ncamadas);
	double delta = z[1] - z[0];
	
	double T_boltz = 1.38064852e-23 * T;

	vector<double> f_eps(ncomp, 0.0);
	vector<double> n_ads(ncomp, 0.0);
	vector<vector<double>> integral(ncomp, vector<double>(Ncamadas, 0.0));

	for (i = 0; i < Ncamadas; i++) {
		for (j = 0; j < ncomp; j++) {
			PSI_z = LEE(z[i], eps_k[j], sigma_ff[j], sigma_ss, rho_atomos);
			PSI_lz = LEE(L - z[i], eps_k[j], sigma_ff[j], sigma_ss, rho_atomos);
			eps = (PSI_z + PSI_lz)/T_boltz;
			f_eps[j] = exp(-eps);
		}


		struct eq EQPROPS = calcEqProps(Pz, x, fb, f_eps, T, Pc, Tc, Zc, c_1, w);
		Pz = EQPROPS.Pz;
		x = EQPROPS.x;



		struct PRMIX RES = calcmix(x, Pz, T, Pc, Tc, Zc, c_1, w);
		//Atualiza o vetor da quantidade adsorvida
		for (j = 0; j < ncomp; j++) {
			integral[j][i] = (RES.dens*x[j] - densbulk*ybulk[j]); // mol/L
		//printf("%f %f %f %f %f %f %f\n", z[i],Pz, RES.dens, x[0], x[1], x[2], x[3]);
		}
	}

	for (j = 0; j < ncomp; j++) {
		n_ads[j] = - A * SIMPS(integral[j], delta, Ncamadas) * 1.0e-7;
	}
	//printf("%f %f %f", n_ads[0], n_ads[1], n_ads[2]);

	return n_ads;
}


vector<double> M_STEELE_NCALC(vector<double> eps_k, double L, double A,
                            size_t Ncamadas,
                            double Pbulk,
							double T,
							vector<double> ybulk,
                            vector<double> Pc,
							vector<double> Tc,
							vector<double> Zc,
                            vector<double> c_1,
							vector<double> w,
							vector<double> sigma_ff,
                   			double sigma_ss){
	size_t i;
	size_t j;
	size_t ncomp = eps_k.size();
	double eps;
	double PSI_z, PSI_lz;
	
	vector<double> x = ybulk;



	struct PRMIX results = calcmix(ybulk, Pbulk, T, Pc, Tc, Zc, c_1, w);
	vector<double> fb = results.fug;
	double densbulk = results.dens;

	double Pz = Pbulk;
	double minint = 0.7;
	double minSigma = minVal(sigma_ff);
	double minimal_space = minint * minSigma;
	vector<double> z = linespace(L / 2.0, minimal_space, Ncamadas);
	double delta = z[1] - z[0];
	//cout << eps_k << " " << L << " " << A << " " << sigma_ff << " " << sigma_ss << endl;
	double T_boltz = 1.38064852e-23 * T;

	vector<double> f_eps(ncomp, 0.0);
	vector<double> n_ads(ncomp, 0.0);
	vector<vector<double>> integral(ncomp, vector<double>(Ncamadas, 0.0));

	for (i = 0; i < Ncamadas; i++) {
		for (j = 0; j < ncomp; j++) {
			PSI_z = STEELE(z[i], eps_k[j], sigma_ff[j], sigma_ss);
			PSI_lz = STEELE(L - z[i], eps_k[j], sigma_ff[j], sigma_ss);
			eps = (PSI_z + PSI_lz)/T_boltz;
			f_eps[j] = exp(-eps);
		}


		struct eq EQPROPS = calcEqProps(Pz, x, fb, f_eps, T, Pc, Tc, Zc, c_1, w);
		Pz = EQPROPS.Pz;
		x = EQPROPS.x;



		struct PRMIX RES = calcmix(x, Pz, T, Pc, Tc, Zc, c_1, w);
		//Atualiza o vetor da quantidade adsorvida
		for (j = 0; j < ncomp; j++) {
			integral[j][i] = (RES.dens*x[j] - densbulk*ybulk[j]); // mol/L
		//printf("%f %f %f %f %f %f %f\n", z[i],Pz, RES.dens, x[0], x[1], x[2], x[3]);
		}
	}

	for (j = 0; j < ncomp; j++) {
		n_ads[j] = - A * SIMPS(integral[j], delta, Ncamadas) * 1.0e-7;
	}
	//printf("%f %f %f", n_ads[0], n_ads[1], n_ads[2]);

	return n_ads;
}


vector<double> M_DRA_NCALC(vector<double> eps0, double z0, double beta,
                            size_t Ncamadas,
                            double Pbulk, double T, vector<double> ybulk,
                            vector<double> Pc, vector<double> Tc, vector<double> Zc,
                            vector<double> c_1, vector<double> w)
{
	size_t i;
	size_t j;
	size_t ncomp = eps0.size();
	double eps;
	double Pz = Pbulk;
	vector<double> x = ybulk;

	struct PRMIX results = calcmix(ybulk, Pbulk, T, Pc, Tc, Zc, c_1, w);
	vector<double> fb = results.fug;
	double densbulk = results.dens;

	double delta = z0 / (double)Ncamadas;

	vector<double> z = linespace(z0, delta, Ncamadas);
	vector<double> f_eps(ncomp, 0.0);
	vector<double> n_ads(ncomp, 0.0);
	vector<vector<double>> integral(ncomp, vector<double>(Ncamadas, 0.0));

	for (i = 0; i < Ncamadas; i++) {
		for (j = 0; j < ncomp; j++) {
			eps = DRA(eps0[j], z0, beta, z[i]);
			f_eps[j] = exp(eps / 8.314 / T);
		}


		struct eq EQPROPS = calcEqProps(Pz, x, fb, f_eps, T, Pc, Tc, Zc, c_1, w);
		Pz = EQPROPS.Pz;
		x = EQPROPS.x;



		struct PRMIX RES = calcmix(x, Pz, T, Pc, Tc, Zc, c_1, w);
		//Atualiza o vetor da quantidade adsorvida
		for (j = 0; j < ncomp; j++) {
			integral[j][i] = (RES.dens*x[j] - densbulk*ybulk[j]); // mol/L
		//printf("%f %f %f %f %f %f %f\n", z[i],Pz, RES.dens, x[0], x[1], x[2], x[3]);
		}
	}

	for (j = 0; j < ncomp; j++) {
		n_ads[j] = SIMPS(integral[j], delta, Ncamadas)/1000.0;
	}
	//printf("%f %f %f", n_ads[0], n_ads[1], n_ads[2]);

	return n_ads;

}


double calcPz_brent_Mix(double Pz,
					vector<double> xnew,
                    vector<double> fb,
                    vector<double> f_eps,
                    double T,
                    vector<double> Pc,
                    vector<double> Tc,
                    vector<double> Zc,
                    vector<double> c_1,
                    vector<double> w)
{
	double sumFb = 0.0;
	struct PRMIX RESmin = calcmix(xnew, Pz, T, Pc, Tc, Zc, c_1, w);
	for (size_t i = 0; i < xnew.size(); i++){
		sumFb +=  fb[i] * f_eps[i] - RESmin.fug[i];
	}
	return sumFb;
}

struct eq calcEqProps(double Pz_old, vector<double> xold, vector<double> fb,
                      vector<double> f_eps, double T, vector<double> Pc,
                      vector<double> Tc, vector<double> Zc, vector<double> c_1,
                      vector<double> w)
{
	double xsum, xsum2;
	double Pz = Pz_old;
	size_t i, k;
	vector<double> xnew = xold;
	struct PRMIX RES = calcmix(xnew, Pz, T, Pc, Tc, Zc, c_1, w);
	xnew = CALCnewX(fb, f_eps, RES.phi, Pz);
	k = 0;
	do {
		// solve for P
		static auto MINFUNC = [&](double x) {
			return calcPz_brent_Mix(x, xnew, fb, f_eps, T, Pc, Tc, Zc, c_1, w);
		};
	    double (*minfun_pointer)(double) = [](double x) { return MINFUNC(x); };
		Pz = brent_zeroin(minfun_pointer, Pz, tol);
		struct PRMIX RES = calcmix(xnew, Pz, T, Pc, Tc, Zc, c_1, w);
		// cout << Pz << endl;
		
		xsum = 0;
		for ( i = 0; i < fb.size(); i++) {
			xsum += fabs(xnew[i]/xold[i]) - 1.0;
			// xsum2 += xnew[i];
		}
		xsum *= xsum;
		// xsum2 -= 1.0;
		if (k > 25){
			break;
		}
		xold = xnew;
		k++;
		// Update composition
		xnew = CALCnewX(fb, f_eps, RES.phi, Pz);
		// cout << k << " " << xsum << endl;
	} while (xsum >= tol);

	struct eq EQ_DATA = { Pz, xnew };
	return EQ_DATA;
}

double MPTAMIX_CALCPZ(vector<double> x, double Pz, vector<double> fb,
                      vector<double> f_eps, double T, vector<double> Pc,
                      vector<double> Tc, vector<double> Zc, vector<double> c_1,
                      vector<double> w)
{
	double fx, fx_prime, fx_prime_prime, Pz_novo;
	size_t iter = 0;
	double condition, Pz_dash, fx_prime_dash;
	do {
		struct PRMIX RES = calcmix(x, Pz, T, Pc, Tc, Zc, c_1, w);

		fx = MPTAfxMIX(fb, f_eps, RES.phi, Pz);
		fx_prime = MPTAdfxMIX(fb, f_eps, RES.phi, Pz);
		fx_prime_prime = MPTAd2fxMIX(fb, f_eps, RES.phi, Pz);

		//Atualiza��es de Pz
		//Pz_novo = Pz - fx / fx_prime;
		//Pz_novo = Pz - fx * fx_prime / (fx_prime * fx_prime - 0.5 * fx * fx_prime_prime); // Halley's method
		//Pz_novo = Pz - fx / fx_prime - fx_prime_prime * fx * fx / (2 * fx_prime * fx_prime * fx_prime); // Schroeder's Method
		
		Pz_dash = Pz - fx / fx_prime;
		fx_prime_dash = MPTAdfxMIX(fb, f_eps, RES.phi, (Pz_dash+Pz) / 2.0);
		Pz_novo = Pz - fx / fx_prime_dash;
		

		condition = fabs(Pz_novo - Pz); // fabs(Pz_novo - Pz); 

		//|| (Pz_novo <= 0)
		if (fx <= tol) {
			return Pz_novo;
		}
			/*
			if (Pz >= Pz_novo) {
				return Pz;
			}
			else {
				return Pz_novo;
			}*/
		Pz = Pz_novo;
		iter++;
	} while (iter < 25);

	return Pz;
}


vector<double> CALCnewX(vector<double> fb,
                        vector<double> f_eps,
                        vector<double> phiz,
                        double Pz)
{
	size_t n = f_eps.size();
	vector<double> xnew(n, 0.0);
	double sum = 0.0;
	size_t i;
	for (i = 0; i < n; i++) {
		xnew[i] = fb[i] * f_eps[i] / (phiz[i] * Pz);
		sum += xnew[i];
	}
	for (i = 0; i < n; i++) {
		xnew[i] = xnew[i]/sum;
	}
	return xnew;
}

double MPTAfxMIX(vector<double> fb,
                 vector<double> f_eps,
                 vector<double> phiz,
                 double Pz)
{
	double sum = 0.0;
	for (size_t i = 0; i < f_eps.size(); i++) {
		sum+= (fb[i] * f_eps[i]) / (phiz[i] * Pz);
	}
	return sum - 1;
}

double MPTAdfxMIX(vector<double> fb,
                  vector<double> f_eps,
                  vector<double> phiz,
                  double Pz)
{
	double sum = 0.0;
	for (size_t i = 0; i < f_eps.size(); i++) {
		sum += (fb[i] * f_eps[i]) / (phiz[i] * Pz * Pz);
	}

	return -sum;
}


double MPTAd2fxMIX(vector<double> fb, vector<double> f_eps, vector<double> phiz, double Pz) {
	double sum = 0.0;
 	for (size_t i = 0; i < f_eps.size(); i++) {
 		sum += (2.0*fb[i] * f_eps[i]) / (phiz[i] * Pz * Pz * Pz);
 	}

 	return sum;
}


