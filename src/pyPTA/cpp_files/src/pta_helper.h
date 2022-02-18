
struct adsorbent
{
	double sigma_ss, rho_atoms;
};

double
MPTAfx(double fb, double f_eps, double phiz, double Pz)
{
	return fb * f_eps / (phiz * Pz) - 1.0;
}

double MPTAdfx(double fb, double f_eps, double phiz, double Pz)
{
	return -fb * f_eps / (phiz * Pz * Pz);
}

double MPTAd2fx(double fb, double f_eps, double phiz, double Pz)
{
	return 2 * fb * f_eps / (phiz * Pz * Pz * Pz);
}
