
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

double lamCDM_ns(){
	return ((double)rand()/(double)RAND_MAX)*0.12 + 0.9;
};

double lamCDM_As(){
	return ((double)rand()/(double)RAND_MAX)*0.2 + 3.0;
};

double SF2_ns(double phistar){
	return 1-(8.00/(phistar*phistar));
};

double SF2_As(double M, double phistar){
	return 23.0259 + log(0.00105542899*pow(10,2.00*M)*pow(phistar,4.00));
};

double SF1_ns(double phistar){
	return 1-(3.00/(phistar*phistar));
};

double SF1_As(double lam, double phistar){
	return 23.0259 + log(0.00844343197*pow(10,lam)*pow(phistar,3.00));
};

double mCM2_As(double M, double sigma, double rdec){
	double a=0.00187631821*rdec*rdec*176.0*pow(10,2.00*(M-sigma));
	double b=0.00105542899*pow(10,2.00*M)*30976.0;
	return 23.0259 + log(a+b);
};

double pCM2_As(double M, double sigma, double rdec){
	double a=0.00187631821*rdec*rdec*176.0*pow(10,2.00*(M-sigma));
	double b=0.00105542899*pow(10,2.00*M)*30976.0;
	return 23.0259 + log(a+b);
};

double dCM2_As(double M, double sigma, double rdec){
	double a=0.00187631821*rdec*rdec*176.0*pow(10,2.00*(M-sigma));
	double b=0.00105542899*pow(10,2.00*M)*30976.0;
	return 23.0259 + log(a+b);
};

double mCM2_As_phi(double M){
	double b=0.00105542899*pow(10,2.00*M)*30976.0;
	return b;
};

double pCM2_As_phi(double M){
	double b=0.00105542899*pow(10,2.00*M)*30976.0;
	return b;
};

double dCM2_As_phi(double M){
	double b=0.00105542899*pow(10,2.00*M)*30976.0;
	return b;
};

double mCM2_As_sigma(double M, double sigma, double rdec){
	double a=0.00187631821*rdec*rdec*176.0*pow(10,2.00*(M-sigma));
	return a;
};

double pCM2_As_sigma(double M, double sigma, double rdec){
	double a=0.00187631821*rdec*rdec*176.0*pow(10,2.00*(M-sigma));
	return a;
};

double dCM2_As_sigma(double M, double sigma, double rdec){
	double a=0.00187631821*rdec*rdec*176.0*pow(10,2.00*(M-sigma));
	return a;
};

double mCM2_ns(double Asphi, double Assig, double M, double m, double phistar){
	double As = Asphi + Assig;
	double a = 1+(Assig/As)*(-(4.00/(phistar*phistar))+(4.00/(phistar*phistar))*(pow(10,(2.00*(m-M)))));
	double b = (Asphi/As)*(-8.00/(phistar*phistar));
	return a+b;
};

double pCM2_ns(double Asphi, double Assig, double M, double m, double phistar){
	double As = Asphi + Assig;
	double a = 1+(Assig/As)*(-(4.00/(phistar*phistar))+(4.00/(phistar*phistar))*(pow(10,(2.00*(m-M)))));
	double b = (Asphi/As)*(-8.00/(phistar*phistar));
	return a+b;
};

double dCM2_ns(double Asphi, double Assig, double M, double m, double phistar){
	double As = Asphi + Assig;
	double a = 1+(Assig/As)*(-(4.00/(phistar*phistar))+(4.00/(phistar*phistar))*(pow(10,(2.00*(m-M)))));
	double b = (Asphi/As)*(-8.00/(phistar*phistar));
	return a+b;
};

double SF2_r(double phistar){
	return (32.00/(phistar*phistar));
};

double SF1_r(double phistar){
	return (8.00/(phistar*phistar));
};

double mCM2_r(double Asphi, double Assig, double phi){
	double As = Asphi + Assig;
	return (Asphi/As)*(32.00/(phi*phi));
};

double pCM2_r(double Asphi, double Assig, double phi){
	double As = Asphi + Assig;
	return (Asphi/As)*(32.00/(phi*phi));
};

double dCM2_r(double Asphi, double Assig, double phi){
	double As = Asphi + Assig;
	return (Asphi/As)*(32.00/(phi*phi));
};

double mCM2_rdec(double m, double sig, double gam){
	return pow(1.0+(8.00*pow(10,(0.5*(gam-m))-(2.00*sig))),-1.00);
};

double pCM2_rdec(double m, double sig, double gam){
	return pow(1.0+(8.00*pow(10,(0.5*(gam-m))-(2.00*sig))),-1.00);
};

double dCM2_rdec(double m, double sig, double gam){
	return pow(1.0+(8.00*pow(10,(0.5*(gam-m))-(2.00*sig))),-1.00);
};

double mCM2_fnl(double Assig, double rdec, double Asphi){
	double As = Asphi + Assig;
	return (5.00/12.00)*pow(Assig/As,2.00)*((3.00/(rdec))-4.00-(2.00*rdec));
}; 

double pCM2_fnl(double Assig, double rdec, double Asphi){
	double As = Asphi + Assig;
	return (5.00/12.00)*pow(Assig/As,2.00)*((3.00/(rdec))-4.00-(2.00*rdec));
};  

double dCM2_fnl(double Assig, double rdec, double Asphi){
	double As = Asphi + Assig;
	return (5.00/12.00)*pow(Assig/As,2.00)*((3.00/(rdec))-4.00-(2.00*rdec));
};  

double mCM2_Nmat(double gam, double m, double sig){
	double a=-2.389 + 0.66666666666*log(pow(10,(4.0*sig)+m-gam));
	if(a<0.00){
		return 0.00;
	}
	else{
		return a;
	}
};  

double pCM2_Nmat(double gam, double m, double sig){
	double a=-2.389 + 0.66666666666*log(pow(10,(4.0*sig)+m-gam));
	if(a<0.00){
		return 0.00;
	}
	else{
		return a;
	}
}; 

double dCM2_Nmat(double gam, double m, double sig){
	double a=-2.389 + 0.66666666666*log(pow(10,(4.0*sig)+m-gam));
	if(a<0.00){
		return 0.00;
	}
	else{
		return a;
	}
};

double mCM2_phistar(double Nmat, double Asphi, double Assig){
	double As = Asphi + Assig;
	return sqrt(4.00*(58.00+0.25*log(Asphi/As)-0.25*Nmat));
};  

double pCM2_phistar(double Nmat, double Asphi, double Assig){
	double As = Asphi + Assig;
	return sqrt(4.00*(58.00+0.25*log(Asphi/As)-0.25*Nmat));
}; 

double dCM2_phistar(double Nmat, double Asphi, double Assig){
	double As = Asphi + Assig;
	return sqrt(4.00*(58.00+0.25*log(Asphi/As)-0.25*Nmat));
};  

double SF2_phistar(double Nmat){
	return sqrt(4.00*(58.00-0.25*Nmat));
};     

double SF1_phistar(double Nmat){
	return sqrt(2.00*(54.00));
}; 

double dist_fnl(double fnl){
	return exp(-0.02*pow((fnl-0.8),2.00));
};

double vary_dist_fnl_1(double fnl,double sdmult){
	return exp(-(0.00428669*pow(sdmult,2.00))*pow((fnl-10.8),2.00));
};

double vary_dist_fnl_2(double fnl,double sdmult){
	return exp(-(0.32*pow(sdmult,2.00))*pow((fnl+1.25),2.00));
};

double mCM2_Gammasigmasampler(double m){
	double a=m;
	return ((double)rand()/(double)RAND_MAX)*(-39.00-a) + a;
};

double pCM2_Gammasigmasampler(double m){
	double a=m;
	return ((double)rand()/(double)RAND_MAX)*(-39.00-a) + a;
};

double dCM2_Gammasigmasampler(double m, double sig){
	double a=(-4.0+m+(4.0*sig));
	return ((double)rand()/(double)RAND_MAX)*(-39.00-a) + a;
};

double SF2_Msampler(){  
	return ((double)rand()/(double)RAND_MAX)*(-2.75) - 4.00;
};

double SF1_lamsampler(){  
	return ((double)rand()/(double)RAND_MAX)*(-6.00) - 7.00;
};

double mCM2_Msampler(){  
	return ((double)rand()/(double)RAND_MAX)*(-11.3333333) - 4.00;
};

double pCM2_Msampler(){  
	return ((double)rand()/(double)RAND_MAX)*(-11.3333333) - 4.00;
};

double dCM2_Msampler(){  
	return ((double)rand()/(double)RAND_MAX)*(-11.3333333) - 4.00;
};

double mCM2_msampler(double M){  
	double a=M-0.30;
	double b=-39.0;
	return ((double)rand()/(double)RAND_MAX)*(b-a) + a;
};

double pCM2_msampler(double M){  
	double a=M-0.30;
	double b=-39.0;
	return ((double)rand()/(double)RAND_MAX)*(b-a) + a;
};

double dCM2_msampler(double M){  
	double a=M-0.30;
	double b=-39.0;
	return ((double)rand()/(double)RAND_MAX)*(b-a) + a;
};

double mCM2_sigmastarsampler_gaussian(double M, double m){
	double u1=((double)rand()/(double)RAND_MAX)*(1.00);
	double u2=((double)rand()/(double)RAND_MAX)*(1.00);
	double rand_norm=sqrt(-2.00*log(u1))*cos(6.28318530718*u2);
	return log10(std::abs(rand_norm*(2.56509966032*pow(10,(2.00*M)-m)*188.0)));
};

double pCM2_sigmastarsampler_gaussian(double M, double m){
	double u1=((double)rand()/(double)RAND_MAX)*(1.00);
	double u2=((double)rand()/(double)RAND_MAX)*(1.00);
	double rand_norm=sqrt(-2.00*log(u1))*cos(6.28318530718*u2);
	return log10(std::abs(rand_norm*(2.56509966032*pow(10,(2.00*M)-m)*188.0)));
};

double dCM2_sigmastarsampler_gaussian(double M, double m){
	double u1=((double)rand()/(double)RAND_MAX)*(1.00);
	double u2=((double)rand()/(double)RAND_MAX)*(1.00);
	double rand_norm=sqrt(-2.00*log(u1))*cos(6.28318530718*u2);
	return log10(std::abs(rand_norm*(2.56509966032*pow(10,(2.00*M)-m)*188.0)));
};

double mCM2_sigmastarsampler(){
	return log10(((double)rand()/(double)RAND_MAX)*0.01);
};

double pCM2_sigmastarsampler(){
	return log10(((double)rand()/(double)RAND_MAX)*0.01);
};

double dCM2_sigmastarsampler(){
	return log10(((double)rand()/(double)RAND_MAX)*0.01);
};

int main()
{
	static const int gridsize=100;
	static const int pres=20;
    static const int smoothing_order=3; 
	double*** r_ns_As_count;
	int maxsam=1000000;  
	double sdmult=3.5;
	double bestfit_ns;
	double bestfit_As;
	double samplens,sampler;
	double SF2_fnlEvidence=0.00;
	double SF1_fnlEvidence=0.00;
	double mCM2_fnlEvidence=0.00;
	double pCM2_fnlEvidence=0.00;
	double dCM2_fnlEvidence=0.00;
	double lamCDM_fnlEvidence=0.00;
	double SF2_fnlmaxlik=0.00;
	double SF1_fnlmaxlik=0.00;
	double mCM2_fnlmaxlik=0.00;
	double pCM2_fnlmaxlik=0.00;
	double dCM2_fnlmaxlik=0.00;
	double lamCDM_fnlmaxlik=0.00;
	double sample_fnlmaxlik=0.00;
	static const double reducedmaxAs=3.4;
	static const double reducedminAs=2.8;
	static const double maxAs=6.91;
	static const double minAs=-2.30;
	static const double maxns=1.02;
	static const double minns=0.90;
	static const double maxr=0.5;
	static const double minr=0.0;
	static const double maxfnl=30.0;
	static const double minfnl=-20.0;
	static const double stepsizeAs=(reducedmaxAs-reducedminAs)/static_cast<double>(gridsize);
	static const double stepsizens=(maxns-minns)/static_cast<double>(gridsize);
	static const double stepsizer=(maxr-minr)/static_cast<double>(gridsize);
	static const double stepsizefnl=(maxfnl-minfnl)/static_cast<double>(gridsize);
	double ns_datum,r_datum,weight_datum,As_datum;
	double a1,a2,a3,a4;
	double b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13;
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13;
	double d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13;
	int ticker=0;
	int date;
	srand(time(NULL));  

	std::cout << "Initialising...";    

	r_ns_As_count=new double**[gridsize];
	for(int r = 0; r < gridsize; ++r) {
		r_ns_As_count[r]=new double*[gridsize];
		for(int ns = 0; ns < gridsize; ++ns) {
			r_ns_As_count[r][ns]=new double[gridsize];
			for(int As = 0; As < gridsize; ++As) { 
				r_ns_As_count[r][ns][As]=0;
			}
		}
	}

	std::cout << "done" << "\n";
	
	std::cout << "Which data? Please enter 2015, 2013 or 0 (BKP): ";
	std::cin >> date;
	std::cout << "\n";
	
	if(date==0){
		std::cout << "Reading 'planck_ns_dataBKP.dat', 'planck_r_dataBKP.dat', 'planck_ln1010Pstar_dataBKP.dat' and 'planck_weight_dataBKP.dat' files...";
	
		std::fstream ns_file("planck_ns_dataBKP.dat", std::ios_base::in);
		std::fstream r_file("planck_r_dataBKP.dat", std::ios_base::in); 
		std::fstream As_file("planck_ln1010Pstar_dataBKP.dat", std::ios_base::in);
		std::fstream weight_file("planck_weight_dataBKP.dat", std::ios_base::in);  

		while(ns_file >> ns_datum && r_file >> r_datum && weight_file >> weight_datum && As_file >> As_datum){
			if(ns_datum<maxns && ns_datum>minns && r_datum<maxr && r_datum>minr && As_datum<reducedmaxAs && As_datum>reducedminAs){
				int i=static_cast<int>(trunc((ns_datum-minns)/stepsizens));
				int j=static_cast<int>(trunc((r_datum-minr)/stepsizer));
				int k=static_cast<int>(trunc((As_datum-reducedminAs)/stepsizeAs));
				r_ns_As_count[j][i][k]+=weight_datum;
			}
		}
		bestfit_ns=0.9663133;
		bestfit_As=3.088923;
	}
	
	if(date==2015){
		std::cout << "Reading 'planck_ns_data.dat', 'planck_r_data.dat', 'planck_ln1010Pstar_data.dat' and 'planck_weight_data.dat' files...";
	
		std::fstream ns_file("planck_ns_data.dat", std::ios_base::in);
		std::fstream r_file("planck_r_data.dat", std::ios_base::in); 
		std::fstream As_file("planck_ln1010Pstar_data.dat", std::ios_base::in);
		std::fstream weight_file("planck_weight_data.dat", std::ios_base::in);  

		while(ns_file >> ns_datum && r_file >> r_datum && weight_file >> weight_datum && As_file >> As_datum){
			if(ns_datum<maxns && ns_datum>minns && r_datum<maxr && r_datum>minr && As_datum<reducedmaxAs && As_datum>reducedminAs){
				int i=static_cast<int>(trunc((ns_datum-minns)/stepsizens));
				int j=static_cast<int>(trunc((r_datum-minr)/stepsizer));
				int k=static_cast<int>(trunc((As_datum-reducedminAs)/stepsizeAs));
				r_ns_As_count[j][i][k]+=weight_datum;
			}
		}
		bestfit_ns=0.9679726;
		bestfit_As=3.089686;
	}
	
	if(date==2013){
		std::cout << "Reading 'planck_ns_data2013.dat', 'planck_r_data2013.dat', 'planck_ln1010Pstar_data2013.dat' and 'planck_weight_data2013.dat' files...";
	
		std::fstream ns_file("planck_ns_data2013.dat", std::ios_base::in);
		std::fstream r_file("planck_r_data2013.dat", std::ios_base::in); 
		std::fstream As_file("planck_ln1010Pstar_data2013.dat", std::ios_base::in);
		std::fstream weight_file("planck_weight_data2013.dat", std::ios_base::in);  

		while(ns_file >> ns_datum && r_file >> r_datum && weight_file >> weight_datum && As_file >> As_datum){
			if(ns_datum<maxns && ns_datum>minns && r_datum<maxr && r_datum>minr && As_datum<reducedmaxAs && As_datum>reducedminAs){
				int i=static_cast<int>(trunc((ns_datum-minns)/stepsizens));
				int j=static_cast<int>(trunc((r_datum-minr)/stepsizer));
				int k=static_cast<int>(trunc((As_datum-reducedminAs)/stepsizeAs));
				r_ns_As_count[j][i][k]+=weight_datum;
			}
		}
		bestfit_ns=0.9624098;
		bestfit_As=3.087056;
	}
	
	if(date!=2013 && date!=2015 && date!=0){
		std::cout << "Input read failure...";
		std::cin.get();
	}
	
	std::cout << "\n" << "Input sample ns: ";
	std::cin >> samplens;
	std::cout << "\n" << "Input sample r: ";
	std::cin >> sampler;
	std::cout << "\n";
	
	std::cout << "done" << "\n";	

	std::cout << "Sampling and binning model parameters..."; 
	
	for(int itr=0;itr<maxsam;itr++){ 
		double SF2_tempNmat=0.00;
		double SF2_tempphistar=SF2_phistar(SF2_tempNmat);
		double SF2_tempMandm1=SF2_Msampler();
		double SF1_tempNmat=0.00;
		double SF1_tempphistar=SF1_phistar(SF1_tempNmat);
		double SF1_templam=SF1_lamsampler();
		double mCM2_tempMandm1=mCM2_Msampler();
		double pCM2_tempMandm1=pCM2_Msampler();
		double dCM2_tempMandm1=dCM2_Msampler();
		double mCM2_tempMandm2=mCM2_msampler(mCM2_tempMandm1);
		double pCM2_tempMandm2=pCM2_msampler(pCM2_tempMandm1);
		double dCM2_tempMandm2=dCM2_msampler(dCM2_tempMandm1);
		double mCM2_tempsigmastar=mCM2_sigmastarsampler();
		double pCM2_tempsigmastar=pCM2_sigmastarsampler();
		double dCM2_tempsigmastar=dCM2_sigmastarsampler();
		double mCM2_tempGammasigma=mCM2_Gammasigmasampler(mCM2_tempMandm2);
		double pCM2_tempGammasigma=pCM2_Gammasigmasampler(pCM2_tempMandm2);
		double dCM2_tempGammasigma=dCM2_Gammasigmasampler(dCM2_tempMandm2,dCM2_tempsigmastar);
		double mCM2_temprdec=mCM2_rdec(mCM2_tempMandm2,mCM2_tempsigmastar,mCM2_tempGammasigma);
		double pCM2_temprdec=pCM2_rdec(pCM2_tempMandm2,pCM2_tempsigmastar,pCM2_tempGammasigma);
		double dCM2_temprdec=dCM2_rdec(dCM2_tempMandm2,dCM2_tempsigmastar,dCM2_tempGammasigma);	
		double mCM2_tempAs=mCM2_As(mCM2_tempMandm1,mCM2_tempsigmastar,mCM2_temprdec);
		double pCM2_tempAs=pCM2_As(pCM2_tempMandm1,pCM2_tempsigmastar,pCM2_temprdec);	
		double dCM2_tempAs=dCM2_As(dCM2_tempMandm1,dCM2_tempsigmastar,dCM2_temprdec);
			
		while(mCM2_temprdec>1.00 || mCM2_tempsigmastar>-2.00 || mCM2_tempAs>maxAs || mCM2_tempAs<minAs){
			mCM2_tempMandm1=mCM2_Msampler();
			mCM2_tempMandm2=mCM2_msampler(mCM2_tempMandm1);
			mCM2_tempsigmastar=mCM2_sigmastarsampler();
			mCM2_tempGammasigma=mCM2_Gammasigmasampler(mCM2_tempMandm2);
			mCM2_temprdec=mCM2_rdec(mCM2_tempMandm2,mCM2_tempsigmastar,mCM2_tempGammasigma);
			mCM2_tempAs=mCM2_As(mCM2_tempMandm1,mCM2_tempsigmastar,mCM2_temprdec);
		}
		
		while(pCM2_temprdec>1.00 || pCM2_tempsigmastar>-2.00 || pCM2_tempAs>maxAs || pCM2_tempAs<minAs || pCM2_As_phi(pCM2_tempMandm1)>0.01*(pCM2_As_phi(pCM2_tempMandm1)+pCM2_As_sigma(pCM2_tempMandm1,pCM2_tempsigmastar,pCM2_temprdec))){
			pCM2_tempMandm1=pCM2_Msampler();
			pCM2_tempMandm2=pCM2_msampler(pCM2_tempMandm1);
			pCM2_tempsigmastar=pCM2_sigmastarsampler();
			pCM2_tempGammasigma=pCM2_Gammasigmasampler(pCM2_tempMandm2);
			pCM2_temprdec=pCM2_rdec(pCM2_tempMandm2,pCM2_tempsigmastar,pCM2_tempGammasigma);
			pCM2_tempAs=pCM2_As(pCM2_tempMandm1,pCM2_tempsigmastar,pCM2_temprdec);
		}
		
		while(dCM2_temprdec>1.00 || dCM2_tempsigmastar>-2.00 || dCM2_temprdec<0.95 || -4.0+(4.0*dCM2_tempsigmastar)+dCM2_tempMandm2 < -39.0 || dCM2_tempAs>maxAs || dCM2_tempAs<minAs || dCM2_As_phi(dCM2_tempMandm1)>0.01*(dCM2_As_phi(dCM2_tempMandm1)+dCM2_As_sigma(dCM2_tempMandm1,dCM2_tempsigmastar,dCM2_temprdec))){
			dCM2_tempMandm1=dCM2_Msampler();
			dCM2_tempMandm2=dCM2_msampler(dCM2_tempMandm1);
			dCM2_tempsigmastar=dCM2_sigmastarsampler();
			dCM2_tempGammasigma=dCM2_Gammasigmasampler(dCM2_tempMandm2,dCM2_tempsigmastar);
			dCM2_temprdec=dCM2_rdec(dCM2_tempMandm2,dCM2_tempsigmastar,dCM2_tempGammasigma);
			dCM2_tempAs=dCM2_As(dCM2_tempMandm1,dCM2_tempsigmastar,dCM2_temprdec);
		}
		
		double mCM2_tempAsphi=mCM2_As_phi(mCM2_tempMandm1);
		double pCM2_tempAsphi=pCM2_As_phi(pCM2_tempMandm1);	
		double dCM2_tempAsphi=dCM2_As_phi(dCM2_tempMandm1);
		double mCM2_tempAssig=mCM2_As_sigma(mCM2_tempMandm1,mCM2_tempsigmastar,mCM2_temprdec);
		double pCM2_tempAssig=pCM2_As_sigma(pCM2_tempMandm1,pCM2_tempsigmastar,pCM2_temprdec);	
		double dCM2_tempAssig=dCM2_As_sigma(dCM2_tempMandm1,dCM2_tempsigmastar,dCM2_temprdec);
		double mCM2_tempNmat=mCM2_Nmat(mCM2_tempGammasigma,mCM2_tempMandm2,mCM2_tempsigmastar);
		double pCM2_tempNmat=pCM2_Nmat(pCM2_tempGammasigma,pCM2_tempMandm2,pCM2_tempsigmastar);
		double dCM2_tempNmat=dCM2_Nmat(dCM2_tempGammasigma,dCM2_tempMandm2,dCM2_tempsigmastar);
		double mCM2_tempphistar=mCM2_phistar(mCM2_tempNmat,mCM2_tempAsphi,mCM2_tempAssig);
		double pCM2_tempphistar=pCM2_phistar(pCM2_tempNmat,pCM2_tempAsphi,pCM2_tempAssig);
		double dCM2_tempphistar=dCM2_phistar(dCM2_tempNmat,dCM2_tempAsphi,dCM2_tempAssig);
		double mCM2_tempns=mCM2_ns(mCM2_tempAsphi,mCM2_tempAssig,mCM2_tempMandm1,mCM2_tempMandm2,mCM2_tempphistar);
		double pCM2_tempns=pCM2_ns(pCM2_tempAsphi,pCM2_tempAssig,pCM2_tempMandm1,pCM2_tempMandm2,pCM2_tempphistar);
		double dCM2_tempns=dCM2_ns(dCM2_tempAsphi,dCM2_tempAssig,dCM2_tempMandm1,dCM2_tempMandm2,dCM2_tempphistar);
		double mCM2_tempr=mCM2_r(mCM2_tempAsphi,mCM2_tempAssig,mCM2_tempphistar);
		double pCM2_tempr=pCM2_r(pCM2_tempAsphi,pCM2_tempAssig,pCM2_tempphistar);
		double dCM2_tempr=dCM2_r(dCM2_tempAsphi,dCM2_tempAssig,dCM2_tempphistar);
		double mCM2_tempfnl=mCM2_fnl(mCM2_tempAssig,mCM2_temprdec,mCM2_tempAsphi);
		double pCM2_tempfnl=pCM2_fnl(pCM2_tempAssig,pCM2_temprdec,pCM2_tempAsphi);
		double dCM2_tempfnl=dCM2_fnl(dCM2_tempAssig,dCM2_temprdec,dCM2_tempAsphi);
		double SF2_tempns=SF2_ns(SF2_tempphistar);
		double SF2_tempr=SF2_r(SF2_tempphistar);
		double SF2_tempAs=SF2_As(SF2_tempMandm1,SF2_tempphistar);
		double SF1_tempns=SF1_ns(SF1_tempphistar);
		double SF1_tempr=SF1_r(SF1_tempphistar);
		double SF1_tempAs=SF1_As(SF1_templam,SF1_tempphistar);		
		double lamCDM_tempAs=lamCDM_As();
		double lamCDM_tempns=lamCDM_ns();
		
		
		if(mCM2_tempns<maxns && mCM2_tempns>minns && mCM2_tempr<maxr && mCM2_tempr>minr && mCM2_tempAs<reducedmaxAs && mCM2_tempAs>reducedminAs && mCM2_tempfnl<maxfnl && mCM2_tempfnl>minfnl){
			int a=static_cast<int>(trunc((mCM2_tempns-minns)/stepsizens));
			int b=static_cast<int>(trunc((mCM2_tempr-minr)/stepsizer));
			int c=static_cast<int>(trunc((mCM2_tempAs-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(mCM2_tempfnl)>0.00){mCM2_fnlEvidence+=sum_r_ns_As*dist_fnl(mCM2_tempfnl);}
			if(sum_r_ns_As*dist_fnl(mCM2_tempfnl)>mCM2_fnlmaxlik){
				mCM2_fnlmaxlik=sum_r_ns_As*dist_fnl(mCM2_tempfnl);
				b1=mCM2_tempns;
				b2=mCM2_tempr;
				b3=mCM2_tempAs;
				b4=mCM2_tempfnl;
				b5=mCM2_tempAsphi;
				b6=mCM2_tempAssig;
				b7=mCM2_tempMandm1;
				b8=mCM2_tempMandm2;
				b9=mCM2_tempsigmastar;
				b10=mCM2_tempphistar;
				b11=mCM2_tempGammasigma;
				b12=mCM2_tempNmat;
				b13=mCM2_temprdec;
			}
		}
		
		if(pCM2_tempns<maxns && pCM2_tempns>minns && pCM2_tempr<maxr && pCM2_tempr>minr && pCM2_tempAs<reducedmaxAs && pCM2_tempAs>reducedminAs && pCM2_tempfnl<maxfnl && pCM2_tempfnl>minfnl){
			int a=static_cast<int>(trunc((pCM2_tempns-minns)/stepsizens));
			int b=static_cast<int>(trunc((pCM2_tempr-minr)/stepsizer));
			int c=static_cast<int>(trunc((pCM2_tempAs-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(pCM2_tempfnl)>0.00){pCM2_fnlEvidence+=sum_r_ns_As*dist_fnl(pCM2_tempfnl);}
			if(sum_r_ns_As*dist_fnl(pCM2_tempfnl)>pCM2_fnlmaxlik){
				pCM2_fnlmaxlik=sum_r_ns_As*dist_fnl(pCM2_tempfnl);
				c1=pCM2_tempns;
				c2=pCM2_tempr;
				c3=pCM2_tempAs;
				c4=pCM2_tempfnl;
				c5=pCM2_tempAsphi;
				c6=pCM2_tempAssig;
				c7=pCM2_tempMandm1;
				c8=pCM2_tempMandm2;
				c9=pCM2_tempsigmastar;
				c10=pCM2_tempphistar;
				c11=pCM2_tempGammasigma;
				c12=pCM2_tempNmat;
				c13=pCM2_temprdec;
			}
		}
		
		if(dCM2_tempns<maxns && dCM2_tempns>minns && dCM2_tempr<maxr && dCM2_tempr>minr && dCM2_tempAs<reducedmaxAs && dCM2_tempAs>reducedminAs && dCM2_tempfnl<maxfnl && dCM2_tempfnl>minfnl){
			int a=static_cast<int>(trunc((dCM2_tempns-minns)/stepsizens));
			int b=static_cast<int>(trunc((dCM2_tempr-minr)/stepsizer));
			int c=static_cast<int>(trunc((dCM2_tempAs-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(dCM2_tempfnl)>0.00){dCM2_fnlEvidence+=sum_r_ns_As*dist_fnl(dCM2_tempfnl);}
			if(sum_r_ns_As*dist_fnl(dCM2_tempfnl)>dCM2_fnlmaxlik){
				dCM2_fnlmaxlik=sum_r_ns_As*dist_fnl(dCM2_tempfnl);
				d1=dCM2_tempns;
				d2=dCM2_tempr;
				d3=dCM2_tempAs;
				d4=dCM2_tempfnl;
				d5=dCM2_tempAsphi;
				d6=dCM2_tempAssig;
				d7=dCM2_tempMandm1;
				d8=dCM2_tempMandm2;
				d9=dCM2_tempsigmastar;
				d10=dCM2_tempphistar;
				d11=dCM2_tempGammasigma;
				d12=dCM2_tempNmat;
				d13=dCM2_temprdec;
			}
		}
		
		if(SF2_tempns<maxns && SF2_tempns>minns && SF2_tempr<maxr && SF2_tempr>minr && SF2_tempAs<reducedmaxAs && SF2_tempAs>reducedminAs){
			int a=static_cast<int>(trunc((SF2_tempns-minns)/stepsizens));
			int b=static_cast<int>(trunc((SF2_tempr-minr)/stepsizer));
			int c=static_cast<int>(trunc((SF2_tempAs-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(0.00)>0.00){SF2_fnlEvidence+=sum_r_ns_As*dist_fnl(0.00);}
			if(sum_r_ns_As*dist_fnl(0.00)>SF2_fnlmaxlik){
				SF2_fnlmaxlik=sum_r_ns_As*dist_fnl(0.00);
				a1=SF2_tempns;
				a2=SF2_tempr;
				a3=SF2_tempAs;
				a4=0.00;
			}
		}

		if(SF1_tempns<maxns && SF1_tempns>minns && SF1_tempr<maxr && SF1_tempr>minr && SF1_tempAs<reducedmaxAs && SF1_tempAs>reducedminAs){
			int a=static_cast<int>(trunc((SF1_tempns-minns)/stepsizens));
			int b=static_cast<int>(trunc((SF1_tempr-minr)/stepsizer));
			int c=static_cast<int>(trunc((SF1_tempAs-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(0.00)>0.00){SF1_fnlEvidence+=sum_r_ns_As*dist_fnl(0.00);}
			if(sum_r_ns_As*dist_fnl(0.00)>SF1_fnlmaxlik){SF1_fnlmaxlik=sum_r_ns_As*dist_fnl(0.00);}
		}
		
		if(lamCDM_tempns<maxns && lamCDM_tempns>minns && lamCDM_tempAs<reducedmaxAs && lamCDM_tempAs>reducedminAs){
			int a=static_cast<int>(trunc((lamCDM_tempns-minns)/stepsizens));
			int b=0;
			int c=static_cast<int>(trunc((lamCDM_tempAs-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(0.00)>0.00){lamCDM_fnlEvidence+=sum_r_ns_As*dist_fnl(0.00);}
			if(sum_r_ns_As*dist_fnl(0.00)>lamCDM_fnlmaxlik){lamCDM_fnlmaxlik=sum_r_ns_As*dist_fnl(0.00);}
		}
		
		if(samplens<maxns && samplens>minns && sampler>minr && sampler<maxr){
			int a=static_cast<int>(trunc((samplens-minns)/stepsizens));
			int b=static_cast<int>(trunc((sampler-minr)/stepsizer));
			int c=static_cast<int>(trunc((bestfit_As-reducedminAs)/stepsizeAs));
			int a_lower=0;
			int b_lower=0;
			int c_lower=0;
			int a_upper=0;
			int b_upper=0;
			int c_upper=0;
			double ha=1.0+exp(std::abs(bestfit_ns-static_cast<double>(a)));
			double hb=1.0+exp(static_cast<double>(b));
			double hc=1.0+exp(std::abs(bestfit_As-static_cast<double>(c)));
			if(a<smoothing_order){a_lower=smoothing_order-a;}
			if(b<smoothing_order){b_lower=smoothing_order-b;}
			if(c<smoothing_order){c_lower=smoothing_order-c;}
			if(a>gridsize-smoothing_order-1){a_upper=-smoothing_order+a+1-gridsize;}
			if(b>gridsize-smoothing_order-1){b_upper=-smoothing_order+b+1-gridsize;}
			if(c>gridsize-smoothing_order-1){c_upper=-smoothing_order+c+1-gridsize;}
			double sum_r_ns_As=0.0;
			double tot=0.0;
			for(int aa=-smoothing_order+a_lower;aa<smoothing_order+1+a_upper;aa++){
				for(int bb=-smoothing_order+b_lower;bb<smoothing_order+1+b_upper;bb++){
					for(int cc=-smoothing_order+c_lower;cc<smoothing_order+1+c_upper;cc++){
						sum_r_ns_As+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb))*r_ns_As_count[b+bb][a+aa][c+cc];
						tot+=sqrt(exp(-pow(a-aa,2.00)/ha)+exp(-pow(c-cc,2.00)/hc)+exp(-pow(b-bb,2.00)/hb));
					}
				}
			}
			sum_r_ns_As*=1.00/tot;
			if(sum_r_ns_As*dist_fnl(0.00)>sample_fnlmaxlik){sample_fnlmaxlik=sum_r_ns_As*dist_fnl(0.00);}
		}
		
		ticker++;
		std::cout << ticker << "\n";
	} 
	
	std::cout << "done" << "\n"; 

	std::cout << "Writing output to 'MCEvidenceSLW_data.dat'...";

	std::string MCEvidence_datafile = "MCEvidenceSLW_data.dat";
	std::ofstream MCEvidence_writeinfile(MCEvidence_datafile.c_str());   

	MCEvidence_writeinfile << std::fixed << std::setprecision(pres) << "n_s stepsize:" << stepsizens << "\r\n"; 
	MCEvidence_writeinfile << "As stepsize:" << stepsizeAs << "\r\n";
	MCEvidence_writeinfile << "r stepsize:" << stepsizer << "\r\n";
	MCEvidence_writeinfile << "f_NL stepsize:" << stepsizefnl << "\r\n"; 
	MCEvidence_writeinfile << "No. of samples:" << maxsam << "\r\n\r\n"; 
	MCEvidence_writeinfile << "SF2 fnl Evidence:" << SF2_fnlEvidence/static_cast<double>(maxsam)<< "\r\n";
	MCEvidence_writeinfile << "SF1 fnl Evidence:" << SF1_fnlEvidence/static_cast<double>(maxsam)<< "\r\n";
	MCEvidence_writeinfile << "pCM2 fnl Evidence:" << pCM2_fnlEvidence/static_cast<double>(maxsam) << "\r\n";
	MCEvidence_writeinfile << "dCM2 fnl Evidence:" << dCM2_fnlEvidence/static_cast<double>(maxsam) << "\r\n";
	MCEvidence_writeinfile << "lamCDM fnl Evidence:" << lamCDM_fnlEvidence/static_cast<double>(maxsam)  << "\r\n";
	MCEvidence_writeinfile << "mCM2 fnl Evidence:" << mCM2_fnlEvidence/static_cast<double>(maxsam) << "\r\n\r\n";
	MCEvidence_writeinfile << "ln(SF2 fnlEvidence / lamCDM fnlEvidence):" << log(SF2_fnlEvidence/(lamCDM_fnlEvidence)) << "\r\n";
	MCEvidence_writeinfile << "ln(SF1 fnlEvidence / lamCDM fnlEvidence):" << log(SF1_fnlEvidence/(lamCDM_fnlEvidence)) << "\r\n";
	MCEvidence_writeinfile << "ln(pCM2 fnlEvidence / lamCDM fnlEvidence):" << log(pCM2_fnlEvidence/(lamCDM_fnlEvidence)) << "\r\n";
	MCEvidence_writeinfile << "ln(dCM2 fnlEvidence / lamCDM fnlEvidence):" << log(dCM2_fnlEvidence/(lamCDM_fnlEvidence)) << "\r\n";
	MCEvidence_writeinfile << "ln(mCM2 fnlEvidence / lamCDM fnlEvidence):" << log(mCM2_fnlEvidence/(lamCDM_fnlEvidence)) << "\r\n\r\n";
	MCEvidence_writeinfile << "ln(pCM2 fnlEvidence / SF2 fnlEvidence):" << log(pCM2_fnlEvidence/(SF2_fnlEvidence)) << "\r\n";
	MCEvidence_writeinfile << "ln(dCM2 fnlEvidence / SF2 fnlEvidence):" << log(dCM2_fnlEvidence/(SF2_fnlEvidence)) << "\r\n";
	MCEvidence_writeinfile << "ln(mCM2 fnlEvidence / SF2 fnlEvidence):" << log(mCM2_fnlEvidence/(SF2_fnlEvidence)) << "\r\n\r\n";
	MCEvidence_writeinfile << "SAMPLE chisq :" << -2.0*log(sample_fnlmaxlik/(lamCDM_fnlmaxlik)) << "\r\n";
	MCEvidence_writeinfile << "SF2: chisq | ns | r | As | fnl | :" << -2.0*log(SF2_fnlmaxlik/(lamCDM_fnlmaxlik)) << "| \t" << a1 << "| \t" << a2 << "| \t" << a3 << "| \t" << a4 <<"\r\n";
	MCEvidence_writeinfile << "SF1: chisq :" << -2.0*log(SF1_fnlmaxlik/(lamCDM_fnlmaxlik)) << "\r\n";
	MCEvidence_writeinfile << "pCM2: chisq | ns | r | As | fnl | Asphi | Assig | M | m | Sig* | Phi* | Gam* | Nmat | rdec | :" << -2.0*log(pCM2_fnlmaxlik/(lamCDM_fnlmaxlik)) << "| \t" << c1 << "| \t" << c2 << "| \t" << c3 << "| \t" << c4 << "| \t" << c5 << "| \t" << c6 << "| \t" << c7 << "| \t" << c8 << "| \t" << c9 << "| \t" << c10 << "| \t" << c11 << "| \t" << c12 << "| \t" << c13 << "\r\n";
	MCEvidence_writeinfile << "dCM2: chisq | ns | r | As | fnl | Asphi | Assig | M | m | Sig* | Phi* | Gam* | Nmat | rdec | :" << -2.0*log(dCM2_fnlmaxlik/(lamCDM_fnlmaxlik)) << "| \t" << d1 << "| \t" << d2 << "| \t" << d3 << "| \t" << d4 << "| \t" << d5 << "| \t" << d6 << "| \t" << d7 << "| \t" << d8 << "| \t" << d9 << "| \t" << d10 << "| \t" << d11 << "| \t" << d12 << "| \t" << d13 << "\r\n";
	MCEvidence_writeinfile << "mCM2: chisq | ns | r | As | fnl | Asphi | Assig | M | m | Sig* | Phi* | Gam* | Nmat | rdec | :" << -2.0*log(mCM2_fnlmaxlik/(lamCDM_fnlmaxlik)) << "| \t" << b1 << "| \t" << b2 << "| \t" << b3 << "| \t" << b4 << "| \t" << b5 << "| \t" << b6 << "| \t" << b7 << "| \t" << b8 << "| \t" << b9 << "| \t" << b10 << "| \t" << b11 << "| \t" << b12 << "| \t" << b13 << "\r\n\r\n";
	
	
	std::cout << "done" << "\n";
	
	MCEvidence_writeinfile.close(); 
	
	std::cout << "All done" << "\n";
	
	std::cout << "EXIT" << "\n";
	std::cin.get();


	return 0;
}
