
//log
//2023-02-06 created 
//2023-02-10 finished but fonction bizzare: bug in initialization ,update or definition
//2023-02-14 debug : pay attention to the demension of array !!!!!!!!!!!!!!!!!!!!!!!
//2023-02-14 afternoon : allow only the atoms on the top can move
//2023-02-16 put 2 atom at the center of the simulation system
//2023-02-20 debug : num_directionb and add stop warning to movea/b
//2023-03-15 separate the energy barrier of the atoms at the edge of one island to break away from this island and throw over it
//2023-03-21 increase the barrier energy of atom on the step edge of islands
//2023-04-13 Schwoebel
//2023-04-17 add estimated physical time
//2023-04-18 increase the energy barrier of same height movement of the atom on the island egde
//2023-04-19 using E_B, E_D, E_M to genaraliser this effect 
//2023-05-25 try to reforme the different energie barriers to a energy landscape general to respect the detailed balance
//2023-06-13 finish the part of different absorbed energies (gamma)
//2023-06-16 ridge-edge
//2023-06-27 change of random number generator that uses the random C++ engine based on the Mersenne Twister 19937 (64 bits) random number generator that should give lowest values of the random number of 1/(2^63-1) (instead of 1/2^32 for gsl that is a bit too large when Gf is small)
//2023-07-20 annealing/growth 
//2023-09-14 corriger la faute dans les probabilité (line 460-514 code précédant), spécialement le problème de gamma
//2023-09-15 change definition of u and w : u nn h<=h0-2 ; w h>=h0
//2023-09-27 keep E_ES and gamma , remove ridge-edge effect (local gamma) to accelate the program
//2023-10-24 using '-gamma2+ES2' takes place the max of '-gamma1+ES1' and '-gamma2+ES2' bacause this is just a reference which can be changed by E_ES (435-478)
//2023-12-12 compte the number of go up and down layer by layer
//2023-12-20 ridge; up-down between dh=2
//2023-01-15 delete EN*

// compile with
// g++ -O3 G-silicene_Gr-SOS-2024-01-15.cpp -lgslcblas -lgsl -lblitz -lm -std=c++14 -o g300_20240115

//run ./gNN_20240115 para coverage_f seed film
// ./g300_20240115 1000 0.5 1 1

    #include <string>
    #include <stdlib.h>
    #include <iostream>
    #include <cstdlib>
    #include <ctime>
    #include <fstream>
	#include <gsl/gsl_rng.h>
    #include <blitz/array.h>
    #include <time.h>
    #include <stdio.h>
    #include <math.h>
    #include <algorithm> 

    // C++ random number generator (C++ RNG)
    #include <random>
    // /////////////
    #define NN 300
//    #define step_length 100
    
    #define nb_layersub 0  //number of substrat layer 
	#define nb_layerf 20    //number max of final layer 
	
	#define ihmax 0 		// initial hmax =0
		
	#define nc 211 //nomber of class:
	
	#define nc_freq 176
	// //////////////
    #define random(x) rand()%(x) 
	
    using namespace std;
    using namespace blitz;
    
    int N2=NN*NN;
    //////////compute the number of movement
    long long int num_move=0 ;
    long long int num_move_1=0 ;
    long long int num_move_2=0 ;
	long long int num_up=0 ;
	long long int num_up_12=0 ;
	long long int num_up_23=0 ;
	long long int num_up2=0 ;
	long long int num_up2_13=0 ;
	long long int num_down=0 ;
	long long int num_down_12=0 ;
	long long int num_down_23=0 ;
	long long int num_down2=0 ;
	long long int num_down2_13=0 ;
    long int num_upfail = 0;
    long int num_downfail = 0;
    
    //estimated evolution time
    double dt=0 ; 
    double dti=0;
    
    clock_t CPU_start, CPU_end; //CPU time
    double CPU_elapsed_time;
    
    int ll=ihmax-nb_layersub + 1 ; //lowest layer
  	int hl=nb_layerf - nb_layersub +ihmax ; //highest layer
  	
    int  modN(int);
    double modNr(double);
    
    int minhmax(void),minhmaxa(void),minhmaxb(void) ;
    int maxhmax(void),maxhmaxa(void),maxhmaxb(void) ;
    
    int num_up2a(int,int,int,int),num_up2b(int,int,int,int),num_down2a(int,int,int,int),num_down2b(int,int,int,int);
    
    int computena(int,int),computenb(int,int); //h-3; down 2   
    int computema(int,int),computemb(int,int); //h-2; down 
    int computeoa(int,int),computeob(int,int); //h-1; diffusion
    int computepa(int,int),computepb(int,int); //h+1; up   
    int computeqa(int,int),computeqb(int,int); //h+2; up2 
    
	int computeclassa(int,int), computeclassb(int,int);
    void movea(int,int), moveb(int,int);
    void move_upa(int,int), move_upb(int,int);
    void move_up2a(int,int), move_up2b(int,int);
    void move_downa(int,int), move_downb(int,int);
    void move_down2a(int,int), move_down2b(int,int);
    void evaporatea(int,int), evaporateb(int,int);
    void deposita(int,int), depositb(int,int);
    void updatea(int,int), updateb(int,int);
    
    int nbcolh(int);
    int nbcoltot(void);
    
    void initialconditions(void);
 	
	Array<int, 2> ha(NN,NN), hb(NN,NN);
    // film height in real space for sublattice alpha (ha) and beta (hb) : ha(i,j) and hb(i,j)
    Array<int, 2> classa(NN,NN), classb(NN,NN);
    // array giving the number of nearest neighbors of point (i,j)
    Array<int, 2> listxa(N2,nc), listxb(N2,nc);
    // we gather the points that are in the same class n for n = 0,1... nc-1
    // for example, that have the same number nppv of ppv,
    // and for an order k=0, 1, 2 ... cardinal(nppv)-1 where cardinal(nppv) < N*N
    // this array gives the coordinate i of the k-th point of the array of points in the classe n (e.g. with nppv=n)
    // where k=order(i,j)
    Array<int, 2> listya(N2,nc), listyb(N2,nc);
    // idem that before but for the coordinate y
    Array<int, 2> ordera(NN,NN), orderb(NN,NN);
    // i.e. the point (i,j) is the order-th point in the list of points that have the same number of ppv
	Array<int, 6> situation2class(4,4,4,4,4,6); // situation2class(n,m,o,p,q,h)
	
	
	int cardinala[nc],cardinalb[nc] ;
    // nb of sites that have same class
    double freq[nc_freq] ;
    // total freqency of each class 
    double freq_up[nc_freq] ;
    double freq_up2[nc_freq] ;
    // freqency only including diffusion and abruption of each class
    double freq_down[nc_freq] ;
    double freq_down2[nc_freq] ;
    
    int class2n[nc] ;
    int class2m[nc] ;
    int class2o[nc] ;
    int class2p[nc] ;
    int class2q[nc] ;
    
    int class2h[nc] ;
    
    int class2freq[nc] ;
    int freq2class[nc_freq] ;
/////////////////////////
int main(int argc, char **argv){
	
	int para;
    double hdep ;
    int seed ;
    int film ; 
    
    sscanf(argv[1],"%d",&para);   // parameter of the simulation
    sscanf(argv[2],"%lg",&hdep);  // deposited thickness
	sscanf(argv[3],"%d",&seed);   // seed used of the simulation
    sscanf(argv[4],"%d",&film);   // film=1 : make film
    
    // READ PARAMETER FILE
    const char *filein = "para-Si_Gr_2024-01-15.ini";
    // input parameters file name
    // Open files to read and write
    string s;                                // data container
    ifstream fin(filein, ios::in);

	int paral ; //use it to check if the parameters of the simulation in the name of programme and parameters file are identical
    double fluxl, Tl, coefl;
    double gamma1l,gamma2l,gamma3l;
    double ES1l,ES2l,ES3l;
    double EN1l,EN2l,EN3l;
    
    double E12l,E23l,E34l;
    double E13l,E24l,E35l;
    
    double ERl;
    
    int testp = 0 ;

    while (testp<1){
    	getline(fin, s, '#');
        paral = atoi(s.c_str());           // parameter read
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        fluxl = atof(s.c_str());           // deposition flux in ML/s
        fin.ignore(80,'\n');

        getline(fin, s, '#');
        Tl = atof(s.c_str());           // temperature in Kelvin
        fin.ignore(80,'\n');
               
        getline(fin, s, '#');
        gamma1l = atof(s.c_str());            // surface energy for h=1
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        gamma2l = atof(s.c_str());            // surface energy for h=2
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        gamma3l = atof(s.c_str());            // surface energy for h=3
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        ES1l = atof(s.c_str());            // binding energy with substrate for h=1
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        ES2l = atof(s.c_str());            // binding energy with substrate for h=2
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        ES3l = atof(s.c_str());            // binding energy with substrate for h=3
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        EN1l = atof(s.c_str());            // binding energy with nn for h=1
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        EN2l = atof(s.c_str());            // binding energy with nn for h=2
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        EN3l = atof(s.c_str());            // binding energy with nn for h=3
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        E12l = atof(s.c_str());            // extra energy barrier between edge(h=1) and ridge(h=2)
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        E23l = atof(s.c_str());            // extra energy barrier between edge(h=2) and ridge(h=3)
        fin.ignore(80,'\n');
		
		getline(fin, s, '#');
        E34l = atof(s.c_str());            // extra energy barrier between edge(h,h>=3) and ridge(h+1,h>=3)
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        E13l = atof(s.c_str());            // extra energy barrier between edge(h=1) and ridge(h=3)
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        E24l = atof(s.c_str());            // extra energy barrier between edge(h=2) and ridge(h+2,h>=3)
        fin.ignore(80,'\n');
		
		getline(fin, s, '#');
        E35l = atof(s.c_str());            // extra energy barrier between edge(h,h>=3) and ridge(h+2,h>=3)
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        ERl = atof(s.c_str());            // extra absorption energy for ridge
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        coefl = atof(s.c_str());            // multiplicative coefficient for diffusion
        fin.ignore(80,'\n');
        
        getline(fin, s, '#');
        fin.ignore(80,'\n');               // separation lign
    
     if (paral == para){
            testp = 1;
            cout<< "Parameters file read" <<endl;
      }      
    }
    const double flux = fluxl ;
    const double T = Tl ;
   	
   	const double gamma1 = gamma1l ;
    const double ES1 = ES1l ;
    const double EN1 = EN1l ;

    const double gamma2 = gamma2l ; 
    const double ES2 = ES2l ;
    const double EN2 = EN2l ;
    
    const double gamma3 = gamma3l ; 
    const double ES3 = ES3l ;
    const double EN3 = EN3l ;
    
    const double ER = ERl ;
    
    const double E12 = E12l ;
    const double E23 = E23l ;
    const double E34 = E34l ;
    
    const double E13 = E13l ;
    const double E24 = E24l ;
    const double E35 = E35l ;
    
    const double coef = coefl ;
	
	
	double hmoy=0.;
	double htest=0.;   // We will get ' .dat' which gives the 'hmax(i,j)' from 'htest'to 'hdep'
	double incr=0.05 ;	//for each 'incr' monolayer during the growth
	double hincr=1./N2/2. ;//each deposit will increase 'hincr' monolayer
	
	double hmovie=0.1;
	double newhmovie;
	double hmoviemoy=0.;
	double increhmovie=hmovie;
	double movieframe=0.01;
	
	double hcheck=0.0 ;
	double incrhcheck=0.001 ;
	
	double timecheck=1800 ;
	double incretimecheck=1800 ;
	
	srand(seed);
	// Generating random numbers with C++11's RNG with Mersenne Twister 19937
    mt19937_64 rng64(seed);
    // uniform distribution in the default [0, 1) range :
    uniform_real_distribution<long double> unif;
    
    long double ran;
    
	// LIST INITIALIZATION
	int n,m,o,p,q,w ;
	int i,j,k,ek,h ; //k [0,nb_layerf-1] is the the index , h [ll,hl] is the height ;k=h-ll
	int c,c_freq ;
	
	////////////////////parameters
	 // Physical parameters
    const double kb = 1.38*1e-23;
    const double ev = 1.6*1e-19;
    const double hbar = 1.05*1e-34;
    
    double rapport;
    rapport = ev/(kb*T);
	
	double nuz ;
    nuz = coef*kb*T/(hbar*2*3.14159); // fréquence caractéristique de vibration atomique nuz=1.0*1e13 typiquement.  cf KMC
                    // with a multiplicative constant coef to adjust diffusion diffusion time scale // add 2021/06/24
    
    double Tdepot = nuz/(flux*2*NN*NN); // useful only for output not for the execution of the program
	double Qf = flux*2*NN*NN/nuz ;   //to be determined the coef corresponding	
	

	
	// Parameters output
	cout << " Growth silicene/Gr(honeycomb) starts !" <<endl ;
    cout << "Parameters " << para << " hdep = " << hdep << " flux(ML/s) = "<< flux << endl;
    cout << "T = " << T << " N = " << NN << endl;
    cout<<" gamma(h=1) = " << gamma1 << " (eV) ; gamma(h=2) = " << gamma2 << " (eV) ; gamma(h=3) = " << gamma3 <<endl ;
    cout<<" E_S(h=1) = " << ES1 << " (eV) ; E_N(h=1) = " << EN1 << " (eV) ; "<<endl ;
    cout<<" E_S(h=2) = " << ES2 << " (eV) ; E_N(h=2) = " << EN2 << " (eV) ; "<<endl ;
    cout<<" E_S(h=3) = " << ES3 << " (eV) ; E_N(h=3) = " << EN3 << " (eV) ; "<<endl ;
    cout<<" E_R = " << ER << " (eV) ; "<<endl ;
    cout<<" E_ES(1<->2) = " << E12 << " (eV) ; E_ES(2<->3) = "<<E23<<" (eV) ; E_ES(3<->4) = "<<E34<<" (eV) ;"<<endl ;
    cout<<" E_ES(1<->3) = " << E13 << " (eV) ; E_ES(2<->4) = "<<E24<<" (eV) ; E_ES(3<->5) = "<<E35<<" (eV) ;"<<endl ;
    cout << " Qf " << Qf << endl ;
    
    //////////compute the position of deposition
	int num_sub = 0;
	int num_flo = 0;
	//int d_atom_1 = 0;
	//int d_atom_2 = 0;
    ///////////////////////////////////////////////////////////////////////// 
    // Start CPU time
	CPU_start = clock();
    // ////////////////
    // PREPARE OUTPUT COORDINATE
    // ////////////////
    // OUTPUT FILES
    string sf1 = "Donnees_G/Para";
    string sf1_movie = "Donnees_G/Para";
    char sf2 [10];
    sprintf (sf2, "%i",para);
    char sf4 [10];
    sprintf (sf4, "%i",NN);
    string sf5 = "_flux_";
    char sf6 [10] ;
    sprintf (sf6, "%f",flux);
    string sf7 = "_T_";
    char sf8 [10] ;
    sprintf (sf8, "%g",T);
    string sf9 = "_seed_";
    char sf10 [10] ;
    sprintf (sf10, "%i",seed);
    string sf11 = "_COV_";   
    
    
    	
    char sf12_movie_begin [10] ;
    sprintf (sf12_movie_begin, "%1.3f",0.); 
    
    string sf12_movie_to="-" ;
    
    char sf12_movie_end [10] ;
    sprintf (sf12_movie_end, "%1.3f",increhmovie);
    
    FILE *fp_movie_a;
    FILE *fp_movie_b;
    
    string sf_movie_a;
    string sf_movie_b;
    
    string sf3_movie_a= "/Movie_a_N_";
	string sf3_movie_b= "/Movie_b_N_";

 
    string sf13 = ".dat";
     
    sf_movie_a= sf1_movie + sf2 + sf3_movie_a + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12_movie_begin + sf12_movie_to + sf12_movie_end + sf13;
	sf_movie_b= sf1_movie + sf2 + sf3_movie_b + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12_movie_begin + sf12_movie_to + sf12_movie_end + sf13;
	
    if (film==1){
    	//
    	size_t sizef_movie_a = sf_movie_a.size() + 1;
    	char * bufferf_movie_a = new char[sizef_movie_a];
    	strncpy(bufferf_movie_a,sf_movie_a.c_str(),sizef_movie_a);
    	fp_movie_a = fopen(bufferf_movie_a,"wb");
    	cout << bufferf_movie_a << endl ;
    	delete [] bufferf_movie_a;
    
    	size_t sizef_movie_b = sf_movie_b.size() + 1;
    	char * bufferf_movie_b = new char[sizef_movie_b];
    	strncpy(bufferf_movie_b,sf_movie_b.c_str(),sizef_movie_b);
    	fp_movie_b = fopen(bufferf_movie_b,"wb");
    	cout << bufferf_movie_b << endl ;
    	delete [] bufferf_movie_b;
    	
    
    }
    
    
    // ////////////////
    // OPEN OUTPUT FILES
    // ////////////////			
    int pointx, pointy ;
    Array<double, 2> xa(NN,NN);
    Array<double, 2> ya(NN,NN); 

    Array<double, 2> xb(NN,NN);
    Array<double, 2> yb(NN,NN); 

    
    // real x,y coordinates of the honeycomb lattice
    for(i=0; i<NN ; i++){
        for(j=0; j<NN; j++){
            // in units of e'_x , e'_y and e'_z
             xa(i,j) = modNr(i + j/2. - 1./4.) ;
             ya(i,j) = pow(3.,1./2.)/2. * j - 1./4./pow(3.,1./2.) ;
             xb(i,j) = modNr(i + j/2. + 1./4.) ;
             yb(i,j) = pow(3.,1./2.)/2. * j + 1./4./pow(3.,1./2.) ;           
        }
	}

	// END real x,y,z coordinates of the honeycomb lattice
	
	int nbatom ; 
	
	//check the proba of deposition
    double min_proba;
    double proba_depo;
    
    min_proba=(double)1./(rng64.max()); 
	
	int heighta[NN][NN]; // We will get hmax of colonne (i,j) of coordinate a1,a2,a3
	int heightb[NN][NN];
		
	//initialization of  freq[] 
	for ( c=0 ; c<nc ; c++ ){
		cardinala[c]=0;
		cardinalb[c]=0;				
	}
	for (c_freq=0; c_freq<nc_freq ;c_freq++){
		freq[c_freq]=0;
		freq_up[c_freq]=0;
		freq_up2[c_freq]=0;
		freq_down[c_freq]=0;
		freq_down2[c_freq]=0;
	}
	//give the frequency corresponding to the intermixing
	
	c=0;
	c_freq=0;
	
	//test point col
	if ((-gamma2+ES2+E12+gamma1+EN1)<0){
	    cout<<" Problem of E* between h=1 and h=2, E12 needs to be larger than "<<-(-gamma2+ES2+gamma1+EN1)<<endl;
	    cout<<" Program is stopped! "<<endl;
	    exit(1);
	}else if ( (ES2+E12+ER)<0){
	    cout<<" Problem of energy barrier of ridge at  h=2, E12 needs to be larger than "<<-ES2-ER<<endl;
	    cout<<" Program is stopped! "<<endl;
	    exit(1);
	}else if ((-gamma3+ES3+E23+gamma2+EN2)<0){
	    cout<<" Problem of E* between h=2 and h=3, E23 needs to be larger than "<<-(-gamma3+ES3+gamma2+EN2)<<endl;
	    cout<<" Program is stopped! "<<endl;
	    exit(1);
	}else if ((ES3+E23)<0){
	    cout<<" Problem of energy barrier of ridge at  h=3, E23 needs to be larger than "<<-ES3<<endl;
	    cout<<" Program is stopped! "<<endl;
	    exit(1);
	}else if ((-gamma3+ES3+E13+gamma1+EN1)<0){
		cout<<" Problem of E* between h=1 and h=3, E13 needs to be larger than "<<-(-gamma3+ES3+gamma1+EN1)<<endl;
	    cout<<" Program is stopped! "<<endl;
	    exit(1);
	}else if ((-gamma3+ES3+E24+gamma2+EN2)<0){
		cout<<" Problem of E* between h=2 and h=4, E13 needs to be larger than "<<-(-gamma3+ES3+gamma2+EN2)<<endl;
	    cout<<" Program is stopped! "<<endl;
	    exit(1);
	}

	
	//give values to some sites of freq[] which have physics meaning 
	for (q=0 ; q<4 ; q++){ //plane
		for (p=0 ; p<=3-q ; p++){
			for (o=0 ; o<=3-q-p ; o++){
				for (n=0 ; n<=3-q-p-o ; n++){
					m = 3-q-p-o-n;
					for (h=0 ;h<6 ; h++){
						c += 1;
						if (c==nc){
							cout<<" The assumed max of class number (nc) is too small , program is stopped !"<< endl;
							exit(1);
						}
						w = p + q ;
						 
    					class2q[c]=q ;
    					class2p[c]=p ;
    					class2o[c]=o ;
    					class2n[c]=n ;
    					class2m[c]=m ; 
    					
     					class2h[c]=h ;
										
						situation2class(m,n,o,p,q,h)=c ;	
			
						if(h>0){				
							c_freq += 1;
							class2freq[c] = c_freq ;
							freq2class[c_freq] = c ;			
														
							if (h==1){
								freq[c_freq] = o*exp(-rapport*(ES1+w*EN1));
								freq_up[c_freq] = p*exp(-rapport*(-gamma2+ES2+E12+gamma1+w*EN1) );
								freq_up2[c_freq] = q*exp(-rapport*(-gamma3+ES3+E13+gamma1+w*EN1) );																	
							}else if (h==2){
								if (n>0){
									freq[c_freq] = o*exp(-rapport*(ER+ES2+w*EN2));
									freq_up[c_freq] = p*exp(-rapport*(ER-gamma3+ES3+E23+gamma2+w*EN2) );
									freq_up2[c_freq] = q*exp(-rapport*(ER-gamma3+ES3+E24+gamma2+w*EN2) );
									freq_down[c_freq] = n*exp(-rapport*(ER+E12+ES2+w*EN2) );
								}else{
									freq[c_freq] = o*exp(-rapport*(ES2+w*EN2));
									freq_up[c_freq] = p*exp(-rapport*(-gamma3+ES3+E23+gamma2+w*EN2) );
									freq_up2[c_freq] = q*exp(-rapport*(-gamma3+ES3+E24+gamma2+w*EN2) );				
									freq_down[c_freq] = n*exp(-rapport*(E12+ES2+w*EN2) );
								}															
							}else if (h==3){
								freq[c_freq] = o*exp(-rapport*(ES3+w*EN3));
								freq_up[c_freq] = p*exp(-rapport*(ES3+E34+w*EN3));
								freq_up2[c_freq] = q*exp(-rapport*(ES3+E35+w*EN3));	
								freq_down[c_freq] = n*exp(-rapport*(ES3+E23+w*EN3) );
								freq_down2[c_freq] = m*exp(-rapport*(ES3+E13+w*EN3) );										       	
							}else if (h==4){
								freq[c_freq] = o*exp(-rapport*(ES3+w*EN3));     	
								freq_up[c_freq] = p*exp(-rapport*(ES3+E34+w*EN3));
								freq_up2[c_freq] = q*exp(-rapport*(ES3+E35+w*EN3));
								freq_down[c_freq] = n*exp(-rapport*(ES3+E34+w*EN3));
								freq_down2[c_freq] = m*exp(-rapport*(ES3+E24+w*EN3));
							}else{
								freq[c_freq] = o*exp(-rapport*(ES3+w*EN3));     	
								freq_up[c_freq] = p*exp(-rapport*(ES3+E34+w*EN3));
								freq_up2[c_freq] = q*exp(-rapport*(ES3+E35+w*EN3));
								freq_down[c_freq] = n*exp(-rapport*(ES3+E34+w*EN3));
								freq_down2[c_freq] = m*exp(-rapport*(ES3+E35+w*EN3));
							}
						}
					}	
				}
			}
		}
	}
				
	cout<<" cmax = "<<c<<endl;
	cout<<" c_freqmax = "<<c_freq<<endl;
	//exit(1);

	//////////////////////////
   	int freq_non_null0=0;
   		for (c=1;c<nc_freq;c++){
   			if (freq[c] != 0){
   				freq_non_null0 +=1 ;
   			}
   		}
 	cout<< "There are " <<freq_non_null0<< "  no void terms of frequency of diffusion "<<endl ;
 	
	/////////////////////////////

	
  //initialize the hmax
  for (i=0;i<NN;i++){
  	for (j=0;j<NN;j++){
  		ha(i,j)=ihmax;
  		hb(i,j)=ihmax;			
  	}
  }
	int ihmaxpu=ihmax + 1 ;
	int hlpu=hl + 1 ;
	int hlmu =hl-1 ;
	
	cout << "hmax initilized" <<endl;


   	initialconditions();
   	cout<<"initialize the tables"<<endl;
   	
   	int freq_non_null=0;
   	
   	for (c=1;c<nc_freq;c++){
   		if (freq[c] != 0){
   			freq_non_null +=1 ;
   		}
   	}
   	int class_c=situation2class(0,0,0,3,0,0);
   	for (c=0;c<nc;c++){
 			if (cardinala[c]>0){
 				cout<< "cardinala["<<c<<"] has "<<cardinala[c]<<" elements"<<endl;
 				cout<<" p = "<<class2p[c]<<" q = "<<class2q[c]<<endl;
 				cout<<" o = "<<class2o[c]<<endl;
 				cout<<" n = "<<class2n[c]<<" m = "<<class2m[c]<<endl;
 			}
 			if (cardinalb[c]>0){
 				cout<<" p = "<<class2p[c]<<" q = "<<class2q[c]<<endl;
 				cout<<" o = "<<class2o[c]<<endl;
 				cout<<" n = "<<class2n[c]<<" m = "<<class2m[c]<<endl;
 			}
 		}
 	if (freq_non_null==freq_non_null0 and cardinala[class_c]==N2 and cardinalb[class_c]==N2){
 		cout << "initialization test done " << endl ; 
 	}else {	
 		cout << "Inizialization makes mistakes, programme is stopped !" << endl;
 		cout<< "critical class is "<<class_c<< " cardinala[cc] = "<<cardinala[class_c]<<" cardinalb[cc] = "<<cardinalb[class_c]<<endl;
 		exit(1);
 	}
//exit(1); 	
///////////////////////////////////////////////////end initialization
/*
for (i=0;i<nc;i++){
	if (cardinal[i]>0){
		cout<<i<<"   "<<cardinal[i]<<endl;
	}
}
//exit(0) ;
*/
	
	int ncmu=nc_freq-1;
	
	//Q[c]=diffa+diffb+downa+downb+upa+upb
	//			Q1    Q2    Q3    Q4  Q5  Q
    double Q[nc_freq] ;
    double Q1[nc_freq];		//a ->
    double Q2[nc_freq];		//a -> + b ->
    double Q3[nc_freq];		//a -> + b -> + a down
    double Q4[nc_freq];		//a -> + b -> + a down + b down
    double Q5[nc_freq];		//a -> + b -> + a down + b down + a down2
    double Q6[nc_freq];		//a -> + b -> + a down + b down + a down2 + b down2 
    double Q7[nc_freq];		//a -> + b -> + a down + b down + a down2 + b down2 + a up
    double Q8[nc_freq];		//a -> + b -> + a down + b down + a down2 + b down2 + a up + b up
    double Q9[nc_freq];		//a -> + b -> + a down + b down + a down2 + b down2 + a up + b up + a up2
    Q[0]=Qf;     //use the fist Q for Qf because this initial first one doesn't have phiscs meaning (=0)
    Q1[0]=0;
	Q2[0]=0;
    Q3[0]=0;
    Q4[0]=0;
    Q5[0]=0;
    Q6[0]=0;
    Q7[0]=0;
    Q8[0]=0;
    Q9[0]=0;
    //////////////////////////////////////////////
    // REJECTION-FREE choice of the site to be moved
    double ra1 ; // random number in (0,Q[nc])
    int ra2 ; // int random  number for finding  the ord of the atom which will move in the class chosen
    double ra3 ;
    
	double del1,del2,del3,del4,del5,del6,del7,del8,del9 ;
    double res ;
    int class_chosen, freq_class_chosen;     
        // //////////////
        // MOVEMENT
        // /////////////
    double dot=1./N2 ;
    double hdeppu=hdep+dot ;    
         
    while (hmoy<hdeppu){
    	// we build the LADDER of RATES       
    	for (c_freq=1;c_freq<nc_freq;c_freq++){
    		c=freq2class[c_freq];
			Q[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq]+freq_down2[c_freq]+freq_up[c_freq]+freq_up2[c_freq])*(cardinala[c]+cardinalb[c]);
			Q1[c_freq]=Q[c_freq-1]+freq[c_freq]*cardinala[c];
			Q2[c_freq]=Q[c_freq-1]+freq[c_freq]*(cardinala[c]+cardinalb[c]);
			Q3[c_freq]=Q[c_freq-1]+freq[c_freq]*(cardinala[c]+cardinalb[c])+freq_down[c_freq]*cardinala[c];
			Q4[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq])*(cardinala[c]+cardinalb[c]);
			Q5[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq])*(cardinala[c]+cardinalb[c])+ freq_down2[c_freq]*cardinala[c] ;
			Q6[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq]+freq_down2[c_freq])*(cardinala[c]+cardinalb[c]);
			Q7[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq]+freq_down2[c_freq])*(cardinala[c]+cardinalb[c])+freq_up[c_freq]*cardinala[c];
    		Q8[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq]+freq_down2[c_freq]+freq_up[c_freq])*(cardinala[c]+cardinalb[c]);
    		Q9[c_freq]=Q[c_freq-1]+(freq[c_freq]+freq_down[c_freq]+freq_down2[c_freq]+freq_up[c_freq])*(cardinala[c]+cardinalb[c])+ freq_up2[c_freq]*cardinala[c];
    	}
		ran= unif(rng64) ;
    	ra1=ran*Q[ncmu];
    	
    	//The estimated exp duration
    	ra3=(double)rand()/RAND_MAX; //random number in (0,1)
		dti=-log(ra3)/Q[ncmu]/nuz;
		proba_depo=(double)Qf/Q[ncmu];
		dt += dti;
        
		if (ra1<=Qf){ // choice of deposition
		  // //////////////
          // Regular CHECKS of the dynamics
          // /////////////
 
         	if (hmoy>=hcheck){  
         		CPU_end = clock();
    			CPU_elapsed_time = ((double) (CPU_end - CPU_start))/CLOCKS_PER_SEC;        	     

				cout << "Running Time = " << CPU_elapsed_time <<" (s) now ! hmoy = "<<hmoy<< endl ;
				cout << "The largest hmax is "<<maxhmax()<<endl;
				cout <<"The smallest hmax is "<<minhmax()<<endl;
				cout<<"There are " << hmoy/hincr <<" time(s) deposition(s) ; " << endl;
				cout<<"There are " << num_move <<" time(s) diffusion(s) ; " << endl;
				cout<<"There are " << num_sub <<" atoms deposited on the substrate "<<endl;
				cout<<"There are " << num_flo <<" atoms deposited on the crystal "<<endl;
				cout<<"There are " << num_up+num_up2 <<" time(s) that one atom moves to a higher position ; " << endl;
				cout<<"There are " << num_down+num_down2 <<"  time(s) that one atom moves to a lower position ; " << endl;
				for (h=1;h<maxhmax()+1;h++){
					cout<<"There is/are "<< nbcolh(h)<<" columns of h = "<<h<<endl;
				}
				cout<<"The estimated physical time is "<<dt<<" s"<<endl; 
				cout<<"The simulation flux (hmoy/physt) = "<<hmoy/dt<<" ML/s ; the parameter flux = "<< flux<<" ML/s"<<endl;
				cout<<"Probability of deposition is "<< proba_depo<<" comparing with the smallest random number "<<min_proba<<endl;
				cout<<endl;				
				
            	hcheck += incrhcheck ;
          	}else if (hmoy< incrhcheck){
          		CPU_end = clock();
    			CPU_elapsed_time = ((double) (CPU_end - CPU_start))/CLOCKS_PER_SEC;
          		cout << "Running Time = " << CPU_elapsed_time <<" (s) now ! hmoy = "<<hmoy<< endl ;
          		cout<<"There are " << num_sub <<" atoms deposited on the substrate "<<endl;
				cout<<"There are " << num_flo <<" atoms deposited on the crystal "<<endl;
          		cout<<"There are " << num_up+num_up2 <<" time(s) that one atom moves to a higher position " << endl;
				cout<<"There are " << num_down+num_down2 <<"  time(s) that one atom moves to a lower position " << endl;
				for (h=1;h<maxhmax()+1;h++){
					cout<<"There is/are "<< nbcolh(h)<<" columns of h = "<<h<<endl;
				}
				cout<<"The estimated physical time is "<<dt<<" s"<<endl;
				cout<<"The simulation flux (hmoy/physt) = "<<hmoy/dt<<" ML/s ; the parameter flux = "<< flux<<" ML/s"<<endl;
          		cout<<"Probability of deposition is "<< proba_depo<<" comparing with the smallest random number "<<min_proba<<endl;
          		cout<<endl;
          	}
          // //////////////
          // END Regular CHECKS of the dynamics
          // /////////////

    		i=random(NN);
    		j=random(NN);
    		if (ha(i,j)>=hlmu or hb(i,j)>=hlmu ){
    			cout << "Warning : hmax is too large :deposition fails ! " << endl;
    			exit (1) ;
    		}else{
    			if (ra1<=Qf/2.){
    				deposita(i,j);
    				if (ha(i,j)==1){
    					num_sub += 1 ;
    				}else{
    					num_flo += 1;
    				}
    			}else{
    				depositb(i,j);
    				if (hb(i,j)==1){
    					num_sub += 1 ;
    				}else{
    					num_flo += 1;
    				}  				
    			}
    			hmoy+=hincr;
    			
//    			cout<< "deposition k = "<< hb(i,j)+1-ll <<endl;  //test
//    			cout << "deposition : hmoy = " << hmoy << endl ;
			}

    		
    }//finish choice of deposition
    
	else { // choice of diffusion 
		for (c=1;c<nc_freq;c++){
			if (ra1<=Q[c]){
				freq_class_chosen=c;
				break;
             }
         }
        // random number in (0,size_class_chosen-1) to find the order of atom which will move in the chosen class
        
        
		res=ra1-Q[freq_class_chosen-1]; 
		del1=Q1[freq_class_chosen]-Q[freq_class_chosen-1];
		del2=Q2[freq_class_chosen]-Q[freq_class_chosen-1];
		del3=Q3[freq_class_chosen]-Q[freq_class_chosen-1];
		del4=Q4[freq_class_chosen]-Q[freq_class_chosen-1];
		del5=Q5[freq_class_chosen]-Q[freq_class_chosen-1];
		del6=Q6[freq_class_chosen]-Q[freq_class_chosen-1];
		del7=Q7[freq_class_chosen]-Q[freq_class_chosen-1];
		del8=Q8[freq_class_chosen]-Q[freq_class_chosen-1];
		del9=Q9[freq_class_chosen]-Q[freq_class_chosen-1];
		
		class_chosen=freq2class[freq_class_chosen];
		
		h=class2h[class_chosen];
				
		if (res<0){
			cout<<"Problems of class chosen of res,program is stopped !"<<endl;
			exit(1);
		}else{
			if (res<=del1){
				ra2 = random(cardinala[class_chosen]);
				movea(ra2,class_chosen);
				num_move += 1;
				if (h==1){
					num_move_1 += 1;
				}else if (h==2){
					num_move_2 += 1; 
				}
			}else if(res<=del2){
				ra2 = random(cardinalb[class_chosen]);
				moveb(ra2,class_chosen);
				num_move += 1;
				if (h==1){
					num_move_1 += 1;
				}else if (h==2){
					num_move_2 += 1; 
				}
			}else if(res<=del3){
				ra2 = random(cardinala[class_chosen]);
				move_downa(ra2,class_chosen);
				num_down += 1;
				if (h==2){
					num_down_12 += 1;
				}else if (h==3){
					num_down_23 += 1; 
				}
			}else if (res<=del4){
				ra2 = random(cardinalb[class_chosen]);
				move_downb(ra2,class_chosen);
				num_down += 1;
				if (h==2){
					num_down_12 += 1;
				}else if (h==3){
					num_down_23 += 1; 
				}
			}else if (res<=del5){
				ra2 = random(cardinala[class_chosen]);
				move_down2a(ra2,class_chosen);
				num_down2 += 1;
				if (h==3){
					num_down2_13 += 1;
				}
			}else if (res<=del6){
				ra2 = random(cardinalb[class_chosen]);
				move_down2b(ra2,class_chosen);
				num_down2 += 1;
				if (h==3){
					num_down2_13 += 1;
				}
			}else if (res<=del7){
				ra2 = random(cardinala[class_chosen]);
				move_upa(ra2,class_chosen);
				num_up += 1;
				if (h==1){
					num_up_12 += 1;
				}else if (h==2){
					num_up_23 += 1; 
				}
			}else if (res<=del8){
				ra2 = random(cardinalb[class_chosen]);
				move_upb(ra2,class_chosen);
				num_up += 1;
				if (h==1){
					num_up_12 += 1;
				}else if (h==2){
					num_up_23 += 1; 
				}
			}else if (res<=del9){
				ra2 = random(cardinala[class_chosen]);
				move_up2a(ra2,class_chosen);
				num_up2 += 1;
				if (h==1){
					num_up2_13 += 1;
				}
			}else{
				ra2 = random(cardinalb[class_chosen]);
				move_up2b(ra2,class_chosen);
				num_up2 += 1;
				if (h==1){
					num_up2_13 += 1;
				}
			}
			
			
			//cout<<"diffusion : num_move = "<<num_move<<endl;
		}

			
	}
	
	//cout<< "diffusion : number of diffusion = "<< num_move << endl;


    // ////////////////
    // OPEN OUTPUT MOVIES '.dat'
    // ////////////////  
    	
    if (film==1 and hmoy>=hmovie){
    	newhmovie=hmovie+increhmovie;
    	if (newhmovie<=hdep){
    		fclose(fp_movie_a);
    		fclose(fp_movie_b);

    		char sf12_movie_begin [10] ;
    		sprintf (sf12_movie_begin, "%1.3f",hmovie); 
    	
    		string sf12_movie_to="-" ;
    
    		char sf12_movie_end [10] ;
    		sprintf (sf12_movie_end, "%1.3f",newhmovie);
    
    		FILE *fp_movie_a;
    		FILE *fp_movie_b;
     
    		string sf_movie_a;
    		string sf_movie_b;
    
    		string sf3_movie_a= "/Movie_a_N_";
    		string sf3_movie_b= "/Movie_b_N_";

    		string sf13 = ".dat";
     
    		sf_movie_a= sf1_movie + sf2 + sf3_movie_a + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12_movie_begin + sf12_movie_to + sf12_movie_end + sf13;
    		sf_movie_b= sf1_movie + sf2 + sf3_movie_b + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12_movie_begin + sf12_movie_to + sf12_movie_end + sf13;
    

  
    			//C2
    		size_t sizef_movie_a = sf_movie_a.size() + 1;
    		char * bufferf_movie_a = new char[sizef_movie_a];
    		strncpy(bufferf_movie_a,sf_movie_a.c_str(),sizef_movie_a);
    		fp_movie_a = fopen(bufferf_movie_a,"wb");
    		cout << bufferf_movie_a << endl ;
    		delete [] bufferf_movie_a;
    
    		size_t sizef_movie_b = sf_movie_b.size() + 1;
    		char * bufferf_movie_b = new char[sizef_movie_b];
    		strncpy(bufferf_movie_b,sf_movie_b.c_str(),sizef_movie_b);
    		fp_movie_b = fopen(bufferf_movie_b,"wb");
    		cout << bufferf_movie_b << endl ;
    		delete [] bufferf_movie_b;
    
    	hmovie = newhmovie ;
		}
    }
    
    // ////////////////
    // OUTPUT writting MOVIES '.dat'
    // ////////////////
    
    if (film==1 and hmoy>=hmoviemoy){
           
           for(i=0; i<NN ; i++){
               for(j=0; j<NN; j++){              		               		
               		heighta[i][j] = ha(i,j) ;
                   	heightb[i][j] = hb(i,j) ;
               }
           }
           
             fwrite(heighta,sizeof(int),NN*NN,fp_movie_a);
             fwrite(heightb,sizeof(int),NN*NN,fp_movie_b);
           
           //cout << "movie output written " << endl ;
           hmoviemoy += movieframe ;
       }
       
    // //////////////
    // OUTPUT writting  .dat .xyz and .ini files
    // /////////////
    
    if (hmoy>=htest){ 

   		   
    // ////////////////
    // OPEN OUTPUT FILES '.dat'
    // ////////////////

    char sf12 [10] ;
    sprintf (sf12, "%1.3f",htest); 
    
    FILE *fp_at;
    FILE *fp_bt;
    
   // OPEN OUTPUT FILES : hmax     // for indirect visualization
    
    string sf_at ;
    string sf_bt ;
    
	string sf3 = "/Data_N_";
    string sf3_at = "/a_N_";
  	string sf3_bt = "/b_N_";
  	  	
    string sf13 = ".dat";

    sf_at= sf1 + sf2 + sf3_at + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12 + sf13;
    sf_bt= sf1 + sf2 + sf3_bt + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12 + sf13;      
           
    size_t sizef_at = sf_at.size() + 1;
    char * bufferf_at = new char[sizef_at];
    strncpy(bufferf_at,sf_at.c_str(),sizef_at);
    fp_at = fopen(bufferf_at,"wb");
    cout << bufferf_at << endl ;
    delete [] bufferf_at;
    
    size_t sizef_bt = sf_bt.size() + 1;
    char * bufferf_bt = new char[sizef_bt];
    strncpy(bufferf_bt,sf_bt.c_str(),sizef_bt);
    fp_bt = fopen(bufferf_bt,"wb");
    cout << bufferf_bt << endl ;
    delete [] bufferf_bt;
    
    // END OPEN OUTPUT .dat FILES 
    	       	// OUTPUT writting a/bt.dat file
				for (i=0; i<NN; i++){
					for (j=0; j<NN; j++ ){
						heighta[i][j]=ha(i,j);
						heightb[i][j]=hb(i,j);
					}
				}
           		
             	fwrite(heighta,sizeof(int),NN*NN,fp_at);
             	fwrite(heightb,sizeof(int),NN*NN,fp_bt);
             	
             	// close output files
     			fclose(fp_at);
     			fclose(fp_bt);
     			
       			// END OUTPUT writting .dat file

    			// OPEN OUTPUT FILE : multicolor xyz format  
   			    /**/
    		
    			FILE *fpxyz2C; // atoms with 2 colors h=1 and h>1
       
    			string sfxyz2C ; 
    
    			string sf3xyz = "/WE_N_";
    			char sf4xyz [10];
    			sprintf (sf4xyz, "%i",NN);
    			string sf5xyz = "_flux_";
    			char sf6xyz [10];
    			sprintf (sf6xyz, "%f",flux);
    			string sf7xyz = "_T_";
    			char sf8xyz [10];
    			sprintf (sf8xyz, "%g",T);
    			string sf9xyz = "_seed_";
    			char sf10xyz [10] ;
    			sprintf (sf10xyz, "%i",seed);
    			string sf11xyz = "_COV_";
    			char sf12xyz [10] ;
    			sprintf (sf12xyz, "%1.3f",htest);
   
    			string sf13xyz2C = "_2C.xyz";
    

   				sfxyz2C = sf1 + sf2 + sf3xyz + sf4xyz + sf5xyz + sf6xyz + sf7xyz +sf8xyz +sf9xyz + sf10xyz + sf11xyz + sf12xyz + sf13xyz2C ;

    			size_t sizefxyz2C = sfxyz2C.size() + 1;
  
    			char * bufferfxyz2C = new char[sizefxyz2C];

    			strncpy(bufferfxyz2C,sfxyz2C.c_str(),sizefxyz2C);

    			fpxyz2C = fopen(bufferfxyz2C,"w");
    			
    			delete [] bufferfxyz2C;
    			/* */
    			// END  OPEN .XYZ FILES

           		// ////////////// //STOP HERE
    			// OUTPUT writting .xyz file
    			// /////////////
    		nbatom=0 ;
    			
    		for (i=0; i<NN; i++){
    			for (j=0; j<NN; j++){
    				nbatom += ha(i,j)+ hb(i,j);
    			}
    		}	

    		fprintf(fpxyz2C,"%d  \n", nbatom);
    		// second line of an xyz file is a comment
    		fprintf(fpxyz2C,"Si \n");  
    		//h_t <= 0
    		for (i=0; i<NN; i++ ){
    			for (j=0; j<NN; j++){
    				if (ha(i,j)>=1){
    					fprintf(fpxyz2C,"Si        %f   %f   %d \n",xa(i,j),ya(i,j),1);
    					for (k=2;k<ha(i,j)+1;k++){
    						fprintf(fpxyz2C,"Al        %f   %f   %d \n",xa(i,j),ya(i,j),k);
    					}
    				}
    				if (hb(i,j)>=1){
    					fprintf(fpxyz2C,"Si        %f   %f   %d \n",xb(i,j),yb(i,j),1);
    					for (k=2;k<hb(i,j)+1;k++){
    						fprintf(fpxyz2C,"Al        %f   %f   %d \n",xb(i,j),yb(i,j),k);
    					}
    				}
    			}
    		}
    			
        		
     //close .xyz files									
     fclose(fpxyz2C);
    
    /* */
    // //////////////
    // END OUTPUT writting .xyz file
    // /////////////


   				// //////////////
    			// OUTPUT writting .ini file
    			// /////////////
    
        	// Compute CPU elapsed time
    			
    		CPU_end = clock();
   			CPU_elapsed_time = ((double) (CPU_end - CPU_start))/CLOCKS_PER_SEC;
   			
   	FILE *fpdata;
   	
	string sfdata ;
    char sf12ini [10] ;
    sprintf (sf12ini, "%1.3f",htest);
	string sf13ini=".ini";
	sfdata=sf1 + sf2 + sf3 + sf4 + sf5 + sf6 + sf7 +sf8 + sf9 + sf10 + sf11 + sf12ini + sf13ini;
	
	size_t sizefdata = sfdata.size() + 1;
    char * bufferfdata = new char[sizefdata];
    strncpy(bufferfdata,sfdata.c_str(),sizefdata);
    fpdata = fopen(bufferfdata,"wb");
    cout << bufferfdata << endl ;
    delete [] bufferfdata; 
    
    //d_atom_1 = cardinala[situation2class(0,3,0,1)]+cardinalb[situation2class(0,3,0,1)];
    //d_atom_1 = cardinala[situation2class(0,3,0,2)]+cardinalb[situation2class(0,3,0,2)];
    
    //write output data
	fprintf(fpdata,"Recall Parameters : \n");
	fprintf(fpdata," Parameters %i  , hdep =  %f (ML) , flux = %f (ML/s) \n",para,htest,flux);
	fprintf(fpdata," T = %f (K), N = %i \n",T,NN);
	fprintf(fpdata,"h<=0 : \n");
	fprintf(fpdata," \n" );
	fprintf(fpdata," gamma (h=1) = %f (eV) \n" , gamma1);
	fprintf(fpdata," gamma (h=2) = %f (eV) \n" , gamma2);
	fprintf(fpdata," gamma (h=3) = %f (eV) \n" , gamma3);
	fprintf(fpdata," \n" );
	fprintf(fpdata," E_S (h=1) = %f (eV) \n" , ES1);
	fprintf(fpdata," E_N (h=1) = %f (eV) \n" , EN1);
	fprintf(fpdata," E_S (h=2) = %f (eV) \n" , ES2);
	fprintf(fpdata," E_N (h=2) = %f (eV) \n" , EN2);
	fprintf(fpdata," E_S (h=3) = %f (eV) \n" , ES3);
	fprintf(fpdata," E_N (h=3) = %f (eV) \n" , EN3);
	
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," E_R = %f (eV) \n" , ER);
	
	fprintf(fpdata," \n" );
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," E_ES (1<->2) = %f (eV) \n" , E12);
	fprintf(fpdata," E_ES (2<->3) = %f (eV) \n" , E23);
	fprintf(fpdata," E_ES (3<->4) = %f (eV) \n" , E34);
	
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," E_ES (1<->3) = %f (eV) \n" , E13);
	fprintf(fpdata," E_ES (2<->4) = %f (eV) \n" , E24);
	fprintf(fpdata," E_ES (3<->5) = %f (eV) \n" , E35);
	
	fprintf(fpdata," \n" );
	
	
	fprintf(fpdata, " The ratio of the frequency of atom moving up from h=1 to h=2 and the frequency of atom moving down from h=2 to h=1 (up/down) is %e \n",exp(-rapport*(max(-gamma1+ES1,-gamma2+ES2)+E12+gamma1+EN1) )/ exp(-rapport*(max(-gamma1+ES1,-gamma2+ES2)+E12+gamma2) ));
	
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," The ratio of the frequency of atom of edge (h=1) leaving the bord and the frequency of moving up is %e \n",exp(-rapport*(ES1+EN1))/exp(-rapport*(max(-gamma1+ES1,-gamma2+ES2)+E12+gamma1+EN1) ));
	
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," The ratio of the frequency of atom of ridge (h=2) leaving the bord and the frequency of moving down is %e \n",exp(-rapport*(ES2))/ exp(-rapport*(max(-gamma1+ES1,-gamma2+ES2)+E12+gamma2) ));
	
	fprintf(fpdata," \n" );
	fprintf(fpdata," \n" );

	//fprintf(fpdata," Tdepot = %f (s) \n", Tdepot);
	fprintf(fpdata," Qf = %e \n",Qf);
	fprintf(fpdata," Estimated_exp_time = %f (s) \n",dt);
	fprintf(fpdata," CPU_elapsed_time =  %f (s) \n",CPU_elapsed_time);
	
	fprintf(fpdata," \n" );


	// output information
	fprintf(fpdata," There are %i atoms deposited on the substrate \n", num_sub);
	fprintf(fpdata," There are %i atoms deposited on the crystal \n", num_flo);
	
	fprintf(fpdata," \n" );
	fprintf(fpdata," There are %i diffusion atoms at h=1 (on substrate) \n", cardinala[situation2class(0,0,3,0,0,1)]+cardinalb[situation2class(0,0,3,0,0,1)]);
	fprintf(fpdata," There are %i diffusion atoms at h=2 \n", cardinala[situation2class(0,0,3,0,0,2)]+cardinalb[situation2class(0,0,3,0,0,2)]);
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," There are %lli time(s) diffusion(s) \n",num_move);
	fprintf(fpdata," There are %lli time(s) diffusion(s) at h=1 \n",num_move_1);
	fprintf(fpdata," There are %lli time(s) diffusion(s) at h=2 \n",num_move_2);
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," There are %lli time(s) the atom moves up  \n",num_up);
	fprintf(fpdata," There are %lli time(s) the atom moves up from h=1 to h=2 \n",num_up_12);
	fprintf(fpdata," There are %lli time(s) the atom moves up from h=2 to h=3 \n",num_up_23);
	fprintf(fpdata," There are %lli time(s) the atom moves up2  \n",num_up2);
	fprintf(fpdata," There are %lli time(s) the atom moves up from h=1 to h=3 \n",num_up2_13);
	
	fprintf(fpdata," \n" );
	fprintf(fpdata," There are %li time(s) the atom fails to move up  \n",num_upfail);
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," There are %lli time(s) the atom moves down \n",num_down);
	fprintf(fpdata," There are %lli time(s) the atom moves down from h=2 to h=1 \n",num_down_12);
	fprintf(fpdata," There are %lli time(s) the atom moves down from h=3 to h=2 \n",num_down_23);
	fprintf(fpdata," There are %lli time(s) the atom moves down2 \n",num_down2);
	fprintf(fpdata," There are %lli time(s) the atom moves down from h=3 to h=1 \n",num_down2_13);
	
	fprintf(fpdata," \n" ); 
	fprintf(fpdata," There are %li time(s) the atom fails to move down \n",num_downfail);
	fprintf(fpdata," \n" );
		
	fprintf(fpdata," The largest hmax is %i \n",maxhmax());
	fprintf(fpdata," The smallest hmax is %i \n",minhmax());
	
	fprintf(fpdata," \n" ); 
	
	for (h=1;h<maxhmax()+1;h++){
		fprintf(fpdata,"There is/are %i columns of h = %i \n",nbcolh(h),h);
	}
	fprintf(fpdata," \n" );
	
	fprintf(fpdata," The ratio of 2D Si is %f \n",(double)(nbcolh(1)-cardinala[situation2class(0,0,3,0,0,1)]-cardinalb[situation2class(0,0,3,0,0,1)])/nbcoltot());
	
	fclose(fpdata);
	
	
	// //////////////
    // END OUTPUT writting .ini file
    // /////////////
    
           	 	cout << "output written : hmoy = " << htest << endl ;

           	 	htest += incr ;
       		}
    	// //////////////
       // END OUTPUT writting  .dat .xyz and .ini files
       // /////////////
       
	
	//time check
		
	// Compute CPU elapsed time
    CPU_end = clock();
    CPU_elapsed_time = ((double) (CPU_end - CPU_start))/CLOCKS_PER_SEC;
	
	if (CPU_elapsed_time>=timecheck){
		proba_depo=(double)Qf/Q[ncmu];
         		
		cout << "Running Time = " << CPU_elapsed_time <<" (s) now ! hmoy = "<<hmoy<< endl ;
		cout << "The largest hmax is "<<maxhmax()<<endl;
		cout <<"The smallest hmax is "<<minhmax()<<endl;
		cout<<"There are " << num_move <<" time(s) diffusion(s) " <<  endl;
		cout<<"There are " << num_sub <<" atoms deposited on the substrate "<<endl;
		cout<<"There are " << num_flo <<" atoms deposited on the crystal "<<endl;
		cout<<"There are " << num_up <<" time(s) that one atom moves to a higher position  " << endl;
		cout<<"There are " << num_down <<"  time(s) that one atom moves to a lower position  " << endl;
		for (h=1;h<maxhmax()+1;h++){
			cout<<"There is/are "<< nbcolh(h)<<" columns of h = "<<h<<endl;
		}
		cout<<"The estimated physical time is "<<dt<<" s"<<endl;
		cout<<"The simulation flux (hmoy/physt) = "<<hmoy/dt<<" ML/s ; the parameter flux = "<< flux<<" ML/s"<<endl;	
		cout<<"Probability of deposition is "<< proba_depo<<" comparing with the smallest random number "<<min_proba<<endl;
		cout<<endl;		
        timecheck += incretimecheck ;
	}	
	
	
	
    }
           
        // //////////////
        // END MOVEMENT
        // /////////////

				
    //
    /* */

    
    // Compute CPU elapsed time
    CPU_end = clock();
    CPU_elapsed_time = ((double) (CPU_end - CPU_start))/CLOCKS_PER_SEC;
    
    // Recall Parameters 
    cout << " Parameters " << para << " ,hdep = " << hdep << " ,flux = "<< flux << endl;
    cout << " T = " << T << " N = " << NN << endl;
    cout<<" E_S(h=1) = " << ES1 << " (eV) ; E_N(h=1) = " << EN1 << " (eV) ; "<<endl ;
    cout<<" E_S(h=2) = " << ES2 << " (eV) ; E_N(h=2) = " << EN2 << " (eV) ; "<<endl ;
    cout<<" E_S(h=3) = " << ES3 << " (eV) ; E_N(h=3) = " << EN3 << " (eV) ; "<<endl ;
    cout<<" E_ES(1<->2) = " << E12 << " (eV) ; E_ES(2<->3) = "<<E23<<" (eV) "<<endl ;
    cout << " Qf =" << Qf << endl ;
    cout<<" The estimated physical time is "<<dt<<" s"<<endl;
    cout <<" Growth done ! CPU_elapsed_time = " << CPU_elapsed_time << endl ;	

	

    return 0;
}



////////////////////////////////////
	
void initialconditions(void){     	
	int i, j ,k,c;
	int N2=NN*NN;
	
	//initialization of cardinal[]  and list(class,order)
	for ( c=0 ; c<nc ; c++ ){	
		
		cardinala[c]=0;
		cardinalb[c]=0;
		
		for (i=0 ; i < N2 ; i++){			
			listxa(i,c) = -1;   // initialization // not in the range 0, 1, ... N-1 of the system
       		listya(i,c) = -1;
       		listxb(i,c) = -1;   // initialization // not in the range 0, 1, ... N-1 of the system
       		listyb(i,c) = -1;
		}
	}
	
	//initialize the class
  	
	for (i=0;i<NN;i++){
   		for (j=0;j<NN;j++){
     			
     		classa(i,j) = computeclassa(i,j);
         	classb(i,j) = computeclassb(i,j);
         
         	cardinala[classa(i,j)] += 1;
         	cardinalb[classb(i,j)] += 1;
       
         	listxa(cardinala[classa(i,j)]-1,classa(i,j)) = i;
         	listya(cardinala[classa(i,j)]-1,classa(i,j)) = j;
         		
         	listxb(cardinalb[classb(i,j)]-1,classb(i,j)) = i;
         	listyb(cardinalb[classb(i,j)]-1,classb(i,j)) = j;

         	ordera(i,j) = cardinala[classa(i,j)]-1;
         	orderb(i,j) = cardinalb[classb(i,j)]-1;
         	
         }      
     }
}
////////////////////////////////////////////////////////
int computeoa(int i, int j){ // nb of Si in the plane of  the site (i,j,k-1)
	int o , im , jm , hd ;
	o=0;
	hd = ha(i,j)-1;

	im = modN(i-1);
  	jm = modN(j-1);
  	
  	if (hb(i,j)==hd){
  		o += 1;
  	}		
  	if (hb(im,j)==hd){
  		o += 1;
  	}
  	if (hb(i,jm)==hd){
  		o += 1;
  	}
  	
	return o;
}

int computeob(int i, int j){ // nb of Si in the plane of  the site (i,j,k-1)
	int o , ip , jp , hd ;
	o = 0;
	hd = hb(i,j)-1;
	
	ip = modN(i+1);
  	jp = modN(j+1);
  	
  	if (ha(i,j)==hd){
  		o += 1;
  	}		
  	if (ha(ip,j)==hd){
  		o += 1;
  	}
  	if (ha(i,jp)==hd){
  		o += 1;
  	}
  	
	return o;
}
int computena(int i, int j){ // nb of Si in the plane of  the site (i,j,k-2)
	int n , im , jm , hd ;
	n = 0;
	hd = ha(i,j)-2;
	
	im = modN(i-1);
  	jm = modN(j-1);

  	if (hb(i,j)==hd){
  		n += 1;
  	}		
  	if (hb(im,j)==hd){
  		n += 1;
  	}
  	if (hb(i,jm)==hd){
  		n += 1;
  	}

	return n;
}

int computenb(int i, int j){ // nb of Si in the plane of  the site (i,j,k-2)
	int n , ip , jp , hd ;
	n = 0;
	hd = hb(i,j)-2;

	ip = modN(i+1);
  	jp = modN(j+1);
  	
  	if (ha(i,j)==hd){
  		n += 1;
  	}		
  	if (ha(ip,j)==hd){
  		n += 1;
  	}
  	if (ha(i,jp)==hd){
  		n += 1;
  	}
  
	return n;
}

int computema(int i, int j){ // nb of Si in the plane of  the site (i,j,k-3)
	int m , im , jm , hd ;
	m = 0;
	hd = ha(i,j)-3;
	
	im = modN(i-1);
  	jm = modN(j-1);

  	if (hb(i,j)<=hd){
  		m += 1;
  	}		
  	if (hb(im,j)<=hd){
  		m += 1;
  	}
  	if (hb(i,jm)<=hd){
  		m += 1;
  	}

	return m;
}

int computemb(int i, int j){ // nb of Si in the plane of  the site (i,j,k-3)
	int m , ip , jp , hd ;
	m = 0;
	hd = hb(i,j)-3;

	ip = modN(i+1);
  	jp = modN(j+1);
  	
  	if (ha(i,j)<=hd){
  		m += 1;
  	}		
  	if (ha(ip,j)<=hd){
  		m += 1;
  	}
  	if (ha(i,jp)<=hd){
  		m += 1;
  	}
  
	return m;
}

int computepa(int i, int j){ // nb of Si in the plane of  the site (i,j,k+1)
	int p , im , jm , hd ;
	p = 0;
	hd = ha(i,j);
	
	im = modN(i-1);
  	jm = modN(j-1);
  	
  	if (hb(i,j)==hd){
  		p += 1;
  	}		
  	if (hb(im,j)==hd){
  		p += 1;
  	}
  	if (hb(i,jm)==hd){
  		p += 1;
  	}
  	
	return p;
}

int computepb(int i, int j){ // nb of Si in the plane of  the site (i,j,k+1)
	int p , ip , jp , hd ;
	p = 0;
	hd = hb(i,j);

	ip = modN(i+1);
  	jp = modN(j+1);
  	
  	if (ha(i,j)==hd){
  		p += 1;
  	}		
  	if (ha(ip,j)==hd){
  		p += 1;
  	}
  	if (ha(i,jp)==hd){
  		p += 1;
  	}

	return p;
}

int computeqa(int i, int j){ // nb of Si in the plane of  the site (i,j,k+1)
	int q , im , jm , hd ;
	q = 0;
	hd = ha(i,j)+1;
	
	im = modN(i-1);
  	jm = modN(j-1);
  	
  	if (hb(i,j)>=hd){
  		q += 1;
  	}		
  	if (hb(im,j)>=hd){
  		q += 1;
  	}
  	if (hb(i,jm)>=hd){
  		q += 1;
  	}
  	
	return q;
}

int computeqb(int i, int j){ // nb of Si in the plane of  the site (i,j,k+1)
	int q , ip , jp , hd ;
	q = 0;
	hd = hb(i,j)+1;

	ip = modN(i+1);
  	jp = modN(j+1);
  	
  	if (ha(i,j)>=hd){
  		q += 1;
  	}		
  	if (ha(ip,j)>=hd){
  		q += 1;
  	}
  	if (ha(i,jp)>=hd){
  		q += 1;
  	}

	return q;
}

//////////////////////////////////////////////////////////////// 
int computeclassa(int i, int j){//give the class of (i,j). 
	int m,n,o,p,q,c,k,ek;
	
	m=computema(i,j);
	n=computena(i,j);
	o=computeoa(i,j);
	p=computepa(i,j);
	q=computeqa(i,j);
	
	k=ha(i,j);
	
	if (k<6 and k>=0){
		ek=k;
	}else{
		ek=5;
	}
	c=situation2class(m,n,o,p,q,ek);

	return c;
}

int computeclassb(int i, int j){//give the class of (i,j,k)
	int m,n,o,p,q,c,k,ek;
	
	m=computemb(i,j);
	n=computenb(i,j);
	o=computeob(i,j);
	p=computepb(i,j);
	q=computeqb(i,j);
	
	k=hb(i,j);
	
	if (k<6 and k>=0){
		ek=k;
	}else{
		ek=5;
	}
	c=situation2class(m,n,o,p,q,ek);

	return c;
}
////////////////////////////////////////////////////////////////	
void updatea(int i, int j){ //update the n,m,p,q of the site (i,j,k)
	int old_class,old_ord,old_num,old_i,old_j, c;
		
	// //raise (i,j,k) from old class and make the last point of that old class to replace the order of (i,j,k)
	
	//find the class and the order of (i,j,k) which will move
	old_class=classa(i,j);
	old_ord=ordera(i,j);
	
	//find the last point of that old class
	old_num=cardinala[old_class];
	
	old_i=listxa(old_num-1,old_class);
	old_j=listya(old_num-1,old_class);
	
	//put it into the old order of (i,j)
    
    listxa(old_ord,old_class) = old_i ;
    listya(old_ord,old_class) = old_j ;

	ordera(old_i,old_j)=old_ord;
	
	//delete the last point of the class
    listxa(old_num-1,old_class) = -1 ;
    listya(old_num-1,old_class) = -1 ;
    
    if (cardinala[old_class]==0){
        cout << "game over (update) , cas "<< old_class << " cardinala[old_class] " <<cardinala[old_class] << endl ;
        exit(1);
    }
    
	cardinala[old_class]-=1;

	// //put (i,j,k) into new class and give it new order

	//find the new class
	c = computeclassa(i,j) ;
	/*
	if (class2k[c]==1){
		cout<<"(a)find an atom whose ek = 1 !"<<endl;
	}
	*/
	//put (i,j,k) into new class
	classa(i,j)= c ;
	cardinala[c]+= 1 ;
	
	// new lists
	listxa(cardinala[c]-1,c) = i ;
	listya(cardinala[c]-1,c) = j ;
	
	//new order
	ordera(i,j)=cardinala[c]-1 ;
	
}

void updateb(int i, int j){ //update the n,m,p,q of the site (i,j,k)
	int old_class,old_ord,old_num,old_i,old_j, c;
		
	// //raise (i,j,k) from old class and make the last point of that old class to replace the order of (i,j,k)
	
	//find the class and the order of (i,j,k) which will move
	old_class=classb(i,j);
	old_ord=orderb(i,j);
	
	//find the last point of that old class
	old_num=cardinalb[old_class];
	
	old_i=listxb(old_num-1,old_class);
	old_j=listyb(old_num-1,old_class);
	
	//put it into the old order of (i,j,k)
    
    listxb(old_ord,old_class) = old_i ;
    listyb(old_ord,old_class) = old_j ;

	orderb(old_i,old_j)=old_ord;
	
	//delete the last point of the class
    listxb(old_num-1,old_class) = -1 ;
    listyb(old_num-1,old_class) = -1 ;
    
    if (cardinalb[old_class]==0){
        cout << "game over (update) , cas "<< old_class << " cardinalb[old_class] " <<cardinalb[old_class] << endl ;
        exit(1);
    }
    
	cardinalb[old_class]-=1;

	// //put (i,j,k) into new class and give it new order

	//find the new class
	c = computeclassb(i,j) ;
	/*
	if (class2k[c]==1){
		cout<<"(b)find an atom whose ek = 1 !"<<endl;
	}
	*/
	//put (i,j,k) into new class
	classb(i,j)= c ;
	cardinalb[c]+= 1 ;
	
	// new lists
	listxb(cardinalb[c]-1,c) = i ;
	listyb(cardinalb[c]-1,c) = j ;
	
	//new order
	orderb(i,j)=cardinalb[c]-1 ;
	
}
//////////	

void evaporatea(int i, int j){
  int im, jm, h;

  im = modN(i-1);
  jm = modN(j-1); 
  h=ha(i,j);
  
  if (h<=0){
  	cout<<"Problem of evaporatea, program is stopped!"<<endl;
  	exit(1);
  }
  	
  ha(i,j) = h-1 ;

  //update  (i,j)
  updatea(i,j);
  
  updateb(i,j);
  updateb(im,j);
  updateb(i,jm);

}

void evaporateb(int i, int j){
  int ip, jp, h;
 
  ip = modN(i+1);
  jp = modN(j+1);
  h=hb(i,j);
  
  if (h<=0){
  	cout<<"Problem of evaporateb, program is stopped!"<<endl;
  	exit(1);
  }
 
  hb(i,j) = h-1 ;

  //update  (i,j,k)
  updateb(i,j);
  
  updatea(i,j);
  updatea(ip,j);
  updatea(i,jp);

}

////////////////////////////////////////////////////////////////
void deposita(int i, int j){
  int im , jm , h;

  im = modN(i-1);
  jm = modN(j-1);
  h=ha(i,j);
  
  if (h>=nb_layerf){
  	cout<<"Problem of deposita, program is stopped!"<<endl;
  	exit(1);
  }
//  cout<<"deposita réussi"<<endl;

  ha(i,j)=h+1 ;
  
  //update  (i,j)
  updatea(i,j);
  
  updateb(i,j);
  updateb(im,j);
  updateb(i,jm);
  
}

void depositb(int i, int j){
  int ip, jp, h;
  
  ip = modN(i+1);
  jp = modN(j+1);
  
  h=hb(i,j);
  
  if (h>=nb_layerf){
  	cout<<"Problem of depositb, program is stopped!"<<endl;
  	exit(1);
  }

  hb(i,j)=h+1 ;
//  cout<<"depositb réussi"<<endl;
  
  //update  (i,j,k)
  updateb(i,j);
  
  updatea(i,j);
  updatea(ip,j);
  updatea(i,jp);

}

/////////////////////////////////////////////////////////////

void movea(int ord, int cas ){ 
	int newi , newj  ;
    int im, jm ;
	int i , j ,k, tk;
	int direction,o ;
	int num_dire=-1 ;
	
	//find which (i,j) will move
	i=listxa(ord,cas);
	j=listya(ord,cas);
	k=ha(i,j);
    if (cardinala[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    
		evaporatea(i,j);
	
       	im = modN(i-1);
    	jm = modN(j-1);
    	
    	o=class2o[cas];
    	
    	if (o > 0){
    		
    		direction = random (o); 
    	
    		tk=hb(i,j);
    		if (tk-k==-1){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=j;
    				goto herea ;
    			}
    		}
    		tk=hb(im,j);
    		if (tk-k==-1){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=im;
    				newj=j;
    				goto herea ;
    			}
    		}
    		tk=hb(i,jm);
    		if (tk-k==-1){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=jm;
    				goto herea ;
    			}
    		}
    		herea :
    		depositb(newi,newj);
    		    	
    	}else{
    		cout<<"Problem of movea, the program is stopped !"<<endl;
    		exit(1);
    	}
	
    }
}

void moveb(int ord, int cas ){ 
	int newi , newj ;
    int ip, jp ;
	int i , j ,k ,tk;
	int direction ,o;
	int num_dire=-1 ;

	//find which (i,j) will move
	i=listxb(ord,cas);
	j=listyb(ord,cas);
	k=hb(i,j);
	
    if (cardinalb[cas]==0){
        	cout << "game over, cas "<< cas << " cardinalb[cas] " << cardinalb[cas] << endl ;
        	exit(1);
    }else{
    
		evaporateb(i,j);
	
    	ip = modN(i+1);
    	jp = modN(j+1);
    	
    	o=class2o[cas];
    	
    	if (o > 0){
    		
    		direction = random (o);  
    		
    		tk=ha(i,j); 	
    		if (tk-k==-1){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=j;
    				goto hereb ;
    			}
    		}
    		
    		tk=ha(ip,j);
    		if (tk-k==-1 ){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=ip;
    				newj=j;
    				goto hereb ;
    			}
    		}
    		
    		tk=ha(i,jp);
    		if (tk-k==-1){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=jp;
    				goto hereb ;
    			}
    		}
    		hereb:
    		deposita(newi,newj);
    	
    	}else{
    		cout<<"Problem of moveb, the program is stopped !"<<endl;
    		exit(1);
    	}
    	
    }

}
//////move up a

void move_upa(int ord, int cas ){ 
	int newi , newj ;
    int im, jm ;
	int i , j ,k, tk;
	int direction,num_direction ;
	int num_dire=-1 ;
	
	//find which (i,j) will move
	i=listxa(ord,cas);
	j=listya(ord,cas);
	k=ha(i,j);
	
    if (cardinala[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    	num_direction=class2p[cas];
    	
		evaporatea(i,j);
	
       	im = modN(i-1);
    	jm = modN(j-1);
    	
    	if (num_direction > 0){
    		direction = random(num_direction);		
    		tk=hb(i,j);
    		if (tk-k==0){
    			num_dire += 1;
				if (num_dire == direction){
    				newi=i;
    				newj=j;
    				goto herea_up ;
    			}
    		}
			tk=hb(im,j);
  			if (tk-k==0){
   				num_dire += 1;
   				 if (num_dire == direction){
    				newi=im;
    				newj=j;
    				goto herea_up ;
				}
    		}
    		tk=hb(i,jm);
    		if (tk-k==0){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=jm;
    				goto herea_up ;
    			}
    		}
    		herea_up :
    		depositb(newi,newj);
    	}else{
    		cout<<"Problem of move_upa, the program is stopped !"<<endl;
    		cout<<" p = "<<num_direction<<endl;
    		exit(1);
    	}
	
    }
}

int num_up2a(int i, int j ,int k, int q){ 
    int im, jm ;
	int tk;
	int num_direction=0 ;
	
    im = modN(i-1);
    jm = modN(j-1);
    		
    tk=hb(i,j);
    if (tk-k==1){
    	num_direction += 1;
    }
    tk=hb(im,j);
    if (tk-k==1){
    	num_direction += 1;
    }
    tk=hb(i,jm);
    if (tk-k==1){
    	num_direction += 1;		
    }
    if (num_direction > q){
    	cout<<"Problem of direction number of upa! Program is stopped!"<<endl;
    	exit(1);
    }
	return num_direction;
}

void move_up2a(int ord, int cas ){ 
	int newi , newj ;
    int im, jm ;
	int i , j ,k, tk;
	int direction,num_direction, q ;
	int num_dire=-1 ;
	
	//find which (i,j) will move
	i=listxa(ord,cas);
	j=listya(ord,cas);
	k=ha(i,j);
	
    if (cardinala[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    	q=class2q[cas];
    	
		evaporatea(i,j);
	
       	im = modN(i-1);
    	jm = modN(j-1);
    	
    	if (q > 0){ 
    		num_direction = num_up2a(i,j,k,q);
    		if (num_direction > 0){
    			direction = random(num_direction);
    		
    			tk=hb(i,j);
    			if (tk-k==1){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=j;
    					goto herea_up ;
    				}
    			}
    			tk=hb(im,j);
    			if (tk-k==1){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=im;
    					newj=j;
    					goto herea_up ;
    				}
    			}
    			tk=hb(i,jm);
    			if (tk-k==1){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=jm;
    					goto herea_up ;
    				}
    			}
    			herea_up :
    			depositb(newi,newj);
    		}else{
    			deposita(i,j);
    			num_upfail += 1;
    			num_up -= 1;
    			if (k==1){
    			    num_up2_13 -= 1;
    			}
    		}
    	}else{
    		cout<<"Problem of move_upa, the program is stopped !"<<endl;
    		cout<<" q = "<<q<<endl;
    		exit(1);
    	}
	
    }
}
///////move up b

void move_upb(int ord, int cas ){ 
	int newi , newj ;
    int ip, jp ;
	int i , j ,k, tk;
	int direction,num_direction ;
	int num_dire=-1 ;
	
	//find which (i,j) will move
	i=listxb(ord,cas);
	j=listyb(ord,cas);
	k=hb(i,j);
	
    if (cardinalb[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    	num_direction=class2p[cas];
    	
		evaporateb(i,j);
	
       	ip = modN(i+1);
    	jp = modN(j+1);
    	
    	if (num_direction > 0){
    		direction = random (num_direction);
    		
    		tk=ha(i,j);
    		if (tk-k==0){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=j;
    				goto hereb_up ;
    			}
    		}
    		tk=ha(ip,j);
    		if (tk-k==0){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=ip;
    				newj=j;
    				goto hereb_up ;
    			}
    		}
    		tk=ha(i,jp);
    		if (tk-k==0){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=jp;
    				goto hereb_up ;
    			}
    		}
    		hereb_up :
    		deposita(newi,newj);
    		
    	}else{
    		cout<<"Problem of move_upb, the program is stopped !"<<endl;
    		cout<<" p = "<<num_direction<<endl;
    		exit(1);
    	}
	
    }
}

int num_up2b(int i, int j ,int k, int q){ 
    int ip, jp ;
	int tk;
	int num_direction=0 ;
	
    ip = modN(i+1);
    jp = modN(j+1);
    		
    tk=ha(i,j);
    if (tk-k==1){
    	num_direction += 1;
    }
    tk=ha(ip,j);
    if (tk-k==1){
    	num_direction += 1;
    }
    tk=ha(i,jp);
    if (tk-k==1){
    	num_direction += 1;		
    }
    if (num_direction > q){
    	cout<<"Problem of direction number of upb! Program is stopped!"<<endl;
    	exit(1);
    }
	return num_direction;
}

void move_up2b(int ord, int cas ){ 
	int newi , newj ;
    int ip, jp ;
	int i , j ,k, tk;
	int direction,num_direction, q ;
	int num_dire=-1 ;
	
	//find which (i,j) will move
	i=listxb(ord,cas);
	j=listyb(ord,cas);
	k=hb(i,j);
	
    if (cardinalb[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    	q=class2q[cas];
    	
		evaporateb(i,j);
	
       	ip = modN(i+1);
    	jp = modN(j+1);
    	
    	if (q > 0){
    		num_direction = num_up2b(i,j,k,q);
    		if (num_direction > 0){
    			direction = random (num_direction);
    		
    			tk=ha(i,j);
    			if (tk-k==1){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=j;
    					goto hereb_up ;
    				}
    			}
    			tk=ha(ip,j);
    			if (tk-k==1){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=ip;
    					newj=j;
    					goto hereb_up ;
    				}
    			}
    			tk=ha(i,jp);
    			if (tk-k==1){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=jp;
    					goto hereb_up ;
    				}
    			}
    			hereb_up :
    			deposita(newi,newj);
    		}else{
    			depositb(i,j);
    			num_upfail += 1;
    			num_up -= 1;
    			if (k==1){
    			    num_up2_13 -= 1;
    			}
    		}
    		
    	}else{
    		cout<<"Problem of move_upb, the program is stopped !"<<endl;
    		cout<<" q = "<<q<<endl;
    		exit(1);
    	}
	
    }
}
//////move down a
void move_downa(int ord, int cas ){ 
	int newi , newj ;
    int im, jm ;
	int i , j ,k, tk;
	int direction, num_direction ;
	int num_dire=-1 ;
	//find which (i,j) will move
	i=listxa(ord,cas);
	j=listya(ord,cas);
	k=ha(i,j);
	
    if (cardinala[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    	num_direction=class2n[cas];

		evaporatea(i,j);
	
       	im = modN(i-1);
    	jm = modN(j-1);
    	

    	if (num_direction > 0){
    		direction = random (num_direction);
    		
    		tk=hb(i,j);
    		if (tk-k==-2){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=j;
    				goto herea_down ;
    			}
    		}
    		tk=hb(im,j);
    		if (tk-k==-2){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=im;
    				newj=j;
    				goto herea_down ;
    			}	
    		}
    		tk=hb(i,jm);
    		if (tk-k==-2){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=jm;
    				goto herea_down ;
    			}
    		}
    		
    		herea_down :   		
    		depositb(newi,newj);	
    	
    	}else{
    		cout<<"Problem of move_downa, the program is stopped !"<<endl;
    		cout<<"n = "<<num_direction<<endl;
    		exit(1);
    	}
	
    }
}

int num_down2a(int i, int j ,int k, int m){ 
    int im, jm ;
	int tk;
	int num_direction=0 ;
	
    im = modN(i-1);
    jm = modN(j-1);
    		
    tk=hb(i,j);
    if (tk-k==-3){
    	num_direction += 1;
    }
    tk=hb(im,j);
    if (tk-k==-3){
    	num_direction += 1;
    }
    tk=hb(i,jm);
    if (tk-k==-3){
    	num_direction += 1;		
    }
    if (num_direction > m){
    	cout<<"Problem of direction number of downa! Program is stopped!"<<endl;
    	exit(1);
    }
	return num_direction;
}
/////////////
void move_down2a(int ord, int cas ){ 
	int newi , newj ;
    int im, jm ;
	int i , j ,k, tk;
	int direction, num_direction, m ;
	int num_dire=-1 ;
	//find which (i,j) will move
	i=listxa(ord,cas);
	j=listya(ord,cas);
	k=ha(i,j);
	
    if (cardinala[cas]==0){
        cout << "game over, cas "<< cas << " cardinala[cas] " << cardinala[cas] << endl ;
        exit(1);
    }else{
    	m=class2m[cas];

		evaporatea(i,j);
	
       	im = modN(i-1);
    	jm = modN(j-1);
    	
    	if (m > 0){ 
    		num_direction = num_down2a(i,j,k,m);
    		if (num_direction > 0){
    			direction = random (num_direction);
    		
    			tk=hb(i,j);
    			if (tk-k==-3){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=j;
    					goto herea_down ;
    				}
    			}
    			tk=hb(im,j);
    			if (tk-k==-3){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=im;
    					newj=j;
    					goto herea_down ;
    				}	
    			}
    			tk=hb(i,jm);
    			if (tk-k==-3){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=jm;
    					goto herea_down ;
    				}
    			}
    		
    		herea_down :   		
    		depositb(newi,newj);	
    		}else{
    			deposita(i,j);
    			num_downfail += 1;
    			num_down -= 1;
    			if (k==3){
    			    num_down2_13 -= 1;
    			}
    		}
    	}else{
    		cout<<"Problem of move_downa, the program is stopped !"<<endl;
    		cout<<"m = "<<m<<endl;
    		exit(1);
    	}
	
    }
}
//////move down b
void move_downb(int ord, int cas ){ 
	int newi , newj  ;
    int ip, jp ;
	int i , j ,k, tk;
	int direction,num_direction ;
	int num_dire=-1 ;
	//find which (i,j) will move
	i=listxb(ord,cas);
	j=listyb(ord,cas);
	k=hb(i,j);
		
    if (cardinalb[cas]==0){
        cout << "game over, cas "<< cas << " cardinalb[cas] " << cardinalb[cas] << endl ;
        exit(1);
    }else{
    	num_direction=class2n[cas];

		evaporateb(i,j);
	
       	ip = modN(i+1);
    	jp = modN(j+1);
    	
    	if (num_direction > 0){
    		direction = random (num_direction);
    		
    		tk=ha(i,j);
    		if (tk-k==-2){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=j;
    				goto hereb_down ;
    			}
    		}
    		tk=ha(ip,j);
    		if (tk-k==-2){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=ip;
    				newj=j;
    				goto hereb_down ;
    			}
    		}
    		tk=ha(i,jp);
    		if (tk-k==-2){
    			num_dire += 1;
    			if (num_dire == direction){
    				newi=i;
    				newj=jp;
    				goto hereb_down ;
    			}
    		}
    		
    		hereb_down :
    		deposita(newi,newj);   		
    			
    	}else{
    		cout<<"Problem of move_downb, the program is stopped !"<<endl;
    		cout<<"n = "<<num_direction<<endl;
    		exit(1);
    	}
	
    }
}


int num_down2b(int i, int j ,int k, int m){ 
    int ip, jp ;
	int tk;
	int num_direction=0 ;
	
    ip = modN(i+1);
    jp = modN(j+1);
    		
    tk=ha(i,j);
    if (tk-k==-3){
    	num_direction += 1;
    }
    tk=ha(ip,j);
    if (tk-k==-3){
    	num_direction += 1;
    }
    tk=ha(i,jp);
    if (tk-k==-3){
    	num_direction += 1;		
    }
    if (num_direction > m){
    	cout<<"Problem of direction number of downa! Program is stopped!"<<endl;
    	exit(1);
    }
	return num_direction;
}
void move_down2b(int ord, int cas ){ 
	int newi , newj  ;
    int ip, jp ;
	int i , j ,k, tk;
	int direction,num_direction,m ;
	int num_dire=-1 ;
	//find which (i,j) will move
	i=listxb(ord,cas);
	j=listyb(ord,cas);
	k=hb(i,j);
		
    if (cardinalb[cas]==0){
        cout << "game over, cas "<< cas << " cardinalb[cas] " << cardinalb[cas] << endl ;
        exit(1);
    }else{
    	m=class2m[cas];

		evaporateb(i,j);
	
       	ip = modN(i+1);
    	jp = modN(j+1);
    	
    	if (m > 0){ 
    		num_direction = num_down2b(i,j,k,m);
    		if (num_direction > 0){
    			direction = random (num_direction);
    		
    			tk=ha(i,j);
    			if (tk-k==-3){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=j;
    					goto hereb_down ;
    				}
    			}
    			tk=ha(ip,j);
    			if (tk-k==-3){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=ip;
    					newj=j;
    					goto hereb_down ;
    				}
    			}
    			tk=ha(i,jp);
    			if (tk-k==-3){
    				num_dire += 1;
    				if (num_dire == direction){
    					newi=i;
    					newj=jp;
    					goto hereb_down ;
    				}
    			}
    		
    		hereb_down :
    		deposita(newi,newj);   		
    		}else{
    			depositb(i,j);
    			num_downfail += 1;
    			num_down -= 1;
    			if (k==3){
    			    num_down2_13 -= 1;
    			}
    		}	
    	}else{
    		cout<<"Problem of move_downb, the program is stopped !"<<endl;
    		cout<<"m = "<<m<<endl;
    		exit(1);
    	}
	
    }
}

////////////////////
int minhmaxa(void){
	int i , j , minh ;
	minh=hl ;
	for (i=0;i<NN;i++){
		for (j=0;j<NN;j++){
			if (ha(i,j)<minh){
				minh=ha(i,j) ;
			}
		}
	}
	return minh ;
}
int minhmaxb(void){
	int i , j , minh ;
	minh=hl ;
	for (i=0;i<NN;i++){
		for (j=0;j<NN;j++){
			if (hb(i,j)<minh){
				minh=hb(i,j) ;
			}
		}
	}
	return minh ;
}
int minhmax(void){
	int mina,minb;
	mina=minhmaxa();
	minb=minhmaxb();
	if (mina<minb){
		return mina;
	}else{
		return minb;
	}
}
////////////////////
int maxhmaxa(void){
	int i , j , maxh ;
	maxh=0 ;
	for (i=0;i<NN;i++){
		for (j=0;j<NN;j++){
			if (ha(i,j)>maxh){
				maxh=ha(i,j) ;
			}
		}
	}
	return maxh ;
}
int maxhmaxb(void){
	int i , j , maxh ;
	maxh=0 ;
	for (i=0;i<NN;i++){
		for (j=0;j<NN;j++){
			if (hb(i,j)>maxh){
				maxh=hb(i,j) ;
			}
		}
	}
	return maxh ;
}
int maxhmax(void){
	int maxa,maxb;
	maxa=maxhmaxa();
	maxb=maxhmaxb();
	if (maxa<maxb){
		return maxb;
	}else{
		return maxa;
	}
}

/////////////////////////////////////////////////////////////cf KMC
int modN(int k){
    
    int resu;
    
    resu = k - NN*floor(1.*k/NN) ;
    
    return resu;
}
//////////////////////////////////////////////////////////////


/////////////////////////cf KMC
double modNr(double x){ //periodic condition
    
    double resu;
    
    resu = x - NN*floor(1.*x/NN) ;
    
    return resu;
}
///////////////////
int nbcolh(int h){
	int i,j;
	int nb=0;
	for (i=0;i<NN;i++){
		for (j=0;j<NN;j++){
			if (ha(i,j)==h){
				nb += 1;
			}
			if (hb(i,j)==h){
				nb += 1;
			}
		}
	}
	return nb;
}
///////////
int nbcoltot(void){
	int i,j;
	int nb=0;
	for (i=0;i<NN;i++){
		for (j=0;j<NN;j++){
			if (ha(i,j)>0){
				nb += 1;
			}
			if (hb(i,j)>0){
				nb += 1;
			}
		}
	}
	return nb;
}
/*
int gradingcondition_k(int i, int k){
	int nk;
	if (i==-1){
		nk=k-1;
	}else if (i==NN){
		nk=k+1;
	}else{
		nk=k;
	}
	
	return nk;
}
*/
	
