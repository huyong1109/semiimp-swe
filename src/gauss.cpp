#include<iostream>
#include<iomanip>
using namespace std;
int main()
{
    int ie = 6;
    int iee = ie-1;

    double a[5][ie], b[ie],x[ie];
    for (int i = 0; i<ie; ++i){
	a[0][i] = 1.0;
	a[1][i] = 3.0;
	a[2][i] = 1.0;
	a[3][i] = 0.0;
	a[4][i] = 0.0;
	b[i] = 5.0;
	x[i] = 0.0;
    }
	//a[4][0] = a[0][0];
	//a[0][0] = 0.0;
	//a[3][0] = a[2][iee];
	//a[2][iee] =0.0;
	//a[3][iee-1] = a[0][iee];
	//a[4][iee-1] = a[2][iee-1];


	cout << "========================="<<endl;
	for (int i = 0; i< ie; ++i){ 
	    for (int j = 0; j< 5; ++j) {
		 cout << setw(10)<<a[j][i] << "\t" ;
	    }
	cout << "|| " <<setw(10)<<b[i] <<"\t"<<setw(10)<< "|| " <<x[i] <<endl;	 
	}

	double y;
	a[3][0]= a[2][iee];
	a[4][0]= a[0][0];
	a[4][iee-1]= a[2][iee-1];
	a[3][iee-1]= a[0][iee];
        for (int i = 1; i < ie-1; ++i) {
	   y = a[0][i]/a[1][i-1];
	   a[0][i] = 0.0;
	   a[1][i] = a[1][i]-a[2][i-1]*y;
	   a[4][i] = a[4][i] -a[4][i-1]*y;
	   b[i]= b[i]- b[i-1]*y;

	   y = a[3][i-1]/a[1][i-1];
	   a[3][i] = a[3][i] -a[2][i-1]*y;
	   a[1][iee]= a[1][iee] - a[4][i-1]*y;
	   b[iee] = b[iee]-b[i-1]*y;
	   
	
	cout << "========================="<<endl;
	for (int i = 0; i< ie; ++i){ 
	    for (int j = 0; j< 5; ++j) {
		 cout <<setw(10)<< a[j][i] << "\t" ;
	    }
	cout << "|| "<<setw(10)<<b[i] <<"\t"<<setw(10)<< "|| " <<x[i] <<endl;	 
	}
        }
    
	a[1][iee]= a[1][iee] - (a[4][iee-1])*(a[3][iee-1])/a[1][iee-1];
	b[iee]= b[iee] - b[iee-1]*(a[3][iee-1])/a[1][iee-1];

	x[iee] = b[iee]/a[1][iee];
	
	cout << "========xn================="<<endl;
	for (int i = 0; i< ie; ++i){ 
	    for (int j = 0; j< 5; ++j) {
		 cout <<setw(10)<< a[j][i] << "\t" ;
	    }
	cout << "|| "<<setw(10)<<b[i] <<"\t"<<setw(10)<< "|| " <<x[i] <<endl;	 
	}

	x[iee-1] = (b[iee-1]-a[4][iee-1]*x[iee])/a[1][iee-1];
	cout << "=========xn-1================"<<endl;
	for (int i = 0; i< ie; ++i){ 
	    for (int j = 0; j< 5; ++j) {
		 cout <<setw(10)<< a[j][i] << "\t" ;
	    }
	cout << "|| "<<setw(10)<<b[i] <<"\t"<<setw(10)<< "|| " <<x[i] <<endl;	 
	}
		
	for (int i = ie-3; i>=0; --i) {
		x[i]= (b[i]-a[2][i]*x[i+1]-a[3][i]*x[iee])/a[1][i];
	cout << "============x_" << i<<"============="<<endl;
	for (int i = 0; i< ie; ++i){ 
	    for (int j = 0; j< 5; ++j) {
		 cout <<setw(10)<< a[j][i] << "\t" ;
	    }
	cout << "|| "<<setw(10)<<b[i] <<"\t"<<setw(10)<< "|| " <<x[i] <<endl;	 
	}
	}
}

