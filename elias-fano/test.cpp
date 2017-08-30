#include "interface.cpp"

struct timespec start, finish;
double elapsed;

using namespace std;

int main(){
	int code;
    CSA* csa;
    csa = (CSA *) malloc(sizeof(CSA));
	uchar myString [] = "abracadabra";
	uchar* tmpBuffer = &myString[0];
	code = build_index(tmpBuffer, 11, NULL, csa);
	ulong numocc;
	ulong temp;
    ulong *occ;
    
	//locate and count tests
	uchar pattern [] = "a";
    locate(csa, &pattern[0], 1, &occ, &numocc);
	if (numocc != 5) cout << "count 1 failed" << endl;
	temp = numocc;
	for (int i=0; i<numocc; i++){
		if (occ[i]==0) temp--;
		if (occ[i]==3) temp--;
		if (occ[i]==5) temp--;
		if (occ[i]==7) temp--;
		if (occ[i]==10) temp--;	
	}
	if (temp != 0 ) cout << "locate 1 failed" << endl;

	uchar pattern2 [] = "abra";
    locate(csa, &pattern2[0], 4, &occ, &numocc);
	if (numocc != 2) cout << "count 2 failed" << endl;
	temp = numocc;
	for (int i=0; i<numocc; i++){
		if (occ[i]==0) temp--;
		if (occ[i]==7) temp--;
	}
	if (temp != 0 ) cout << "locate 2 failed" << endl;

	uchar pattern3 [] = "r";
    locate(csa, &pattern3[0], 1, &occ, &numocc);
	if (numocc != 2) cout << "count 3 failed" << endl;
	temp = numocc;
	for (int i=0; i<numocc; i++){
		if (occ[i]==2) temp--;
		if (occ[i]==9) temp--;
	}
	if (temp != 0 ) cout << "locate 3 failed" << endl;

	uchar pattern4 [] = "abracadabra";
    locate(csa, &pattern4[0], 11, &occ, &numocc);	
	if (numocc != 1) cout << "count 4 failed" << endl;
	temp = numocc;
	for (int i=0; i<numocc; i++){
		if (occ[i]==0) temp--;
	}
	if (temp != 0 ) cout << "locate 4 failed" << endl;

	uchar pattern5 [] = "z";
    locate(csa, &pattern5[0], 1, &occ, &numocc);
	if (numocc != 0) cout << "count 5 failed" << endl;

	string ext;
	for (int i = 0; i<11; i++){
		for(int j = i; j<11; j++){
			ext = extract(csa, i, j); 
			for(int k = i; k<=j; k++){
				if (ext[k-i]!=myString[k]){
					cout << myString[k] << endl;
					cout << "extract failed" << endl;
					return 1;
				}
			}
		}		
	}
	return 0;
}

