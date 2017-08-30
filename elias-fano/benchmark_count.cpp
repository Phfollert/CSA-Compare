
#include "interface.cpp"
#include "iostream"

struct timespec start, finish;
using namespace std;

#define times 10

double diff(timespec start, timespec end, ulong numocc) //time in seconds
{
    double el = (end.tv_sec - start.tv_sec)*1000000000;
    el += (end.tv_nsec - start.tv_nsec);
    return el/numocc;
}

double error(vector<double> values, double mean){
  int size = times;
  double error = 0;
  for(int i = 0; i<size; i++){
    error += pow(values[i]-mean,2);
  }
  error = sqrt(error/(size-1))/sqrt(size);
  return error;
}

int main(){
    CSA* csa;
    double size;
    unsigned char* buffer;
    unsigned char pattern;
    int fd,r,i,j;
    ulong numocc;
    ulong *occ;
	ofstream result1,result2,result3,result4,result5,result6,error1,error2,error3,error4,error5,error6;
  	result1.open("dna50mb_count.txt");
  	result2.open("dna100mb_count.txt");
  	result3.open("dna200mb_count.txt");
  	result4.open("english50mb_count.txt");
  	result5.open("english100mb_count.txt");
  	result6.open("english200mb_count.txt");
  	error1.open("dna50mb_count_error.txt");
  	error2.open("dna100mb_count_error.txt");
  	error3.open("dna200mb_count_error.txt");
  	error4.open("english50mb_count_error.txt");
  	error5.open("english100mb_count_error.txt");
  	error6.open("english200mb_count_error.txt");
	double elapsed,m1,m2,m3,m4,m5,m6,e1,e2,e3,e4,e5,e6;
	vector<double> v1(times);
	vector<double> v2(times);
	vector<double> v3(times);
	vector<double> v4(times);
	vector<double> v5(times);
	vector<double> v6(times);
    //dna 50
    csa = (CSA *) malloc(sizeof(CSA));
    buffer = (unsigned char*)malloc(52428800);
    fd = open ("../dna.50MB", O_RDONLY);
    r = read(fd, buffer, 52428800);
    close(fd);
    build_index(buffer, 52428800, NULL, csa);
    index_size(csa, &size);
    cout << "50MB dna" << endl;
    cout << "size in MB" << endl << size << endl;
    for (i=2000; i<2100; i++){
		m1 = 0;
        for(j=0; j<times; j++){
			clock_gettime(CLOCK_MONOTONIC, &start);
            count(csa, buffer, i+1, &numocc);			
	        clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = diff(start,finish,numocc);
			m1 += elapsed;
			v1[j] = elapsed;
        }
		m1 = m1/times;
		e1 = error(v1,m1);
		result1 << m1;		
    	error1 << e1;
		if (i!=2099){
			result1 << ',';
			error1 << ',';		
		}
    }
    free(buffer);
    free(csa);

    //dna 100
    csa = (CSA *) malloc(sizeof(CSA));
    buffer = (unsigned char*)malloc(52428800*2);
    fd = open ("../dna.100MB", O_RDONLY);
    r = read(fd, buffer, 52428800*2);
    close(fd);
    build_index(buffer, 52428800*2, NULL, csa);
    index_size(csa, &size);
    cout << "100MB dna" << endl;
    cout << "size in MB" << endl << size << endl;
    for (i=2000; i<2100; i++){
		m2 = 0;
        for(j=0; j<times; j++){
			clock_gettime(CLOCK_MONOTONIC, &start);
            count(csa, buffer, i+1, &numocc);
			clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = diff(start,finish,numocc);
			m2 += elapsed;
			v2[j] = elapsed;
        }
		m2 = m2/times;
		e2 = error(v2,m2);
		result2 << m2;		
    	error2 << e2;
		if (i!=2099){
			result2 << ',';
			error2 << ',';		
		}
    }
    free(buffer);
    free(csa);

    //dna 200
    csa = (CSA *) malloc(sizeof(CSA));
    buffer = (unsigned char*)malloc(52428800*4);
    fd = open ("../dna.200MB", O_RDONLY);
    r = read(fd, buffer, 52428800*4);
    close(fd);
    build_index(buffer, 52428800*4, NULL, csa);
    index_size(csa, &size);
    cout << "200MB dna" << endl;
    cout << "size in MB" << endl << size << endl;
    for (i=2000; i<2100; i++){
		m3 = 0;
        for(j=0; j<times; j++){			
        	clock_gettime(CLOCK_MONOTONIC, &start);
            count(csa, buffer, i+1, &numocc);
	        clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = diff(start,finish,numocc);
			m3 += elapsed;
			v3[j] = elapsed;
        }
		m3 = m3/times;
		e3 = error(v3,m3);
		result3 << m3;		
    	error3 << e3;
		if (i!=2099){
			result3 << ',';
			error3 << ',';		
		}
    }
    free(buffer);
    free(csa);

    //english 50
    csa = (CSA *) malloc(sizeof(CSA));
    buffer = (unsigned char*)malloc(52428800);
    fd = open ("../english.50MB", O_RDONLY);
    r = read(fd, buffer, 52428800);
    close(fd);
    build_index(buffer, 52428800, NULL, csa);
    index_size(csa, &size);
    cout << "50MB english" << endl;
    cout << "size in MB" << endl << size << endl;
    for (i=2000; i<2100; i++){
		m4 = 0;
        for(j=0; j<times; j++){
        	clock_gettime(CLOCK_MONOTONIC, &start);
            count(csa, buffer, i+1, &numocc);			
        	clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = diff(start,finish,numocc);
			m4 += elapsed;
			v4[j] = elapsed; 
		}
		m4 = m4/times;
		e4 = error(v4,m4);
		result4 << m4;		
    	error4 << e4;
		if (i!=2099){
			result4 << ',';
			error4 << ',';		
		}
    }
    free(buffer);
    free(csa);

    //english 100
    csa = (CSA *) malloc(sizeof(CSA));
    buffer = (unsigned char*)malloc(52428800*2);
    fd = open ("../english.100MB", O_RDONLY);
    r = read(fd, buffer, 52428800*2);
    close(fd);
    build_index(buffer, 52428800*2, NULL, csa);
    index_size(csa, &size);
    cout << "100MB english" << endl;
    cout << "size in MB" << endl << size << endl;
    for (i=2000; i<2100; i++){
		m5 = 0;
        for(j=0; j<times; j++){
        	clock_gettime(CLOCK_MONOTONIC, &start);
            count(csa, buffer, i+1, &numocc);
        	clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = diff(start,finish,numocc);
			m5 += elapsed;
			v5[j] = elapsed;
        }
		m5 = m5/times;
		e5 = error(v5,m5);
		result5 << m5;		
    	error5 << e5;
		if (i!=2099){
			result5 << ',';
			error5 << ',';		
		}
    }
    free(buffer);
    free(csa);

    //english 200
    csa = (CSA *) malloc(sizeof(CSA));
    buffer = (unsigned char*)malloc(52428800*4);
    fd = open ("../english.200MB", O_RDONLY);
    r = read(fd, buffer, 52428800*4);
    close(fd);
    build_index(buffer, 52428800*4, NULL, csa);
    index_size(csa, &size);
    cout << "200MB english" << endl;
    cout << "size in MB" << endl << size << endl;
    for (i=2000; i<2100; i++){
		m6 = 0;
        for(j=0; j<times; j++){			
        	clock_gettime(CLOCK_MONOTONIC, &start);
            count(csa, buffer, i+1, &numocc);
        	clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = diff(start,finish,numocc);
			m6 += elapsed;
			v6[j] = elapsed;
        }
		m6 = m6/times;
		e6 = error(v6,m6);
		result6 << m6;		
    	error6 << e6;
		if (i!=2099){
			result6 << ',';
			error6 << ',';		
		}
    }
    free(buffer);
    free(csa);

	
	result1.close();
	result2.close();
	result3.close();
	result4.close();
	result5.close();
	result6.close();
	error1.close();
	error2.close();
	error3.close();
	error4.close();    
	error5.close();
	error6.close();
	return 0;
}

