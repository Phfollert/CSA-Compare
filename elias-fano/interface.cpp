#include <stdlib.h>
#include "comparray4.cpp"
#include "qsufsort.cpp"
#include <string.h>

#ifndef uchar
#define uchar unsigned char
#endif
#ifndef ulong
#define ulong unsigned long
#endif


using namespace sdsl;

/* Three function to variables to manage parameters */
static bool is_delimeter(char *delimiters, char c) {
    int i = 0, len_delimiters = strlen(delimiters);
    bool is = false;
    for (i = 0; i < len_delimiters; i++)
        if (c == delimiters[i]) is = true;
    return is;
}

static void parse_parameters(char *options, int *num_parameters, char ***parameters, char *delimiters) {
    int i = 0, j = 0, temp = 0, num = 0, len_options = strlen(options);
    char *options_temp;
    while (i < len_options) {
        while ((i < len_options) && is_delimeter(delimiters, options[i])) i++;
        temp = i;
        while ((i < len_options) && !is_delimeter(delimiters, options[i])) i++;
        if (i != temp) num++;
    }
    (*parameters) = (char **) malloc(num * sizeof(char *));
    i = 0;
    while (i < len_options) {
        while ((i < len_options) && is_delimeter(delimiters, options[i])) i++;
        temp = i;
        while ((i < len_options) && !is_delimeter(delimiters, options[i])) i++;
        if (i != temp) {
            (*parameters)[j] = (char *) malloc((i - temp + 1) * sizeof(char));
            options_temp = options + temp;
            strncpy((*parameters)[j], options_temp, i - temp);
            ((*parameters)[j])[i - temp] = '\0';
            j++;
        }
    }
    *num_parameters = num;
}

static void free_parameters(int num_parameters, char ***parameters) {
    int i = 0;
    for (i = 0; i < num_parameters; i++)
        free((*parameters)[i]);
    free((*parameters));
}


/*////////////////////
//Building the Index//
////////////////////*/

int build_index(uchar *text, ulong length, char *build_options, CSA *csa) {
    char delimiters[] = " =;";
    int j, num_parameters;
    char **parameters;
    int rankb_w = 32*8;
    int free_text = false; /* don't free text by default */
    if (build_options != NULL) {
        parse_parameters(build_options, &num_parameters, &parameters, delimiters);
        for (j = 0; j < num_parameters; j++) {
            if ((strcmp(parameters[j], "samplerate") == 0) && (j < num_parameters - 1)) {
                rankb_w = atoi(parameters[j + 1])*8;
                j++;
            } else if (strcmp(parameters[j], "free_text") == 0)
                free_text = true;
        }
        free_parameters(num_parameters, &parameters);
    }

    int n = length;

    /* make the SA */
    int i, *x, *p;
    int k, l;
    p = (int *) malloc((n + 1) * sizeof *p);
    x = (int *) malloc((n + 1) * sizeof *x);
    if (!p || !x) {
        return 1;
    }
    for (i = 0; i < n; i++) {
        x[i] = text[i];
    }
    l = 0;
    k = UCHAR_MAX + 1;

    suffixsort(x, p, n, k, l);
    free(x);
    p[0] = n;
    /* End Make SA */
    for (i = 0; i <= n; ++i) p[i]++; //1 to n+1 instead of 0 to n

    create_csa(n, p, text, rankb_w, csa);
    free(p);
    if (free_text) free(text);
    return 0;
}

void free_index(CSA *csa){
    free(csa->SA);
    free(csa->ISA);
    free(csa->C2);
    free(csa->S);
	delete[] csa->ranks;
	delete[] csa->selects;
	delete[] csa->vectors;
    free(csa);
}

int index_size(CSA* csa, double *size) {
    ulong ans = sizeof(csa);
    int m = csa->m;
    int sample = csa->sample;
    int n = csa->n;
    ans += size_in_bytes(csa->rankB) + size_in_bytes(csa->Bv)+size_in_bytes(csa->Dv)+size_in_bytes(csa->rankD);
    ans += (m+1)*(sizeof(unsigned char)) + (m+2)*sizeof(int); //S & C2
    ans += 2*(n/sample+1)*sizeof(int); // SA & ISA
    for (int i = 0; i<=m; i++){
        ans += size_in_bytes(csa->vectors[i])+size_in_bytes(csa->ranks[i])+size_in_bytes(csa->selects[i]);
    }
    *size = 1.0*ans/(1024*1024);
    return 0;
}
/*////////////////////
//Querying the Index//
////////////////////*/
int count(CSA *csa, uchar *pattern, ulong length, ulong *numocc) {
    int l, r, len;
    len = csa_bsearch(pattern, length, csa, &l, &r);
    *numocc = r - l + 1;
    return 0;
}

int locate(CSA *csa, uchar *pattern, ulong length, ulong **occ, ulong *numocc) {
    int l, r, len;
    csa_bsearch(pattern, length, csa, &l, &r);
    *numocc = r - l + 1;
    (*occ) = csa_locate(csa, l, r);
    return 0;
}

/*///////////////////////
//Accessing the index//
///////////////////////*/
char *extract(CSA *csa, ulong from, ulong to) {
    ulong n = csa->n;
    int sample = csa->sample;
    int pos, k;
    if (to >= n) to = n;
    ulong j;
    //cout << "from " << from << endl;
    char* ans = (char*) malloc (sizeof(char)*(to-from+1));
    j = (from + 1) / sample * sample;

    //cout << "from inicial " << from << endl;
    pos = csa->ISA[j/sample];
    //cout << "pos inicial " << pos << endl;

    if(j == 0) j = 1;

    for (k = j; k < from+1; k++) {
        pos = csa_psi(csa, pos);
        //cout << "pos nueva " << pos << endl;
    }



    //cout << "where" << endl;
    for (k = from; k <= to; k++) {
        ans[k - from] = csa->S[csa->rankD(pos+1)-1];
        //cout << "pos " << pos << endl;
        //cout << ans[k-from] << endl;
        pos = csa_psi(csa, pos);
    }

    return ans;
}
