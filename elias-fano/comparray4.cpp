#include <stdio.h>
#include <stdlib.h>
#include "comparray4.h"
#include "iostream"
#include <algorithm>

#ifndef min
#define min(x, y) ((x)<(y)?(x):(y))
#endif

#define dprintf

using namespace sdsl;
using namespace std;

//Creates psi and assigns it to p, it also creates elias-fano bitvectors.
void psisort(int *p, unsigned char *s, int n, CSA* csa, bit_vector* bvs) {
    int i;
    int m = csa->m;
    int x, c;
    for (i = 0; i <= n; i++) {
        x = p[i] - 1;
        c = s[x];
        bvs[csa->map[c]][i] = 1;
        //cout << "rip" << endl;
    }

    csa ->vectors = new sd_vector<>[m+1];
    for( int i = 0; i<=m; i++){
        csa -> vectors[i] = sd_vector<>(bvs[i]);
        csa -> ranks[i] = rank_support_sd<1>(&csa->vectors[i]);
        csa -> selects[i] = select_support_sd<1>(&csa->vectors[i]);
    }
}

void create_csa(int n, int *p, unsigned char *s, int rankb_w, CSA *csa) {
    int i, v, b, x, b2, d, w, m, size;
    int C[SIGMA];
    unsigned char last;
    csa -> sample = rankb_w;
    csa -> n = n;
    int C2[SIGMA+2];
    bit_vector Dv(n+1,0);
    for (i = 0; i < SIGMA; ++i) {
        C[i] = 0;
    }

    Dv[0] = 1;
    last = s[p[1]-1];
    Dv[1] = 1;
    for (i=2; i<=n; i++){
        if(last!=s[p[i]-1]) {
            last = s[p[i]-1];
            Dv[i]=1;
        }
    }
    //cout << Dv << endl;
    m = 0;
    for (i = 0; i < n; ++i) {
        if (C[s[i]]++ == 0){
            m++;
        }
    }

    unsigned char* S = (unsigned char *)malloc((m+1)*sizeof(unsigned char));

    csa->Dv = rrr_vector<63>(Dv);
    csa->rankD = rrr_vector<63>::rank_1_type(&(csa->Dv));

    S[0] = 0;
    for (m = 0, v = 1, i = 0; i < SIGMA; i++) {
        if (C[i] > 0) {
            m++;
            S[m] = i;
            C2[m] = v;
            v += C[i];
        }
    }
    C2[m + 1] = v;

    csa->C2 = (int *) malloc ((m+2)*sizeof(int));

    for(int i =1; i<=m+1; i++){
        csa->C2[i]=C2[i];
    }

    for (v = 0, i = 0; i < SIGMA; i++) {
        v = v + C[i];
        C[i] = v;
    }

    bit_vector bvs[m+1];

    csa -> ranks = new rank_support_sd<1>[m+1];
    csa -> selects = new select_support_sd<1>[m+1];

    for(i = 0; i<m; i++){
        for(int k=S[i];k<S[i+1]; k++){
            csa->map[k] = i;
        }
        bvs[i] = bit_vector(n+1,0);
    }

    for (int k = S[m]; k<=SIGMA; k++){
        csa->map[k] = m;
    }

    bvs[m] = bit_vector(n+1,0);

    csa->m = m;


    csa->SA = (int *) malloc((n/rankb_w+1) * sizeof(int));
    csa->ISA = (int *) malloc((n/rankb_w+1) * sizeof(int));

    psisort(p, s - 1, n, csa, bvs);

    csa->ISA[0] = csa->selects[0](1);
    csa->SA[0] = n;

    if (csa->ISA == NULL) {
        perror("ISA malloc failled\n");
        exit(1);
    }
    int cont = 1;

    bit_vector B(n+1, 0);

    for (i = 0; i <= n; i++) {
        if (p[i] % rankb_w == 0) {
            csa->ISA[p[i]/rankb_w] = i;
            B[i] = 1;
            if (p[i] != n+1) csa->SA[cont++] = p[i]-1;
        }
    }

    B[0] = 1;
    csa->Bv = rrr_vector<63>(B);
    csa->rankB = rrr_vector<63>::rank_1_type(&(csa->Bv));
    csa->S = S;
}

inline
int csa_psi(CSA *csa, int i) {
    //cout << i << endl;
    int pos = csa->rankD(i+1)-1;
    //cout << pos << endl;
    int lessthan = csa->C2[pos];
    //cout << i << " " << lessthan << endl;
    int ret = csa->selects[pos](i-lessthan+1);
    //cout << ret << endl;
    //cout << "no" << endl;
    return ret;
}

ulong *csa_locate(CSA *csa,int l, int r)
{
    ulong *I;
    int k,i, j;
    int sample;
    int *sa;

    sample = csa->sample;
    sa = csa->SA;

    //cout << endl;
    //cout << csa->Bv << endl;

    I = (ulong 	*)malloc((r-l+1)*sizeof(*I));

    for (i = l; i<=r; i++){
        k = 0;
        j = i;
        while (csa->Bv[j] == 0) {
            j = csa_psi(csa, j);
            k++;
        }
        I[i-l] = sa[csa->rankB(j+1)-1] - k;
    }
    return I;
}

/* backward search */
int csa_bsearch(unsigned char *key, int keylen, CSA *csa, int *li, int *ri) {
    int m = csa->m;
    int r,l,ltmp,rtmp, pos;
    unsigned char c = key[keylen-1];
    //cout << "c: " << c << endl;
    /*for (int i = 0; i < SIGMA+1; i++){
        cout << csa->map[i] << " ";
    }
    cout << endl;*/

    l = csa->C2[csa->map[c-1]+1];
    r = csa->C2[csa->map[c]+1]-1;
    /*for (int i = 0; i<SIGMA; i++){
        cout << csa->C[i] << " ";
    }

    cout << endl;*/

    for (int k = keylen-2; k>=0; k--) {
        if (l > r) {
            cout << "pattern doesn't exist" << endl;
            *li = 1;
            *ri = 0;
            return 0;
        }
        c = key[k];
        ltmp = csa->ranks[csa->map[c]](l);
        rtmp = csa->ranks[csa->map[c]](r+1)-1;
        l = ltmp + csa->C2[csa->map[c-1]+1];
        r = rtmp + csa->C2[csa->map[c-1]+1];
    }
    *ri=r;
    *li=l;
}

