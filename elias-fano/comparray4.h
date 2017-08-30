/* comparray4.cpp
   Copyright (C) 2005, K. Sadakane, all rights reserved.

   This file contains an implementation of CSA.
   For more information, see 

   K. Sadakane. Compressed text databases with efficient query
     algorithms based on the compressed suffix array.
     In Proceedings 11th Annual International Symposium on Algorithms
     and Computation (ISAAC)}, LNCS v. 1969, pages 410--421, 2000.

   K. Sadakane. Succinct representations of lcp information and 
     improvements in the compressed suffix arrays.
     In Proceedings 13th Annual ACM-SIAM Symposium on Discrete
     Algorithms (SODA), 2002.

   K. Sadakane. New text indexing functionalities of the compressed
     suffix arrays. Journal of Algorithms, 48(2):294--313, 2003.


   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>

#define SIGMA 256

using namespace sdsl;

typedef struct csa {
    int m, n;
    int sample;
    unsigned char map[SIGMA];
    int *SA, *ISA, *C2;
    unsigned char *S;
    rrr_vector<63> Dv, Bv;
    rrr_vector<63>::rank_1_type rankD, rankB;
    sd_vector<>* vectors;
    rank_support_sd<1>* ranks;
    select_support_sd<1>* selects;
} CSA;

