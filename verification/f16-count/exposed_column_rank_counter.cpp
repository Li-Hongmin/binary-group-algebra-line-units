// Exact conditional-rank enumeration for Section 8 of the manuscript.
// C++17, standard library only, one CPU, O(q^2) storage; q in {4,16}.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using U = std::uint64_t;
using Wide = __uint128_t;
struct Fraction {
    U n=1, d=1;
    Fraction(U a=1, U b=1): n(a), d(b) {
        U g=std::gcd(n,d); n/=g; d/=g;
    }
    Fraction times(Fraction x) const { return Fraction(n*x.n,d*x.d); }
    bool operator==(Fraction x) const { return Wide(n)*x.d==Wide(x.n)*d; }
};
std::ostream& operator<<(std::ostream& o,Fraction x) {
    return o << "{\"numerator\":" << x.n << ",\"denominator\":" << x.d << "}";
}
template<class T> void array_json(const std::vector<T>& x) {
    std::cout << '[';
    for (std::size_t i=0;i<x.size();++i) { if(i) std::cout << ','; std::cout << x[i]; }
    std::cout << ']';
}
void require(bool p,const char* message) { if(!p) throw std::runtime_error(message); }
struct Field {
    int q, modulus;
    int mul[16][16]={}, inv[16]={};
    Field(int q_,int mod):q(q_),modulus(mod) {
        for(int a=0;a<q;++a) for(int b=0;b<q;++b) {
            int x=a,y=b,v=0;
            while(y) { if(y&1) v^=x; y>>=1; x<<=1; if(x&q) x^=modulus; }
            mul[a][b]=v;
        }
        for(int a=1;a<q;++a) {
            for(int b=1;b<q;++b) if(mul[a][b]==1) inv[a]=b;
            require(inv[a]!=0,"field inverse missing");
        }
    }
};
struct Term { int first, second; };
using Columns=std::vector<std::vector<Term>>;
// Original nonzero-g line sum; no affine formula is used to generate terms.
Columns line_terms(const Field& f,int d1,int d2,int relation,const std::vector<int>& exposed) {
    Columns out;
    for(int x:exposed) {
        std::vector<Term> col;
        for(int a=0;a<2;++a) if(relation&(1<<a))
            for(int g=1;g<f.q;++g) if((a^g)==x)
                col.push_back({a^f.mul[d1][g],a^f.mul[d2][g]});
        out.push_back(col);
    }
    return out;
}
void verify_affine_indices(const Field& f,int d1,int d2) {
    for(int a=0;a<2;++a) for(int x=0;x<f.q;++x) {
        auto t=line_terms(f,d1,d2,1<<a,{x});
        if(a==x) require(t[0].empty(),"nonzero diagonal");
        else {
            require(t[0].size()==1,"line coefficient not unique");
            require(t[0][0].first==(f.mul[1^d1][a]^f.mul[d1][x]),"first affine index mismatch");
            require(t[0][0].second==(f.mul[1^d2][a]^f.mul[d2][x]),"second affine index mismatch");
        }
    }
}
std::vector<int> exposed_set(const Field& f,int d1,int d2) {
    bool edge[16][16]={};
    int s=(d1!=1)+(d2!=1);
    for(int d:{d1,d2}) if(d!=1) {
        int k=f.mul[f.inv[d]][1^d];
        for(int x=0;x<f.q;++x) for(int a=0;a<2;++a) for(int b=0;b<2;++b)
            if(a!=b) edge[x][x^f.mul[k][a^b]]=true;
    }
    bool available[16]={};
    for(int x=2;x<f.q;++x) available[x]=true;
    std::vector<int> out;
    for(int x=2;x<f.q;++x) if(available[x]) {
        out.push_back(x); available[x]=false;
        for(int y=0;y<f.q;++y) if(edge[x][y]) available[y]=false;
    }
    for(int x:out) for(int y:out) require(!edge[x][y],"exposed set has forbidden edge");
    require(int(out.size())>=(f.q-2)/(2*s+1),"exposed set smaller than guaranteed");
    return out;
}
// Reference from independent local truth tables. Role-labelled variables in
// different columns must be disjoint before multiplying any probabilities.
Fraction truth_table_reference(const Columns& cols) {
    U used=0;
    Fraction out;
    for(const auto& col:cols) {
        std::vector<int> vars;
        for(auto t:col) { vars.push_back(t.first); vars.push_back(16+t.second); }
        std::sort(vars.begin(),vars.end());
        vars.erase(std::unique(vars.begin(),vars.end()),vars.end());
        U local=0;
        for(int v:vars) local|=U(1)<<v;
        require((used&local)==0,"cross-column variables collide"); used|=local;
        U zeros=0,total=U(1)<<vars.size();
        for(U b=0;b<total;++b) {
            int parity=0;
            for(auto t:col) {
                auto i=std::lower_bound(vars.begin(),vars.end(),t.first)-vars.begin();
                auto j=std::lower_bound(vars.begin(),vars.end(),16+t.second)-vars.begin();
                parity^=((b>>i)&1)&((b>>j)&1);
            }
            zeros+=parity==0;
        }
        out=out.times(Fraction(zeros,total));
    }
    return out;
}
int conditional_rank(const Columns& cols,unsigned a) {
    unsigned basis[16]={};
    int rank=0;
    for(const auto& col:cols) {
        unsigned row=0;
        for(auto t:col) if((a>>t.first)&1U) row^=1U<<t.second;
        while(row) {
            int p=31-__builtin_clz(row);
            if(basis[p]) row^=basis[p];
            else { basis[p]=row; ++rank; break; }
        }
    }
    return rank;
}
U direct_pairs(const Columns& cols,int q) {
    U zeros=0;
    for(unsigned a=0;a<(1U<<q);++a) for(unsigned b=0;b<(1U<<q);++b) {
        bool zero=true;
        for(const auto& col:cols) {
            int parity=0;
            for(auto t:col) parity^=((a>>t.first)&1U)&((b>>t.second)&1U);
            if(parity) { zero=false; break; }
        }
        zeros+=zero;
    }
    return zeros;
}
int main(int argc,char** argv) {
    int q=16, max_slopes=1000000;
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        if(arg=="--q" && i+1<argc) q=std::stoi(argv[++i]);
        else if(arg=="--max-slopes" && i+1<argc) max_slopes=std::stoi(argv[++i]);
        else { std::cerr << "usage: --q {4|16} [--max-slopes N]\n"; return 2; }
    }
    if(q!=4 && q!=16) { std::cerr << "unsupported q\n"; return 2; }
    auto start=std::chrono::steady_clock::now();
    Field f(q,q==16 ? 0x13 : 0x7);
    const int expected_slopes=(q-1)*(q-1)-1, expected_cases=3*expected_slopes;
    std::vector<std::array<int,3>> expected, completed;
    for(int d1=1;d1<q;++d1) for(int d2=1;d2<q;++d2) if(d1!=1||d2!=1)
        for(int c=1;c<4;++c) expected.push_back({d1,d2,c});
    int evaluated=0,validated=0,counterexamples=0,unknown=0,completed_slopes=0;
    bool stop=false;
    std::string failure;
    for(int d1=1;d1<q&&!stop;++d1) for(int d2=1;d2<q&&!stop;++d2) {
        if(d1==1&&d2==1) continue;
        if(completed_slopes>=max_slopes) { stop=true; break; }
        int s=(d1!=1)+(d2!=1),t=2-s;
        std::vector<int> exposed;
        try { verify_affine_indices(f,d1,d2); exposed=exposed_set(f,d1,d2); }
        catch(const std::exception& e) { failure=e.what(); stop=true; break; }
        int this_slope=0;
        for(int c=1;c<4;++c) {
            Columns cols=line_terms(f,d1,d2,c,exposed);
            Fraction reference;
            try { reference=truth_table_reference(cols); }
            catch(const std::exception& e) { failure=e.what(); stop=true; break; }
            std::vector<U> histogram(static_cast<std::size_t>(q+1),0);
            for(unsigned a=0;a<(1U<<q);++a) ++histogram[static_cast<std::size_t>(conditional_rank(cols,a))];
            U zeros=0;
            for(int r=0;r<=q;++r) zeros+=histogram[static_cast<std::size_t>(r)]*(U(1)<<(q-r));
            U total=U(1)<<(2*q);
            Fraction probability(zeros,total),bound;
            for(std::size_t k=0;k<exposed.size();++k) bound=bound.times(Fraction(3,4));
            bool inequality=Wide(zeros)*bound.d<=Wide(bound.n)*total;
            bool agreement=probability==reference;
            bool brute_agreement=q!=4||direct_pairs(cols,q)==zeros;
            std::string verdict;
            failure.clear();
            if(!agreement||!brute_agreement) {
                verdict="UNKNOWN"; ++unknown; stop=true;
                failure=!agreement?"rank/truth-table disagreement":"rank/direct-pair disagreement";
            } else if(!inequality) {
                verdict="COUNTEREXAMPLE"; ++counterexamples; stop=true;
                failure="exact probability exceeds bound";
            } else { verdict="VALIDATED"; ++validated; }
            ++evaluated; ++this_slope; completed.push_back({d1,d2,c});
            std::cout << "{\"kind\":\"case\",\"q\":"<<q<<",\"modulus\":"<<f.modulus
              <<",\"slope\":["<<d1<<','<<d2<<",1],\"E\":[0,1],\"relation\":["<<(c&1)<<','<<((c>>1)&1)
              <<"],\"s\":"<<s<<",\"t\":"<<t<<",\"exposed\":";
            array_json(exposed);
            std::cout << ",\"conditional_arrays\":"<<(U(1)<<q)<<",\"rank_histogram\":";
            array_json(histogram);
            std::cout << ",\"zero_pairs\":"<<zeros<<",\"total_pairs\":"<<total
              <<",\"probability\":"<<probability<<",\"reference_probability\":"<<reference
              <<",\"bound\":"<<bound<<",\"verdict\":\""<<verdict<<"\",\"failure\":";
            if(failure.empty()) std::cout << "null"; else std::cout << '"'<<failure<<'"';
            std::cout << "}\n" << std::flush;
            if(stop) break;
        }
        if(this_slope==3) ++completed_slopes;
    }
    int unresolved=expected_cases-evaluated+unknown;
    std::string verdict=counterexamples?"COUNTEREXAMPLE":(unresolved?"UNKNOWN":"VALIDATED");
    auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout << "{\"kind\":\"summary\",\"q\":"<<q<<",\"expected_slopes\":"<<expected_slopes
      <<",\"expected_cases\":"<<expected_cases<<",\"completed_slopes\":"<<completed_slopes
      <<",\"evaluated_cases\":"<<evaluated<<",\"validated_cases\":"<<validated
      <<",\"counterexample_cases\":"<<counterexamples<<",\"unknown_evaluated_cases\":"<<unknown
      <<",\"unresolved_remaining\":"<<unresolved<<",\"missing_cases\":[";
    bool first=true;
    for(auto key:expected) if(std::find(completed.begin(),completed.end(),key)==completed.end()) {
        if(!first) std::cout<<','; first=false;
        std::cout << "{\"slope\":["<<key[0]<<','<<key[1]<<",1],\"relation\":["<<(key[2]&1)<<','<<((key[2]>>1)&1)<<"]}";
    }
    std::cout << "],\"verdict\":\""<<verdict<<"\",\"elapsed_seconds\":"<<elapsed<<",\"failure\":";
    if(failure.empty()) std::cout<<"null"; else std::cout<<'"'<<failure<<'"';
    std::cout << "}\n";
    return verdict=="VALIDATED"?0:1;
}
