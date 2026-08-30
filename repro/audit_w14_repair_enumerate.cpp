#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

struct Sig {
    uint16_t r[14];
    uint8_t d;
    bool operator==(Sig const& o) const {
        if(d!=o.d) return false;
        for(int i=0;i<d;i++) if(r[i]!=o.r[i]) return false;
        return true;
    }
};

static inline vector<uint16_t> kernel_basis(vector<uint16_t> rows, int n){
    int R=0; vector<int> piv;
    for(int c=0;c<n;c++){
        int p=-1; for(int i=R;i<(int)rows.size();i++) if((rows[i]>>c)&1){p=i;break;}
        if(p<0) continue;
        swap(rows[R],rows[p]);
        for(int i=0;i<(int)rows.size();i++) if(i!=R && ((rows[i]>>c)&1)) rows[i]^=rows[R];
        piv.push_back(c); R++; if(R==(int)rows.size()) break;
    }
    vector<char> isp(n,0); for(int c:piv) isp[c]=1;
    vector<uint16_t> out;
    for(int f=0;f<n;f++) if(!isp[f]){
        uint16_t x=(uint16_t)(1u<<f);
        for(int i=0;i<R;i++) if((rows[i]>>f)&1) x ^= (uint16_t)(1u<<piv[i]);
        out.push_back(x);
    }
    return out;
}

static inline int popc(uint16_t x){ return __builtin_popcount((unsigned)x); }

static inline Sig signature(int n, const vector<uint16_t>& nbr, const vector<pair<int,int>>& edges,
                            const uint8_t* b, int extra_u=-1, int extra_v=-1){
    vector<uint16_t> rows; rows.reserve(n+2);
    uint16_t Zmask=0, Ymask=0;
    for(int v=0;v<n;v++){
        if(b[v]==0) rows.push_back(nbr[v]);
        else if(b[v]==1){ rows.push_back((uint16_t)(nbr[v] ^ (1u<<v))); Ymask |= (1u<<v); }
        else { rows.push_back((uint16_t)(1u<<v)); Zmask |= (1u<<v); }
    }
    if(extra_u>=0) rows.push_back((uint16_t)(1u<<extra_u));
    if(extra_v>=0) rows.push_back((uint16_t)(1u<<extra_v));
    auto bas=kernel_basis(rows,n);
    Sig s{}; s.d=(uint8_t)bas.size();
    uint16_t nonZ=(uint16_t)(((1u<<n)-1u) ^ Zmask);
    for(int i=0;i<(int)bas.size();i++){
        uint16_t a=bas[i], z=0;
        uint16_t aa=a;
        while(aa){ uint16_t lb=aa & (uint16_t)(-aa); int v=__builtin_ctz((unsigned)lb); z ^= nbr[v]; aa ^= lb; }
        uint16_t q=(uint16_t)((a & nonZ) | (z & Zmask));
        int ep=0; for(auto [u,v]:edges) ep ^= (((a>>u)&1) & ((a>>v)&1));
        int yc=popc((uint16_t)(a & Ymask));
        if(yc&1){ cerr<<"odd ycount\n"; abort(); }
        int sign=ep ^ ((yc/2)&1);
        s.r[i]=(uint16_t)(q | (sign<<n));
    }
    int R=0;
    for(int c=0;c<n;c++){
        int p=-1; for(int i=R;i<s.d;i++) if((s.r[i]>>c)&1){p=i;break;}
        if(p<0) continue;
        swap(s.r[R],s.r[p]);
        for(int i=0;i<s.d;i++) if(i!=R && ((s.r[i]>>c)&1)) s.r[i]^=s.r[R];
        R++;
    }
    if(R!=s.d){ cerr<<"dependent q images\n"; abort(); }
    return s;
}

struct Mask {
    static const int W=22;
    uint64_t a[W];
    bool operator==(Mask const&o) const { return memcmp(a,o.a,sizeof(a))==0; }
};
struct MH { size_t operator()(Mask const&m) const noexcept { uint64_t h=0x9e3779b97f4a7c15ULL; for(int i=0;i<Mask::W;i++){ uint64_t x=m.a[i]; h^=x+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);} return (size_t)h; } };

static vector<pair<int,int>> wheel_edges(int n){
    vector<pair<int,int>> E; for(int i=1;i<n;i++) E.push_back({0,i}); int m=n-1;
    for(int i=1;i<n;i++){ int u=i, j=(i % m)+1; if(u>j) swap(u,j); E.push_back({u,j}); }
    sort(E.begin(),E.end()); E.erase(unique(E.begin(),E.end()),E.end()); return E;
}
static vector<uint16_t> nbrs(int n, const vector<pair<int,int>>& E){ vector<uint16_t>N(n); for(auto [u,v]:E){N[u]|=1u<<v;N[v]|=1u<<u;} return N; }

static inline void trits(uint64_t x,int n,uint8_t*b){ for(int i=n-1;i>=0;i--){b[i]=x%3;x/=3;} }

static Mask collision_mask(int n, const vector<pair<int,int>>& E, const vector<uint16_t>& Nt,
                           const vector<vector<pair<int,int>>>& Ec, const vector<vector<uint16_t>>& Nc,
                           const uint8_t*b){
    int m=E.size(), H=1+2*m; vector<Sig> sig(H);
    sig[0]=signature(n,Nt,E,b);
    for(int j=0;j<m;j++){
        sig[1+2*j]=signature(n,Nc[j],Ec[j],b);
        sig[2+2*j]=signature(n,Nt,E,b,E[j].first,E[j].second);
    }
    Mask M{}; int k=0;
    for(int i=0;i<H;i++) for(int j=i+1;j<H;j++,k++) if(sig[i]==sig[j]) M.a[k>>6] |= 1ULL<<(k&63);
    return M;
}

static bool disjoint(const Mask&a,const Mask&b,int words){ for(int i=0;i<words;i++) if(a.a[i]&b.a[i]) return false; return true; }

int main(int argc,char**argv){
    int n= argc>1?atoi(argv[1]):10; int nth=argc>2?atoi(argv[2]):omp_get_max_threads(); omp_set_num_threads(nth);
    auto E=wheel_edges(n); auto Nt=nbrs(n,E); int m=E.size(), H=1+2*m, np=H*(H-1)/2, words=(np+63)/64;
    vector<vector<pair<int,int>>> Ec(m); vector<vector<uint16_t>> Nc(m);
    for(int j=0;j<m;j++){ Ec[j]=E; Ec[j].erase(Ec[j].begin()+j); Nc[j]=nbrs(n,Ec[j]); }
    uint64_t total=1; for(int i=0;i<n;i++) total*=3;
    cerr<<"W"<<n<<" contexts="<<total<<" edges="<<m<<" H="<<H<<" pairs="<<np<<" threads="<<nth<<"\n";
    double t0=omp_get_wtime();
    vector<unordered_set<Mask,MH>> local(nth);
    #pragma omp parallel
    {
        int tid=omp_get_thread_num(); auto &S=local[tid]; S.reserve((size_t)(total/nth/4+1000)); uint8_t b[16];
        #pragma omp for schedule(static)
        for(uint64_t x=0;x<total;x++){
            trits(x,n,b); Mask M=collision_mask(n,E,Nt,Ec,Nc,b); S.insert(M);
        }
    }
    cerr<<"enum secs="<<(omp_get_wtime()-t0)<<" local sizes"; size_t sum=0; for(auto&s:local){cerr<<" "<<s.size();sum+=s.size();} cerr<<" sum="<<sum<<"\n";
    vector<Mask> V; V.reserve(sum);
    for(auto &s:local){ for(auto const&m:s) V.push_back(m); s.clear(); s.rehash(0); }
    local.clear(); local.shrink_to_fit();
    auto lessmask=[](Mask const&a,Mask const&b){ for(int i=0;i<Mask::W;i++){ if(a.a[i]<b.a[i]) return true; if(a.a[i]>b.a[i]) return false; } return false; };
    sort(V.begin(),V.end(),lessmask);
    V.erase(unique(V.begin(),V.end()),V.end());
    cerr<<"distinct collision masks="<<V.size()<<" merge secs="<<(omp_get_wtime()-t0)<<"\n";
    { map<int,size_t> hist; for(auto const &mm:V){int bc=0;for(int w=0;w<words;w++)bc+=__builtin_popcountll(mm.a[w]);hist[bc]++;} cerr<<"bitcount hist head:";int z=0;for(auto [bc,c]:hist){cerr<<" "<<bc<<":"<<c;if(++z>=15)break;}cerr<<"\n"; }
    if(n==14){
        ofstream fo("w14_masks.bin", ios::binary); uint64_t NN=V.size(); uint32_t ww=words, npp=np; fo.write((char*)&NN,8); fo.write((char*)&ww,4); fo.write((char*)&npp,4); for(auto const&m:V)fo.write((char*)m.a,words*8); fo.close(); cerr<<"dumped w14_masks.bin bytes="<<(16+NN*words*8)<<"\n"; return 0;
    }
    vector<int> freq(np,0);
    size_t total_setbits=0;
    for(auto const &mm:V) for(int w=0;w<words;w++){ uint64_t x=mm.a[w]; while(x){int z=__builtin_ctzll(x);freq[w*64+z]++;total_setbits++;x&=x-1;} }
    int K = argc>3?atoi(argv[3]):8;
    vector<int> rankbits(np); iota(rankbits.begin(),rankbits.end(),0);
    sort(rankbits.begin(),rankbits.end(),[&](int a,int b){ if(freq[a]!=freq[b]) return freq[a]>freq[b]; return a<b; });
    struct Key { uint64_t lo=0,hi=0; bool operator==(Key const&o)const{return lo==o.lo&&hi==o.hi;} };
    struct KeyHash { size_t operator()(Key const&k)const noexcept{return k.lo^(k.hi+0x9e3779b97f4a7c15ULL+(k.lo<<6)+(k.lo>>2));} };
    struct Group { array<uint16_t,8> key{}; vector<int> ids; };
    unordered_map<Key,int,KeyHash> gid; gid.reserve(V.size()/8+100);
    vector<Group> groups;
    auto key_for=[&](Mask const &mm){
        array<uint16_t,8> a{}; int got=0;
        for(int p:rankbits){ if((mm.a[p>>6]>>(p&63))&1ULL){ a[got++]=(uint16_t)p; if(got==K)break; } }
        if(got<K){ cerr<<"mask has fewer than K bits\n"; abort(); }
        sort(a.begin(),a.begin()+K); return a;
    };
    auto packkey=[&](array<uint16_t,8> const&a){ Key k{}; int sh=0; for(int i=0;i<K;i++){ uint64_t val=a[i]; if(sh<64){ k.lo |= val<<sh; if(sh>53) k.hi |= val>>(64-sh); } else k.hi |= val<<(sh-64); sh+=11; } return k; };
    for(int i=0;i<(int)V.size();i++){
        auto a=key_for(V[i]); Key k=packkey(a); auto it=gid.find(k); int g;
        if(it==gid.end()){g=groups.size();gid.emplace(k,g);Group G;G.key=a;groups.push_back(std::move(G));} else g=it->second;
        groups[g].ids.push_back(i);
    }
    cerr<<"grouped K="<<K<<" groups="<<groups.size()<<" avg="<<(double)V.size()/groups.size()<<" total_setbits="<<total_setbits<<"\n";
    size_t N=V.size(), BW=(N+63)/64;
    vector<uint64_t> has((size_t)np*BW,0);
    for(size_t i=0;i<N;i++) for(int w=0;w<words;w++){uint64_t x=V[i].a[w];while(x){int z=__builtin_ctzll(x);int p=w*64+z;has[(size_t)p*BW+(i>>6)]|=1ULL<<(i&63);x&=x-1;}}
    cerr<<"index built MB="<<(has.size()*8.0/1048576.0)<<" secs="<<(omp_get_wtime()-t0)<<"\n";
    atomic<long long> found(-1); atomic<unsigned long long> cand_total(0), checks(0);
    #pragma omp parallel
    {
        vector<uint64_t> cand(BW);
        vector<int> cids; cids.reserve(100000);
        #pragma omp for schedule(dynamic,1)
        for(long long gg=0;gg<(long long)groups.size();gg++){
            if(found.load(memory_order_relaxed)>=0) continue;
            auto const&G=groups[gg];
            int p0=G.key[0]; for(size_t q=0;q<BW;q++) cand[q]=~has[(size_t)p0*BW+q];
            if(N%64)cand.back()&=(1ULL<<(N%64))-1;
            for(int t=1;t<K;t++){int p=G.key[t];for(size_t q=0;q<BW;q++)cand[q]&=~has[(size_t)p*BW+q];}
            cids.clear(); for(size_t q=0;q<BW;q++){uint64_t x=cand[q];while(x){int z=__builtin_ctzll(x);cids.push_back((int)(q*64+z));x&=x-1;}}
            cand_total.fetch_add(cids.size(),memory_order_relaxed);
            for(int ai:G.ids){
                auto const&A=V[ai];
                for(int bj:cids){checks.fetch_add(1,memory_order_relaxed); if(disjoint(A,V[bj],words)){found.store(((long long)ai<<32)|(unsigned)bj,memory_order_relaxed);break;}}
                if(found.load(memory_order_relaxed)>=0)break;
            }
        }
    }
    if(found.load()>=0){long long q=found.load();cout<<"FOUND_DISJOINT "<<(q>>32)<<" "<<(q&0xffffffffLL)<<" groups="<<groups.size()<<" checks="<<checks.load()<<"\n";return 2;}
    cout<<"NO_TWO_COVER distinct="<<V.size()<<" groups="<<groups.size()<<" candidate_sum="<<cand_total.load()<<" checks="<<checks.load()<<" sec="<<(omp_get_wtime()-t0)<<"\n";
    return 0;
}
