#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
struct Mask{uint64_t a[22];};
struct Key{array<uint16_t,16>a{}; bool operator==(Key const&o)const{return a==o.a;}};
struct KH{size_t operator()(Key const&x)const noexcept{uint64_t h=0;for(int i=0;i<16;i++)h^=(uint64_t)x.a[i]+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
static inline bool disjoint(const Mask&a,const Mask&b,int words){for(int w=0;w<words;w++)if(a.a[w]&b.a[w])return false;return true;}
int main(){
 const int K=16; ifstream f("w14_masks.bin",ios::binary);uint64_t N;uint32_t words,np;f.read((char*)&N,8);f.read((char*)&words,4);f.read((char*)&np,4);vector<Mask>V(N);for(auto&m:V){memset(m.a,0,sizeof m.a);f.read((char*)m.a,words*8);}cerr<<"loaded N="<<N<<"\n";
 vector<int>freq(np);for(auto&m:V)for(int w=0;w<(int)words;w++){uint64_t x=m.a[w];while(x){int z=__builtin_ctzll(x);freq[w*64+z]++;x&=x-1;}}
 vector<int>rank(np);iota(rank.begin(),rank.end(),0);sort(rank.begin(),rank.end(),[&](int a,int b){if(freq[a]!=freq[b])return freq[a]>freq[b];return a<b;});
 struct Group{Key key;vector<int>ids;};unordered_map<Key,int,KH>gm;gm.reserve(N/8);vector<Group>G;G.reserve(40000);
 for(int i=0;i<(int)N;i++){Key key;int got=0;for(int p:rank)if((V[i].a[p>>6]>>(p&63))&1ULL){key.a[got++]=p;if(got==K)break;}sort(key.a.begin(),key.a.end());auto it=gm.find(key);int g;if(it==gm.end()){g=G.size();gm.emplace(key,g);G.push_back({key,{}});}else g=it->second;G[g].ids.push_back(i);}cerr<<"groups="<<G.size()<<"\n";
 size_t BW=(N+63)/64;vector<uint64_t>has((size_t)np*BW);for(size_t i=0;i<N;i++)for(int w=0;w<(int)words;w++){uint64_t x=V[i].a[w];while(x){int z=__builtin_ctzll(x),p=w*64+z;has[(size_t)p*BW+(i>>6)]|=1ULL<<(i&63);x&=x-1;}}cerr<<"index ready MB="<<has.size()*8.0/1048576<<"\n";
 atomic<long long> found(-1);atomic<unsigned long long>checks(0),candtot(0);double t0=omp_get_wtime();
 #pragma omp parallel
 { vector<uint64_t> cand(BW);vector<int>ids;ids.reserve(512);
   #pragma omp for schedule(dynamic,1)
   for(long long gg=0;gg<(long long)G.size();gg++){
    if(found.load(memory_order_relaxed)>=0)continue;auto const&gr=G[gg];int p0=gr.key.a[0];for(size_t q=0;q<BW;q++)cand[q]=~has[(size_t)p0*BW+q];if(N%64)cand.back()&=(1ULL<<(N%64))-1;for(int t=1;t<K;t++){int p=gr.key.a[t];for(size_t q=0;q<BW;q++)cand[q]&=~has[(size_t)p*BW+q];}
    ids.clear();for(size_t q=0;q<BW;q++){uint64_t x=cand[q];while(x){int z=__builtin_ctzll(x);ids.push_back(q*64+z);x&=x-1;}}candtot.fetch_add(ids.size(),memory_order_relaxed);
    for(int ai:gr.ids){for(int bj:ids){checks.fetch_add(1,memory_order_relaxed);if(disjoint(V[ai],V[bj],words)){found.store(((long long)ai<<32)|(unsigned)bj,memory_order_relaxed);break;}}if(found.load(memory_order_relaxed)>=0)break;}
   }
 }
 cerr<<"search sec="<<omp_get_wtime()-t0<<" checks="<<checks.load()<<" cand_sum="<<candtot.load()<<"\n";
 if(found.load()>=0){auto q=found.load();cout<<"FOUND_DISJOINT "<<(q>>32)<<" "<<(q&0xffffffffLL)<<"\n";return 2;}cout<<"NO_TWO_COVER W14 distinct_collision_masks="<<N<<" groups="<<G.size()<<" checks="<<checks.load()<<"\n";return 0;
}
