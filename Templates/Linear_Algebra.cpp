namespace Linear_Algebra {
	using i64 = long long;
    using d64 = long double;
    using i128 = __int128;
    using pii = std::pair<int, int>;
    using pII = std::pair<std::vector<int>, std::vector<int>>;
    using pll = std::pair<i64, i64>;
    using pLL = std::pair<std::vector<i64>, std::vector<i64>>;
	constexpr int maxn = 1e6 + 10;
	constexpr int mod = 998244353;
	constexpr int N = 2e5 + 10000;
	constexpr double eps = 1e-9;
    constexpr d64 EPS = 1e-18;

    struct matrix { // 1-based index(start from 1 not 0)
        static const int MOD = 1e9 + 7;
	    std::vector<std::vector<int>> G; int N;
	    matrix(int n = 0) {N = n; G.assign(n + 1, std::vector<int>(n + 1, 0));}
	    matrix(const matrix &other) {
		    N = other.N;
		    G = other.G;
		}
	    inline void clear() {for(auto &row : G) fill(row.begin(), row.end(), 0);}
	    inline int size() const {return (int)G.size() - 1;}
	    friend matrix operator * (const matrix &a, const matrix &b) {
	        int n = a.N; matrix c(n);
	        for(int i = 1; i <= n; i++) for(int k = 1; k <= n; k++) 
                if(a.G[i][k]) for(int j = 1; j <= n; j++) 
					c.G[i][j] = (c.G[i][j] + 1ll * a.G[i][k] * b.G[k][j] % MOD) % MOD;
	        return c;
	    }
	    friend bool operator == (const matrix &a, const matrix &b) {
	    	int n = a.N, m = b.N;
		    if(n != m) return false;
		    for(int i = 1; i <= n; i++) for(int j = 1; j <= n; j++) {
		        if(a.G[i][j] != b.G[i][j]) return false;
		    }
		    return true;
		}
		friend bool operator != (const matrix &a, const matrix &b) {
		    return !(a == b);
		}
		friend matrix operator + (const matrix &a, const matrix &b) {
	        int n = a.N; matrix c(n);
	        for(int i = 1; i <= n; i++) for(int j = 1; j <= n; j++)
	            c.G[i][j] = (a.G[i][j] + b.G[i][j]) % MOD;
	        return c;
	    }
	    friend matrix operator - (const matrix &a, const matrix &b) {
	        int n = a.N; matrix c(n);
	        for(int i = 1; i <= n; i++) for(int j = 1; j <= n; j++) {
                c.G[i][j] = (a.G[i][j] - b.G[i][j]) % MOD;
                if(c.G[i][j] < 0) c.G[i][j] += MOD;
            }
	        return c;
	    }
	    matrix& operator = (const matrix &other) {
	        if(this != &other) G = other.G;
	        return *this;
	    }
	    friend std::istream & operator >> (std::istream &in, matrix &m) {
	        int n = m.N;
	        for(int i = 1; i <= n; i++) {
	            for(int j = 1; j <= n; j++) {
	                in >> m.G[i][j];
	                m.G[i][j] %= MOD;
	            }
	        }
	        return in;
	    }
	    friend std::ostream & operator << (std::ostream &out, const matrix &m) {
	        int n = m.N;
	        for(int i = 1; i <= n; i++) for(int j = 1; j <= n; j++) {
	            out << m.G[i][j] << (j == n ? '\n' : ' ');
	        }
	        return out;
	    }
		inline matrix mpow(i64 p) {
	        int n = N; matrix res(n);
	        for(int i = 1; i <= n; i++) res.G[i][i] = 1;
	        matrix base = *this;
	        while(p != 0) {
	            if(p & 1) res = res * base;
				base = base * base;
	            p >>= 1;
	        }
	        return *this = res;
	    }
	    matrix inv() const {
	        int n = N; matrix A = *this, I(n), U(0);
	        for(int i = 1; i <= n; i++) I.G[i][i] = 1;
	        for(int i = 1; i <= n; i++) {
	            int pivot = i;
	            for(int j = i; j <= n; j++) if(A.G[j][i] != 0) {
					pivot = j; break;
				}
	            if(A.G[pivot][i] == 0) {
	            	std::cout << "No Solution\n";
	            	return U;
				}
	            if(pivot != i) {
	                std::swap(A.G[i], A.G[pivot]);
	                std::swap(I.G[i], I.G[pivot]);
	            }
	            i64 inv_pivot = qpow(A.G[i][i], MOD - 2, MOD);
	            for(int j = 1; j <= n; j++) {
	                A.G[i][j] = 1LL * A.G[i][j] * inv_pivot % MOD;
	                I.G[i][j] = 1LL * I.G[i][j] * inv_pivot % MOD;
	            }
	            for(int r = 1; r <= n; r++) {
	                if(r == i || A.G[r][i] == 0) continue;
	                i64 factor = A.G[r][i];
	                for(int c = 1; c <= n; c++) {
	                    A.G[r][c] = (A.G[r][c] - 1LL * factor * A.G[i][c]) % MOD;
	                    if(A.G[r][c] < 0) A.G[r][c] += MOD; 
	                    I.G[r][c] = (I.G[r][c] - 1LL * factor * I.G[i][c]) % MOD;
	                    if(I.G[r][c] < 0) I.G[r][c] += MOD; 
	                }
	            }
	        }
	        return I;
	    }
	};
	std::vector<double> gauss(matrix A, std::vector<int> b) {
	    int n = A.N;
	    std::vector<std::vector<double>> M(n, std::vector<double>(n));
	    std::vector<double> B(n);
	    for(int i = 0; i < n; i++) B[i] = (double)b[i + 1];
	    for(int i = 0; i < n; i++)
	        for(int j = 0; j < n; j++)
	            M[i][j] = (double)A.G[i + 1][j + 1];
	    int rank = 0;
	    for(int i = 0; i < n; i++) {
	        int pivot = -1;
	        for(int j = rank; j < n; j++) {
	            if(fabs(M[j][i]) > EPS) {
	                pivot = j;
	                break;
	            }
	        }
	        if(pivot == -1) continue;
	        std::swap(M[rank], M[pivot]);
	        std::swap(B[rank], B[pivot]);
	        double div = M[rank][i];
	        for(int k = i; k < n; k++) M[rank][k] /= div;
	        B[rank] /= div;
	        for(int j = 0; j < n; j++) {
	            if(j == rank) continue;
	            if(fabs(M[j][i]) > EPS) {
	                double factor = M[j][i];
	                for(int k = i; k < n; k++) {
	                    M[j][k] -= factor * M[rank][k];
	                }
	                B[j] -= factor * B[rank];
	            }
	        }
	        rank++;
	    }
	    for(int i = rank; i < n; i++)
	        if(fabs(B[i]) > EPS) return {-1}; // no answer
	    if(rank < n) return {0}; // inf answers
	    std::vector<double> x(n + 1, 1);
	    for(int i = 1; i <= n; i++) x[i] = B[i - 1];
	    return x;
	}
	i64 mxor(std::vector<i64> a) {
		int n = a.size(), k = 0;
		auto Gauss = [&]() -> void {
			for(int i = 63; i >= 0; i--) {
				for(int j = k; j < n; j++) {
					if((a[j] >> i) & 1) {
						std::swap(a[j], a[k]);
						break;
					}
				}
				if(((a[k] >> i) & 1) == 0) continue;
				for(int j = 0; j < n; j++)
					if(j != k && ((a[j] >> i) & 1)) 
						a[j] ^= a[k];
				++k; if(k == n) break;
			}
		}; Gauss(); i64 ans = 0ll;
		for(int i = 0; i < k; i++) ans ^= a[i];
		return ans;
 	}
 	struct LinearBasis {
	    static const int bit = 63;
	    i64 base[bit]; std::vector<i64> v;
	    int zero = 0, rank = 0;
	    LinearBasis() {memset(base, 0, sizeof(base)); v.clear();}
	    bool ask(i64 x) {
	    	for(int i = bit - 1; i >= 0; i--)
	    		if(x & (1ll << i)) x ^= base[i];
	    	return x == 0;
		}
	    void insert(i64 x) {
	        if(x == 0) {return zero = 1, void();}
	        for(int i = bit - 1; i >= 0; i--) {
	            if(!(x >> i & 1)) continue;
	            if(!base[i]) {
	                base[i] = x;
	                rank++;
	                return;
	            }
	            x ^= base[i];
	        }
	        if(x == 0) zero = 1;
	    }
        void merge(LinearBasis &b) {
            for(int i = bit - 1; i >= 0; i--) 
                if(b.base[i]) insert(b.base[i]);
        }
	    i64 max_xor() {
	        i64 res = 0;
	        for(int i = bit - 1; i >= 0; i--)
	            if((res ^ base[i]) > res) res ^= base[i];
	        return res;
	    }
	    i64 min_xor() {
	        if(zero) return 0;
	        for(int i = 0; i < bit; i++)
	            if(base[i]) return base[i];
	        return 0;
	    }
	    void rebuild() {
	        i64 tmp[bit]; v.clear();
			memcpy(tmp, base, sizeof(base));
	        for(int i = bit - 1; i >= 0; i--)
	            for(int j = i - 1; j >= 0; j--)
	                if(tmp[i] & (1ll << j)) tmp[i] ^= tmp[j];
	        for(int i = 0; i < bit; i++) if(tmp[i]) v.push_back(tmp[i]);
	    }
	    i64 kth_xor(i64 k) {
            k -= zero; if(k == 0) return 0ll;
	        int r = v.size(); // rebuild first
	        i64 total = 1LL << r, res = 0ll;
	        if(k >= total) return -1;
	        for(int i = bit - 1; i >= 0; i--)
	            if(k & (1ll << i)) res ^= v[i];
	        return res;
	    }
	};
}