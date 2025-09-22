namespace Number_Theory {
    using std::gcd;
    using std::lcm;
    using std::min;
    using std::max;
	using i64 = long long;
    using d64 = long double;
    using i128 = __int128;
    using pii = std::pair<int, int>;
    using pII = std::pair<std::vector<int>, std::vector<int>>;
    using pll = std::pair<i64, i64>;
    using pLL = std::pair<std::vector<i64>, std::vector<i64>>;
	constexpr int maxn = 1e6 + 10;
	constexpr int mod = 998244353;
	constexpr int MOD = 1e9 + 7;
	constexpr int N = 2e5 + 10000;
	constexpr double eps = 1e-9;
    constexpr d64 EPS = 1e-18;
    std::vector<int> phi, prime, ind;

	i64 qpow(i64 a, i64 b, i64 p) {
		i64 ans = 1ll % p;
		for(; b != 0; b >>= 1) {
			if(b & 1) ans = (i128)ans * a % p;
			a = (i128)a * a % p;
		} 
        return ans;
	} 
	i64 qmul(i64 x, i64 y, i64 p) {
	    i64 ans = 0; x %= p, y %= p;
	    for(; y != 0; y >>= 1) {
	        if(y & 1) ans = (ans + x) % p;
	        x = (x + x) % p;
	    }
	    return ans;
	}
	i64 qmod(i64 x, i64 p) {
		if(x >= p) return x - p;
		else return x;
	}
    i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
        if(!b) return x = 1, y = 0, a;
        i64 g = exgcd(b, a % b, y, x);
        y = y - (a / b) * x; return g;
    }
    i64 congruence(i64 a, i64 b, i64 m) { // ax\equivb(mod m)
	    i64 x, y, d = exgcd(a, m, x, y);
	    if(b % d) return -1;
		else m = m / d;
	    return ((i128)x * (b / d) % m + m) % m;
	}
	i64 invmod(i64 a, i64 p) {
		return congruence(a, 1ll, p);
	}
    i64 CRT(std::vector<i64> &r, std::vector<i64> &a) {
	    i128 n = 1, ans = 0;
	    int k = a.size();
	    for(auto x : r) n *= x;
	    for(int i = 0; i < k; i++) {
	        i128 m = n / r[i]; i64 R, y;
	        exgcd(m % r[i], r[i], R, y);
	        ans = (ans + (i128)a[i] * m % n * R % n) % n;
	    }
	    i64 res = (ans % n + n) % n;
	    return res;
	}
    i64 exCRT(std::vector<i64> &r, std::vector<i64> &a) {
	    i128 A = a[0]; i64 R = r[0];
	    for(int i = 1; i < (int)r.size(); i++) {
	        i64 p = a[i], q = r[i];
	        i64 x = congruence(R, (p - A) % q, q);
	        if(x == -1) return -1;
	        A += (i128)R * x;
	        R = std::lcm(R, q);
	        A %= R; if(A < 0) A += R;
	    }
	    return (A % R + R) % R;
	}
	i64 getphi(int n) { // phi(n)
		i64 ans = 1;
		for(int i = 2; i * i <= n; ++i) {
			if(n % i != 0) continue;
			ans *= i - 1, n /= i;
			while(!(n % i)) ans *= i, n /= i;
		} 
		if(n > 1) return ans * (n - 1);
		else return ans;
	}
	i64 exEuler(i64 a, std::string b, i64 m) { 
        // a^b mod m = a^(b mod phim + phim) mod m (b >= phim)
		i64 f = 0, phim = getphi(m), B = 0;  // B = b mod m
		for(auto &c : b) {
			B = 10ll * B + c - '0';
			if(B >= phim) f = 1, B %= phim; // f == 1 -> b >= phim
		}
		B += f ? phim : 0;
		return qpow(a, B, m);
	}
	void EulerSieve(int n) {
        phi.clear(), prime.clear();
		std::vector<int> vis(n + 2, 0);
        phi.resize(n + 2), prime.resize(n + 1);
		int cnt = 0; phi[1] = 1;
		for(int i = 2; i <= n + 1; i++) {
			if(!vis[i]) prime[++cnt] = i, phi[i] = i - 1;
			for(int j = 1; i * prime[j] <= n + 1 && j <= cnt; j++) {
				vis[i * prime[j]] = true;
				if(i % prime[j] == 0) {
					phi[i * prime[j]] = phi[i] * prime[j];
					break;
				} else phi[i * prime[j]] = phi[i] * phi[prime[j]];
			}
		}
	}
	int MPR(int n, int p, int l = 1, int r = 50) { // minimum primitive root
		int q = p; // EulerSieve first, p = phi[n]
    	std::vector<int> pri;
    	for(int i = 2; i * i <= q; i++) {
    		if(q % i == 0) pri.push_back(i);
    		while(q % i == 0) q /= i;
		} if(q > 1) pri.push_back(q);
		auto check = [&](int i) -> bool {
			if(MATH::qpow(i, p, n) != 1) return 0;
			for(auto j : pri) 
				if(MATH::qpow(i, p / j, n) == 1) 
					return 0;
			return 1;
		}; 
		for(int i = l; i <= r; i++) // n^0.25
			if(check(i)) return i;
		return 0;
	}
	i64 BSGS(i64 a, i64 b, i64 p) { // a^x \equiv b (mod p) x = min
		a %= p, b %= p;
		if(b == 1) return 0ll;
		i64 m = ceil(sqrt(p));
		std::unordered_map<i64, i64> hash;
		for(i64 j = 0; j < m; j++) {
			hash[b] = j; // non-decreasing
			b = (i128)b * a % p;
		} // x = m * i - j, (a^m)^i \equiv ba^j (mod p)
		i64 M = qpow(a, m, p), A = 1ll; // M = a^m, A = (a^m)^i
		for(i64 i = 1; i <= m; i++) {
			A = (i128)A * M % p;
			if(hash.count(A)) {
				i64 ans = i * m - hash[A];
				if(ans >= 0) return ans;
			}
		}
		return -1;
	}
	i64 exBSGS(i64 a, i64 b, i64 p) { // gcd(a,p)!=1
		a %= p, b %= p;
		if(b == 1 || p == 1) return 0;
		i64 k = 0, d, A = 1ll;
		while((d = gcd(a, p)) != 1) {
			if(b % d) return -1;
			++k, b /= d, p /= d;
			A = (i128)A * (a / d) % p; // (a^k)/D
			if(A == b) return k;
		}
		b = (i128)b * invmod(A, p) % p;
		i64 ans = BSGS(a, b, p);
		if(ans == -1) return -1;
		else return ans + k;
	}
	void IndSieve(i64 g, i64 p) { // g has to be p's primitive root
		i64 B = sqrt(p * sqrt(p) / log(p)); // Block size
		std::unordered_map<i64, i64> hash;
		i64 w = invmod(qpow(g, B, p), p);
		i64 n = sqrt(p) + 1, Phi = p - 1; // if p is prime
		prime.clear(); prime.resize(n + 2);
		ind.clear(); ind.resize(n + 2);
		std::vector<i64> vis(n + 2, 0);
		for(i64 i = 0ll, j = 1ll; i < B; i++) {
			hash[j] = i, j = (i128)j * g % p;
		}
		auto bsgs = [&](i64 a) -> i64 { // BSGS
			for(int i = 0; i <= p / B; i++) {
				if(hash.count(a)) return i * B + hash[a];
				a = (i128)a * w % p;
			}
			return -1;
		};
		ind[1] = 0, ind[0] = bsgs(p - 1);
		for(int i = 2, cnt = 0; i <= n; i++) {
			if(!vis[i]) prime[++cnt] = i, ind[i] = bsgs(i);
			for(int j = 1; j <= cnt && i * prime[j] <= n; j++) {
				i64 k = prime[j]; vis[i * k] = 1;
				ind[i * k] = qmod(ind[i] + ind[k], Phi);
				if(i % k == 0) break;
			}
		}
	}
	i64 Ind(i64 a, i64 p) { // Discrete Log
		i64 lim = sqrt(p) + 1;
		if(a <= lim) return ind[a];
		i64 b = p / a, c = p % a, Phi = p - 1;
		if(c < a - c) return qmod(qmod(ind[0] + Ind(c, p), Phi) - ind[b] + Phi, Phi);
		else return qmod(Ind(a - c, p) - ind[b + 1] + Phi, Phi);
	}
	i64 conv(i64 a, i64 h, i64 p) { // convert to base h:ind_h(a)
		// M = getphi(p), t = ind_g(a), k = ind_g(h)
        i64 M = p - 1; // prime || p^(k-1)(p-1)
        i64 t = Ind(a, p);
        i64 k = Ind(h, p);
	    return congruence(k, t, M);
	}
    i64 DSM(i64 n, i64 k, int flag = 0) { 
        // Divisor Summatory Method
        // sum(i*floor(n/i),i=1~n,k=n)=sum(f(i),i=1~n), f(i)=sum(divisors of i)
		// f(x)+f(x+1)+...f(y)=DSM(y,y)-DSM(x-1,x-1)=(1~y)-(1~(x-1))
		i64 l, r, t, ans = 0ll;
		if(flag == 0) for(l = 1; l <= n; l = r + 1) {
			t = k / l; if(!t) break;
			r = min(n, k / t);
			ans += t * (l + r) * (r - l + 1) / 2ll;
		}
		if(flag == 1) for(l = 1; l <= n; l = r + 1) { // sum(floor(k/i),i=1~n)
			t = k / l; if(!t) break;
			r = min(n, k / t);
			ans += t * (r - l + 1);
		}
        /*
        // sum(C*floor(n1/i)*floor(n2/i)...*floor(nk/i),i=1~min({nj},j=1~k))
        for(l = 1; l <= min{n_1,n_2...n_k}; l = r + 1) {
        	i = 1~k, t_i = n_i / l; if(!t_i) break;
        	i = 1~k, r = min({n_i / t_i});
        	ans += prod(t_i,i=1~k) * (sum of coefficients)
        } 
        */
		return ans;
	}
    double fixn(double tar, int n = 2) {
        d64 p = std::pow(10.0L, n);
        d64 v = (d64)tar * p + ((tar > 0) ? 0.5L : -0.5L);
        i64 iv = (i64)(v);
        d64 res = (d64)iv / p;
        return (double)res;
    }
}