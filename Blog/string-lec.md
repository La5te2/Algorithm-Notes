### 1. KMP

##### 基本定义

先从最简单的***KMP***开始。（~~一点也不简单~~）

> 模式串匹配：给定一个文本串S和需要在文本串中搜索的模式串T，$|S|>>|T|$ ，查询在S中T出现的次数，位置等信息。

暴力匹配的时间复杂度大约为 $O(|S||T|)$，而 ***KMP*** 时间复杂度大约为 $O(|S|+|T|)$。

> 朴素的单模式串匹配枚举每一个文本串元素，然后从这一位开始不断向后比较，每次比较失败之后都要从头开始重新比对，而KMP每次失配之后，不会从头重新开始枚举，而是根据已经得知的数据，从某个特定的位置开始匹配。

当某一位失配时，将前一位跳跃到之前匹配过的某一位——使用失配数组来完成这一跳跃。

例子：（以下字符串首位**下标从1开始**）

```
T:ABABC
S:ABABABCAA
//MOVE 2 POSITIONS
T:  ABABC
S:ABABABCAA
```

预处理时应当考虑当**模式串T**的第 $i$ 位失配时需跳转到的新位置。假设有两个指针$s$ 和 $t$，其中 $s\in S$, $t\in T$。***KMP*** 需要保证文本串 $S$ 中的指针 $s$ 永远不后退（因为在文本串中，之前匹配过的所有字符已经没有用了——都是匹配完成或者已经失配的）。

这意味着当模式串的第 $i$ 位失配时，前 $i-1$ 位完成匹配。将模式串 $T$ 的指针 $t$ 移动到 $t'$ 的位置 ，使得 $T(1,2,...,t'-1)$ 和 $S(s-t’+1,...,s-1)$ 逐位相等，从而保证指针 $s$ 不用回退。且 $t'$ 的位置仅与模式串本身有关。

所以 ***KMP*** 数组（即用于确定失配后变化位置的数组，针对模式串而非文本串）记录的是：

> 移位法则：在**模式串$T$**中，对于每一位 $T(i)$ ,其 ***KMP*** 数组应当是记录一个位置 $j$ ，满足 $T(i)=T(j)(j\leq i)$，且 $T(k)=T(k+i-j)(k=1,2,...j-1;j\neq 1)$
>
> //定义真子串：字符串中最长相等真前缀和真后缀

形式化理解：对于模式串 $T$ 而言的每个前缀 $T'$，用 ***KMP*** 数组记录**到 $T'$ 为止的模式串前缀**的**真前缀**和**真后缀**的最大相同的位置。

当需要用 ***KMP*** 数组换位时，当且仅当未完全匹配，失配后仅部分前缀需要快速移动，因此 ***KMP*** 中不包括模式串本身。（包括模式串则匹配完成，与失配的定义矛盾；若跳过整个 $T$，则无意义）

令 $kmp[i]$ 表示当匹配到模式串 $T$ 的第 $i$ 位之后失配该跳转到的位置（可以跳过的字符个数）。

对于模式串的第一位和第二位而言，只能回跳到第一位，因为是 ***KMP*** 是要将真前缀跳跃到与它相同的真后缀上去（通常也可以反着理解），所以当 $i =1,2$ 时,相同的真前缀只会是 $T(1)$ 这一个字符。

***KMP*** 的代码需要实现两个部分：$kmp$ 数组的生成和模式串与文本串的匹配。其中，生成部分等价于模式串与模式串的匹配，从而得出 $kmp$ 数组。（$kmp$ 数组元素的值相当于 $fail$ 指针指向的位置）
$$
t'=kmp[i-1]+1
$$

```cpp
void KMP(string s, string t) {
   //std::cin >> s >> t;
   s = " " + s, t = " " + t; 
   int n = t.size() - 1;
   int m = s.size() - 1;
   std::vector<int> kmp(n + 1);
   for(int i = 2, j = 0; i <= n; i++) {
      while(t[i] != t[j + 1] && j) j = kmp[j];
      if(t[j + 1] == t[i]) ++ j;
      kmp[i] = j;
   }
   for(int i = 1, j = 0; i <= m; i++) {
      while(s[i] != t[j + 1] && j) j = kmp[j];//失配回跳
      if(t[j + 1] == s[i]) ++ j;
      if(j == n) {
         //std::cout << i - n + 1 << '\n';//position
         j = kmp[j];//回跳，继续匹配
      }
   }
   //for(int i = 1; i <= n; i++) std::cout << kmp[i] << " \n"[i == n];
   /*
   		T:ababc
   		K:00120
   */
}
```



##### 最小循环子串

> 模式串 $T$ 的最小循环节称为最小循环子串，其长度为 $L = |T| - kmp[|T|]$。
>
> 最小循环子串若存在，任意一段可以表示为：

$$
T(K,K+1,...L+K-1),K=0,1,...
$$

计算一个字符串由多少个相同子串连接而成：

```cpp
string s;
cin >> s;
s = " " + s;
int n = s.size() - 1;
std::vector<int> f(n + 1);
for(int i = 2, j = 0; i <= n; i++) {
	if(s[i] != s[j + 1] && j) j = f[j];
	if(s[i] == s[j + 1]) ++ j;
	f[i] = j;  
}	
if(n % (n - f[n]) == 0) //%=0 说明最小循环子串存在
	cout << n / (n - f[n]) << '\n';
else cout << 1 << '\n';
```



##### 失配树

给定字符串 $S$ ，定义 $pre_k = S(1,2,...k),suf_k=S(|S|-k+1,...|S|)$，其中 $1\leq k \leq |S|$。

定义 $Border(S)$（简写为 $B(S)$）为对于 $\forall i\in [1,|S|]$，满足 $pre_i=suf_i$ 的 $pre_i$ 的集合。 $B(S)$ 中的每个元素称为 $S$ 的 “$border$”。

> 定义 $S$ 的一个 $border$ 为既是 $S$ 的真前缀又是 $S$ 的真后缀的字符串。

现有 $m$ 组询问，每组询问给定 $i,j$，求 $S$ 的 $pre_i$ 和 $pre_j$ 的最长公共 $border$ 的长度。

由 ***KMP*** 定义得，$B(S)=\{S_{1...kmp[|S|]},S_{1,...kmp[kmp[|S|]]},S_{1,...kmp[kmp[kmp[|S|]]],...}\}$，且 $B(S)$ 一定是有限集合。注意到所求 $border$ 一定是两个前缀的公共后缀和公共前缀，因此在求出两串的 `lcp` 以后只需要在其中任意一个串上找到其最长的长度不超过 `lcp` 长度的 $border$ ，那么该串即为两串的最长公共 $border$。

首先对原串运行 ***KMP***，考虑 $S$ 的每个前缀，将前缀和其最长 $border$ 连边，构成`fail`树，在树上运行 $LCA$，得到两个给定前缀的最长公共 $border$。注意，因为 $border$ 不能是整个串，两个点必须都跳至少一次，即当两个前缀在`fail`树上是祖先—后代关系时，结果是祖先的父亲而非祖先自身。



练习题：

- [x] [POI2006 OKR-Periods of Words](https://www.luogu.com.cn/problem/P3435)
- [BOI2009 Radio Transmission](https://www.luogu.com.cn/problem/P4391)
- [EC Final 2022\] Binary String](https://www.luogu.com.cn/problem/P9717)
- [P10634 BZOJ2372 music](https://www.luogu.com.cn/problem/P10634)
- [POI2012\] PRE-Prefixuffix](https://www.luogu.com.cn/problem/P3546)
- [USACO15FEB\] Censoring G](https://www.luogu.com.cn/problem/P3121)
- [NOI2014 动物园](https://www.luogu.com.cn/problem/P2375)
- [SCOI2016 围棋](https://www.luogu.com.cn/problem/P3290)
- [POI2007 OSI-Axes of Symmetry](https://www.luogu.com.cn/problem/P3454)
- [POI2005\] SZA-Template](https://www.luogu.com.cn/problem/P3426)
- [P4173 残缺的字符串](https://www.luogu.com.cn/problem/P4173)
- [Cnoi2021 符文破译](https://www.luogu.com.cn/problem/P8112)
- [BOI2009 Radio Transmission](https://www.luogu.com.cn/problem/P4391)
- [CERC2018 The ABCD Murderer](https://www.luogu.com.cn/problem/P7456)
- [COCI2011-2012#4 KRIPTOGRAM](https://www.luogu.com.cn/problem/P8085)
- [蓝桥杯 2024 国 A 重复的串](https://www.luogu.com.cn/problem/P10581)



### 2. Z函数



