# Terasure库用户手册


---

>Terasure算法库实现了在2^8的伽罗华域内的reed solomon编码，MSR编码，MBR编码，Hitchhiker编码
使用了AVX2向量指令集，每次计算32Byte的数据，大大地提高了运算效率

## 公用头文件介绍
### arch.h
支持AVX2向量指令集运算的头文件

```
typedef __m256i encode_t 
//将类型__m256i定义为encode_t类型
```
```
void init_arch()
//初始化函数，每次使用前都需要调用
```
```
__attribute__((always_inline))
static inline encode_t xor_region(encode_t input1, encode_t input2)
//实现了两个向量之间的异或操作(伽罗华域内的加减法都是异或操作)
```
```
__attribute__((always_inline))
static inline encode_t multiply_region(encode_t input, uint8_t x)
//实现了一个向量与一个整数之间的乘法(伽罗华域内的乘法实现)
```
```
__attribute__((always_inline))
static inline encode_t zero()
//给一个AVX2向量赋值为0
```
```
__attribute__((always_inline))
static void inline store(encode_t *dst, encode_t val)
//给一个向量赋值到内存空间，等同于*dst = val
```
```
__attribute__((always_inline))
static void inline prefetch(const encode_t *src)
//预取一个内存地址处的数据
```
>使用```__attribute__((always_inline))```是为了提升性能



### gf.h
预处理在2^{8}的伽罗华域中乘除法的计算
```
void gf_init()
//打表预处理乘除法的计算
```
```
static inline uint8_t gf_div(uint8_t a, uint8_t b)
static inline uint8_t gf_mul(uint8_t a, uint8_t b)
//伽罗华域中两个整数的乘除法计算
```


## reed solomon编码
rs编码是最基本的纠删码，运行效率非常高，cauchy矩阵的结构有systematic的特性，vandermonde矩阵可以加密所有数据。

使用时只有broken_id和broken_num需要用户显式指定
### 接口介绍(include)
#### reed_solomon.h
纠删码算法的函数声明

```
struct rs_conf_t {
	int n;
	int k;
	int r;

	void *(*allocator)(size_t);

	void (*deallocator)(void *);

	uint8_t **matrix; //matrix of coeficients for encoding
	int matrix_type; //0 cauchy, 1 vandermonde, 2 cauchy-vandermonde
};
typedef struct rs_conf_t rs_conf;
//reed solomon纠删码算法的基本结构
//n是编码后的节点数量
//k是编码前的数据组数
//r是冗余的节点数量，用于恢复
//allocator和deallocator是动态内存分配函数，一般是malloc和free
//matrix是编码时的系数矩阵
//matrix_type是矩阵类型，0是cauchy矩阵，1是vandermonde矩阵，2是cauchy-vandermonde矩阵
```

```
struct rs_regenerate_context_t {
	uint8_t **inv_matrix;
	int *broken_id;
	int broken_num;
	int *recover_id;
};
typedef struct rs_regenerate_context_t rs_regenerate_context;
//这是做数据恢复时的结构
//inv_matrix是用于恢复的编码矩阵的逆矩阵
//broken_id是损坏的节点的id
//broken_num是损坏的节点个数
//recover_id是用于恢复的节点id
```
```
void rs_init(rs_conf *conf, int n, int k, void *(allocator)(size_t), void (*deallocator)(void *));
//用于初始化纠删码算法的结构，给n，k，r赋值，设置动态内存分配函数

void vandermonde_matrix(rs_conf *conf);
//给纠删码算法结构里的系数矩阵初始化为vanderomonde矩阵，矩阵维度为n*k，因为vandermonde矩阵没有systematic的特性

void cauchy_matrix(rs_conf *conf);
//给纠删码算法结构里的系数矩阵初始化为cauchy矩阵，矩阵维度为r*k，因为cauchy矩阵有systematic的特性，编码后的前k个节点存储的是源数据，计算量较小

void rs_encode(int len, rs_conf *conf, uint8_t **data, uint8_t **output, uint8_t **memory_pre_allocated);
//reed solomon的编码过程，最终编码后的结果会被存储在output中
//len是每个节点的数据大小
//conf是纠删码算法的结构
//data是编码前的数据，维度为k*1*len
//output存储编码后的所有节点的数据，维度为n*1*len
//memory_pre_allocated是使用cauchy矩阵时暂存计算结果的矩阵，使用vandermonde矩阵的时候设为NULL，维度为r*1*len

void rs_decode(int len, rs_conf *conf, rs_regenerate_context *context,uint8_t **data, uint8_t **output);
//reed solomon的解码过程，如果有节点损坏，可用剩余数据，恢复出源数据，把源数据存储在output中
//len和conf如上述
//context是记录当前节点损坏情况，以及恢复时用的逆矩阵
//data是编码后的节点的数据，其中损坏的节点的数据可以设为任意值，维度为n*1*len
//output存储解码后恢复出的源数据，维度为k*1*len

void rs_regenerate(int len, rs_regenerate_context *context, rs_conf *conf, uint8_t **data, uint8_t **output);
//reed solomon的恢复损坏节点的过程，直接把损坏的节点存储的数据给恢复出来，存于output中
//len，conf，context，data和解码一样
//output存储损坏的节点存储的编码后的数据

void rs_regenerate_context_init(rs_conf *conf, rs_regenerate_context *context,int *broken_id, int broken_num);
//初始化context的内容
//broken_id为损坏的节点的id
//broken_num为损坏的节点的数量，不能超过r个

void rs_free_regenerate_context(rs_conf *conf, rs_regenerate_context *context);
//释放context中动态分配的内存空间

void rs_free_conf(rs_conf *conf);
//释放conf中动态分配的内存空间

```

### 函数实现(src)
实现头文件中声明的函数，以及一些其他的功能函数

#### reed_solomon.c
(头文件中解释过的函数就不再解释了)
这些函数都不需要用户显式调用
```
__attribute__((always_inline))
static inline void super_fast_cal(int len, const int K, const int R, uint8_t **data, uint8_t **output, uint8_t **Matrix)
//用于编码计算的函数
//len为每个节点存储的数据长度
//K和R为系数矩阵的维度，系数矩阵为R*K
//data为编码前的数据
//output位编码后的数据
//Matrix为系数矩阵
```
```
__attribute__((always_inline))
static inline void super_fast_reg(int len, const int K, const int R, uint8_t **data, uint8_t **output, uint8_t **Matrix, int *recover_id)
//用于恢复和解码的计算函数
//前面一些参数和上述相同
//recover_id表明用于恢复和解码的节点id，在计算过程中仅仅是节点的index和编码时不同
```
```
bool check_broken(int x, int *broken_id, int broken_num)
//检查节点x是否是损坏的
```
```
void invert_matrix(uint8_t **in_mat, uint8_t **out_mat, int n){
//计算in_mat的逆矩阵，输出到out_mat，n是矩阵维度
```
>以上函数都不需要直接被用户调用，都在头文件中声明的接口函数中被调用了

### 测试(test)
用于测试编解码正确性和效率的文件

#### rs_test_cauchy.c
对于要使用AVX2向量指令集的数组，都需要对齐，否则会出现段错误，如下：
```
for(int i = 0; i < k; i++){
		posix_memalign((void **)&(data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
		posix_memalign((void **)&(decode_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
		for(int j = 0; j < DATA_SIZE; j++){
			data[i][j] = (uint8_t)(rand() & 0xff);
			decode_data[i][j] = 0;
		}
	}

	for(int i = 0; i < r; i++)
		posix_memalign((void **)&(memory_pre_allocated[i]), 64, sizeof(uint8_t) * DATA_SIZE);

	for(int i = 0; i < n; i++)
		posix_memalign((void **)&(encode_data[i]), 64, sizeof(uint8_t) * DATA_SIZE);
```
```
void test_decode_correctness(rs_conf conf, uint8_t **data, uint8_t **encode_data, uint8_t **memory_pre_allocated, uint8_t **decode_data)
//测试解码的正确性

void test_regenerate_correctness(rs_conf conf, uint8_t **data, uint8_t **encode_data, uint8_t **memory_pre_allocated, uint8_t **decode_data)
//测试恢复损坏数据的正确性
```
两个函数内，都随机生成了4个broken_id，测试一百轮，都得到了正确结果。

#### rs_test_vandermonde.c
vandermonde矩阵和上述cauchy矩阵测试相同。也能够得到正确结果。

#### rs_benchmark_cauchy.c
```
void test_encode_speed(rs_conf conf, uint8_t **data, uint8_t **encode_data, uint8_t **memory_pre_allocated)
//测试编码的性能

void test_decode_speed(rs_conf conf, uint8_t **encode_data, uint8_t **decode_data){
//测试解码的性能
```
对于性能的计算，使用公式如下：
```
printf("regenerate Throughput: %.2fGB/s\n", test_round * (double) DATA_SIZE * n / ((end - start) / (double) CLOCKS_PER_SEC) / 1024 / 1024 / 1024);
//测试test_round轮
//DATA_SIZE是每个节点的数据量
```

### 测试结果
每个节点数据大小为64MB，跑100轮取平均
#### 个人台式机评测结果
```
/**********Cauchy Performance************/
Total Clock Time: 5.78s
Encode Throughput: 15.13GB/s
Total Clock Time: 6.11s
Regenerate Throughput: 14.33GB/s
```

#### 服务器评测结果
```
/**********Cauchy Performance************/
Total Clock Time: 8.59s
Encode Throughput: 10.19GB/s
Total Clock Time: 8.89s
Regenerate Throughput: 9.84GB/s
```

## MBR编码
MBR编码可以最小化带宽使用，当一个节点损坏的时候，最小化带宽使用来再生这个节点，当不多于n-k个节点损坏的时候，可以直接恢复出原数据

实现时的二维数组都用一维数组来表示

使用时只有broken_id和broken_num需要用户显式指定
### 接口介绍(include)
#### mbr.h
```
struct mbr_conf_t {
	int n;
	int k;
	int d;

	void *(*allocator)(size_t);

	void (*deallocator)(void *);

	uint8_t **matrix; //matrix of coeficients for encoding
};
typedef struct mbr_conf_t mbr_conf;
//n为编码后的节点数量
//d为每个节点存的数据量
//k为recover所有数据需要的最小数据量
//一般有k<=d == n-1

struct mbr_recover_context_t{
	uint8_t **inv_matrix;
	uint8_t **matrix_half;
	int *broken_id;
	int broken_num;
	int *recover_id;
};
typedef struct mbr_recover_context_t mbr_recover_context;
//inv_matrix和matrix_half是计算时使用到的矩阵，不需要用户显式设置
//broken_id是顺坏的数据节点
//broken_num是损坏的节点个数，最多可以损坏n-k个
//recover_id是用于恢复的k个节点

struct mbr_regenerate_context_t {//regenerate single node
	uint8_t **inv_matrix;
	int broken_id;
	uint8_t *broken_vector;
};
typedef struct mbr_regenerate_context_t mbr_regenerate_context;
//inv_matrix是再生时使用到的矩阵
//broken_id是需要再生的节点
//broken_vector是再生节点的系数矩阵行
```

```
void mbr_init(mbr_conf *conf, int n, int k, int d, void *(allocator)(size_t), void (*deallocator)(void *));
//初始化mbr_conf结构的函数

void mbr_encode(int len, mbr_conf *conf, uint8_t **data, uint8_t **output, uint8_t **memory_pre_allocated);
//mbr编码函数，变量使用和rs编码类似
//data维度为d*d*len
//output维度为n*d*len
//memory_pre_allocated维度为(n-k)*d*len

void mbr_recover(int len, mbr_conf *conf, mbr_recover_context *context, uint8_t **data, uint8_t **output);
//mbr恢复数据的函数
//data维度n*d*len
//output维度d*d*len

void mbr_recover_context_init(mbr_conf *conf, mbr_recover_context *context, int *broken_id, int broken_num);
//mbr恢复数据的context_init函数

void mbr_regenerate(int len, mbr_conf *conf, mbr_regenerate_context *context, uint8_t **data, uint8_t **output);
//mbr再生一个节点数据的函数
//data维度为n*d*len
//output维度为1*d*len

void mbr_regenerate_context_init(mbr_conf *conf, mbr_regenerate_context *context,int broken_id);
//mbr再生节点的context_init

void mbr_free_regenerate_context(mbr_conf *conf, mbr_regenerate_context *context);
void mbr_free_conf(mbr_conf *conf);
void mbr_free_recover_context(mbr_conf *conf, mbr_recover_context *context);
//动态分配的内存释放函数
```
### 测试结果
``n=14, d=13, k=10, DATA_SIZE=(1<<22)*d*d``
#### 个人台式机评测结果
```
Total Clock Time: 8.00s
Encode Throughput: 8.89GB/s
Total Clock Time: 32.08s
Recover Throughput: 2.06GB/s
Total Clock Time: 6.06s
Regenerate Throughput: 0.84GB/s
```

#### 服务器评测结果
```
Total Clock Time: 34.27s
Encode Throughput: 2.07GB/s
Total Clock Time: 46.38s
Recover Throughput: 1.42GB/s
Total Clock Time: 8.49s
Regenerate Throughput: 0.60GB/s
```

## MSR编码
MSR编码能够最小化存储容量
### 接口介绍(include)
#### msr.h
```
struct msr_conf_t {
    int n;
    int k;
    int r;

    int alpha;
    int beta;
    int groups;

    void *(*allocate)(size_t);

    void (*deallocate)(void *);

    uint8_t *node_companion;
    uint8_t *z_companion;
    uint8_t *theta;

    int nodes_round_up;
};
typedef struct msr_conf_t msr_conf;
//存储MSR编码使用的变量

struct msr_encode_context_t {
    size_t encoding_buf_size;

    uint8_t *matrix;

    int *survived;
    int survive_cnt;


    int *erased;
    int erase_cnt;

    bool *is_erased;
    int *erase_id;

    int *sigmas;
    int sigma_max;

};
typedef struct msr_encode_context_t msr_encode_context;
//存储MSR encode时的变量

struct msr_regenerate_context_t {
    size_t regenerate_buf_size;

    uint8_t *matrix;
    uint8_t *u_matrix;

    int broken;

    int *z_num;
    int *z_pos;
    int *z_comp_pos;


};
typedef struct msr_regenerate_context_t msr_regenerate_context;
//存储MSR再生时的变量
```
```
void msr_fill_encode_context(msr_encode_context *context, const msr_conf *conf, uint8_t **data);
//初始化msr encode context的函数

void msr_fill_regenerate_context(msr_regenerate_context *context, const msr_conf *conf, int broken);
//初始化msr regenerate context的函数

void msr_encode(int len, const msr_encode_context *context, const msr_conf *conf, uint8_t *buf, uint8_t **data,
                uint8_t **output);
//msr encode的函数

void msr_regenerate(int len, const msr_regenerate_context *context, const msr_conf *conf, uint8_t *buf, uint8_t **data,
                    uint8_t *output);
//msr regenerate一个节点的函数

void msr_get_regenerate_offset(int len, const msr_regenerate_context *context, const msr_conf *conf, int *offsets);

int msr_init(msr_conf *conf, int n, int k, void *(*allocate)(size_t), void (*deallocate)(void *));
//初始化msr conf的函数

void msr_free_conf(msr_conf *conf);
void msr_free_encode_context(const msr_conf *conf,msr_encode_context *context);
void msr_free_regenerate_context(const msr_conf *conf,msr_regenerate_context *context);
//三个free函数
```

### 测试结果
#### 个人台式机评测结果
```
n:14 r:4
Total Clock Time: 0.34s
Encode Throughput: 5478.82MB/s
Total Clock Time: 0.40s
Decode Throughput: 4657.37MB/s
Total Clock Time: 0.11s
Regenerate Throughput: 1235.36MB/s
```
#### 服务器评测结果
```
n:14 r:4
Total Clock Time: 0.50s
Encode Throughput: 3758.10MB/s
Total Clock Time: 0.84s
Decode Throughput: 2236.96MB/s
Total Clock Time: 0.26s
Regenerate Throughput: 516.22MB/s
```

## Hitchhiker
### 测试结果
#### 个人台式机评测结果
```
n:14 r:4
Total Clock Time: 0.22s
Encode Throughput: 8453.02MB/s
Total Clock Time: 0.11s
Regenerate Throughput: 1206.96MB/s
```
#### 服务器评测结果
```
n:14 r:4
Total Clock Time: 0.31s
Encode Throughput: 6061.45MB/s
Total Clock Time: 0.19s
Regenerate Throughput: 706.41MB/s
```
## 机器配置
### 个人台式机配置
CPU配置如下：
```
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                8
On-line CPU(s) list:   0-7
Thread(s) per core:    2
Core(s) per socket:    4
Socket(s):             1
NUMA node(s):          1
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 158
Model name:            Intel(R) Core(TM) i7-7700K CPU @ 4.20GHz
Stepping:              9
CPU MHz:               4399.951
CPU max MHz:           4500.0000
CPU min MHz:           800.0000
BogoMIPS:              8400.00
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              8192K
NUMA node0 CPU(s):     0-7
```
内存情况:
```
              total        used        free      shared  buff/cache   available
Mem:       32861844     4306892    25276836      794340     3278116    27228680
Swap:      33472508           0    33472508
```
### 服务器配置
CPU如下：
```
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                24
On-line CPU(s) list:   0-23
Thread(s) per core:    2
Core(s) per socket:    6
Socket(s):             2
NUMA node(s):          2
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 79
Model name:            Intel(R) Xeon(R) CPU E5-2643 v4 @ 3.40GHz
Stepping:              1
CPU MHz:               1744.824
BogoMIPS:              6808.11
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              20480K
```
内存情况：
```
              total        used        free      shared  buff/cache   available
Mem:      131807172    24388436    26256000      295284    81162736   105717360
Swap:      32899068     5417040    27482028
```

## 参考资料

 1. Alexandros G. Dimakis, P. Brighten Godfrey, Yunnan Wu, Martin Wainwright and Kannan Ramchandran,          "Network Coding for Distributed Storage Systems".
 2. K. V. Rashmi, Nihar B. Shah, and P. Vijay Kumar, "Optimal Exact-Regenerating Codes for Distributed        Storage at the MSR and MBR Points via a Product-Matrix Construction".
 3. James S. Plank, Kevin M. Greenan, Ethan L. Miller, "Screaming Fast Galois Field Arithmetic Using Intel     SIMD Instructions".
 4. KV Rashmi, Preetum Nakkiran, Jingyan Wang, Nihar B. Shah, and Kannan Ramchandran, "Having Your Cake and Eating It Too: Jointly Optimal Erasure Codes for I/O, Storage, and Network-bandwidth".
 5. K. V. Rashmi, Nihar B. Shah, P. Vijay Kumar, Kannan Ramchandran, "Explicit Construction of Optimal Exact Regenerating Codes for Distributed Storage".
 6. Nihar B. Shah, "On Minimizing Data-Read and Download for Storage-Node Recovery".
 7. M. Nikhil Krishnan and P. Vijay Kumar, "On MBR codes with replication".
 8. Sian-Jheng Lin, Wei-Ho Chung, "Novel Repair-by-Transfer Codes and Systematic Exact-MBR Codes with Lower Complexities and Smaller Field Sizes".