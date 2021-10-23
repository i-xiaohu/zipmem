还有些东西是创建FM-index的，我需要添加上去，就命名程序名为fmidx

虽然qpz情况很大程度上缓解，但是浪费时间的情况还是很严重，我明明还有很多事情要做，也有很多技能需要培养，可现在很多时间浪费在看视频上。

我单独写一个check文件，用来检查种子的一致性，各种情况的统计，比如种子的数量，或者种子的正确率等等。
再检查AS，NM，POS的一致性。
双端问题没有时间再做了，论文抓紧
另外给supplementary seeding设置几个等级
既然BWA-MEM是ground truth，完全可以根据Zipmem效果不好的种子丢失情况，来判断优化的方向。

如果是双端，那么需要一次智能匹配，这个阶段不妨交给双端的seeding，如果这不耗时的话。如果耗时，可以交给后续的extend，这样就不算我们的时间消耗了。
对了，测试一下，直接单端压缩reads1，和同时压缩reads1和reads2，出来的reads1效果，是不是一样的。
如果不是，那么测试这两个谁更适合我想到的双端排序。

我那个toy有点乱，其实可以稍微整理一下。如果hfastq输出进度条有问题，那就直接像whisper一样，按点刷新。

包括CS-generating阶段，总感觉它有提升的空间。

介绍完seeding之后，加上一句话，无论哪种seeding，都会在新型数据压缩结果下面产生大量的重复。
我比较担心，CS-seeding其实有进步空间。尤其在supplementary seeding这方面；
在人类数据上，check函数跑的有点慢；我需要优化一下，并且给点及时的反馈。

李恒在BWA-SW里讲过关于real dataset比对准确度的评价方式：时间统一使用CPU time，AS作为评价比对位置正确的唯一指标；而我只是衡量了AS的一致率，高低差异没有体现出来。比对质量的差异同样没有体现出来（Li,H 的方法是把both low mapping舍弃，而我只保留both high mapping，这是我的问题所在）。

论文里忘记声明所有的程序都是默认参数执行

## 测试实验

我暂时用以数据集为主导整理实验结果，运行服务器为comput14；20核心；Intel(R) Xeon(R) CPU E5-2680 v2 @ 2.80GHz；测试线程数为16.

### spring_SRR1562082_1

单纯查看种子的情况，完整匹配占绝大多数，错误数也在容忍范围之内。
Exactly matched seeds   97.39
Mismatched number of noisy seeds        1.24

#### Zipmem

```plaintext
gencs_time      11.7    1.5     7.8X
csseed_time     12.4    0.9     14.4X
supseed_time    16.5    1.8     9.3X
input_reads     5825771
read_len        101
Compression ration:     20.86
Mismatches number:      0.09
Reused seeds:           5.73
Supplementary seeds:    0.02
```

#### 比对一致性

```plaintext
MAX_rss:     3.8 GB
CPU_time:    00:00:52  i.e.  52.9 sec
Wall_clock:  00:00:16  i.e.  16.4 sec
Job_CPU:     3.23

Primary: 5825771
POS:    99.95 % = 5822646 / #pri
AS:     98.85 % = 5759008 / #pri
NM:     99.90 % = 5820047 / #pri
```

### spring_ERP001775_s5_1

#### Zipmem

目前的时间消耗非常离谱，我还在调查。不明白CS-gen为什么耗费了那么长时间，Supplementary seeding反而时间占比很低，跟我的预期差距很大。

```plaintext
MAX_rss:     13.7 GB
CPU_time:    153:55:08  i.e.  554108.2 sec
Wall_clock:  13:01:32  i.e.  46892.1 sec
Job_CPU:     11.82
```

在人类数据上，可以看到包含错误的种子数其实还是变多了的；这足够证明，与CS的差异，很多都是排序错误，而并非序列本身的错误。

```plaintext
reads number            1075379007
seeds number            27.48
Exactly matched seeds   72.66
Mismatched number of noisy seeds        1.84
```

调查CPU利用率总是跟不上的原因。把zipmem输出seeds的方式改为gzwrite输出，在SMP3上测试，发现读入花费3-5s，处理花费3-5S，且CPU利用率比较高，但是输出耗费时间是处理时间的数倍，大约几十秒，是绝对的瓶颈。改为使用fprintf标准输出后，平均每次输出大约花费4S，比gzwrite函数函数要快5-40倍；看来无法通过gzip的方式减少硬盘写入时间，压缩时间就已经严重超载了。

关键就在于zipmem太快了，处理时长为3-5S，这几乎和读入一样快，所以处理时间必然小于IO时间，必然不能实现重叠；如果写入是瓶颈还可以多个线程去写，但是读入是瓶颈，就很难办了。
解决办法有三个：
1）加快所有数据集的读入速度。推荐方法有：~~多线程读入~~、读入原fastq文件、~~优化读入函数~~。
2）更换数据集，让处理时间长于IO时间
3）把seeding和extending联合起来；

读入原fastq文件已经比读入gzip文件快一些了，但是读入时间好像还是很长，我不输出任何文件，CPU利用率也只有13.16（具体结果如下），那么最好的方案应该就是第三个了；

需要在补充材料中解释一下，我虽然提供了zipmem和BWA-MEM seeding，但我记录的时间来源于将seeding和extending结合在一起的时间，而如果真的是先seeding，再extending，总时间消耗会更长，因为包括了seeds的IO时间。

```plaintext
注：运行机器为空闲状态的SMP3
Input_time      2956.1
Process_time    52677.7 3434.2  15.3X
output_time     393.6
gencs_time      2731.8  322.5   8.5X
  Allocating    198.2   99.3
  Assemble      2495.4  202.6
  Sumup         37.9    20.4
csseed_time     23678.9 1467.5  16.1X
supseed_time    26266.9 1644.2  16.0X
input_reads     1075379007
read_len        101
Compression ration      14.83
Mismatches number       0.26
Reused seeds            27.48
Supplementary seeds     4.38

MAX_rss:     13.9 GB
CPU_time:    14:52:09  i.e.  53529.2 sec  
Wall_clock:  01:07:47  i.e.  4067.5 sec
Job_CPU:     13.16
```

#### BWA-MEM seeding

```plaintext
MAX_rss:     10.1 GB
CPU_time:    39:03:25  i.e.  140605.5 sec
Wall_clock:  02:10:03  i.e.  7803.5 sec
Job_CPU:     18.02
```

#### 比对一致性

但让人不满意的是，zipmem的种子延伸消耗的内存特别夸张，我需要程序运行角度盯防一下，期望是不影响最终结果的极少数情况。

```plaintext
Extending zipmem seeds.
MAX_rss:     51.2 GB
CPU_time:    27:54:49  i.e.  100489.7 sec
Wall_clock:  02:13:11  i.e.  7991.7 sec
Job_CPU:     12.57

Extending BWA-MEM seeds
MAX_rss:     27.7 GB
CPU_time:    27:23:39  i.e.  98619.8 sec
Wall_clock:  02:10:59  i.e.  7859.3 sec
Job_CPU:     12.55
```

人类效果如下，如果过滤掉低质量的，效果应该会更好

```plaintext
Primary: 1075379007
POS:    99.35 % = 1068383868 / #pri
AS:     97.92 % = 1052975186 / #pri
NM:     99.58 % = 1070866663 / #pri
```
