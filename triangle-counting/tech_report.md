# rapids@hkust 技术报告

## 基本算法介绍

Main Logic, 详细请查看[triangle_counting_main.cpp](triangle_counting_main.cpp)文件: 

* Load Bin File (Direct IO).
* Swap Elements in the edge list, and In-Place Parallel Sort.
* Convert Edge List to CSR (Fusion了去重和去self-loops和deg-histogram以及CSR转化).
* Graph Reordering (Degree-Descending).
* Triangle Counting.

### 注意点

转为Degree-Oriented Directed Graph(DODG)的图之后, 每个triangle只count一次。
Triangle Counting部分划分了两个Range, 第一个Range为degree大的点, 第二的Range为其他点。

* 第一个部分采用了Sparse-Bitmap和Dense Bitmap intersection
* 第二个部分采用了向量化的Merge on the sorted arrays

Edge List to CSR转化时候内存不足，所以输出sorted edge list到文件先, 然后mmap进来, 利用操作系统的page-swap来使得程序可以正确运行。
sorted edge list的访问是sequential的, 所以操作系统的page-prefetch可以正确得预测。

## 并行化设计

* 预处理阶段: 
在[pre_processing.h](pre_processing.h) (`EdgeListHistogram`和`Reorder` (Graph Reordering))
和[pre_processing_dodg.h](pre_processing_dodg.h) (预处理Edge List到Degree Oriented Directed Graph(DODG), degree-descending) 
大量使用了自研的深度优化的并行primitives, 详细请查看 [util/local_buffer.h](util/local_buffer.h), [util/primitives.h](util/primitives.h) 

Triangle Counting部分划分了两个Range, 第一个Range为degree大的点, 第二的Range为其他点

* 整体的并行采用了对点的计算划分的方式。
* 第一个部分使用了[util/libpopcnt.h](util/libpopcnt.h)进行向量化处理(利用了`AVX2`指令集)。
* 第二个部分采用了向量化的Merge on the sorted arrays。

## 算法优化

所有部分都是并行的，并且Triangle-Counting部分使用了向量化(利用了`AVX2`指令集)。

## 详细算法设计与实现

* Load Bin File (Direct IO).

```cpp
        // Load Bin File (DIO).
        auto file_name = string_option->value(0);
        auto file_fd = open(file_name.c_str(), O_RDONLY | O_DIRECT, S_IRUSR | S_IWUSR);
        Edge *edge_lst = (Edge *) memalign(PAGE_SIZE, size + IO_REQ_SIZE);
        auto max_omp_threads = omp_get_max_threads();
        Timer io_timer;
        size_t read_size = 0;
#pragma omp parallel num_threads(IO_QUEUE_DEPTH)
        {
#pragma omp for schedule(dynamic, 1) reduction(+:read_size)
            for (size_t i = 0; i < size; i += IO_REQ_SIZE) {
                auto it_beg = i;
                auto *chars = reinterpret_cast<uint8_t *>(edge_lst);
                auto ret = pread(file_fd, chars + it_beg, IO_REQ_SIZE, it_beg);
                if (ret != IO_REQ_SIZE) {
                    log_error("Err, %zu, %zu, %zu, %d", i, it_beg, IO_REQ_SIZE, ret);
                } else {
                    read_size += ret;
                }
            }
#pragma omp single
            log_info("%zu, %zu", read_size, size);
        }
        log_info("IO Time: %.6lfs, DIO-QPS: %.6lf GB/s", io_timer.elapsed(), size / io_timer.elapsed() / pow(1024, 3));
        log_info("Load File Time: %.9lfs", global_timer.elapsed());
```

* Swap Elements in the edge list, and In-Place Parallel Sort.

```cpp
Timer sort_timer;
        int32_t max_node_id = 0;
#pragma omp parallel for reduction(max: max_node_id) schedule(dynamic, 32*1024)
        for (size_t i = 0u; i < num_edges; i++) {
            if (edge_lst[i].first > edge_lst[i].second) {
                swap(edge_lst[i].first, edge_lst[i].second);
            }
            max_node_id = max(max_node_id, max(edge_lst[i].first, edge_lst[i].second));
        }
        log_info("Populate File Time: %.9lfs", global_timer.elapsed());
        // In-Place Parallel Sort.
        ips4o::parallel::sort(edge_lst, edge_lst + num_edges, [](Edge l, Edge r) {
            if (l.first == r.first) {
                return l.second < r.second;
            }
            return l.first < r.first;
        });
        log_info("Sort Time: %.9lfs", sort_timer.elapsed());
        auto num_vertices = static_cast<uint32_t >(max_node_id) + 1;
        log_info("Pre-Process Edge List Time: %.9lf s", global_timer.elapsed());
```

* Convert Edge List to CSR (Fusion了去重和去self-loops和deg-histogram以及CSR转化). `ConvertEdgeListToDODGCSR`请查看[pre_processing_dodg.h](pre_processing_dodg.h)

```zsh
graph_t g{.n=num_vertices, .m = 0, .adj=nullptr, .row_ptrs=nullptr};
        uint32_t *deg_lst;
        g.adj = nullptr;

        ConvertEdgeListToDODGCSR(num_edges, edge_lst, num_vertices, deg_lst, g.row_ptrs, g.adj,
                                 max_omp_threads, [&](size_t it) {
                    return !(edge_lst[it].first == edge_lst[it].second
                             || (it > 0 && edge_lst[it - 1] == edge_lst[it]));
                });
        g.m = g.row_ptrs[num_vertices];
        log_info("Undirected Graph G = (|V|, |E|): %lld, %lld", g.n, g.m);
        log_info("Mem Usage: %s KB", FormatWithCommas(getValue()).c_str());
```

* Graph Reordering (Degree-Descending)。详细请查看[pre_processing.h](pre_processing.h)

```cpp
  // 3rd: Reordering.
        vector<int32_t> new_dict;
        vector<int32_t> old_dict;
        munmap(edge_lst, size);

        auto *tmp_mem_blocks = (int32_t *) malloc(size / 2);
        auto *org = g.adj;
        ReorderDegDescendingDODG(g, new_dict, old_dict, tmp_mem_blocks, deg_lst);
        free(org);
        free(deg_lst);
```

* Triangle Counting. 详细请查看 [triangle_counting.h](triangle_counting.h)。

```cpp
        log_info("Mem Usage: %s KB", FormatWithCommas(getValue()).c_str());
        size_t tc_cnt = 0;
        tc_cnt = CountTriBMPAndMergeWithPackDODG(g, max_omp_threads);
        log_info("Mem Usage: %s KB", FormatWithCommas(getValue()).c_str());
        log_info("There are %zu triangles in the input graph.", tc_cnt);
        printf("There are %zu triangles in the input graph.\n", tc_cnt);    
```

## 实验结果与分析

* `./tc -f ~/datasets/s28.e15.kron.edgelist.bin`

```
2019-11-25 08:20:19 INFO  (et: 0.000050 s) (func: main)  triangle_counting_main.cpp:42: File size: 32212254720
2019-11-25 08:20:19 INFO  (et: 0.000160 s) (func: main)  triangle_counting_main.cpp:43: #of Edges: 4026531840
2019-11-25 08:22:46 INFO  (et: 146.223736 s) (func: main)  triangle_counting_main.cpp:75: 32212254720, 32212254720
2019-11-25 08:22:46 INFO  (et: 146.224009 s) (func: main)  triangle_counting_main.cpp:77: IO Time: 146.223792s, DIO-QPS: 0.205165 GB/s
2019-11-25 08:22:46 INFO  (et: 146.224067 s) (func: main)  triangle_counting_main.cpp:81: Load File Time: 146.224024750s
2019-11-25 08:22:49 INFO  (et: 149.948342 s) (func: main)  triangle_counting_main.cpp:92: Populate File Time: 149.948299180s
2019-11-25 08:24:10 INFO  (et: 230.916343 s) (func: main)  triangle_counting_main.cpp:100: Sort Time: 84.692248286s
2019-11-25 08:24:10 INFO  (et: 230.916409 s) (func: main)  triangle_counting_main.cpp:102: Pre-Process Edge List Time: 230.916366628 s
2019-11-25 08:24:11 INFO  (et: 231.378374 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:35: [ConvertEdgeListToDODGCSR]: InitTime: 0.461940972 s
2019-11-25 08:26:03 INFO  (et: 343.240148 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:52: [ConvertEdgeListToDODGCSR]: Histogram Time: 112.323715825 s
2019-11-25 08:26:03 DEBUG (et: 343.605481 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:60: 3973862397
2019-11-25 08:26:05 INFO  (et: 345.166313 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:68: Mem Usage: 61,339,504 KB
2019-11-25 08:26:41 INFO  (et: 382.073342 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:71: Mem Usage: 61,274,684 KB
2019-11-25 08:26:43 INFO  (et: 383.318484 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:82: Allocate Inside (adj_lst)...
2019-11-25 08:26:43 INFO  (et: 383.318561 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:85: [ConvertEdgeListToDODGCSR]: PrefixSum Time: 152.402128395 s
2019-11-25 08:28:59 INFO  (et: 519.503705 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:104: [ConvertEdgeListToDODGCSR]: Total Conversion Time: 288.587272254 s
2019-11-25 08:28:59 INFO  (et: 519.503765 s) (func: main)  triangle_counting_main.cpp:118: Undirected Graph G = (|V|, |E|): 268435456, 3973862397
2019-11-25 08:28:59 INFO  (et: 519.503905 s) (func: main)  triangle_counting_main.cpp:119: Mem Usage: 50,129,628 KB
2019-11-25 08:29:00 INFO  (et: 520.651031 s) (func: BucketSortSmallBuckets)  primitives.h:210: [BucketSortSmallBuckets]: Mem Size Buckets: 1195391488, Bucket#: 2334749
2019-11-25 08:29:01 INFO  (et: 521.427602 s) (func: ReorderDegDescendingDODG)  pre_processing_dodg.h:143: Deg-descending time:  1.539978042 s
2019-11-25 08:29:04 INFO  (et: 524.973543 s) (func: Reorder)  pre_processing.h:57: [Reorder]: Finish PrefixSum Time: 3.545900899 s
2019-11-25 08:29:43 INFO  (et: 564.085329 s) (func: Reorder)  pre_processing.h:76: [Reorder]: Finish Reorder Time: 39.112 s
2019-11-25 08:29:45 INFO  (et: 565.144850 s) (func: main)  triangle_counting_main.cpp:133: Mem Usage: 19,741,036 KB
2019-11-25 08:29:46 INFO  (et: 566.774529 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:62: Stop Deg at [0, 1818]
2019-11-25 08:29:46 INFO  (et: 566.775100 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:64: finish init row_ptrs_end, max d: 4652, time: 1.630203395s
2019-11-25 08:29:47 INFO  (et: 567.149394 s) (func: PackWords)  triangle_counting.h:32: Finish Indexing: 2.004496779s
2019-11-25 08:49:56 INFO  (et: 1776.829108 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:126: Forward cost: 1211.684 s, Mem Usage: 19854352 KB
2019-11-25 08:49:56 INFO  (et: 1776.829154 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:127: Triangle Cnt: 199196078202
2019-11-25 08:49:56 INFO  (et: 1776.836080 s) (func: main)  triangle_counting_main.cpp:136: Mem Usage: 19,854,352 KB
2019-11-25 08:49:56 INFO  (et: 1776.836107 s) (func: main)  triangle_counting_main.cpp:137: There are 199196078202 triangles in the input graph.
There are 199196078202 triangles in the input graph.
```

* `./tc -f ~/datasets/s29.e10.kron.edgelist.bin`

```
2019-11-25 08:50:57 INFO  (et: 0.000037 s) (func: main)  triangle_counting_main.cpp:42: File size: 42949672960
2019-11-25 08:50:57 INFO  (et: 0.000119 s) (func: main)  triangle_counting_main.cpp:43: #of Edges: 5368709120
2019-11-25 08:54:13 INFO  (et: 195.302496 s) (func: main)  triangle_counting_main.cpp:75: 42949672960, 42949672960
2019-11-25 08:54:13 INFO  (et: 195.302799 s) (func: main)  triangle_counting_main.cpp:77: IO Time: 195.302638s, DIO-QPS: 0.204810 GB/s
2019-11-25 08:54:13 INFO  (et: 195.302850 s) (func: main)  triangle_counting_main.cpp:81: Load File Time: 195.302818199s
2019-11-25 08:54:18 INFO  (et: 200.265308 s) (func: main)  triangle_counting_main.cpp:92: Populate File Time: 200.265275975s
2019-11-25 08:56:06 INFO  (et: 308.553620 s) (func: main)  triangle_counting_main.cpp:100: Sort Time: 113.250744824s
2019-11-25 08:56:06 INFO  (et: 308.553687 s) (func: main)  triangle_counting_main.cpp:102: Pre-Process Edge List Time: 308.553655583 s
2019-11-25 08:56:07 INFO  (et: 309.457397 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:35: [ConvertEdgeListToDODGCSR]: InitTime: 0.903684799 s
2019-11-25 08:58:48 INFO  (et: 470.728103 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:52: [ConvertEdgeListToDODGCSR]: Histogram Time: 162.174390994 s
2019-11-25 08:58:49 DEBUG (et: 471.371367 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:60: 5325948976
2019-11-25 08:58:51 INFO  (et: 473.586071 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:68: Mem Usage: 61,735,428 KB
2019-11-25 08:59:45 INFO  (et: 527.864921 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:71: Mem Usage: 61,688,800 KB
2019-11-25 08:59:47 INFO  (et: 529.554735 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:82: Allocate Inside (adj_lst)...
2019-11-25 08:59:47 INFO  (et: 529.554803 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:85: [ConvertEdgeListToDODGCSR]: PrefixSum Time: 221.001091691 s
2019-11-25 09:03:16 INFO  (et: 738.164126 s) (func: ConvertEdgeListToDODGCSR)  pre_processing_dodg.h:104: [ConvertEdgeListToDODGCSR]: Total Conversion Time: 429.610414191 s
2019-11-25 09:03:16 INFO  (et: 738.164187 s) (func: main)  triangle_counting_main.cpp:118: Undirected Graph G = (|V|, |E|): 536870909, 5325948976
2019-11-25 09:03:16 INFO  (et: 738.164326 s) (func: main)  triangle_counting_main.cpp:119: Mem Usage: 55,572,564 KB
2019-11-25 09:03:18 INFO  (et: 740.110395 s) (func: BucketSortSmallBuckets)  primitives.h:210: [BucketSortSmallBuckets]: Mem Size Buckets: 1339237888, Bucket#: 2615699
2019-11-25 09:03:19 INFO  (et: 741.507095 s) (func: ReorderDegDescendingDODG)  pre_processing_dodg.h:143: Deg-descending time:  2.935480386 s
2019-11-25 09:03:26 INFO  (et: 748.816384 s) (func: Reorder)  pre_processing.h:57: [Reorder]: Finish PrefixSum Time: 7.309232582 s
2019-11-25 09:04:20 INFO  (et: 802.833119 s) (func: Reorder)  pre_processing.h:76: [Reorder]: Finish Reorder Time: 54.017 s
2019-11-25 09:04:22 INFO  (et: 804.423657 s) (func: main)  triangle_counting_main.cpp:133: Mem Usage: 29,219,024 KB
2019-11-25 09:04:25 INFO  (et: 807.236215 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:62: Stop Deg at [0, 1584]
2019-11-25 09:04:25 INFO  (et: 807.236300 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:64: finish init row_ptrs_end, max d: 4252, time: 2.812600182s
2019-11-25 09:04:25 INFO  (et: 807.597704 s) (func: PackWords)  triangle_counting.h:32: Finish Indexing: 3.174003793s
2019-11-25 09:29:16 INFO  (et: 2298.745532 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:126: Forward cost: 1494.322 s, Mem Usage: 29330520 KB
2019-11-25 09:29:16 INFO  (et: 2298.745583 s) (func: CountTriBMPAndMergeWithPackDODG)  triangle_counting.h:127: Triangle Cnt: 164015236714
2019-11-25 09:29:16 INFO  (et: 2298.752393 s) (func: main)  triangle_counting_main.cpp:136: Mem Usage: 29,330,520 KB
2019-11-25 09:29:16 INFO  (et: 2298.752419 s) (func: main)  triangle_counting_main.cpp:137: There are 164015236714 triangles in the input graph.
There are 164015236714 triangles in the input graph.
```

## 程序代码模块说明

### 核心并行预处理和Triangle-Counting代码

文件 | 功能
--- | ---
[pre_processing.h](pre_processing.h) | `EdgeListHistogram`和`Reorder` (Graph Reordering)
[pre_processing_dodg.h](pre_processing_dodg.h) | 预处理Edge List到Degree Oriented Directed Graph(DODG), degree-descending
[triangle_counting.h](triangle_counting.h) | Hybrid (1) Sparse-Bitmap against Bitmap For Dense Range and (2) Vectorized Merge on Sorted Arrays
[triangle_counting_main.cpp](triangle_counting_main.cpp) | Main Logic

### 核心并行相关基础代码

文件/目录 | 功能
--- | ---
[ips4o](ips4o) | header-only in-place parallel sort (3rd-party dependency), 用于sort edge list 
[util/local_buffer.h](util/local_buffer.h), [util/primitives.h](util/primitives.h) | 自研的OpenMP-based parallel primitives, 在graph pre-processing中大量使用
[util/libpopcnt.h](util/libpopcnt.h) | 向量化的pop the counts of 1-bits in the array elements (3rd-party dependency)
[util/search_util.h](util/search_util.h), [util/set_inter_cnt_utils.h](util/set_inter_cnt_utils.h) | 自研的向量化的set-intersection方法
[util/graph.h](util/graph.h) | Compressed Sparse Row (CSR) graph的表示

### 其他Util文件

文件 | 功能
--- | ---
[util/popl.h](util/popl.h) | 命令行解析
[util/log.h](util/log.h), [util/log.cpp](util/log.cpp)  | 调试日志输出
[util/pretty_print.h](util/pretty_print.h) | 输出stl到流
[util/timer.h](util/timer.h) | 计时
[util/util.h](util/util.h) | 采集内存信息

## 程序代码编译说明

* 需要支持c++14标准的`g++`, 链接了pthread库, 使用了OpenMP, micro-architecture设为了`native` (本地编译), `-O3`优化

```bash
make
```

## 程序代码运行使用说明

```bash
./tc -f input_file_path
```