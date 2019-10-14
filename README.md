# TriangleCounting

Triangle Counting 

## GPU23 Statistics

* LJ

```
2019-10-14 09:52:58 INFO  (ts: 1571017978.369168 s, et: 0.000458 s)  triangle_counting_main.cpp:35: File size: 551950184
2019-10-14 09:52:58 INFO  (ts: 1571017978.369430 s, et: 0.000719 s)  triangle_counting_main.cpp:36: #of Edges: 68993773
2019-10-14 09:52:58 INFO  (ts: 1571017978.374253 s, et: 0.005542 s)  triangle_counting_main.cpp:50: Populate Mem Time: 0.005448280s
2019-10-14 09:52:58 INFO  (ts: 1571017978.410091 s, et: 0.041381 s)  triangle_counting_main.cpp:59: Load File Time: 0.041286879s
2019-10-14 09:52:58 INFO  (ts: 1571017978.410395 s, et: 0.041685 s)  triangle_counting_main.cpp:62: [main] Mem: 540624 KB
2019-10-14 09:52:58 INFO  (ts: 1571017978.461428 s, et: 0.092718 s)  primitives.h:123: [Histogram]: Local Comp Time: 0.050986211s
2019-10-14 09:52:58 INFO  (ts: 1571017978.482728 s, et: 0.114018 s)  primitives.h:165: [BucketSort]: Histogram, Time: 0.072286072s
2019-10-14 09:52:58 INFO  (ts: 1571017978.498967 s, et: 0.130257 s)  primitives.h:173: [BucketSort]: Before Scatter, Time: 0.088524448s
2019-10-14 09:52:58 INFO  (ts: 1571017978.724555 s, et: 0.355845 s)  primitives.h:185: [BucketSort]: Before Sort, Time: 0.314113486s
2019-10-14 09:52:58 INFO  (ts: 1571017978.821269 s, et: 0.452558 s)  pre_processing.h:48: Finish Sort, 0.410826820s
2019-10-14 09:52:58 INFO  (ts: 1571017978.889418 s, et: 0.520708 s)  pre_processing.h:61: New # of edges: 42851237, Elapsed: 0.478976451s
2019-10-14 09:52:58 INFO  (ts: 1571017978.889594 s, et: 0.520883 s)  triangle_counting_main.cpp:64: [main] Mem: 1174080 KB
2019-10-14 09:52:58 INFO  (ts: 1571017978.889626 s, et: 0.520916 s)  triangle_counting_main.cpp:67: Load Edge List Time: 0.520823817 s
2019-10-14 09:52:58 INFO  (ts: 1571017978.889652 s, et: 0.520942 s)  triangle_counting_main.cpp:74: Undirected Graph G = (|V|, |E|): 4847571, 42851237
2019-10-14 09:52:58 INFO  (ts: 1571017978.898393 s, et: 0.529683 s)  pre_processing.h:82: [ConvertEdgeListToCSR]: InitTime: 0.008714014 s
2019-10-14 09:52:58 INFO  (ts: 1571017978.927682 s, et: 0.558972 s)  pre_processing.h:101: [ConvertEdgeListToCSR]: Histogram Time: 0.038005929 s
2019-10-14 09:52:58 INFO  (ts: 1571017978.956137 s, et: 0.587427 s)  pre_processing.h:110: [ConvertEdgeListToCSR]: Histogram Time: 0.066458966 s
2019-10-14 09:52:58 INFO  (ts: 1571017978.972863 s, et: 0.604153 s)  pre_processing.h:125: [ConvertEdgeListToCSR]: PrefixSum Time: 0.083185233 s
2019-10-14 09:52:59 INFO  (ts: 1571017979.212722 s, et: 0.844012 s)  pre_processing.h:140: [ConvertEdgeListToCSR]: Total Conversion Time: 0.323045887 s
2019-10-14 09:52:59 INFO  (ts: 1571017979.222553 s, et: 0.853843 s)  primitives.h:204: [BucketSortSmallBuckets]: Mem Size Buckets: 335104320, Bucket#: 20334
2019-10-14 09:52:59 INFO  (ts: 1571017979.266346 s, et: 0.897636 s)  tc_utils.h:212: Deg-descending time:  0.053531586 s
2019-10-14 09:52:59 INFO  (ts: 1571017979.297355 s, et: 0.928645 s)  tc_utils.h:138: [Reorder]: Finish PrefixSum Time: 0.030907616 s
2019-10-14 09:52:59 INFO  (ts: 1571017979.419592 s, et: 1.050882 s)  tc_utils.h:157: [Reorder]: Finish Reorder Time: 0.122 s
2019-10-14 09:52:59 INFO  (ts: 1571017979.420641 s, et: 1.051931 s)  triangle_counting_main.cpp:85: [main] Mem: 1351636 KB
2019-10-14 09:52:59 INFO  (ts: 1571017979.422455 s, et: 1.053745 s)  triangle_counting_main.cpp:87: [main] Mem: 812620 KB
2019-10-14 09:52:59 INFO  (ts: 1571017979.445952 s, et: 1.077242 s)  triangle_counting.h:39: finish init row_ptrs_end, max d: 685, time: 0.023464221s
2019-10-14 09:52:59 INFO  (ts: 1571017979.458670 s, et: 1.089960 s)  triangle_counting.h:59: Finish Indexing: 0.036187656s
2019-10-14 09:52:59 INFO  (ts: 1571017979.877863 s, et: 1.509153 s)  triangle_counting.h:123: Forward cost: 0.455 s, Mem Usage: 859872 KB
2019-10-14 09:52:59 INFO  (ts: 1571017979.877933 s, et: 1.509223 s)  triangle_counting.h:124: Triangle Cnt: 285730264
2019-10-14 09:52:59 INFO  (ts: 1571017979.893038 s, et: 1.524327 s)  triangle_counting_main.cpp:102: There are 285730264 triangles in the input graph.
There are 285730264 triangles in the input graph.
```

* scale-free 24

```
2019-10-14 09:51:50 INFO  (ts: 1571017910.239402 s, et: 0.000466 s)  triangle_counting_main.cpp:35: File size: 2147483648
2019-10-14 09:51:50 INFO  (ts: 1571017910.239644 s, et: 0.000708 s)  triangle_counting_main.cpp:36: #of Edges: 268435456
2019-10-14 09:51:50 INFO  (ts: 1571017910.244978 s, et: 0.006042 s)  triangle_counting_main.cpp:50: Populate Mem Time: 0.005950599s
2019-10-14 09:51:50 INFO  (ts: 1571017910.369982 s, et: 0.131045 s)  triangle_counting_main.cpp:59: Load File Time: 0.130953605s
2019-10-14 09:51:50 INFO  (ts: 1571017910.370192 s, et: 0.131256 s)  triangle_counting_main.cpp:62: [main] Mem: 2099024 KB
2019-10-14 09:51:50 INFO  (ts: 1571017910.800830 s, et: 0.561894 s)  primitives.h:123: [Histogram]: Local Comp Time: 0.430607920s
2019-10-14 09:51:50 INFO  (ts: 1571017910.877869 s, et: 0.638933 s)  primitives.h:165: [BucketSort]: Histogram, Time: 0.507646563s
2019-10-14 09:51:50 INFO  (ts: 1571017910.898543 s, et: 0.659608 s)  primitives.h:173: [BucketSort]: Before Scatter, Time: 0.528319571s
2019-10-14 09:51:52 INFO  (ts: 1571017912.152406 s, et: 1.913470 s)  primitives.h:185: [BucketSort]: Before Sort, Time: 1.782184679s
2019-10-14 09:51:52 INFO  (ts: 1571017912.658172 s, et: 2.419236 s)  pre_processing.h:48: Finish Sort, 2.287949928s
2019-10-14 09:51:52 INFO  (ts: 1571017912.865714 s, et: 2.626778 s)  pre_processing.h:61: New # of edges: 260379850, Elapsed: 2.495492544s
2019-10-14 09:51:52 INFO  (ts: 1571017912.865863 s, et: 2.626927 s)  triangle_counting_main.cpp:64: [main] Mem: 4770384 KB
2019-10-14 09:51:52 INFO  (ts: 1571017912.865886 s, et: 2.626950 s)  triangle_counting_main.cpp:67: Load Edge List Time: 2.626859096 s
2019-10-14 09:51:52 INFO  (ts: 1571017912.865900 s, et: 2.626964 s)  triangle_counting_main.cpp:74: Undirected Graph G = (|V|, |E|): 16777216, 260379850
2019-10-14 09:51:52 INFO  (ts: 1571017912.878070 s, et: 2.639134 s)  pre_processing.h:82: [ConvertEdgeListToCSR]: InitTime: 0.012150170 s
2019-10-14 09:51:53 INFO  (ts: 1571017913.235555 s, et: 2.996619 s)  pre_processing.h:101: [ConvertEdgeListToCSR]: Histogram Time: 0.369639190 s
2019-10-14 09:51:53 INFO  (ts: 1571017913.316113 s, et: 3.077177 s)  pre_processing.h:110: [ConvertEdgeListToCSR]: Histogram Time: 0.450197406 s
2019-10-14 09:51:53 INFO  (ts: 1571017913.336516 s, et: 3.097580 s)  pre_processing.h:125: [ConvertEdgeListToCSR]: PrefixSum Time: 0.470598717 s
2019-10-14 09:51:54 INFO  (ts: 1571017914.867881 s, et: 4.628945 s)  pre_processing.h:140: [ConvertEdgeListToCSR]: Total Conversion Time: 2.001964456 s
2019-10-14 09:51:54 INFO  (ts: 1571017914.891653 s, et: 4.652717 s)  primitives.h:204: [BucketSortSmallBuckets]: Mem Size Buckets: 1039362560, Bucket#: 406001
2019-10-14 09:51:55 INFO  (ts: 1571017915.027087 s, et: 4.788151 s)  tc_utils.h:212: Deg-descending time:  0.159038445 s
2019-10-14 09:51:55 INFO  (ts: 1571017915.132506 s, et: 4.893570 s)  tc_utils.h:138: [Reorder]: Finish PrefixSum Time: 0.105344981 s
2019-10-14 09:51:56 INFO  (ts: 1571017916.139366 s, et: 5.900430 s)  tc_utils.h:157: [Reorder]: Finish Reorder Time: 1.007 s
2019-10-14 09:51:56 INFO  (ts: 1571017916.140623 s, et: 5.901687 s)  triangle_counting_main.cpp:85: [main] Mem: 5114264 KB
2019-10-14 09:51:56 INFO  (ts: 1571017916.147351 s, et: 5.908415 s)  triangle_counting_main.cpp:87: [main] Mem: 3017108 KB
2019-10-14 09:51:56 INFO  (ts: 1571017916.246365 s, et: 6.007429 s)  triangle_counting.h:39: finish init row_ptrs_end, max d: 1773, time: 0.098964865s
2019-10-14 09:51:56 INFO  (ts: 1571017916.287457 s, et: 6.048521 s)  triangle_counting.h:59: Finish Indexing: 0.140054370s
2019-10-14 09:52:06 INFO  (ts: 1571017926.241663 s, et: 16.002727 s)  triangle_counting.h:123: Forward cost: 10.094 s, Mem Usage: 3017120 KB
2019-10-14 09:52:06 INFO  (ts: 1571017926.241742 s, et: 16.002806 s)  triangle_counting.h:124: Triangle Cnt: 10286638314
2019-10-14 09:52:06 INFO  (ts: 1571017926.264234 s, et: 16.025298 s)  triangle_counting_main.cpp:102: There are 10286638314 triangles in the input graph.
There are 10286638314 triangles in the input graph.
```