# TriangleCounting

Triangle Counting 

## GPU23 Statistics

* LJ

```
2019-10-12 19:45:04 INFO  (ts: 1570880704.109176 s, et: 0.000479 s)  triangle_counting_main.cpp:38: File size: 551950184
2019-10-12 19:45:04 INFO  (ts: 1570880704.109349 s, et: 0.000652 s)  triangle_counting_main.cpp:39: #of Edges: 68993773
2019-10-12 19:45:04 INFO  (ts: 1570880704.397944 s, et: 0.289246 s)  triangle_counting_main.cpp:45: Load File Time: 0.289153531s
2019-10-12 19:45:04 INFO  (ts: 1570880704.542928 s, et: 0.434231 s)  primitives.h:142: Before Scatter, Time: 0.144884439s
2019-10-12 19:45:04 INFO  (ts: 1570880704.750749 s, et: 0.642052 s)  primitives.h:154: Before Sort, Time: 0.352706202s
2019-10-12 19:45:05 INFO  (ts: 1570880705.077249 s, et: 0.968551 s)  pre_processing.h:50: Finish Sort, 0.679206171s
2019-10-12 19:45:05 INFO  (ts: 1570880705.141157 s, et: 1.032459 s)  pre_processing.h:71: New # of edges: 42851237, Elapsed: 0.743113643s
2019-10-12 19:45:05 INFO  (ts: 1570880705.141384 s, et: 1.032686 s)  triangle_counting_main.cpp:49: Mem: 1703264 KB
2019-10-12 19:45:05 INFO  (ts: 1570880705.181305 s, et: 1.072608 s)  triangle_counting_main.cpp:51: Mem: 1164460 KB
2019-10-12 19:45:05 INFO  (ts: 1570880705.181348 s, et: 1.072650 s)  triangle_counting_main.cpp:54: load edge list bin time: 1.072558114 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.181367 s, et: 1.072669 s)  triangle_counting_main.cpp:61: Undirected Graph G = (|V|, |E|): 4847571, 42851237
2019-10-12 19:45:05 INFO  (ts: 1570880705.295944 s, et: 1.187247 s)  pre_processing.h:121: before csr transform time: 0.114558055 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.533120 s, et: 1.424423 s)  pre_processing.h:147: edge list to csr time: 0.351734480 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.542910 s, et: 1.434213 s)  primitives.h:172: Mem Size Buckets: 335104320, Bucket#: 20334
2019-10-12 19:45:05 INFO  (ts: 1570880705.584995 s, et: 1.476297 s)  tc_utils.h:212: Deg-descending time:  0.051814155 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.612565 s, et: 1.503867 s)  tc_utils.h:137: init ordering structures time: 0.027498378 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.749218 s, et: 1.640521 s)  tc_utils.h:153: parallel transform and sort: 0.137 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.752506 s, et: 1.643809 s)  tc_utils.h:157: parallel transform after copying: 0.140 s
2019-10-12 19:45:05 INFO  (ts: 1570880705.753348 s, et: 1.644650 s)  triangle_counting_main.cpp:72: Mem: 1620816 KB
2019-10-12 19:45:05 INFO  (ts: 1570880705.754493 s, et: 1.645795 s)  triangle_counting_main.cpp:74: Mem: 1286040 KB
2019-10-12 19:45:05 INFO  (ts: 1570880705.777662 s, et: 1.668965 s)  triangle_counting.h:39: finish init row_ptrs_end, max d: 685, time: 0.023147091s
2019-10-12 19:45:05 INFO  (ts: 1570880705.787509 s, et: 1.678811 s)  triangle_counting.h:59: Finish Indexing: 0.032993409s
2019-10-12 19:45:06 INFO  (ts: 1570880706.178646 s, et: 2.069949 s)  triangle_counting.h:123: Forward cost: 0.424 s, Mem Usage: 1346980 KB
2019-10-12 19:45:06 INFO  (ts: 1570880706.178725 s, et: 2.070028 s)  triangle_counting.h:124: Triangle Cnt: 285730264
2019-10-12 19:45:06 INFO  (ts: 1570880706.194810 s, et: 2.086112 s)  triangle_counting_main.cpp:89: There are 285730264 triangles in the input graph.
There are 285730264 triangles in the input graph.
```

* scale-free 24

```
2019-10-12 19:45:17 INFO  (ts: 1570880717.297371 s, et: 0.000406 s)  triangle_counting_main.cpp:38: File size: 2147483648
2019-10-12 19:45:17 INFO  (ts: 1570880717.297558 s, et: 0.000593 s)  triangle_counting_main.cpp:39: #of Edges: 268435456
2019-10-12 19:45:18 INFO  (ts: 1570880718.338713 s, et: 1.041748 s)  triangle_counting_main.cpp:45: Load File Time: 1.041671807s
2019-10-12 19:45:19 INFO  (ts: 1570880719.097543 s, et: 1.800578 s)  primitives.h:142: Before Scatter, Time: 0.758751896s
2019-10-12 19:45:20 INFO  (ts: 1570880720.300542 s, et: 3.003578 s)  primitives.h:154: Before Sort, Time: 1.961752646s
2019-10-12 19:45:21 INFO  (ts: 1570880721.359645 s, et: 4.062680 s)  pre_processing.h:50: Finish Sort, 3.020855848s
2019-10-12 19:45:21 INFO  (ts: 1570880721.567403 s, et: 4.270438 s)  pre_processing.h:71: New # of edges: 260379850, Elapsed: 3.228613263s
2019-10-12 19:45:21 INFO  (ts: 1570880721.567565 s, et: 4.270600 s)  triangle_counting_main.cpp:49: Mem: 7280608 KB
2019-10-12 19:45:21 INFO  (ts: 1570880721.693603 s, et: 4.396638 s)  triangle_counting_main.cpp:51: Mem: 5183460 KB
2019-10-12 19:45:21 INFO  (ts: 1570880721.693655 s, et: 4.396690 s)  triangle_counting_main.cpp:54: load edge list bin time: 4.396614722 s
2019-10-12 19:45:21 INFO  (ts: 1570880721.693671 s, et: 4.396705 s)  triangle_counting_main.cpp:61: Undirected Graph G = (|V|, |E|): 16777216, 260379850
2019-10-12 19:45:22 INFO  (ts: 1570880722.604084 s, et: 5.307119 s)  pre_processing.h:121: before csr transform time: 0.910398409 s
2019-10-12 19:45:24 INFO  (ts: 1570880724.358538 s, et: 7.061573 s)  pre_processing.h:147: edge list to csr time: 2.664851507 s
2019-10-12 19:45:24 INFO  (ts: 1570880724.393935 s, et: 7.096970 s)  primitives.h:172: Mem Size Buckets: 1039362560, Bucket#: 406001
2019-10-12 19:45:24 INFO  (ts: 1570880724.528179 s, et: 7.231214 s)  tc_utils.h:212: Deg-descending time:  0.169494234 s
2019-10-12 19:45:24 INFO  (ts: 1570880724.633725 s, et: 7.336760 s)  tc_utils.h:137: init ordering structures time: 0.105483913 s
2019-10-12 19:45:25 INFO  (ts: 1570880725.699383 s, et: 8.402418 s)  tc_utils.h:153: parallel transform and sort: 1.066 s
2019-10-12 19:45:25 INFO  (ts: 1570880725.708659 s, et: 8.411694 s)  tc_utils.h:157: parallel transform after copying: 1.075 s
2019-10-12 19:45:25 INFO  (ts: 1570880725.709359 s, et: 8.412394 s)  triangle_counting_main.cpp:72: Mem: 7911668 KB
2019-10-12 19:45:25 INFO  (ts: 1570880725.715637 s, et: 8.418672 s)  triangle_counting_main.cpp:74: Mem: 5877448 KB
2019-10-12 19:45:25 INFO  (ts: 1570880725.816571 s, et: 8.519606 s)  triangle_counting.h:39: finish init row_ptrs_end, max d: 1773, time: 0.100895222s
2019-10-12 19:45:25 INFO  (ts: 1570880725.853292 s, et: 8.556328 s)  triangle_counting.h:59: Finish Indexing: 0.137616642s
2019-10-12 19:45:35 INFO  (ts: 1570880735.742032 s, et: 18.445067 s)  triangle_counting.h:123: Forward cost: 10.026 s, Mem Usage: 5883072 KB
2019-10-12 19:45:35 INFO  (ts: 1570880735.742115 s, et: 18.445150 s)  triangle_counting.h:124: Triangle Cnt: 10286638314
2019-10-12 19:45:35 INFO  (ts: 1570880735.758235 s, et: 18.461270 s)  triangle_counting_main.cpp:89: There are 10286638314 triangles in the input graph.
There are 10286638314 triangles in the input graph.
```