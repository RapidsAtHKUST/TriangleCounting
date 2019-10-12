```
==101615== 23,712 bytes in 39 blocks are possibly lost in loss record 7 of 14
==101615==    at 0x4C2B955: calloc (vg_replace_malloc.c:711)
==101615==    by 0x40126B4: _dl_allocate_tls (in /usr/lib64/ld-2.17.so)
==101615==    by 0x58837AB: pthread_create@@GLIBC_2.2.5 (in /usr/lib64/libpthread-2.17.so)
==101615==    by 0x545599F: ??? (in /usr/lib64/libgomp.so.1.0.0)
==101615==    by 0x41230A: int RemoveDuplicates<int, unsigned long>(std::pair<int, int>*&, unsigned long&, std::pair<int, int>*&) (pre_processing.h:17)
==101615==    by 0x4033F1: main (triangle_counting_main.cpp:48)
==101615== 
==101615== 19,390,288 bytes in 1 blocks are possibly lost in loss record 10 of 14
==101615==    at 0x4C29BC3: malloc (vg_replace_malloc.c:299)
==101615==    by 0x40F847: void ConvertEdgeListToCSR<int>(unsigned int, std::pair<int, int>*, unsigned int, unsigned int*&, unsigned int*&, int*&, int) (pre_processing.h:80)
==101615==    by 0x4034FB: main (triangle_counting_main.cpp:64)
==101615== 
==101615== 19,390,288 bytes in 1 blocks are possibly lost in loss record 11 of 14
==101615==    at 0x4C29BC3: malloc (vg_replace_malloc.c:299)
==101615==    by 0x40F853: void ConvertEdgeListToCSR<int>(unsigned int, std::pair<int, int>*, unsigned int, unsigned int*&, unsigned int*&, int*&, int) (pre_processing.h:81)
==101615==    by 0x4034FB: main (triangle_counting_main.cpp:64)
==101615== 
==101615== 275,975,092 bytes in 1 blocks are possibly lost in loss record 12 of 14
==101615==    at 0x4C29BC3: malloc (vg_replace_malloc.c:299)
==101615==    by 0x41238B: int RemoveDuplicates<int, unsigned long>(std::pair<int, int>*&, unsigned long&, std::pair<int, int>*&) (pre_processing.h:53)
==101615==    by 0x4033F1: main (triangle_counting_main.cpp:48)
==101615== 
==101615== 551,950,184 bytes in 1 blocks are possibly lost in loss record 13 of 14
==101615==    at 0x4C29BC3: malloc (vg_replace_malloc.c:299)
==101615==    by 0x4033D2: main (triangle_counting_main.cpp:46)
==101615== 
==101615== 551,950,184 bytes in 1 blocks are possibly lost in loss record 14 of 14
==101615==    at 0x4C29BC3: malloc (vg_replace_malloc.c:299)
==101615==    by 0x41239D: int RemoveDuplicates<int, unsigned long>(std::pair<int, int>*&, unsigned long&, std::pair<int, int>*&) (pre_processing.h:54)
==101615==    by 0x4033F1: main (triangle_counting_main.cpp:48)
==101615== 
==101615== LEAK SUMMARY:
==101615==    definitely lost: 0 bytes in 0 blocks
==101615==    indirectly lost: 0 bytes in 0 blocks
==101615==      possibly lost: 1,418,679,748 bytes in 44 blocks
==101615==    still reachable: 331,992 bytes in 160 blocks
==101615==         suppressed: 0 bytes in 0 blocks
==101615== Reachable blocks (those to which a pointer was found) are not shown.
==101615== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==101615== 
==101615== ERROR SUMMARY: 6 errors from 6 contexts (suppressed: 0 from 0)
==101615== ERROR SUMMARY: 6 errors from 6 contexts (suppressed: 0 from 0)

```