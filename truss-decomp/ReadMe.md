## Parallel Truss Decomposition

### Build

```zsh
cmake ~/workspace/yche/git-repos/OutOfCoreSCAN/truss-decomp 
```

### Files (Important)

Files | Comment
--- | ---
[util](util) | graph, log, stat, intersection, md5, print, timer, util, serialization
[extern_variables.cpp](extern_variables.cpp) | extern variables
[local_buffer.h](local_buffer.h) | cache-aware local buffer
[iter_helper.h](iter_helper.h), [iter_helper.cpp](iter_helper.cpp) | level-iteration logics
[parallel_all_edge_cnc.h](parallel_all_edge_cnc.h) | triangle-counting utils
[pkt_support_update_utils.h](pkt_support_update_utils.h) | support-update utils
[radix_hash_map.h](radix_hash_map.h) | radix-partitioning-based map
[set_utils.h](set_utils.h) | galloping-based set-intersection

### File Organizations

Folder | Commment
--- | ---
[reordering_utils](reordering_utils) | reordering utils
[playground](playground) | test some language features
[cmake](cmake) | `findxxx` cmake files