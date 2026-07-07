# ClusterLuck

Aggressive PACE 2026 exact-track solver package.

Build:

```bash
make
```

Run:

```bash
./solve < instance.nw
```

The entrypoint intentionally TLEs for two-tree instances with more than 500 leaves:

```bash
#p 2 n, n > 500
```

Otherwise it dispatches to the two-tree RS solver for `t=2`, and to the empirical/stable many-tree core for `t>2`.

## Update note for branch `new`

This branch contains the corrected submit-ready packaging for ClusterLuck. The previous pushed folder was incomplete and did not include the intended combined submit wrapper/build target and final submit binary.

This branch includes:

- `solver_submit`: Linux static submit executable.
- `combined_wrapper.cpp`: wrapper combining the many-tree and two-tree routes into one executable.
- Makefile target for building `solver_submit`.
- Bundled many-tree source.
- Bundled h56/HiGHS two-tree route source.

ClusterLuck is the aggressive/practical version.
