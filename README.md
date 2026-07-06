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
