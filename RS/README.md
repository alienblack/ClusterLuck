# RS heuristic submit

PACE heuristic solver for rooted agreement forests.

Build the static Linux submit binary with:

```sh
make build/rs_submit
cp build/rs_submit main_submit
```

The repository vendors the Whidden rSPR code used by
`src/whidden/whidden_bridge.cpp` under `third_party/whidden`.
