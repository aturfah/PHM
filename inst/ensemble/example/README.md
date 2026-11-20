# example/

Worked examples, each in its own folder with a ready-made `config.R` and a
step-by-step README.

| Folder | Dataset | Shows |
|---|---|---|
| [`double_ring/`](double_ring/) | The package's `double_ring` data: 5000 points on two concentric rings (2D, simulated) | One analysis run three ways: in serial, on one machine with 5 replicates at a time, and on a Slurm cluster. |

Run an example from a copy of its folder, not in place, so results stay out of
the package or repo. Each example's README shows how.
