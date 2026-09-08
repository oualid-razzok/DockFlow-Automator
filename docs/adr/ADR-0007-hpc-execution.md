# ADR 0007: HPC execution model (slurm templates + resumable chunking)

- **Status**: accepted (work-order item 22)
- **Date**: 2026-09-08

## Context

Real virtual screening runs on clusters: thousands of ligands, hours
of CPU.  DockFlow is a single-node Python application; building a full
workflow manager (dependencies, retries, scheduling policy) would
duplicate slurm/nextflow and cannot be tested without a scheduler.

## Decision

1. DockFlow stays single-node; HPC support means **generating correct
   job files, not managing a queue**.  `scripts/generate_sbatch.py`
   emits a slurm array-job template (one array task = one ligand chunk
   = one `dockflow dock --resume` call) with CPU/memory/exhaustiveness
   placeholders filled from the project config.
2. The load-bearing enabler is **checkpoint/resume**: every run
   directory carries `progress.json` (per-ligand completion), so an
   interrupted array task restarts exactly where it died.  This is
   tested in CI; the slurm template itself is linted, not executed.
3. `dockflow info --check-hpc` reports whether the current environment
   looks like a batch node (env vars `SLURM_JOB_ID` etc.) and warns
   about `cpu x parallel` oversubscription inside a job.
4. Multi-node MPI-style docking is REJECTED: Vina's parallelism is
   shared-memory only; pretending otherwise would be dishonest.

## Consequences

- Cluster users bring their own scheduler knowledge; DockFlow's
  contract is "resumable + honest about cores".
- The template targets slurm (the dominant academic scheduler); other
  schedulers get the same chunked-resume pattern documented in
  USER_GUIDE.
