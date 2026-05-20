
Beforehand, make sure  you had exported the required environment variables.

```sh
export CPLEX_STUDIO_DIR="/opt/ibm/ILOG/CPLEX_Studio201/"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/ibm/ILOG/CPLEX_Studio201/opl/bin/x86-64_linux/
export JULIA_COPY_STACKS=1
```

## Docker IP with CPLEX

Self-contained **builder + runner** for this repo (no dependency on a parent COPAlgorithms tree). From **this directory** (repository root):

```bash
export CPLEX_HOST_PATH="/path/to/CPLEX_Studio2211"   # host Studio install; optional if already at default path
docker compose -f docker-compose.ip.yml build
docker compose -f docker-compose.ip.yml run --rm cbrp-original-ip-builder
docker compose -f docker-compose.ip.yml run --rm cbrp-original-ip-runner
```

| Variable | Default | Role |
|----------|---------|------|
| `IP_ARTIFACTS_DIR` | `./build/ip-artifacts` | Host folder mounted at `/artifacts` (staging writes `cbrp-original-ip/CBRP_original` and `.julia` here). |
| `CPLEX_HOST_PATH` | `/opt/ibm/ILOG/CPLEX_Studio2211` | Host CPLEX Studio root for the bind mount. |
| `CPLEX_ROOT_DIR` | `/opt/ibm/ILOG/CPLEX_Studio2211` | Mount point inside the container (must match `CPLEX_ROOT_DIR` passed to the staging script). |

Images: `cbrp-original-ip-builder:latest`, `cbrp-original-ip-runner:latest`. Dockerfiles: [`docker/Dockerfile.ip-builder`](docker/Dockerfile.ip-builder), [`docker/Dockerfile.ip-runner`](docker/Dockerfile.ip-runner). Staging: [`scripts/stage_ip_artifacts.sh`](scripts/stage_ip_artifacts.sh).

Inside the runner, `WORKDIR` is the staged Julia project. Example solves (same as below):

```bash
julia --threads=1 --project=. src/run.jl data/campinas-sparse/1.sbrp --out solutions/complete_smoke --intersection-cuts --ip
```

**Carlos sparse digraph + Path-CBRP MILP:** use `--instance-type carlos`, `--ip`, `--no-cbrp-metric-closure`, and `--path-cbrp-mip`. That keeps the street digraph (no Floyd–Warshall metric closure) and runs the arc-indexed Path-CBRP model with a global travel + service time bound and compact arc MTZ. Do not combine `--path-cbrp-mip` with `--brkga`. CPLEX time cap for the Path MILP: `--time-limit` (seconds; `0` defaults to 3600).

**Path-CBRP subtour elimination:** `--subcycle-separation` (`first`|`best`|`all`|`none`, default `all`) selects which violated SECs to add per separation round. `--subcycle-separation-engine` (`root`|`callback`, default `root`) chooses the host:

| Engine | Behavior |
|--------|----------|
| `root` | Pre-MIP LP cutting-plane loop (max-flow on relaxed `x`,`y`), then binary MIP solve |
| `callback` | Single CPLEX generic callback: Path SECs as **user cuts** at each `RELAXATION` (fractional LP) and **lazy constraints** at integer `CANDIDATE` points (requires `--path-cbrp-mip`; can be much slower on large Carlos instances) |

**Optional no MTZ:** `--no-path-cbrp-mtz` skips compact arc MTZ constraints (keeps `w` and depot bound). Requires `--path-cbrp-mip` and `--subcycle-separation-engine callback` with `--subcycle-separation` not `none` (otherwise the run errors: connectivity must come from SEC separation).

Inequalities: \(\sum_{a \in \delta^{+}(S)} x_a \ge y_{b,i} + y_{b',j} - 1\). Logs include `maxFlowCuts` (total), `maxFlowUserCuts`, `maxFlowLazyCuts`, `maxFlowCutsTime`, `subcycleSeparationEngine`, and `pathCbrpMtzEnabled`. Use `none` to skip SECs (faster, MTZ-only when MTZ is enabled).

**Callback debug (stdout):** with `engine=callback`, each batch of submitted cuts prints as `[PathSEC] RELAXATION user: +N cuts (cum. user=…, lazy=…)` before CPLEX’s next `User` / `UserPurge2` line. Disable with `PATH_CBRP_SEC_CALLBACK_LOG=0`.

Example:

```sh
julia --threads=1 --project=. src/run.jl data/carlos/notified-alto-santo/notified-alto-santo-1000-2021.txt --instance-type carlos --ip --no-cbrp-metric-closure --path-cbrp-mip --out solutions/pcbrp_smoke --time-limit 60
```

**Path-CBRP MIP warm start** (`--path-cbrp-warm-sol PATH`, `--path-cbrp-warm-sol-format`): optional hints before the Path MILP solve. Formats:

| Format | Flag | File |
|--------|------|------|
| `brkga-sol` | `--path-cbrp-warm-sol-format brkga-sol` | `writeSolution` output (`PREFIX_brkga.sol`); expanded via metric closure |
| `path-cbrp-mtz` | `--path-cbrp-warm-sol-format path-cbrp-mtz` | Article reference: `X: i j` arcs, `Y: node block_id` (`y[node][block]=1`) |
| `auto` (default) | omit or `auto` | Detects `X:` lines → `path-cbrp-mtz`, else `brkga-sol` |

**BRKGA `.sol` (two commands):**

```sh
julia --threads=1 --project=. src/run.jl INSTANCE.txt --instance-type carlos --brkga --out solutions/my_prefix --vehicle-time-limit 120
julia --threads=1 --project=. src/run.jl INSTANCE.txt --instance-type carlos --ip --no-cbrp-metric-closure --path-cbrp-mip \
  --path-cbrp-warm-sol solutions/my_prefix_brkga.sol --path-cbrp-warm-sol-format brkga-sol \
  --out solutions/pcbrp_warm --vehicle-time-limit 120 --time-limit 3600
```

**Article `path-cbrp-mtz` reference (one command):**

```sh
julia --threads=1 --project=. src/run.jl data/carlos/notified-alto-santo/notified-alto-santo-1000-2016.txt \
  --instance-type carlos --ip --no-cbrp-metric-closure --path-cbrp-mip \
  --path-cbrp-warm-sol solutions/results-cbrp-article/path-cbrp-mtz/notified-alto-santo-1000-2016.txt \
  --path-cbrp-warm-sol-format path-cbrp-mtz \
  --out solutions/pcbrp_warm --vehicle-time-limit 120 --time-limit 3600
```

Use the **same** `INSTANCE.txt`, `--vehicle-time-limit`, and `--drop-zero-profit-blocks` (if any) as when the warm file was produced. For **`path-cbrp-mtz`**, every `X:` arc must exist on the sparse digraph produced by the Carlos reader (depot dummy arcs included); otherwise warm start is skipped with a warning. Path `--time-limit` caps CPLEX only; BRKGA uses `conf/config.conf` plus its own `--time-limit` on the BRKGA run. Logs include `warmStartUsed` when the warm pipeline succeeds (MIP hints or hard-fix below).

**Julia-only hard-fix warm start:** there is no CLI flag. When calling `runPathCbrpMipModel` from Julia (REPL, script, or tests), set `app["path_cbrp_fix_warm_start"] = true` (or the string `"true"`) so Path-CBRP adds explicit equality constraints `fix_warm_x[k]` and `fix_warm_y[k]` on the binary `x` and `y` variables (instead of MIP starts). CPLEX must then satisfy MTZ, time, and linking at those values—useful to tell apart an inconsistent warm pattern from model issues. Continuous `w` (MTZ potentials) are left free. On success, logs set `warmStartFixed` to `"true"` (otherwise `"false"`).

**Tests:** from this directory: `julia --project=. -e 'using Test; include("test/runtests.jl")'` (CPLEX Path-CBRP smoke may skip if CPLEX fails).

**Path-CBRP time-budget sweep (experiment):** repeatedly solves the Path-CBRP MILP on the same Carlos sparse instance while increasing `data.T` (minutes) from `--min-T` by `--time-step` until `num_serviced_blocks >= num_positive_profit_blocks` (blocks with total profit &gt; 0 from the Carlos reader; zero-profit blocks are never chosen in the max-profit MILP), CPLEX throws, or caps `--max-T` (default 400) / `--max-iterations` (default 100) apply. The CSV includes `num_positive_profit_blocks`; success terminal is `all_positive_profit_blocks`. One CSV row per solve (`--out-csv` required). CPLEX per-solve wall cap: `--time-limit` (seconds; `0` → 3600). This script does not use `src/run.jl` or `logs/log`.

```sh
julia --threads=1 --project=. experiments/path_cbrp_time_sweep.jl \
  data/carlos/notified-alto-santo/notified-alto-santo-1000-2021.txt \
  --out-csv /tmp/path_sweep.csv --min-T 10 --time-step 5 --max-T 400 --max-iterations 100 --time-limit 30
```

You can get the parameters list by typing:

```sh
julia --threads=auto --project=. src/run.jl

positional arguments:
  instance              Instance file path

optional arguments:
  --cluster-size-profits
                        true if you want to consider the blocks'
                        profits as the size of the block, and false
                        otherwise
  --unitary-profits     true if you want to consider the blocks'
                        profits as 1, and false otherwise
  --ip                  true if you want to run the I.P. model, and
                        false otherwise
  --brkga               true if you want to run the BRKGA, and false
                        otherwise
  --brkga-conf BRKGA-CONF
                        BRKGA config file directory (default:
                        "conf/config.conf")
  --vehicle-time-limit VEHICLE-TIME-LIMIT
                        Vehicle time limit in minutes (default: "120")
  --instance-type INSTANCE-TYPE
                        Instance type (matheus|carlos) (default:
                        "matheus")
  --nosolve             Not solve flag
  --out OUT             Path to write the solution found
  --batch BATCH         Batch file path
  --intersection-cuts   Intersection cuts for the complete model
  --subcycle-separation SUBCYCLE-SEPARATION
                        Subcycle separation strategy: first|best|all|none
                        (default: "all")
  --subcycle-separation-engine ENGINE
                        Path-CBRP SEC host: root|callback (default: "root")
  --y-integer           Fix the variable y, for the complete model,
                        when running the separation algorithm
  --z-integer           Fix the variable z, for the complete model,
                        when running the separation algorithm
  --w-integer           Fix the variable w, for the complete model,
                        when running the separation algorithm
  -h, --help            show this help message and exit
```

And to run the IP model for a single instance:

```sh
julia --threads=auto --project=. src/run.jl data/campinas-random/1.sbrp --out solutions/complete_campinas-random1 --intersection-cuts --ip

Version identifier: 20.1.0.0 | 2020-11-10 | 9bedb6d68
CPXPARAM_TimeLimit                               3600
Tried aggregator 1 time.
LP Presolve eliminated 32 rows and 1 columns.
Reduced LP has 1499 rows, 1027 columns, and 9937 nonzeros.
Presolve time = 0.00 sec. (3.12 ticks)
Initializing dual steep norms . . .

Iteration log . . .
Iteration:     1   Dual objective     =            18.000000
Version identifier: 20.1.0.0 | 2020-11-10 | 9bedb6d68
CPXPARAM_TimeLimit                               3600
Probing fixed 0 vars, tightened 30 bounds.
Probing time = 0.14 sec. (0.72 ticks)
Clique table members: 92.
MIP emphasis: balance optimality and feasibility.
MIP search method: dynamic search.
Parallel mode: deterministic, using up to 20 threads.
Root relaxation solution time = 0.02 sec. (1.41 ticks)

        Nodes                                         Cuts/
   Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap

      0     0       18.0000    16                     18.0000        0
*     0+    0                           18.0000       18.0000             0.00%
      0     0        cutoff             18.0000       18.0000        0    0.00%
Elapsed time = 0.32 sec. (13.98 ticks, tree = 0.01 MB, solutions = 1)

Root node processing (before b&c):
  Real time             =    0.32 sec. (14.04 ticks)
Parallel b&c, 20 threads:
  Real time             =    0.00 sec. (0.00 ticks)
  Sync time (average)   =    0.00 sec.
  Wait time (average)   =    0.00 sec.
                          ------------
Total (root+branch&cut) =    0.32 sec. (14.04 ticks)
Version identifier: 20.1.0.0 | 2020-11-10 | 9bedb6d68
CPXPARAM_TimeLimit                               3600
Warning:  No solution found from 1 MIP starts.
Retaining values of one MIP start for possible repair.
Found incumbent of value 18.000000 after 0.04 sec. (22.71 ticks)

        Nodes                                         Cuts/
   Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap

*     0+    1                           18.0000       18.0000             0.00%

Root node processing (before b&c):
  Real time             =    0.04 sec. (22.76 ticks)
Parallel b&c, 20 threads:
  Real time             =    0.00 sec. (0.00 ticks)
  Sync time (average)   =    0.00 sec.
  Wait time (average)   =    0.00 sec.
                          ------------
Total (root+branch&cut) =    0.04 sec. (22.76 ticks)
instance,|V|,|A|,|B|,T,model,initialLP,yLP,yLPTime,lazyCuts,cost,solverTime,relativeGAP,nodeCount,meters,tourMinutes,blocksMeters,numVisitedBlocks,intersectionCutsTime,intersectionCuts1,intersectionCuts2
1,31,930,5,120.0,IP,18.0,18.0,0.348148434,0,18.00,0.059470715,0.0,0,5027.0,14.673,1585.0,5,0.307596037,0,0
```

The last two lines are the summarized results.

For running the genetic algorithm:

```sh
julia --threads=auto --project=. src/run.jl data/campinas-random/1.sbrp --out solutions/brkga_campinas-random1 --brkga

------------------------------------------------------
> Experiment started at 2024-09-04T13:11:17.225
> Configuration: conf/config.conf
> Algorithm Parameters:

>  - population_size 2000
>  - elite_percentage 0.3
>  - mutants_percentage 0.15
>  - num_elite_parents 2
>  - total_parents 3
>  - bias_type LOGINVERSE
>  - num_independent_populations 3
>  - pr_number_pairs 0
>  - pr_minimum_distance 0.15
>  - pr_type PERMUTATION
>  - pr_selection BESTSOLUTION
>  - alpha_block_size 1.0
>  - pr_percentage 1.0
>  - exchange_interval 200
>  - num_exchange_indivuduals 2
>  - reset_interval 600
> Seed: 12345
> Stop rule: GENERATIONS
> Stop argument: 1000000
> Maximum time (s): 100.0
> Number of parallel threads for decoding: 20
------------------------------------------------------

[13:11:17.266] Generating initial tour...
Initial profit: 18.0

[13:11:17.481] Building BRKGA data...
New population size: 50

[13:11:17.539] Initializing BRKGA data...

[13:11:17.916] Warming up...

[13:11:18.332] Evolving...
* Iteration | Cost | CurrentTime
Performing path relink at 200...
- No improvement found | Elapsed time: 0.02
Performing path relink at 400...
- No improvement found | Elapsed time: 0.02
Performing path relink at 600...
- No improvement found | Elapsed time: 0.01
Performing path relink at 800...
- No improvement found | Elapsed time: 0.01
Performing path relink at 1000...
- No improvement found | Elapsed time: 0.02
Performing path relink at 1200...
- No improvement found | Elapsed time: 0.02
Performing path relink at 1400...
- No improvement found | Elapsed time: 0.02
Performing path relink at 1600...
- No improvement found | Elapsed time: 0.01
Performing path relink at 1800...
- No improvement found | Elapsed time: 0.02
[13:12:58.332] End of optimization

Total number of iterations: 1954
Last update iteration: 0
Total optimization time: 100.00
Last update time: 0.00
Large number of iterations between improvements: 0
Total path relink time: 0.15
Total path relink calls: 9
Number of homogenities: 0
Improvements in the elite set: 0
Best individual improvements: 0
Last index 5
BRKGA cost:18.0
Blocks cost:18.0
Calculated cost:18.0
instance,|V|,|A|,|B|,T,model,cost,solverTime,meters,tourMinutes,blocksMeters,numVisitedBlocks
1,31,930,5,120.0,BRKGA,18.0,100.00000500679016,3122.0,11.815500000000002,1585.0,5
```

For running several instances at once, you can use the `--batch` flag.

```sh
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random/intersection_cuts_relaxed_y/complete.batch > log_campinas_random_complete
```

In this case, the file `log_campinas_random_complete` contains all the logs, so, if you want to extract only the summaries, you can use the `grep` command.

```sh
grep -A1 "instance" log_campinas_random_complete

instance,|V|,|A|,|B|,T,model,initialLP,yLP,yLPTime,lazyCuts,cost,solverTime,relativeGAP,nodeCount,meters,tourMinutes,blocksMeters,numVisitedBlocks,intersectionCutsTime,intersectionCuts1,intersectionCuts2
1,31,930,5,120.0,IP,18.0,18.0,0.513351861,0,18.00,0.034353803,0.0,0,5027.0,14.673,1585.0,5,0.339560982,0,0
--
instance,|V|,|A|,|B|,T,model,initialLP,yLP,yLPTime,lazyCuts,cost,solverTime,relativeGAP,nodeCount,meters,tourMinutes,blocksMeters,numVisitedBlocks,intersectionCutsTime,intersectionCuts1,intersectionCuts2
2,58,3306,5,120.0,IP,21.0,21.0,0.403875981,0,21.00,0.714356294,0.0,0,4062.0,17.7165,2583.0,5,0.001910256,0,0
--
instance,|V|,|A|,|B|,T,model,initialLP,yLP,yLPTime,lazyCuts,cost,solverTime,relativeGAP,nodeCount,meters,tourMinutes,blocksMeters,numVisitedBlocks,intersectionCutsTime,intersectionCuts1,intersectionCuts2
3,22,462,5,120.0,IP,27.0,27.0,1.518863074,0,27.00,0.082502341,0.0,0,4466.0,15.357000000000003,1924.0,5,0.024349873,0,0
--
instance,|V|,|A|,|B|,T,model,initialLP,yLP,yLPTime,lazyCuts,cost,solverTime,relativeGAP,nodeCount,meters,tourMinutes,blocksMeters,numVisitedBlocks,intersectionCutsTime,intersectionCuts1,intersectionCuts2
4,25,600,5,120.0,IP,14.0,14.0,1.306037827,0,14.00,0.017199381,0.0,0,2200.0,14.122500000000002,2405.0,5,0.009738489,0,0
--
instance,|V|,|A|,|B|,T,model,initialLP,yLP,yLPTime,lazyCuts,cost,solverTime,relativeGAP,nodeCount,meters,tourMinutes,blocksMeters,numVisitedBlocks,intersectionCutsTime,intersectionCuts1,intersectionCuts2
5,50,2450,5,120.0,IP,25.0,25.0,0.364455516,0,25.00,0.078959983,0.0,0,5433.0,15.903000000000002,1723.0,5,0.015641141,0,0
```

For replicating all the paper experiments, simply run the bash script `sh.sh`.
