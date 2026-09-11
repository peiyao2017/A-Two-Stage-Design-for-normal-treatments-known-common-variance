Optimal Parameter Search for the Two-Stage Design
================

# Overview

This folder contains the R functions used to search for optimal design
parameters for the proposed two-stage treatment-selection procedure with
known common variance and known reference mean.

The main R file is:

`R functions for optimal parameters search.R`

Before running any simulation or optimization, first run this file so
that all required R functions are defined in the current R session.

``` r
source("C:/Users/n01654025/Downloads/R functions/R functions for optimal parameters search.R")
```

If the working directory is already the `R functions` folder, use:

``` r
source("R functions for optimal parameters search.R")
```

The file defines three helper functions and two main optimization
functions.

# Functions

## `stage1_power()`

This helper function computes the Stage-I power under the least
favorable configuration.

``` r
stage1_power(n1, gamma, delta1, delta2, k, sigma)
```

Arguments:

- `n1`: Stage-I sample size for each experimental treatment.
- `gamma`: Stage-I screening cutoff increment.
- `delta1`: Lower efficacy margin.
- `delta2`: Target efficacy margin.
- `k`: Number of experimental treatments.
- `sigma`: Known common standard deviation.

## `stage2_sample_size()`

This helper function computes the Stage-II sample size required to
achieve a specified Stage-II power.

``` r
stage2_sample_size(power2, alpha, delta2, sigma)
```

Arguments:

- `power2`: Desired Stage-II power, $1-\beta_2$.
- `alpha`: One-sided Stage-II significance level.
- `delta2`: Target efficacy margin.
- `sigma`: Known common standard deviation.

## `sampling_quantities()`

This helper function computes the sampling characteristics of a
candidate design.

``` r
sampling_quantities(n1, n2, gamma, delta1, delta2, k, sigma)
```

Arguments:

- `n1`: Stage-I sample size per experimental treatment.
- `n2`: Stage-II sample size for the selected treatment.
- `gamma`: Stage-I screening cutoff increment.
- `delta1`: Lower efficacy margin.
- `delta2`: Target efficacy margin.
- `k`: Number of experimental treatments.
- `sigma`: Known common standard deviation.

The function returns the probability of terminating after Stage I under
the global null, expected sample sizes under the global null and LFC,
the sampling-cost criterion, and the maximum possible sample size.

## `optimize_design1()`

Use this function when only the overall power is specified.

``` r
optimize_design1(
  n1_range, gamma_range, overall_power, alpha,
  delta1, delta2, k, sigma
)
```

Arguments:

- `n1_range`: Candidate integer values of the Stage-I sample size.
- `gamma_range`: Candidate values of the Stage-I screening cutoff
  increment.
- `overall_power`: Target overall power, $1-\beta$.
- `alpha`: One-sided Stage-II significance level.
- `delta1`: Lower efficacy margin.
- `delta2`: Target efficacy margin.
- `k`: Number of experimental treatments.
- `sigma`: Known common standard deviation.

For each candidate `(n1, gamma)` pair, the function computes Stage-I
power and retains pairs for which Stage-I power is strictly larger than
the target overall power. It then determines the required Stage-II
power, computes `n2`, and selects the design minimizing the
sampling-cost criterion.

The returned vector contains `n1`, `n2`, `gamma`, `P_terminate_null`,
`EN_null`, `EN_LFC`, `cost`, `Nmax`, `power_stage1`, and `power_stage2`.

## `optimize_design2()`

Use this function when both the overall power and Stage-II power are
specified.

``` r
optimize_design2(
  n1_range, gamma_range, overall_power, stage2_power,
  alpha, delta1, delta2, k, sigma
)
```

Arguments:

- `n1_range`: Candidate integer values of the Stage-I sample size.
- `gamma_range`: Candidate values of the Stage-I screening cutoff
  increment.
- `overall_power`: Target overall power, $1-\beta$.
- `stage2_power`: Prespecified Stage-II power, $1-\beta_2$.
- `alpha`: One-sided Stage-II significance level.
- `delta1`: Lower efficacy margin.
- `delta2`: Target efficacy margin.
- `k`: Number of experimental treatments.
- `sigma`: Known common standard deviation.

The Stage-II power must be strictly larger than the overall power. The
function first computes `n2`, then searches over `(n1, gamma)` and
retains designs satisfying the Stage-I power requirement. Among feasible
designs, it selects the one minimizing the sampling-cost criterion.

The returned vector contains the same quantities as
`optimize_design1()`.

# Illustrative Example 1: Table 1

The first example reproduces the simulation used for Table 1.

The simulation settings are:

- $k=2,3,4,5$;
- overall power $1-\beta=0.70,0.75,0.80,0.85,0.90$;
- $\delta_2=0.20,0.30$;
- $\delta_1=0$;
- $\alpha=0.05$;
- $\sigma=1$;
- `n1_range = 10:500`;
- `gamma_range = seq(0, 0.20, by = 0.005)`.

``` r
k_values <- 2:5
overall_power_values <- seq(0.70, 0.90, by = 0.05)
delta2_values <- c(0.20, 0.30)

delta1 <- 0
alpha <- 0.05
sigma <- 1

n1_range <- 10:500
gamma_range <- seq(0, 0.20, by = 0.005)

table_results <- list()
count <- 0

for (k in k_values) {
  for (overall_power in overall_power_values) {
    for (delta2 in delta2_values) {

      fit <- optimize_design1(
        n1_range = n1_range,
        gamma_range = gamma_range,
        overall_power = overall_power,
        alpha = alpha,
        delta1 = delta1,
        delta2 = delta2,
        k = k,
        sigma = sigma
      )

      count <- count + 1

      table_results[[count]] <- data.frame(
        k = k,
        overall_power = overall_power,
        delta2 = delta2,
        n1 = as.integer(fit["n1"]),
        n2 = as.integer(fit["n2"]),
        gamma = as.numeric(fit["gamma"]),
        E0_N = as.numeric(fit["EN_null"]),
        cost = as.numeric(fit["cost"]),
        Nmax = as.integer(fit["Nmax"]),
        one_minus_pi0 = as.numeric(fit["P_terminate_null"]),
        stage1_power = as.numeric(fit["power_stage1"]),
        stage2_power = as.numeric(fit["power_stage2"])
      )
    }
  }
}

design_table1 <- do.call(rbind, table_results)
design_table1 <- unique(design_table1)
rownames(design_table1) <- NULL

if (nrow(design_table1) != 40) {
  warning(paste("Expected 40 rows, but obtained", nrow(design_table1), "rows."))
}
```

## Print Table 1

``` r
design_table1
```

    ##    k overall_power delta2  n1  n2 gamma      E0_N      cost Nmax one_minus_pi0
    ## 1  2          0.70    0.2  47 196 0.010  235.4972  257.7689  290     0.2780756
    ## 2  2          0.70    0.3  20  88 0.005  105.2080  114.5095  128     0.2589994
    ## 3  2          0.75    0.2  57 213 0.000  273.7500  296.8857  327     0.2500000
    ## 4  2          0.75    0.3  26  96 0.020  119.9427  131.9814  148     0.2922634
    ## 5  2          0.80    0.2  73 239 0.015  312.4423  344.9694  385     0.3035887
    ## 6  2          0.80    0.3  32 105 0.010  140.3282  153.2800  169     0.2730643
    ## 7  2          0.85    0.2  94 266 0.020  365.4790  406.6339  454     0.3327856
    ## 8  2          0.85    0.3  41 120 0.030  162.1640  180.6329  202     0.3319668
    ## 9  2          0.90    0.2 122 309 0.025  438.4815  493.2368  553     0.3706100
    ## 10 2          0.90    0.3  55 138 0.045  193.1046  219.2771  248     0.3977925
    ## 11 3          0.70    0.2  64 223 0.035  364.3182  385.7800  415     0.2272728
    ## 12 3          0.70    0.3  28 102 0.060  161.1497  171.5447  186     0.2436302
    ## 13 3          0.75    0.2  77 246 0.035  418.1937  444.0992  477     0.2390500
    ## 14 3          0.75    0.3  34 108 0.040  187.5687  197.5580  210     0.2076974
    ## 15 3          0.80    0.2  95 271 0.040  480.9966  515.0777  556     0.2767654
    ## 16 3          0.80    0.3  42 120 0.055  214.6540  228.9499  246     0.2612170
    ## 17 3          0.85    0.2 120 293 0.040  565.1242  606.4479  653     0.2999172
    ## 18 3          0.85    0.3  53 133 0.065  249.8165  269.5612  292     0.3171689
    ## 19 3          0.90    0.2 152 342 0.045  675.3443  734.2549  798     0.3586424
    ## 20 3          0.90    0.3  69 148 0.070  299.8660  326.3590  355     0.3725269
    ## 21 4          0.70    0.2  76 246 0.050  500.8584  521.9196  550     0.1997625
    ## 22 4          0.70    0.3  33 111 0.065  223.7201  232.0383  243     0.1736930
    ## 23 4          0.75    0.2  90 269 0.045  576.3060  599.8523  629     0.1958886
    ## 24 4          0.75    0.3  41 117 0.075  255.3196  266.7559  281     0.2194910
    ## 25 4          0.80    0.2 110 293 0.050  662.6502  694.9189  733     0.2401019
    ## 26 4          0.80    0.3  48 133 0.070  295.5195  309.0669  325     0.2216578
    ## 27 4          0.85    0.2 136 318 0.050  776.5006  816.8684  862     0.2688661
    ## 28 4          0.85    0.3  61 139 0.075  345.4406  363.1932  383     0.2702115
    ## 29 4          0.90    0.2 172 355 0.050  934.2240  986.8153 1043     0.3064112
    ## 30 4          0.90    0.3  76 160 0.075  415.1360  438.7493  464     0.3053999
    ## 31 5          0.70    0.2  85 268 0.060  644.6714  665.4865  693     0.1803307
    ## 32 5          0.70    0.3  38 115 0.075  288.5147  295.7518  305     0.1433502
    ## 33 5          0.75    0.2 101 284 0.050  743.8232  764.2634  789     0.1590731
    ## 34 5          0.75    0.3  45 128 0.085  328.9591  339.7264  353     0.1878193
    ## 35 5          0.80    0.2 122 308 0.055  854.9160  884.0920  918     0.2048183
    ## 36 5          0.80    0.3  54 136 0.075  381.5953  392.9525  406     0.1794461
    ## 37 5          0.85    0.2 148 339 0.055  999.4691 1037.1691 1079     0.2346045
    ## 38 5          0.85    0.3  66 147 0.075  446.7676  461.1832  477     0.2056627
    ## 39 5          0.90    0.2 186 374 0.055 1200.5096 1250.6497 1304     0.2767124
    ## 40 5          0.90    0.3  82 169 0.080  534.5477  556.1010  579     0.2630315
    ##    stage1_power stage2_power
    ## 1     0.7993268    0.8757370
    ## 2     0.7965112    0.8788326
    ## 3     0.8346105    0.8986227
    ## 4     0.8314306    0.9020597
    ## 5     0.8641323    0.9257842
    ## 6     0.8662606    0.9235096
    ## 7     0.8976394    0.9469282
    ## 8     0.8950701    0.9496463
    ## 9     0.9285189    0.9692856
    ## 10    0.9279932    0.9698347
    ## 11    0.7690986    0.9101564
    ## 12    0.7634355    0.9169079
    ## 13    0.8047157    0.9320062
    ## 14    0.8070496    0.9293110
    ## 15    0.8419353    0.9501918
    ## 16    0.8425253    0.9495264
    ## 17    0.8833376    0.9622595
    ## 18    0.8807952    0.9650371
    ## 19    0.9183852    0.9799810
    ## 20    0.9207454    0.9774689
    ## 21    0.7511203    0.9319413
    ## 22    0.7486991    0.9349550
    ## 23    0.7903148    0.9489890
    ## 24    0.7936675    0.9449801
    ## 25    0.8313597    0.9622790
    ## 26    0.8288291    0.9652170
    ## 27    0.8740104    0.9725285
    ## 28    0.8757530    0.9705933
    ## 29    0.9154386    0.9831353
    ## 30    0.9145018    0.9841424
    ## 31    0.7381741    0.9482857
    ## 32    0.7430920    0.9420098
    ## 33    0.7831911    0.9576207
    ## 34    0.7815619    0.9596169
    ## 35    0.8258143    0.9687408
    ## 36    0.8264421    0.9680049
    ## 37    0.8681095    0.9791392
    ## 38    0.8702532    0.9767272
    ## 39    0.9119572    0.9868884
    ## 40    0.9109889    0.9879374

``` r
knitr::kable(
  design_table1,
  digits = c(0, 2, 2, 0, 0, 3, 2, 2, 0, 2, 2, 2),
  caption = "Optimal designs when only the overall power is prespecified."
)
```

| k | overall_power | delta2 | n1 | n2 | gamma | E0_N | cost | Nmax | one_minus_pi0 | stage1_power | stage2_power |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.70 | 0.2 | 47 | 196 | 0.010 | 235.50 | 257.77 | 290 | 0.28 | 0.80 | 0.88 |
| 2 | 0.70 | 0.3 | 20 | 88 | 0.005 | 105.21 | 114.51 | 128 | 0.26 | 0.80 | 0.88 |
| 2 | 0.75 | 0.2 | 57 | 213 | 0.000 | 273.75 | 296.89 | 327 | 0.25 | 0.83 | 0.90 |
| 2 | 0.75 | 0.3 | 26 | 96 | 0.020 | 119.94 | 131.98 | 148 | 0.29 | 0.83 | 0.90 |
| 2 | 0.80 | 0.2 | 73 | 239 | 0.015 | 312.44 | 344.97 | 385 | 0.30 | 0.86 | 0.93 |
| 2 | 0.80 | 0.3 | 32 | 105 | 0.010 | 140.33 | 153.28 | 169 | 0.27 | 0.87 | 0.92 |
| 2 | 0.85 | 0.2 | 94 | 266 | 0.020 | 365.48 | 406.63 | 454 | 0.33 | 0.90 | 0.95 |
| 2 | 0.85 | 0.3 | 41 | 120 | 0.030 | 162.16 | 180.63 | 202 | 0.33 | 0.90 | 0.95 |
| 2 | 0.90 | 0.2 | 122 | 309 | 0.025 | 438.48 | 493.24 | 553 | 0.37 | 0.93 | 0.97 |
| 2 | 0.90 | 0.3 | 55 | 138 | 0.045 | 193.10 | 219.28 | 248 | 0.40 | 0.93 | 0.97 |
| 3 | 0.70 | 0.2 | 64 | 223 | 0.035 | 364.32 | 385.78 | 415 | 0.23 | 0.77 | 0.91 |
| 3 | 0.70 | 0.3 | 28 | 102 | 0.060 | 161.15 | 171.54 | 186 | 0.24 | 0.76 | 0.92 |
| 3 | 0.75 | 0.2 | 77 | 246 | 0.035 | 418.19 | 444.10 | 477 | 0.24 | 0.80 | 0.93 |
| 3 | 0.75 | 0.3 | 34 | 108 | 0.040 | 187.57 | 197.56 | 210 | 0.21 | 0.81 | 0.93 |
| 3 | 0.80 | 0.2 | 95 | 271 | 0.040 | 481.00 | 515.08 | 556 | 0.28 | 0.84 | 0.95 |
| 3 | 0.80 | 0.3 | 42 | 120 | 0.055 | 214.65 | 228.95 | 246 | 0.26 | 0.84 | 0.95 |
| 3 | 0.85 | 0.2 | 120 | 293 | 0.040 | 565.12 | 606.45 | 653 | 0.30 | 0.88 | 0.96 |
| 3 | 0.85 | 0.3 | 53 | 133 | 0.065 | 249.82 | 269.56 | 292 | 0.32 | 0.88 | 0.97 |
| 3 | 0.90 | 0.2 | 152 | 342 | 0.045 | 675.34 | 734.25 | 798 | 0.36 | 0.92 | 0.98 |
| 3 | 0.90 | 0.3 | 69 | 148 | 0.070 | 299.87 | 326.36 | 355 | 0.37 | 0.92 | 0.98 |
| 4 | 0.70 | 0.2 | 76 | 246 | 0.050 | 500.86 | 521.92 | 550 | 0.20 | 0.75 | 0.93 |
| 4 | 0.70 | 0.3 | 33 | 111 | 0.065 | 223.72 | 232.04 | 243 | 0.17 | 0.75 | 0.93 |
| 4 | 0.75 | 0.2 | 90 | 269 | 0.045 | 576.31 | 599.85 | 629 | 0.20 | 0.79 | 0.95 |
| 4 | 0.75 | 0.3 | 41 | 117 | 0.075 | 255.32 | 266.76 | 281 | 0.22 | 0.79 | 0.94 |
| 4 | 0.80 | 0.2 | 110 | 293 | 0.050 | 662.65 | 694.92 | 733 | 0.24 | 0.83 | 0.96 |
| 4 | 0.80 | 0.3 | 48 | 133 | 0.070 | 295.52 | 309.07 | 325 | 0.22 | 0.83 | 0.97 |
| 4 | 0.85 | 0.2 | 136 | 318 | 0.050 | 776.50 | 816.87 | 862 | 0.27 | 0.87 | 0.97 |
| 4 | 0.85 | 0.3 | 61 | 139 | 0.075 | 345.44 | 363.19 | 383 | 0.27 | 0.88 | 0.97 |
| 4 | 0.90 | 0.2 | 172 | 355 | 0.050 | 934.22 | 986.82 | 1043 | 0.31 | 0.92 | 0.98 |
| 4 | 0.90 | 0.3 | 76 | 160 | 0.075 | 415.14 | 438.75 | 464 | 0.31 | 0.91 | 0.98 |
| 5 | 0.70 | 0.2 | 85 | 268 | 0.060 | 644.67 | 665.49 | 693 | 0.18 | 0.74 | 0.95 |
| 5 | 0.70 | 0.3 | 38 | 115 | 0.075 | 288.51 | 295.75 | 305 | 0.14 | 0.74 | 0.94 |
| 5 | 0.75 | 0.2 | 101 | 284 | 0.050 | 743.82 | 764.26 | 789 | 0.16 | 0.78 | 0.96 |
| 5 | 0.75 | 0.3 | 45 | 128 | 0.085 | 328.96 | 339.73 | 353 | 0.19 | 0.78 | 0.96 |
| 5 | 0.80 | 0.2 | 122 | 308 | 0.055 | 854.92 | 884.09 | 918 | 0.20 | 0.83 | 0.97 |
| 5 | 0.80 | 0.3 | 54 | 136 | 0.075 | 381.60 | 392.95 | 406 | 0.18 | 0.83 | 0.97 |
| 5 | 0.85 | 0.2 | 148 | 339 | 0.055 | 999.47 | 1037.17 | 1079 | 0.23 | 0.87 | 0.98 |
| 5 | 0.85 | 0.3 | 66 | 147 | 0.075 | 446.77 | 461.18 | 477 | 0.21 | 0.87 | 0.98 |
| 5 | 0.90 | 0.2 | 186 | 374 | 0.055 | 1200.51 | 1250.65 | 1304 | 0.28 | 0.91 | 0.99 |
| 5 | 0.90 | 0.3 | 82 | 169 | 0.080 | 534.55 | 556.10 | 579 | 0.26 | 0.91 | 0.99 |

Optimal designs when only the overall power is prespecified.

# Illustrative Example 2: Table 2

The second example uses both prespecified overall power and Stage-II
power.

The simulation settings are:

- $k=2,3,4,5$;
- overall power $1-\beta=0.70,0.80$;
- Stage-II power $1-\beta_2=0.80,0.90$, subject to Stage-II power being
  strictly larger than overall power;
- $\delta_2=0.20,0.30$;
- $\delta_1=0$;
- $\alpha=0.05$;
- $\sigma=1$;
- `n1_range = 10:500`;
- `gamma_range = seq(0, 0.20, by = 0.005)`.

For overall power 0.70, Stage-II powers 0.80 and 0.90 are used. For
overall power 0.80, only Stage-II power 0.90 is used.

``` r
k_values <- 2:5
overall_power_values <- c(0.70, 0.80)
delta2_values <- c(0.20, 0.30)

delta1 <- 0
alpha <- 0.05
sigma <- 1

n1_range <- 10:500
gamma_range <- seq(0, 0.20, by = 0.005)

table2_results <- list()
count <- 0

for (k in k_values) {
  for (overall_power in overall_power_values) {

    if (overall_power == 0.70) {
      stage2_power_values <- c(0.80, 0.90)
    }

    if (overall_power == 0.80) {
      stage2_power_values <- c(0.90)
    }

    for (stage2_power in stage2_power_values) {
      for (delta2 in delta2_values) {

        fit <- optimize_design2(
          n1_range = n1_range,
          gamma_range = gamma_range,
          overall_power = overall_power,
          stage2_power = stage2_power,
          alpha = alpha,
          delta1 = delta1,
          delta2 = delta2,
          k = k,
          sigma = sigma
        )

        count <- count + 1

        table2_results[[count]] <- data.frame(
          k = k,
          overall_power = overall_power,
          stage2_power = stage2_power,
          delta2 = delta2,
          n1 = as.integer(fit["n1"]),
          n2 = as.integer(fit["n2"]),
          gamma = as.numeric(fit["gamma"]),
          E0_N = as.numeric(fit["EN_null"]),
          cost = as.numeric(fit["cost"]),
          Nmax = as.integer(fit["Nmax"]),
          one_minus_pi0 = as.numeric(fit["P_terminate_null"]),
          stage1_power = as.numeric(fit["power_stage1"])
        )
      }
    }
  }
}

design_table2 <- do.call(rbind, table2_results)
design_table2 <- unique(design_table2)
rownames(design_table2) <- NULL

if (nrow(design_table2) != 24) {
  warning(paste("Expected 24 rows, but obtained", nrow(design_table2), "rows."))
}
```

## Print Table 2

``` r
design_table2
```

    ##    k overall_power stage2_power delta2  n1  n2 gamma     E0_N      cost Nmax
    ## 1  2           0.7          0.8    0.2  76 155 0.000 268.2500  286.0511  307
    ## 2  2           0.7          0.8    0.3  34  69 0.010 118.1085  126.7343  137
    ## 3  2           0.7          0.9    0.2  45 215 0.030 232.7371  260.9498  305
    ## 4  2           0.7          0.9    0.3  20  96 0.045 103.7338  116.3311  136
    ## 5  2           0.8          0.9    0.2  83 215 0.000 327.2500  352.2856  381
    ## 6  2           0.8          0.9    0.3  38  96 0.015 144.3335  157.1496  172
    ## 7  3           0.7          0.8    0.2 110 155 0.020 454.2744  468.8593  485
    ## 8  3           0.7          0.8    0.3  48  69 0.005 203.6400  208.1334  213
    ## 9  3           0.7          0.9    0.2  67 215 0.035 366.5367  387.6998  416
    ## 10 3           0.7          0.9    0.3  29  96 0.035 164.7733  172.6692  183
    ## 11 3           0.8          0.9    0.2 119 215 0.015 533.2244  551.8645  572
    ## 12 3           0.8          0.9    0.3  54  96 0.045 234.0461  245.4433  258
    ## 13 4           0.7          0.8    0.2 132 155 0.020 664.1071  673.2447  683
    ## 14 4           0.7          0.8    0.3  59  69 0.040 294.7601  299.6911  305
    ## 15 4           0.7          0.9    0.2  85 215 0.045 513.9851  532.1188  555
    ## 16 4           0.7          0.9    0.3  37  96 0.045 230.8942  236.7956  244
    ## 17 4           0.8          0.9    0.2 145 215 0.040 747.6692  770.4014  795
    ## 18 4           0.8          0.9    0.3  64  96 0.050 334.2845  342.8348  352
    ## 19 5           0.7          0.8    0.2 149 155 0.020 888.3010  894.0132  900
    ## 20 5           0.7          0.8    0.3  67  69 0.060 393.3382  398.4775  404
    ## 21 5           0.7          0.9    0.2  98 215 0.040 679.2885  691.0314  705
    ## 22 5           0.7          0.9    0.3  43  96 0.035 304.0922  307.3056  311
    ## 23 5           0.8          0.9    0.2 161 215 0.035 990.6427 1004.9246 1020
    ## 24 5           0.8          0.9    0.3  72  96 0.065 438.7562  447.0976  456
    ##    one_minus_pi0 stage1_power
    ## 1     0.25000000    0.8768011
    ## 2     0.27378947    0.8750089
    ## 3     0.33610651    0.7782160
    ## 4     0.33610651    0.7782160
    ## 5     0.25000000    0.8889688
    ## 6     0.28819309    0.8891936
    ## 7     0.19822993    0.8751039
    ## 8     0.13565175    0.8750752
    ## 9     0.23006194    0.7779416
    ## 10    0.18986162    0.7779462
    ## 11    0.18035164    0.8891656
    ## 12    0.24952027    0.8889872
    ## 13    0.12188959    0.8750287
    ## 14    0.14840490    0.8752312
    ## 15    0.19076708    0.7778080
    ## 16    0.13651825    0.7779028
    ## 17    0.22014330    0.8889925
    ## 18    0.18453697    0.8891466
    ## 19    0.07547740    0.8750228
    ## 20    0.15451915    0.8753979
    ## 21    0.11958856    0.7784703
    ## 22    0.07195643    0.7779249
    ## 23    0.13654581    0.8891381
    ## 24    0.17962306    0.8891386

``` r
knitr::kable(
  design_table2,
  digits = c(0, 2, 2, 2, 0, 0, 3, 2, 2, 0, 2, 2),
  caption = "Optimal designs when both overall power and Stage-II power are prespecified."
)
```

| k | overall_power | stage2_power | delta2 | n1 | n2 | gamma | E0_N | cost | Nmax | one_minus_pi0 | stage1_power |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.7 | 0.8 | 0.2 | 76 | 155 | 0.000 | 268.25 | 286.05 | 307 | 0.25 | 0.88 |
| 2 | 0.7 | 0.8 | 0.3 | 34 | 69 | 0.010 | 118.11 | 126.73 | 137 | 0.27 | 0.88 |
| 2 | 0.7 | 0.9 | 0.2 | 45 | 215 | 0.030 | 232.74 | 260.95 | 305 | 0.34 | 0.78 |
| 2 | 0.7 | 0.9 | 0.3 | 20 | 96 | 0.045 | 103.73 | 116.33 | 136 | 0.34 | 0.78 |
| 2 | 0.8 | 0.9 | 0.2 | 83 | 215 | 0.000 | 327.25 | 352.29 | 381 | 0.25 | 0.89 |
| 2 | 0.8 | 0.9 | 0.3 | 38 | 96 | 0.015 | 144.33 | 157.15 | 172 | 0.29 | 0.89 |
| 3 | 0.7 | 0.8 | 0.2 | 110 | 155 | 0.020 | 454.27 | 468.86 | 485 | 0.20 | 0.88 |
| 3 | 0.7 | 0.8 | 0.3 | 48 | 69 | 0.005 | 203.64 | 208.13 | 213 | 0.14 | 0.88 |
| 3 | 0.7 | 0.9 | 0.2 | 67 | 215 | 0.035 | 366.54 | 387.70 | 416 | 0.23 | 0.78 |
| 3 | 0.7 | 0.9 | 0.3 | 29 | 96 | 0.035 | 164.77 | 172.67 | 183 | 0.19 | 0.78 |
| 3 | 0.8 | 0.9 | 0.2 | 119 | 215 | 0.015 | 533.22 | 551.86 | 572 | 0.18 | 0.89 |
| 3 | 0.8 | 0.9 | 0.3 | 54 | 96 | 0.045 | 234.05 | 245.44 | 258 | 0.25 | 0.89 |
| 4 | 0.7 | 0.8 | 0.2 | 132 | 155 | 0.020 | 664.11 | 673.24 | 683 | 0.12 | 0.88 |
| 4 | 0.7 | 0.8 | 0.3 | 59 | 69 | 0.040 | 294.76 | 299.69 | 305 | 0.15 | 0.88 |
| 4 | 0.7 | 0.9 | 0.2 | 85 | 215 | 0.045 | 513.99 | 532.12 | 555 | 0.19 | 0.78 |
| 4 | 0.7 | 0.9 | 0.3 | 37 | 96 | 0.045 | 230.89 | 236.80 | 244 | 0.14 | 0.78 |
| 4 | 0.8 | 0.9 | 0.2 | 145 | 215 | 0.040 | 747.67 | 770.40 | 795 | 0.22 | 0.89 |
| 4 | 0.8 | 0.9 | 0.3 | 64 | 96 | 0.050 | 334.28 | 342.83 | 352 | 0.18 | 0.89 |
| 5 | 0.7 | 0.8 | 0.2 | 149 | 155 | 0.020 | 888.30 | 894.01 | 900 | 0.08 | 0.88 |
| 5 | 0.7 | 0.8 | 0.3 | 67 | 69 | 0.060 | 393.34 | 398.48 | 404 | 0.15 | 0.88 |
| 5 | 0.7 | 0.9 | 0.2 | 98 | 215 | 0.040 | 679.29 | 691.03 | 705 | 0.12 | 0.78 |
| 5 | 0.7 | 0.9 | 0.3 | 43 | 96 | 0.035 | 304.09 | 307.31 | 311 | 0.07 | 0.78 |
| 5 | 0.8 | 0.9 | 0.2 | 161 | 215 | 0.035 | 990.64 | 1004.92 | 1020 | 0.14 | 0.89 |
| 5 | 0.8 | 0.9 | 0.3 | 72 | 96 | 0.065 | 438.76 | 447.10 | 456 | 0.18 | 0.89 |

Optimal designs when both overall power and Stage-II power are
prespecified.

# Notes

The optimization uses a grid search over `n1_range` and `gamma_range`. A
finer grid can provide a more detailed search but will require more
computation time.

If no feasible design is found, enlarge the search ranges.

To generate a GitHub README, knit this file using the `github_document`
output format. This creates a `README.md` file that can be uploaded
together with the R source file.
