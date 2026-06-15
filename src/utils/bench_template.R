# mini benchmark ----
library(bench)

k <- 10000

# Benchmarking code with time and memory tracking
bench_results <- bench::mark(
  threads_1 = fastKmeans(coords, k = k, iter.max = 15, threads = 1),
  threads_3 = fastKmeans(coords, k = k, iter.max = 15, threads = 3),
  threads_8 = fastKmeans(coords, k = k, iter.max = 15, threads = 8),
  iterations = 3,   # Number of times to run each (keep it low for heavy computations)
  check = FALSE,# CRITICAL: k-means results might vary slightly between runs/threads, this prevents bench from throwing an error
  memory = TRUE
)
print(bench_results)


bench_results2 <- bench::mark(
  threads_16 = fastKmeans(coords, k = k, iter.max = 15, threads = 16),
  threads_24 = fastKmeans(coords, k = k, iter.max = 15, threads = 24),
  threads_32 = fastKmeans(coords, k = k, iter.max = 15, threads = 32),
  iterations = 3,   # Number of times to run each (keep it low for heavy computations)
  check = FALSE,# CRITICAL: k-means results might vary slightly between runs/threads, this prevents bench from throwing an error
  memory = TRUE
)
print(bench_results2)


bench_results3 <- bench::mark(
  threads_10 = fastKmeans(coords, k = k, iter.max = 15, threads = 10),
  threads_13 = fastKmeans(coords, k = k, iter.max = 15, threads = 13),
  threads_16 = fastKmeans(coords, k = k, iter.max = 15, threads = 16),
  iterations = 3,   # Number of times to run each (keep it low for heavy computations)
  check = FALSE # CRITICAL: k-means results might vary slightly between runs/threads, this prevents bench from throwing an error
)
print(bench_results3)

saveRDS(bench_results2, file = "bench_results.rds")