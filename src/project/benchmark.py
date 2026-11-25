
import time
import statistics


def benchmark_function(func, *args, repeats=5, warmup=1, **kwargs):
    """
    Benchmark a function and return execution statistics.
    """

    # Warm-up runs
    for _ in range(warmup):
        func(*args, **kwargs)

    execution_times = []

    for _ in range(repeats):
        start = time.perf_counter_ns()

        func(*args, **kwargs)

        end = time.perf_counter_ns()

        execution_times.append(end - start)

    return {
        "runs": repeats,
        "average_ms": statistics.mean(execution_times) / 1_000_000,
        "median_ms": statistics.median(execution_times) / 1_000_000,
        "min_ms": min(execution_times) / 1_000_000,
        "max_ms": max(execution_times) / 1_000_000,
    }
