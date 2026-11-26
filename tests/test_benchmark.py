
from src.project.benchmark import benchmark_function


def sample_function():
    return sum(range(1000))


def test_benchmark_function():
    result = benchmark_function(
        sample_function,
        repeats=3,
        warmup=1
    )

    assert result["runs"] == 3
    assert result["average_ms"] >= 0
    assert result["median_ms"] >= 0
    assert result["min_ms"] >= 0
    assert result["max_ms"] >= result["min_ms"]
