# CELERES

This repository contains the source code for CELERES, a linkable ring signature
scheme, together with protocol tests, benchmarks, and a signature-size script.

CELERES comes in two flavors, selected at build time via `CELERES_TAG_MODE`:

- `self_orthogonal` - the tag is a compressed self-orthogonal representative.
- `rref` - the tag is a full packed RREF matrix.

Three security levels are provided: `l1`, `l3` and `l5`.

## Third-Party Code

The implementation contains copied or adapted code/ideas from the following
repositories:
- [LESS official repository](https://github.com/less-sig/LESS) - the code comes without any license.
- [Calamari and Falafel](https://github.com/less-sig/LESS) - the code comes without any license.
- [SPECK](https://github.com/SPECK-Signature/SPECK/tree/main) - the code comes without any license.

## Build

Reference build:

```sh
cmake -S . -B build_ref -DIMPL=ref -DCMAKE_BUILD_TYPE=Release
cmake --build build_ref --parallel
```

AVX2 build:

```sh
cmake -S . -B build_avx2 -DIMPL=avx2 -DCMAKE_BUILD_TYPE=Release
cmake --build build_avx2 --parallel
```

Example with ring size 4 and the `rref` tag flavor:

```sh
cmake -S . -B build_avx2_r4_rref -DIMPL=avx2 -DRING_SIZE=4 \
      -DCELERES_TAG_MODE=rref -DCMAKE_BUILD_TYPE=Release
cmake --build build_avx2_r4_rref --parallel
```

## CMake Options

- `IMPL=ref|avx2`
- `RING_SIZE=<r>`
- `CELERES_TAG_MODE=self_orthogonal|rref`
- `BUILD_TESTS=ON|OFF`
- `BUILD_BENCHMARKS=ON|OFF`
- `BENCH_ITERS=<n>`

## Tests

```sh
ctest --test-dir build_avx2 --output-on-failure
```

Run one family:

```sh
ctest --test-dir build_avx2 -R "test_celeres_l" --output-on-failure
ctest --test-dir build_avx2 -R "test_celeres_link_l" --output-on-failure
```

## Benchmarks

```sh
./build_avx2/bench_celeres_l1
```

Replace `l1` with `l3` or `l5` for the other parameter sets.

## Signature Sizes

```sh
python3 scripts/signature_sizes.py --ring-sizes 2,4,64
```
