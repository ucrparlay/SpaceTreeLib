# PSI: A library for Parallel Spatial Indexes

[![build and test](https://img.shields.io/github/actions/workflow/status/ucrparlay/SpaceTreeLib/ci.yml?branch=main&logo=github&label=build%20%26%20test)](https://github.com/ucrparlay/SpaceTreeLib/actions/workflows/ci.yml)
[![coverage](https://img.shields.io/codecov/c/github/ucrparlay/SpaceTreeLib/main?logo=codecov&logoColor=white)](https://codecov.io/gh/ucrparlay/SpaceTreeLib)
[![docker](https://img.shields.io/github/actions/workflow/status/ucrparlay/SpaceTreeLib/docker-build.yml?branch=main&logo=docker&logoColor=white&label=docker)](https://github.com/ucrparlay/SpaceTreeLib/actions/workflows/docker-build.yml)
[![C++20, header-only](https://img.shields.io/badge/C%2B%2B-20%20%C2%B7%20header--only-00599C?logo=cplusplus&logoColor=white)](CMakeLists.txt)
[![license](https://img.shields.io/github/license/ucrparlay/SpaceTreeLib?color=blue)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1145%2F3774934.3786412-1f6feb)](https://doi.org/10.1145/3774934.3786412)

PSI is a high-performance parallel library for a collection of spatial indexes -- kd-trees, quad/oct-trees, and 1-D trees over a space-filling curve -- which are:
- 🚀 Blazing fast, and scales to billions of input size,
- 🔀 Highly parallel, achieves almost linear speedup for hundreds of cores,
- 🎳 Supports tree construction, batch updates (with handling of imbalance), and a wide range of queries,
- 📝 Provide strong theoretical guarantees for the running time, parallelism, and I/O (cache) efficiency,
- 🛠️ Easy to adapt, integrate, and deploy.

## Docs
- [Quick Start](docs/QUICK_START.md)
- [User Manual](docs/MANUAL.md)
- [Docker](docs/DOCKER_QUICK_REFERENCE.md)
- [Artifact Evaluation](docs/ARTIFACT_EVALUATION.md)
- [How to Pick the Tree 🧐](docs/TREE_ANSWER.md)
- [Contributing](CONTRIBUTING.md)
- [Third-party code](THIRD_PARTY.md)

## Citation
If you use our code, please cite our papers:
```tex
@inproceedings{men2026dynamic,
    author = {Men, Ziyang and Huang, Bo and Gu, Yan and Sun, Yihan},
    title = {Parallel Dynamic Spatial Indexes},
    year = {2026},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    booktitle = {Proceedings of the 31st ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming},
    location = {Sydney, Australia},
    series = {PPoPP '26},
    doi = {10.1145/3774934.3786412}
}

@article{men2025parallel,
  title={Parallel kd-tree with Batch Updates},
  author={Men, Ziyang and Shen, Zheqi and Gu, Yan and Sun, Yihan},
  journal={Proceedings of the ACM on Management of Data},
  volume={3},
  number={1},
  pages={1--26},
  year={2025},
  publisher={ACM New York, NY, USA}
}
```

## License

PSI is MIT licensed — see [LICENSE](LICENSE).

It also redistributes third-party code under other terms; see
[THIRD_PARTY.md](THIRD_PARTY.md). In particular, the Hilbert curve
implementation requires this notice:

> Hilbert Curve implementation copyright 1998, Rice University

Have a good day! ☀️

