# PSI: A library for Parallel Spatial Indexes
PSI is a high-performance parallel library for a collection of spatial indexes, e.g., Kd-trees, Quad/Oct-trees, and R-trees, which are:
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

@article{men2026parallel,
  title={Parallel Dynamic Spatial Indexes},
  author={Men, Ziyang and Huang, Bo and Gu, Yan and Sun, Yihan},
  journal={arXiv preprint arXiv:2601.05347},
  year={2026}
}
```

Have a good day! ☀️

