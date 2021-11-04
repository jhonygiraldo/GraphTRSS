# Graph Time-varying Signal Reconstruction via Sobolev Smoothness
Authors: [Jhony H Giraldo](https://sites.google.com/view/jhonygiraldo), [Arif Mahmood](https://itu.edu.pk/faculty-itu/dr-arif-mahmood/), [Dorina Thanou](https://people.epfl.ch/dorina.thanou?lang=en), and [Thierry Bouwmans](https://sites.google.com/site/thierrybouwmans/)
- - - -
![Pipeline](https://github.com/jhonygiraldo/GraphTRSS/blob/main/doc/pipeline.png)
- - - -
**Abstract**: Graph Signal Processing (GSP) is an emerging research field that tries to extend the concepts of digital signal processing to graphs. GSP has multiple applications in sensor networks, machine learning, and image processing to name a few. The sampling and reconstruction of static graph signals have played a central role in GSP, where the bandlimitedness of the underlying signals is a common prior assumption. However, many real-world graph signals are time-varying in nature, and some works in the literature have focused on the smoothness of the temporal differences rather than on the bandlimitedness of the graph signals. In this work, we also assume that the temporal differences of graph signals are smooth, and then we introduce a novel algorithm based on the extension of a Sobolev smoothness function from static to time-varying graph signals for reconstruction from samples. We dub our algorithm as Graph Time-varying signal Reconstruction via Sobolev Smoothness (GraphTRSS). We explore some theoretical aspects of the converge rate of GraphTRSS by studying the condition number of the Hessian associated with our optimization problem. Our algorithm has the advantage of converging faster than other methods that are based on Laplacian operators, without requiring expensive eigenvalue decompositions or matrix inversions. GraphTRSS outperforms several methods for time-varying graph signal reconstruction on two COVID-19 datasets. In addition, our algorithm also performs well on two environmental datasets for the reconstruction of particulate matter and sea surface temperature.
- - - -
## Getting started

This repository is the code released for our paper Graph Time-varying Signal Reconstruction via Sobolev Smoothness. [Journal Paper](https://doi.org/).

#### Installation

Clone this repository
```bash
git clone https://github.com/jhonygiraldo/GraphTRSS  
```
Inside the repository, initialize and clone the submodules
```bash
git submodule init
git submodule update
```
Install gspboox.
- - - -
## Citation

If you use our code, please cite

        @article{giraldo2022graph,
          title={Graph Time-varying Signal Reconstruction via Sobolev Smoothness},
          author={Giraldo, Jhony H and Mahmood, Arif and Garcia-Garcia, Belmar and Thanou, Dorina and Bouwmans, Thierry},
          journal={IEEE Transactions on Signal and Information Processing over Networks},
          year={2022}
        }
        
        @inproceedings{giraldo2020minimization,
          title={On the minimization of Sobolev norms of time-varying graph signals: Estimation of new Coronavirus disease 2019 cases},
          author={Giraldo, Jhony H and Bouwmans, Thierry},
          booktitle={IEEE International Workshop on Machine Learning for Signal Processing},
          pages={1--6},
          year={2020},
          organization={IEEE}
        }

- - - -
## References

- This repository builds upon, thus borrows code from [gspbox](https://github.com/epfl-lts2/gspbox) and from the paper "Qiu, K., Mao, X., Shen, X., Wang, X., Li, T., & Gu, Y. (2017). Time-varying graph signal reconstruction. IEEE Journal of Selected Topics in Signal Processing, 11(6), 870-88". Therefore please consider citing the following papers if you use this repository:
- Qiu, K., Mao, X., Shen, X., Wang, X., Li, T., & Gu, Y. (2017). Time-varying graph signal reconstruction. IEEE Journal of Selected Topics in Signal Processing, 11(6), 870-88
- Perraudin, N., Paratte, J., Shuman, D., Martin, L., Kalofolias, V., Vandergheynst, P., & Hammond, D. K. (2014). GSPBOX: A toolbox for signal processing on graphs. arXiv preprint arXiv:1408.5781.

- We also use code from [Blue-Noise](https://github.com/jhonygiraldo/Blue-Noise-Sampling-on-Graphs) and from the paper "Anis, A., Gadde, A., & Ortega, A. (2016). Efficient sampling set selection for bandlimited graph signals using graph spectral proxies. IEEE Transactions on Signal Processing, 64(14), 3775-3789." for the deterministic sampling methods. Therefore, if you use the codes of deterministic sampling please consider citing the following papers:
- Parada-Mayorga, A., Lau, D. L., Giraldo, J. H., & Arce, G. R. (2019). Blue-noise sampling on graphs. IEEE Transactions on Signal and Information Processing over Networks, 5(3), 554-569.
- Anis, A., Gadde, A., & Ortega, A. (2016). Efficient sampling set selection for bandlimited graph signals using graph spectral proxies. IEEE Transactions on Signal Processing, 64(14), 3775-3789.
