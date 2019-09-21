# scClassify

<img src="man/figures/scClassifySticker.png" align="right" width="300"/>

Single cell classification via cell-type hierarchies based on ensemble learning and sample size estimation.


[![](https://travis-ci.com/SydneyBioX/scClassify.svg?token=Mgp3nhHdKGyRbLPymfkS&branch=master)]()

## Installation


Install `hopach` package from Bioconductor:

```r
BiocManager::install("hopach")
```


```r
library(devtools)
devtools::install_github("SydneyBioX/scClassify")
```

## Pretrained models

Currently available pre-trained scClassify models (in `scClassifyTrainModel` class)

|          Tissue        | Organism | Training Data | Accession |  Summary                                | Download `.rds`  |
| :--: | :--: | :--: | :--: | :--: | :--: |
| Primary visual cortex  |   mouse  |  Tasic (2018) | [GSE115746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115746)    | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Tasic2018_mouseNeuronal.html)|                  |
| Primary visual cortex  |   mouse  |  Tasic (2016) | [GSE71585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585)    | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Tasic2016_mouseNeuronal.html)|                  |
| Visual cortex          |   mouse  |  Hrvatin      | [GSE102827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102827)  |   [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Hrvatin_mouseNeuronal.html)        |                  |
|    Lung                |   mouse  |  Cohen        |   [GSE119228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119228)        |    [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Cohen_mouseLung.html)        |                  |
|    Kidney              |   mouse  |  Park         |  [GSE107585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107585)         | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Park_mouseKidney.html)          |                  |
|    Liver               |   human  |  MacParland   | [GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469)          | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/MacParlandres_humanLiver.html)   |                  |
|    Liver               |   human  |  Aizarani     |  [GSE124395](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124395)         | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Aizarani_humanLiver.html)        |                  |
|    Pancreas            |   human  |  Xin          |           |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Xin_humanPancreas.html)          |                  |
|    Pancreas            |   human  |  Wang         |           |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Wang_humanPancreas.html)          |                  |
|    Pancreas            |   human  |  Lawlor       |           |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Lawlor_humanPancreas.html)           |                  |
|    Pancreas            |   human  |  Segerstolpe  |           |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Segerstolpe_humanPancreas.html)           |                  |
|    Pancreas            |   human  |  Muraro       |           | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Muraro_humanPancreas.html)           |                  |
|    Pancreas            |   human  |  Baron        |           |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Baron_humanPancreas.html)          |                  |
|    Pancreas            |   human  |  joint        |           |           |                  |
|    Melanoma            |   human  |  Li           | [GSE123139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139)          |   [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Li_humanMelanoma.html)         |                  |
|    PBMC                |   human  |  Ding (joint) |           |           |                  |
|    Tabula Muris        |   mouse  |  Tabula Muris |           |   [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/TabulaMuris.html)        |                  |

