# scClassify: hierarchical classification of cells

<img src="man/figures/scClassifySticker.png" align="right" width="200"/>

Single cell classification via cell-type hierarchies based on ensemble learning and sample size estimation.


[![](https://travis-ci.org/SydneyBioX/scClassify.svg?branch=master)]()

## Installation


Install Bioconductor packages `S4Vectors`, `hopach` and `limma` packages using `BiocManager`:

```r
# install.packages("BiocManager")
BiocManager::install(c("S4Vectors", "hopach", "limma"))
```

Then install `scClassify` using `devtools`:

```r
library(devtools)
devtools::install_github("SydneyBioX/scClassify")
```

## Vignette and Shiny app

You can find the vignette at this webiste: https://sydneybiox.github.io/scClassify/index.html.

Also, you can find our interactive shiny application (beta) at this webiste:
http://shiny.maths.usyd.edu.au/scClassify.



## Pretrained models

Currently available pre-trained scClassify models (in `scClassifyTrainModel` class)

|          Tissue        | Organism | Training Data | Accession |  Summary                                | Download `.rds`  |
| :--: | :--: | :--: | :--: | :--: | :--: |
| Primary visual cortex  |   mouse  |  Tasic (2018) | [GSE115746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115746)    | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Tasic2018_mouseNeuronal.html)|    [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainTasicV2resClass.rds)    |
| Primary visual cortex  |   mouse  |  Tasic (2016) | [GSE71585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585)    | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Tasic2016_mouseNeuronal.html)|  [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainTasicV1resClass.rds)   |
| Visual cortex          |   mouse  |  Hrvatin      | [GSE102827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102827)  |   [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Hrvatin_mouseNeuronal.html)        |  [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainHrvatinresClass.rds)  |
|    Lung                |   mouse  |  Cohen        |   [GSE119228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119228)        |    [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Cohen_mouseLung.html)        |  [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainCohenresClass.rds)    |
|    Kidney              |   mouse  |  Park         |  [GSE107585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107585)         | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Park_mouseKidney.html)          |   [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainParkresClass.rds)     |
|    Liver               |   human  |  MacParland   | [GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469)          | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/MacParlandres_humanLiver.html)   | [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainMacParlandresClass.rds)  |
|    Liver               |   human  |  Aizarani     |  [GSE124395](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124395)         | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Aizarani_humanLiver.html)        |  [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainAizaraniresClass.rds) |
|    Pancreas            |   human  |  Xin          |   [GSE81608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81608)         |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Xin_humanPancreas.html)          |  [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainXinClass.rds)                   |
|    Pancreas            |   human  |  Wang         |   [GSE83139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83139)        |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Wang_humanPancreas.html)          |     [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainWangClass.rds)                |
|    Pancreas            |   human  |  Lawlor       |   [GSE86469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469)         |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Lawlor_humanPancreas.html)           |     [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainLawlorClass.rds)                |
|    Pancreas            |   human  |  Segerstolpe  |   [E-MTAB-5061](https://www.ebi.ac.uk/arrayexpress/E-MTAB-5061)        |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Segerstolpe_humanPancreas.html)           |  [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainSegerstolpeClass.rds)                  |
|    Pancreas            |   human  |  Muraro       |    [GSE85241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241)         | [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Muraro_humanPancreas.html)           |     [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainMuraroClass.rds)               |
|    Pancreas            |   human  |  Baron        |   [GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133)         |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Baron_humanPancreas.html)          |   [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainBaronClass.rds)               |
|    Pancreas            |   human  |  joint        |    -     |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Joint_humanPancreas.html)          |        [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/jointPancreasClass.rds)          |
|    Melanoma            |   human  |  Li           | [GSE123139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123139)          |   [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Li_humanMelanoma.html)         |    [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainLiresClass.rds)     |
|    PBMC                |   human  |  Ding (joint) |    -   |  [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/Joint_humanPBMC.html)   |           [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/jointPBMCClass.rds)        |
|    Tabula Muris        |   mouse  |  Tabula Muris | [GSE109774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109774)           |   [link](https://SydneyBioX.github.io/scClassify/articles/webOnly/TabulaMuris.html)        |    [link](http://www.maths.usyd.edu.au/u/yingxinl/wwwnb/scClassify/trainTMresClass.rds)              |



## Contact us

If you have any enquiries, espcially about performing `scClassify` to classify your cells or to build your own models, please contact <yingxin.lin@sydney.edu.au> or <bioinformatics@maths.usyd.edu.au>.


## Reference

**scClassify: hierarchical classification of cells**

Yingxin Lin, Yue Cao, Hani J Kim, Agus Salim, Terence P. Speed, Dave Lin, Pengyi Yang, Jean Yee Hwa Yang

bioRxiv 776948; doi: https://doi.org/10.1101/776948


