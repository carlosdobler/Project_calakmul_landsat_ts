# Notes for Calakmul Landsat Time Series Analysis

## Important references

The paper from DeVries et al. (2016) "Characterizing forest change..." is useful to: a) copy part of the methods (especially pre-processing), and b) see how surface reflectance images are transformed via tasseled cap (Crist's "A Tasseled Cap equivalent transformation for reflectance factor data" (1985) is referenced for the coefficients). BFAST here is used to monitor change, and as such, the study is designed to detect a single disturbance - in other words, they use BFAST in its "online" mode.

The difference between offline (my method) and online modes of change detection is described in Zhu's review "Change detection using Landsat time series" (2017).

The Distrubance Index, derived from tasseled cap bands, is described in Healey et al's "Comparison of Tasseled Cap-based Landsat data structures..." (2005). In the paper they argue that tasseled cap-derived data is better to detect forest disturbances than original Landsat data. Accuracy did not var greatly among tasseled cap data structures.


## Google Cloud

### R and RStudio Server resources
https://grantmcdermott.com/2017/05/30/rstudio-server-compute-engine/ (setup through gcloud, then launches RStudio Server)

https://cloudyr.github.io/googleComputeEngineR/index.html (everything done from RStudio (desktop) to then launch RStudio Server)

### Add storage
https://cloud.google.com/compute/docs/disks/add-persistent-disk (add, format, and mount instructions)

### Pricing
[Virtual machines](https://cloud.google.com/compute/vm-instance-pricing)  
[Persistent disks](https://cloud.google.com/compute/disks-image-pricing)




