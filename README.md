# TError

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3590703.svg)](https://doi.org/10.5281/zenodo.3590703)

**TError** is a *Matlab* package designed to quantify the uncertainty of eruption source parameters calculated from tephra deposits. Inputs of the code are a range of field-based, model-based and empirical parameters (i.e., clast diameter, crosswind and downwind ranges, thickness measurement, area of isopach contours, bulk deposit density, empirical constants and wind speed).

**TError** is a contraction of Tephra and Error. Any bad joke is not intentional.

## Repository content

The repository contains the following files:

File | Description
----- | ------
`TError_sensitivity.m` | Systematic sensitivity analysis mode of TError
`TError_propagation.m` | Stochastic error propagation mode of TError
`TError.pdf` | User manual
`LICENSE.md` | License associated with the code
`README.md` | This file
`isopach_example.txt` | Example of a file containing isopach infomation
`dep/` | Dependencies to the code

## Usage
For both `TError_sensitivity.m` and `TError_propagation.m`, edit the header of the files and run the script in *Matlab*.

### Output
Upon completion, a folder named `OUTPUT` containing sub-folders named after the run name are created.


## Additional documentation
Instructions are provided in the the user manual `TError.pdf`. Updates are presented [here](https://e5k.github.io/pages/terror).


## Citation
Please cite **TError** as:
> Biass, S., Bagheri, G., Aeberhard, W., Bonadonna, C., 2014. TError: towards a better quantification of the uncertainty propagated during the characterization of tephra deposits. Stat. Volcanol. 1, 1–27. https://doi.org/10.5038/2163-338X.1.2

> Biass S, Bagheri G, Aeberhard W, Bonadonna C, 2014. TError v1.0. [doi:10.5281/zenodo.3590703](https://doi.org/10.5281/zenodo.3590702)

## License
GBF is a free and open source software releaser under GPL 3. See the documentation or the file `LICENSE.md` for further information.

## Related scripts
The method to estimate the plume height from the method of Carey and Sparks (1986) can be found in [this code](https://github.com/e5k/CareySparks86_Matlab) and documented [here](https://e5k.github.io/pages/cs86).

Volume calculations can be found in [this code](https://github.com/e5k/TephraFits) and documented [here](https://e5k.github.io/pages/tephrafits).

## Dependencies
This script uses *Matlab* dependencies `nhist.m` and `linspecer.m` by [Jonathan C Lansey](https://www.mathworks.com/matlabcentral/profile/authors/302713-jonathan-c-lansey). The calculation of the mass eruption rate with the method of Degruyter and Bonadonna (2012) is made possible with the script found in:
> Degruyter, W., Bonadonna, C., 2012. Improving on mass flow rate estimates of volcanic eruptions. Geophys Res Lett 39. https://doi.org/10.1029/2012GL052566