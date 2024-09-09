# repoTemplate - Bioanalytics Repository Template
This repository contains a template for a reproducible research project. The idea of reproducible research is that all the components of your work, such as the directory of input and output files, and evidently the structure of the analysis scripts, are prespecified. In this way, your work can be easily checked by adviser, mentors, collaborators, others working in your area, and the results can be easily regenerated.


## Folders
Brief description of the contents of each folder. 

- `input/`
  - This folder should contain all the necessary inputs of your project.
  - The later includes the files that contain the input data (could be organized in subfolders). The files should be given appropriate names that correspond to their usage. The term input data refers to the the constant/reference data necessary to run the project, as well as, example data used to produce the example output files inside the `output` folder.
  - This folder should also include a `README.md` file describing the format and type of the input data, as well as, the script that utilizes them.
- `output/`
  - In this folder should be placed all the outputs of your project, produced using the `input` data. The outputs should be named appropriately.
  - If the project produces more than a few files, they should be organized in subfolders. In case the output files are produced gradually by different scripts, you should consider organizing the subfolders according to the implemented scripts. The names of the subfolders should be the names of the corresponding scripts.
  - This folder should also include a `README.md` file describing the outputs and listing the values of the input parameters used to produce the example output files.
- `scripts/`
  - This folder should contain all the code, organized in scripts. The code should be split in steps and the names of the scripts should follow the [tidyverse style guide](https://style.tidyverse.org/files.html).
  - It should also include a `README.md` file describing the individual scripts. Their individual inputs, outputs and how they are connected.
- `docs/`
  - This folder is optional. It should include markdown documents describing the workflow of the project, along with examples for different user scenarios.



## Getting Started

### Dependencies

- Describe any prerequisites, libraries, OS version, etc., needed before installing program.

### Installing

- How/where to download your program
- Any modifications needed to be made to files/folders

### Executing program

- How to run the program using the example input data to reproduce the example outputs.
- Step-by-step commands

```
code blocks for commands
```



## Inputs

* Describe the format and specifications of the input files/data.
* List of the input parameters, their type, accepted and default values.



## Outputs

* Describe the output files of the project.

  

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/BiodataAnalysisGroup/bioanalytics/blob/main/LICENSE) file for details

