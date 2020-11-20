# The CDC _Cyclospora cayetanensis_ complete genotyping workflow

This repository contains the workflow developed by Joel Barratt at the Centers for Disease Control and Prevention (CDC) for genotyping Cyclospora cayetanensis. This workflow is comprised of three modules that perform three major tasks required by CDC to identify genetically related clusters of Cyclospora cayetaneisis infections. 

 1. **MODULE 1** - Assign haplotypes to each specimen in accordance with the CDC typing method discussed [here](https://www.cambridge.org/core/journals/epidemiology-and-infection/article/evaluation-of-an-ensemblebased-distance-statistic-for-clustering-mlst-datasets-using-epidemiologically-defined-clusters-of-cyclosporiasis/F13FB2483E3DE1CF5C33706B4A1A7182).

 2. **MODULE 2** - Examine the genotype information generated by MODULE 1 and assess the relationship between each possible pair of specimens using this information. This second module is based on an updated version of the CDCs Eukaryotyping ensemble discribed [here](https://github.com/Joel-Barratt/Eukaryotyping).

 3. **MODULE 3** - Will predict the most appropriate number of clusters in the population under analysis using a set of reference specimens of known genetic linkage.  

_Please cite this [manuscript](https://www.cambridge.org/core/journals/parasitology/article/genotyping-genetically-heterogeneous-cyclospora-cayetanensis-infections-to-complement-epidemiological-case-linkage/0C51FBFFB172DF50357C1D171E9B8657) and this [manuscript](https://www.cambridge.org/core/journals/epidemiology-and-infection/article/evaluation-of-an-ensemblebased-distance-statistic-for-clustering-mlst-datasets-using-epidemiologically-defined-clusters-of-cyclosporiasis/F13FB2483E3DE1CF5C33706B4A1A7182#fndtn-information) if you use this workflow:_

```
1. Barratt, JLN, S Park, FS Nascimento, J Hofstetter, M Plucinski, S Casillas, RS Bradbury, MJ Arrowood, Y Qvarnstrom, E Talundzic (2019) Genotyping genetically heterogeneous Cyclospora cayetanensis infections to complement epidemiological case linkage. Parasitology:1–9 doi:10.1017/S0031182019000581


2. Nascimento, FS, JLN Barratt, K Houghton, M Plucinski, J Kelley, S Casillas, C Bennett, C Snider, R Tuladhar, J Zhang, B Clemons, S Madison-Antenucci, A Russell, E Cebelinski, J Haan, T Robinson, MJ Arrowood, E Talundzic, RS Bradbury, and Y Qvarnstrom (2020) Evaluation of an ensemble-based distance statistic for clustering MLST datasets using epidemiologically defined clusters of cyclosporiasis. Epidemiology & Infection: 148, e172, 1–10. https://doi.org/10.1017/
S0950268820001697
```

## Getting started

>These instructions will help you set up and run this code on your local machine for development and testing purposes. See deployment for notes for information on how to deploy the project on a live system.

This code was developed and tested using a Mac running OSX Catalina 10.15.3. Subsequent instructions are provided _only_ for installing it on an OSX system.

First create a local copy of this repository:

`git clone git@github.com:Joel-Barratt/Complete-Cyclospora-typing-workflow.git` 

### Prerequisites for OSX Catalina

>Prerequisites for installation of this code

#### Xcode Command Line Tools

Install Xcode

```bash
xcode-select --install
```
Check Xcode is included in your $PATH (e.g., /Library/Developer/CommandLineTools)

```bash
xcode-select -p
```

#### Local R package

Go to [CRAN](https://cran.r-project.org/bin/macosx/), download and install `R-3.6.3.pkg (notarized, for Catalina)`

Check that R is correctly installed

```bash
R -h   # this should return a help print out with options for R  
```

#### Other prerequisites

There are several other essential prerequisites and dependencies that must be installed in order to run this workflow correctly. This section will be updated in the future. Please stand by for updates.


### Running this code

>While in the Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_BETA_1.001 directory with all the files from the cloned github run (example):

```bash
bash MODULE_1_hap_caller.sh \
-D TEST_DATA \
-C  /path/to/where/you/unzipped/this/software \
-L 10 \
-R 500 \
-T 10 \
-H Y
```
> This will identify the haplotypes in 2 test samples. A produce a genotype file should be produced for two specimens _C_TEST1_20_ and _C_TEST2_20_ in the Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_BETA_1.001/HAPLOTYPE_CALLER_CYCLO_V2/SPECIMEN_GENOTYPES/ directory. A new haplotype data sheet will also be generated in the Cyclospora_typing_workflow_MacOS_High_Sierra_BETA_1.001/haplotype_sheets directory. To understand the arguments supplied while in the  Complete_Cyclospora_typing_workflow_MacOS_High_Sierra_BETA_1.001 directory with all the files from the cloned github run:

```bash
bash MODULE_1_hap_caller.sh -h
```
>This should provide a brief help listing all possible arguments.

## Additional Information on Project Participation

For additional detailed information on how to run these workflows and how to participate in this project please contact the authors of this [manuscript](https://www.cambridge.org/core/journals/parasitology/article/genotyping-genetically-heterogeneous-cyclospora-cayetanensis-infections-to-complement-epidemiological-case-linkage/0C51FBFFB172DF50357C1D171E9B8657) and this [manuscript](https://www.cambridge.org/core/journals/epidemiology-and-infection/article/evaluation-of-an-ensemblebased-distance-statistic-for-clustering-mlst-datasets-using-epidemiologically-defined-clusters-of-cyclosporiasis/F13FB2483E3DE1CF5C33706B4A1A7182#fndtn-information) directly to ensure your work remains compatible with the CDCs Cyclospora cayetanensis genotyping procedures. 


## Deployment

<!-- need to update once on SciComp and CDCgov github -->
This section will be updated in the future.


## Contributing

 <!-- need to update @Joel -->

This section will be updated in the future.


## Acknowledgments

* Sincere thanks to the scientists and programmers who developed the various pre-requisite tools used in this workflow.


## Public Domain

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This soruce code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.


## Privacy
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
Surveillance Platform [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/privacy.html](http://www.cdc.gov/privacy.html).

## Contributing
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page are subject to the [Presidential Records Act](http://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and our [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
