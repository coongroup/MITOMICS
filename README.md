# MITOMICS
Code base for MITOMICS project

### Directory Structure
```bash
.
├── LICENSE
├── README.md
├── reference
│   └── color_codes.txt
└── src
    ├── analysis
    │   ├── IJM
    │   │   ├── README.md
    │   │   ├── dynamic_range_boxplots.R
    │   │   ├── multi_omics_heatmap.R
    │   └── YS
    │       ├── README.md
    │       ├── functions.R
    │       ├── make.R
    │       ├── plan.R
    │       └── templates
    │           └── embedding_plot.html
    └── processing
        ├── DRB
        │   └── ImputeMissingH3kValues
        │       └── ImputeMissingH3kValues
        │           ├── App.config
        │           ├── Biomolecule.cs
        │           ├── Condition.cs
        │           ├── Global.cs
        │           ├── ImputeMissingH3kValues.csproj
        │           ├── LfqValues.cs
        │           ├── bin
        │           │   └── Debug
        │           │       ├── ImputeMissingH3kValues.exe
        │           │       ├── ImputeMissingH3kValues.exe.config
        │           │       ├── ImputeMissingH3kValues.pdb
        │           │       ├── LumenWorks.Framework.IO.dll
        │           │       ├── LumenWorks.Framework.IO.xml
        │           │       ├── MathNet.Numerics.dll
        │           │       └── MathNet.Numerics.xml
        │           └── obj
        │               └── Debug
        │                   └── DesignTimeResolveAssemblyReferencesInput.cache
        ├── IJM
        │   ├── lipidomics
        │   │   ├── 1.0-lipidomics_processing.py
        │   │   ├── 2.0-lipidomics_processing.py
        │   │   ├── 3.0-lipidomics_processing.py
        │   │   ├── 4.0-lipidomics_processing.py
        │   │   ├── README.md
        │   ├── metabolomics
        │   │   ├── 1.0-metabolomics_processing.py
        │   │   ├── 2.0-metabolomics_processing.py
        │   │   ├── README.md
        │   └── proteomics
        │       ├── 1.0-proteomics_processing.py
        │       ├── 2.0-proteomics_processing.py
        │       ├── 3.0-proteomics_processing.R
        │       ├── 4.0-proteomics_processing.py
        │       ├── README.md
        └── KAO
            ├── GC_runOrderCorrection_R_20200810.R
            ├── README.R
            └── Time_stamp.R

```

### Data availability

All associated mass spectrometry RAW files and search results were deposited into MassIVE data repository and can be accessed using the following link ftp://MSV000086685@massive.ucsd.edu (currently requires reviewer credentials)