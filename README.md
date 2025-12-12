# Pan-repository conjugated metabolome

Life’s chemistry is far richer than current biochemical maps reveal. While metabolomics has catalogued tens of thousands of small molecules, deeper layers of biochemical diversity remain largely hidden. Within this uncharted space are **conjugated metabolites** - molecules formed when two or more metabolites are covalently joined through amidation, esterification, ether formation, alkylamination, or related linkages. Such conjugates act as microbial signals, detoxification intermediates, or endogenous regulators, yet their diversity is poorly characterized. Using a reverse spectral search strategy, a data centric approach that compares pairs of MS/MS spectra and filters by substructure masses, we inferred paired substructures for 24,227,439 spectra from public metabolomics repositories.

<img src="gitfigs/workflow.svg" width="850"/>

<table width="100%" border="0">
<tr>
  <td width="60%" align="center" valign="center">
    <img src="gitfigs/revcos.svg" />
  </td>
  <td width="40%" align="center">
    <img src="gitfigs/class.svg" />
  </td>
  <!-- <td width="30%" align="center">
    <img src="gitfigs/distributions.png" />
  </td> -->
</tr>
</table>


## Data
- MS/MS reference libraries
  - GNPS ([ALL_GNPS_NO_PROPOGATED.msp](https://external.gnps2.org/gnpslibrary), downloaded on Nov 11, 2024)
  - MassBank ([MassBank_NIST.msp](https://github.com/MassBank/MassBank-data/releases/tag/2024.06), 2024.06 release)
  - MoNA ([MoNA-export-LC-MS-MS_Spectra.msp](https://mona.fiehnlab.ucdavis.edu/downloads), downloaded on Oct 17, 2024)
  - NIST20 (Commercially available)
- Datasets
  - [MSV000086131](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000086131) (Mammalian feces)
  - [MSV000082261](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000082261) (Diabetes study)
  - [MSV000082433](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000082433) (Human feces, losartan conjugates)
  - [MSV000083306](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000083306) (Tomato seedling extracts)
  - [MSV000098638](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000098638) (Microbial cultures)
  - [MSV000099690](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000099690) (Chemical synthesis & RT matching)
- Delta masses
  - as summarized in [this ref](https://www.biorxiv.org/content/10.1101/2024.10.07.617109v2)
  - LIPID MAPS ([LMSD](https://www.lipidmaps.org/databases/lmsd/download), 2024-10-23)


## Web app
For users to quickly explore the conjugated metabolome results, we built a [web app](https://conjugated-metabolome.streamlit.app).


## Citation
> S. Xing. Navigation of the pan-repository conjugated metabolome. https://github.com/Philipbear/conjugated_metabolome


## License
This work is licensed under the Apache License 2.0.
