# Navigating the conjugated metabolome

We mined 1.32 billion MS/MS spectra across public metabolomics repositories using reverse spectral searching and delta-mass infererence to map **metabolite conjugations**. We generated structural hypotheses for 24,227,439 spectra, encompassing 217,291 substructure pairs with dual spectral support and 3,412,720 additional candidates with single-match support. These results provide a pan-repository map of potential conjugation chemistry, establish a resource for structural discovery, and offer a framework to further explore the scope and diversity of the conjugated metabolome.

<img src="gitfigs/workflow.svg" width="850"/>

<table width="100%" border="0">
<tr>
  <td width="65%" align="center" valign="center">
    <img src="gitfigs/revcos.svg" />
  </td>
  <td width="35%" align="center">
    <img src="gitfigs/class.svg" />
  </td>
  <!-- <td width="30%" align="center">
    <img src="gitfigs/distributions.png" />
  </td> -->
</tr>
</table>

## Reverse spectral search
`Reverse spectral search` is an MS/MS similarity framework originally proposed for spectral identification and later extended to improve robustness to chimeric spectra. Here, we repurpose reverse spectral searching for substructure annotation. 
Two implementations of reverse cosine similarity were provided:
- [matchms](https://github.com/matchms/matchms)-based reverse cosine: [revcos.py](https://github.com/Philipbear/conjugated_metabolome/blob/main/bin/revcos_matchms/revcos.py)
- [Flash](https://github.com/YuanyueLi/FlashEntropySearch) reverse cosine: [flash_cos.py](https://github.com/Philipbear/conjugated_metabolome/blob/main/bin/main/flash_cos.py)


## Data availability
- Conjugate search results of 149.9 million clustered MS/MS: [Zenodo link](https://zenodo.org/records/17245769)
- GNPS conjugated metabolome libraries: [GNPS external library link](https://external.gnps2.org/gnpslibrary)
- Web app to explore conjugated metabolome results: [Web app link](https://conjugated-metabolome.streamlit.app)
- Datasets
  - [MSV000086131](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000086131) (Mammalian feces)
  - [MSV000082261](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000082261) (Diabetes study)
  - [MSV000082433](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000082433) (Human feces, losartan conjugates)
  - [MSV000083306](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000083306) (Tomato seedling extracts)
  - [MSV000098638](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000098638) (Microbial cultures)
  - [MSV000099690](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000099690) (Chemical synthesis & RT matching)
  - [MSV000096359](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000096359) (Human urine cohort for drug conjugates)
- MS/MS reference libraries
  - GNPS ([ALL_GNPS_NO_PROPOGATED.msp](https://external.gnps2.org/gnpslibrary), downloaded on Nov 11, 2024)
  - MassBank ([MassBank_NIST.msp](https://github.com/MassBank/MassBank-data/releases/tag/2024.06), 2024.06 release)
  - MoNA ([MoNA-export-LC-MS-MS_Spectra.msp](https://mona.fiehnlab.ucdavis.edu/downloads), downloaded on Oct 17, 2024)
  - NIST20 (Commercially available)


## Citation
> S. Xing et al. Navigating the conjugated metabolome. To be preprinted. https://github.com/Philipbear/conjugated_metabolome


## License
This work is licensed under the Apache License 2.0.
