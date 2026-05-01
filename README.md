# Navigating the conjugated metabolome

We mined 1.32 billion MS/MS spectra across public metabolomics repositories using reverse spectral searching and delta-mass inference to map **metabolite conjugations**. We generated structural hypotheses for 24,227,439 spectra, encompassing 217,291 substructure pairs with dual spectral support and 3,412,720 additional candidates with single-match support. These results provide a pan-repository map of potential conjugation chemistry, establish a resource for structural discovery, and offer a framework to further explore the scope and diversity of the conjugated metabolome.

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

In both implementations, both the similarity score and the number of matched peaks will be output. Codes were tested using Python 3.10 on macOS (14.6, M2 Max) and Linux (Ubuntu 20.04). Packages of `numpy` and `numba` are required. For the `matchms`-based implementation, each spectral comparison takes several milliseconds.


## Interpreting predictions
Reverse cosine matching identifies substructure-level relationships between a query spectrum and a reference. We recommend the following before treating any prediction as a candidate identification:
1. **Treat every match as substructure-level rather than exact.** A reverse cosine match to a specific reference compound (e.g., chenodeoxycholic acid) should be read as evidence for the broader substructure class to which that compound belongs (e.g., a dihydroxy bile acid scaffold), since closely related isomers often produce nearly indistinguishable MS/MS spectra. Stereochemistry, double-bond position, and the exact regiochemistry of conjugation are not resolved by this approach.
2. **Consult the compound class-specific FDRs** reported in Supplementary Table 2 of the manuscript before interpreting an individual annotation. For chemical classes with high empirical FDRs, additional filtering or orthogonal validation is essential. In such cases, we recommend raising the score and matched-peak thresholds beyond the defaults used here, and limiting interpretation to predictions that also satisfy the chemical-context filters described below.
3. **Apply chemical and biological context filters** to refine candidate sets. Effective filters include, but are not limited to: requiring co-occurrence of a predicted conjugate with its proposed parent compound within the same LC-MS file (as we used for the drug conjugate analysis); restricting candidates to those carrying diagnostic fragment ions or neutral losses for a moiety of interest using MassQL queries (as demonstrated for the steroid-phosphoethanolamine conjugates); and requiring detection across multiple independent datasets to filter out study-specific artefacts.
4. **Distinguish in vivo conjugates from in-source or in-droplet artefacts.** Some conjugates may form post-chromatography within the electrospray source or arise from co-elution of unrelated species. Retention time matching against a synthetic standard is the most direct way to distinguish genuine biological conjugates from such artefacts; consistent detection across orthogonal chromatographic methods provides a useful preliminary check before committing to chemical synthesis.
5. **Always validate with synthetic standards before drawing biological conclusions.** The synthetic accessibility of conjugation chemistry is one of the principal advantages of this resource: amide, ester, and related linkages are typically straightforward to prepare in milligram quantities, enabling MS/MS and retention time matching at level 1 confidence45 for any prediction of biological interest.


## Data availability
- Conjugate search results of 149.9 million clustered MS/MS: [Zenodo link](https://zenodo.org/records/17245769)
- GNPS conjugated metabolome libraries: [Zenodo link](https://zenodo.org/records/19837386)
- Web app to explore conjugated metabolome results: [Web app link](https://conjugated-metabolome.gnps2.org)
- Datasets
  - [MSV000086131](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000086131) (Mammalian feces)
  - [MSV000082261](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000082261) (Diabetes study)
  - [MSV000082433](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000082433) (Human feces, losartan conjugates)
  - [MSV000083306](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000083306) (Tomato seedling extracts)
  - [MSV000098638](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000098638) (Microbial cultures)
  - [MSV000096359](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000096359) (Human urine cohort for drug conjugates)
  - [MSV000097935](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000097935) (Human serum samples for drug conjugates)
  - [MSV000099690](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000099690) (Chemical synthesis & RT matching)
- MS/MS reference libraries
  - GNPS ([ALL_GNPS_NO_PROPOGATED.msp](https://external.gnps2.org/gnpslibrary), downloaded on Nov 11, 2024)
  - MassBank ([MassBank_NIST.msp](https://github.com/MassBank/MassBank-data/releases/tag/2024.06), 2024.06 release)
  - MoNA ([MoNA-export-LC-MS-MS_Spectra.msp](https://mona.fiehnlab.ucdavis.edu/downloads), downloaded on Oct 17, 2024)
  - NIST20 (Commercially available)


## Citation
> Shipei Xing, Abubaker Patan, Julius Agongo, Harsha Gouda, Vincent Charron-Lamoureux, Yasin El Abiead, Zhewen Hu, Haoqi Nina Zhao, Ipsita Mohanty, Jasmine Zemlin, Wilhan Donizete Gonçalves Nunes, Lindsey A Burnett, Mingxun Wang, Dionicio Siegel, Pieter C Dorrestein. **Navigating the conjugated metabolome.** Preprint at https://doi.org/10.64898/2026.02.06.704496 (2026).


## License
This work is licensed under the Apache License 2.0.
