# Data acquisition and access

Use the exact-source and study-entry registers to retrieve original objects. Check the recorded SHA-256 before normalizing or unpacking. A study accession is an upstream entry; it supplies no exact identity for a project derivative on its own. Some moving links carry a recorded snapshot hash. The two S32 map entries distinguish the original LF digest and the locked CRLF digest. [Exact source links and hashes](../data/INTEGRATED_SOURCE_URLS_AND_SHA256.csv) · [Historical source identities](../data/PUBLIC_SOURCE_IDENTITIES.csv) · [Upstream study register](../data/UPSTREAM_STUDY_FILES.csv)

Example to obtain a file whose register provides a direct URL and SHA-256:

~~~powershell
curl.exe --fail --location "<official-file-URL>" --output "<local-filename>"
Get-FileHash -LiteralPath "<local-filename>" -Algorithm SHA256
~~~

## By route

| Route | Public source and use | Exact materials and access conditions |
| --- | --- | --- |
| S31 | Obtain UKBEC CEL through the GEO series; use the specified GTEx v11 TPM/annotations and HPA/HBTB v25.1 transcript archive with Ensembl 109. Use the historical digest on moving download URLs. | The portable entries are included. Obtain the seven exact derived inputs and UKBEC saved predictions from the data custodian for the applicable use. GTEx controlled data, if needed, follow dbGaP. |
| S32 Bo | PRJNA905082 provides raw sequence reads. Public fixed software references are linked by commit in the source register. | The three originally supplied processed count/VSD/specimen files remain under supplier conditions. The gene axis is included. The remaining five B5 members and unfurnished fold-level references, rules and sample bindings retain their individual access conditions. |
| S32 GEO / PsychENCODE | Obtain the original two GEO count/annotation releases. Follow the Sousa et al. PsychENCODE adult count and metadata release and the recorded Ensembl/XSAnno version. | Original source downloads can be compared to the register. The research runners are included. The exact research target matrices, masks, sample/label eligibility and unmatched research assets retain their corresponding access conditions. |
| S32 AHBA | Download the two RNA-seq archives H0351_2001_rnaseq.zip and H0351_2002_rnaseq.zip from the Allen RNA-sequencing section and compare the SHA-256 values in the source register. | The transform and split target matrix and index are included. Reuse the fixed public model and state. Obtain the precise Bo donor-equal Pearson centroids under their recorded access conditions. [Target matrix and index](../ahba/README.md) |
| S33 ACT | Use GSE104687 as the study accession. All 21 direct TPM URLs and hashes of the selected profiles matched previously through Allen file metadata. The project roster contains the exact selected links. Obtain the public fixed RANK model from its recorded Git commit. Request controlled original reads from NIAGADS DSS, accession NG00059, and access approved files through DARM. | Obtain the linked 21-specimen selection through the applicable data stewards. The configuration, original rank background, code and saved regional evidence are in the ACT directory. [ACT access](../act/README.md) |
| B17 adult GSE25219 | Locate the supplementary CEL files on GEO GSE25219 and use the gene-level GPL5175 sample pages. After you obtain the approved sample selection, fetch those sample files individually or use the original GSE25219_RAW.tar series archive. | Obtain the exact 46 feature vectors, source/target axes, original B13/B17 comparison state, labels, native cache/reference and run lock from authorized project/data custodians. A GEO download alone supplies the raw study, while the B17 scorer starts at the frozen PLM features. |

For named-recipient replay, request the exact route, needed files and scientific purpose from the project corresponding author and any relevant provider. Use their approved access method. The public tables document local replay outcomes; the provider's terms govern redistribution of submitted data and derivatives.

[NCBI GEO data-use notice](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) · [GTEx terms](https://gtexportal.org/home/license) · [HPA licence](https://www.proteinatlas.org/about/licence) · [Bo RNA-seq BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA905082) · [Sousa et al. primary study](https://doi.org/10.1126/science.aan3456) · [Allen AHBA data](https://human.brain-map.org/static/download?rw=t) · [Allen ACT data](https://aging.brain-map.org/download/index) · [Allen terms](https://alleninstitute.org/legal/terms-of-use) · [NIAGADS controlled ACT reads](https://dss.niagads.org/datasets/ng00059/)
