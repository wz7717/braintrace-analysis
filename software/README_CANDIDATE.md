# Using the frozen software in this candidate

The [November 20 software entry](README.md) and [version binding](../PAPER_VERSION_BINDING.json) identify the current macaque model. The original integration inventory below describes its historical release.

The software files derive from the reviewed deployed commit 2bf9563. Use the pinned Python 3.11.9 environment and the app dependency and package instructions. Set up the compatible scientific environment, then assemble the files into a separate runtime directory. [Streamlit deployment instructions](docs/STREAMLIT_CLOUD.md) · [Frozen input contract](docs/INPUT_FORMATS.md)

From this software directory on a host with the prescribed dependencies:

~~~text
python -B tools/prepare_runtime.py --assets runtime_assets --out PATH_TO_NEW_RUNTIME
~~~

Streamlit can use the stated source layout and prepare its runtime itself. Select the macaque route and run the supplied GSE128537 / 10_10 and GSE179328 / LID_100660 examples. They start with the exact prepared expression, masks and mapping from the scoped web edition. The fixed reference choices are beam_015 and beam_023. [Repository and model documentation](README.md)

The human CEL route uses five linked payloads from the original release, covering the RAW and RANK models, the frozen PLM reference and two provenance indices. Their checksums and original-release pointers are in [the human input reference list](EXISTING_HUMAN_INPUT_REFERENCES.csv). Follow the conditions of that release and place the matching files at the indicated paths for a local CEL run. Also set up the documented R/PLM environment and obtain the original public CEL. The older software README describes the full deployed set.

For GSM755211 obtain the original CEL/GZ from GEO and verify the recorded SHA-256. The existing Demo already binds the full frozen route. [GEO CEL source and checksum](GEO_CEL_SOURCE.json) · [Existing human fixed assets and usage terms](https://github.com/wz7717/braintrace-web/tree/9059a83c9d8c153982c75e6f5815fe932f89a9b4/human_fixed_assets/HUMAN_FIXED_ASSET_SET_v1) · [Demo](https://braintrace-web.streamlit.app/)
