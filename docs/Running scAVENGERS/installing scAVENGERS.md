# Installing scAVENGERS
## Installation
The required softwares and packages are provided as a conda environment. Conda (https://www.anaconda.com/products/distribution) helps users to manage multiple packages easily. Running below will clone scAVENGERS repository create an environment named 'scavengers.'
```
git clone https://github.com/kaistcbfg/scAVENGERS
conda env create -f scAVENGERS/envs/environment.yaml
```
## Requirements
The program is successfully tested under the environment below. Some of the requirements are automatically installed when `scAVENGERS pipeline` is executed.

### External tools
|name|version|remarks|
|---|---|---|
|freebayes|1.3.5|
|strelka2|2.9.10|Automatically downloaded by the shell script when the pipeline is excuted.|
|vartrix|1.1.22|
|bcftools|1.12|
|souporcell|2.0|Automatically downloaded and compiled by the shell script when the pipeline is executed. Troublet program will only be used among the whole pipeline.|
|snakemake|7.3.8|

### Python and its packages
|name|version|remarks|
|---|---|---|
|python|3.8.12|
|numba|0.54.1|
|numpy|1.20.3|
|pandas|1.4.0|
|scikit-learn|1.0.2|
|scipy|1.7.3|
|pystan|2.19.1.1|Provided in another environment used for ambient variant detection|

### Rust
|name|version|remarks|
|---|---|---|
|rust|1.59.0|

These dependencies are provided in a conda environment below.
```
name: scavengers
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=1_gnu
  - aioeasywebdav=2.4.0=py38h578d9bd_1001
  - aiohttp=3.8.1=py38h0a891b7_1
  - aiosignal=1.2.0=pyhd8ed1ab_0
  - amply=0.1.4=py_0
  - appdirs=1.4.4=pyh9f0ad1d_0
  - async-timeout=4.0.2=pyhd8ed1ab_0
  - attmap=0.13.2=pyhd8ed1ab_0
  - attrs=21.4.0=pyhd8ed1ab_0
  - backports=1.0=py_2
  - backports.functools_lru_cache=1.6.4=pyhd8ed1ab_0
  - bc=1.07.1=h7f98852_0
  - bcftools=1.12=h45bccc9_1
  - bcrypt=3.2.0=py38h497a2fe_2
  - binutils_impl_linux-64=2.36.1=h193b22a_2
  - boto3=1.21.38=pyhd8ed1ab_0
  - botocore=1.24.38=pyhd8ed1ab_0
  - brotlipy=0.7.0=py38h0a891b7_1004
  - bzip2=1.0.8=h7f98852_4
  - c-ares=1.18.1=h7f98852_0
  - ca-certificates=2022.6.15=ha878542_0
  - cachetools=5.0.0=pyhd8ed1ab_0
  - certifi=2022.6.15=py38h578d9bd_0
  - cffi=1.15.0=py38h3931269_0
  - charset-normalizer=2.0.12=pyhd8ed1ab_0
  - coin-or-cbc=2.10.7=h3786ebc_0
  - coin-or-cgl=0.60.3=he7e83c3_2
  - coin-or-clp=1.17.6=h256e9bb_3
  - coin-or-osi=0.108.6=h3b589db_2
  - coin-or-utils=2.11.6=h573740c_0
  - coincbc=2.10.7=0_metapackage
  - colorama=0.4.4=pyh9f0ad1d_0
  - conda=4.13.0=py38h578d9bd_1
  - conda-package-handling=1.8.1=py38h0a891b7_1
  - configargparse=1.5.3=pyhd8ed1ab_0
  - connection_pool=0.0.3=pyhd3deb0d_0
  - cryptography=36.0.2=py38h2b5fc30_1
  - datrie=0.8.2=py38h497a2fe_3
  - decorator=5.1.1=pyhd8ed1ab_0
  - defusedxml=0.7.1=pyhd8ed1ab_0
  - docutils=0.18.1=py38h578d9bd_1
  - dropbox=11.29.0=pyhd8ed1ab_0
  - filechunkio=1.8=py_2
  - filelock=3.6.0=pyhd8ed1ab_0
  - freebayes=1.3.5=py38ha193a2f_3
  - frozenlist=1.3.0=py38h0a891b7_1
  - ftputil=5.0.3=pyhd8ed1ab_0
  - gcc_impl_linux-64=11.2.0=h82a94d6_15
  - gitdb=4.0.9=pyhd8ed1ab_0
  - gitpython=3.1.27=pyhd8ed1ab_0
  - google-api-core=2.5.0=pyhd8ed1ab_0
  - google-api-python-client=2.43.0=pyhd8ed1ab_0
  - google-auth=2.6.3=pyh6c4a22f_0
  - google-auth-httplib2=0.1.0=pyhd8ed1ab_0
  - google-cloud-core=2.2.2=pyh6c4a22f_0
  - google-cloud-storage=2.1.0=pyh6c4a22f_0
  - google-crc32c=1.1.2=py38h8838a9a_2
  - google-resumable-media=2.1.0=pyh6c4a22f_0
  - googleapis-common-protos=1.56.0=py38h578d9bd_0
  - grpcio=1.45.0=py38ha0cdfde_0
  - gsl=2.6=he838d99_2
  - htslib=1.12=h9093b5e_1
  - httplib2=0.20.4=pyhd8ed1ab_0
  - icu=70.1=h27087fc_0
  - idna=3.3=pyhd8ed1ab_0
  - importlib-metadata=4.11.3=py38h578d9bd_1
  - importlib_resources=5.6.0=pyhd8ed1ab_1
  - iniconfig=1.1.1=pyh9f0ad1d_0
  - jinja2=3.1.1=pyhd8ed1ab_0
  - jmespath=1.0.0=pyhd8ed1ab_0
  - joblib=1.1.0=pyhd8ed1ab_0
  - jsonschema=4.4.0=pyhd8ed1ab_0
  - jupyter_core=4.9.2=py38h578d9bd_0
  - kernel-headers_linux-64=2.6.32=he073ed8_15
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.19.3=h3790be6_0
  - ld_impl_linux-64=2.36.1=hea4e1c9_2
  - libarchive=3.5.2=hb890918_2
  - libblas=3.9.0=12_linux64_openblas
  - libcblas=3.9.0=12_linux64_openblas
  - libcrc32c=1.1.2=h9c3ff4c_0
  - libcurl=7.83.1=h7bff187_0
  - libdeflate=1.7=h7f98852_5
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=h516909a_1
  - libffi=3.4.2=h7f98852_5
  - libgcc-devel_linux-64=11.2.0=h0952999_15
  - libgcc-ng=12.1.0=h8d9b700_16
  - libgfortran-ng=11.2.0=h69a702a_11
  - libgfortran5=11.2.0=h5c6108e_11
  - libgomp=12.1.0=h8d9b700_16
  - libiconv=1.16=h516909a_0
  - liblapack=3.9.0=12_linux64_openblas
  - liblapacke=3.9.0=12_linux64_openblas
  - libllvm11=11.1.0=hf817b99_2
  - libmamba=0.24.0=hd8a31e3_1
  - libmambapy=0.24.0=py38h923e62a_1
  - libnghttp2=1.47.0=h727a467_0
  - libnsl=2.0.0=h7f98852_0
  - libopenblas=0.3.18=pthreads_h8fe5266_0
  - libprotobuf=3.20.0=h6239696_0
  - libsanitizer=11.2.0=he4da1e4_15
  - libsodium=1.0.18=h36c2ea0_1
  - libsolv=0.7.22=h6239696_0
  - libssh2=1.10.0=ha56f1ee_2
  - libstdcxx-ng=12.1.0=ha89aaad_16
  - libxml2=2.9.14=h22db469_0
  - libzlib=1.2.11=h36c2ea0_1013
  - llvmlite=0.37.0=py38h4630a5e_1
  - logmuse=0.2.6=pyh8c360ce_0
  - lz4-c=1.9.3=h9c3ff4c_1
  - lzo=2.10=h516909a_1000
  - mamba=0.24.0=py38h1abaa86_1
  - markupsafe=2.1.1=py38h0a891b7_1
  - multidict=6.0.2=py38h0a891b7_1
  - nbformat=5.3.0=pyhd8ed1ab_0
  - ncurses=6.2=h58526e2_4
  - nomkl=1.0=h5ca1d4c_0
  - numba=0.54.1=py38h4bf6c61_0
  - numpy=1.20.3=py38h9894fe3_1
  - oauth2client=4.1.3=py_0
  - openssl=1.1.1p=h166bdaf_0
  - packaging=21.3=pyhd8ed1ab_0
  - pandas=1.4.0=py38h43a58ef_0
  - parallel=20211222=ha770c72_0
  - paramiko=2.10.3=pyhd8ed1ab_0
  - peppy=0.31.2=pyhd8ed1ab_2
  - perl=5.32.1=1_h7f98852_perl5
  - pip=21.3.1=pyhd8ed1ab_0
  - plac=1.3.5=pyhd8ed1ab_0
  - pluggy=1.0.0=py38h578d9bd_3
  - ply=3.11=py_1
  - prettytable=3.2.0=pyhd8ed1ab_0
  - protobuf=3.20.0=py38hfa26641_4
  - psutil=5.9.0=py38h0a891b7_1
  - pulp=2.6.0=py38h578d9bd_1
  - py=1.11.0=pyh6c4a22f_0
  - pyasn1=0.4.8=py_0
  - pyasn1-modules=0.2.7=py_0
  - pybind11-abi=4=hd8ed1ab_3
  - pycosat=0.6.3=py38h0a891b7_1010
  - pycparser=2.21=pyhd8ed1ab_0
  - pygments=2.11.2=pyhd8ed1ab_0
  - pynacl=1.5.0=py38h0a891b7_1
  - pyopenssl=22.0.0=pyhd8ed1ab_0
  - pyparsing=3.0.8=pyhd8ed1ab_0
  - pyrsistent=0.18.1=py38h0a891b7_1
  - pysftp=0.2.9=py_1
  - pysocks=1.7.1=py38h578d9bd_5
  - pytest=7.1.1=py38h578d9bd_1
  - python=3.8.12=hb7a2778_2_cpython
  - python-dateutil=2.8.2=pyhd8ed1ab_0
  - python-fastjsonschema=2.15.3=pyhd8ed1ab_0
  - python-irodsclient=1.1.3=pyhd8ed1ab_0
  - python_abi=3.8=2_cp38
  - pytz=2021.3=pyhd8ed1ab_0
  - pyu2f=0.1.5=pyhd8ed1ab_0
  - pyyaml=6.0=py38h0a891b7_4
  - ratelimiter=1.2.0=py_1002
  - readline=8.1=h46c0cb4_0
  - reproc=14.2.4=h295c915_1
  - reproc-cpp=14.2.4=h295c915_1
  - requests=2.27.1=pyhd8ed1ab_0
  - retry=0.9.2=py_0
  - rsa=4.8=pyhd8ed1ab_0
  - ruamel_yaml=0.15.100=py38h27cfd23_0
  - rust=1.59.0=h30be4a0_0
  - rust-std-x86_64-unknown-linux-gnu=1.59.0=hc1431ca_0
  - s3transfer=0.5.2=pyhd8ed1ab_0
  - samtools=1.3.1=0
  - scikit-learn=1.0.2=py38h1561384_0
  - scipy=1.7.3=py38h56a6a73_0
  - setuptools=60.9.3=py38h578d9bd_0
  - six=1.16.0=pyh6c4a22f_0
  - slacker=0.14.0=py_0
  - smart_open=5.2.1=pyhd8ed1ab_0
  - smmap=3.0.5=pyh44b312d_0
  - snakemake=7.3.8=hdfd78af_0
  - snakemake-minimal=7.3.8=pyhdfd78af_0
  - sqlite=3.37.0=h9cd32fc_0
  - stone=3.3.1=pyhd8ed1ab_0
  - stopit=1.1.2=py_0
  - sysroot_linux-64=2.12=he073ed8_15
  - tabix=1.11=hdfd78af_0
  - tabixpp=1.1.0=he28291e_6
  - tabulate=0.8.9=pyhd8ed1ab_0
  - threadpoolctl=3.0.0=pyh8a188c0_0
  - tk=8.6.11=h27826a3_1
  - tomli=2.0.1=pyhd8ed1ab_0
  - toposort=1.7=pyhd8ed1ab_0
  - tqdm=4.62.3=pyhd8ed1ab_0
  - traitlets=5.1.1=pyhd8ed1ab_0
  - typing-extensions=4.1.1=hd8ed1ab_0
  - typing_extensions=4.1.1=pyha770c72_0
  - ubiquerg=0.6.1=pyh9f0ad1d_0
  - uritemplate=4.1.1=pyhd8ed1ab_0
  - urllib3=1.26.9=pyhd8ed1ab_0
  - vartrix=1.1.22=h27ee8bf_0
  - vcflib=1.0.2=h3198e80_5
  - veracitools=0.1.3=py_0
  - wcwidth=0.2.5=pyh9f0ad1d_2
  - wheel=0.37.1=pyhd8ed1ab_0
  - wrapt=1.14.0=py38h0a891b7_1
  - xz=5.2.5=h516909a_1
  - yaml=0.2.5=h7f98852_2
  - yaml-cpp=0.7.0=h27087fc_1
  - yarl=1.7.2=py38h0a891b7_2
  - yte=1.2.1=py38h578d9bd_0
  - zipp=3.8.0=pyhd8ed1ab_0
  - zlib=1.2.11=h36c2ea0_1013
  - zstd=1.5.2=h8a70e8d_1
prefix: /home/rotate/anaconda3/envs/scavengers
```
2. After you cloned the repository, you can run scAVENGERS by executing the command below.
```
$SCAVENGERS_DIRECTORY/scAVENGERS [pipeline|cluster] [options]
```