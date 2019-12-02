# Get Conda if you don't have it
set -ex
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

OPT_DIR=${DIR}/opt

mkdir -p OPT_DIR

CONDA=Miniconda3-latest-Linux-x86_64.sh
PYTHON_DIR=${OPT_DIR}/miniconda3

if [[ ! -d ${PYTHON_DIR} ]]; then
wget -q https://repo.continuum.io/miniconda/${CONDA}\
    && sh ${CONDA} -b -p ${PYTHON_DIR}\
    && rm -f ${CONDA}
fi

# Build an environment according to the yml
$PYTHON_DIR/condabin/conda env create -f STR_environment.yml


