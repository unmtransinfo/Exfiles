#!/usr/bin/env bash

ARGS=2
E_BADARGS=85

GTEX_IRODS_DIR="/renci/irods/gtex/"
DOCKER_DIR="/home/dockeruser/venv"
CONTAINER_NAME="dc_jupyter"

# The paths where the data files reside
export GTEX_ANNOTATIONS="${GTEX_IRODS_DIR}/Annotations"
export GTEX_RNASEQ="${GTEX_IRODS_DIR}/RNA-seq"

if [ $# -ne "$ARGS" ]
then
  echo "Usage: `basename $0` <username> <password>"
  exit $E_BADARGS
fi

# The script to start the irods Jupyter container
./helium run jupyter -venv -U $1 -P $2

# Add additional setup to Irods container in order to run UNM notebook
echo 'Creating directories and installing dependencies...'
docker exec -ti ${CONTAINER_NAME} sudo /home/dockeruser/venv/bin/pip3 install --upgrade pip
docker exec -ti ${CONTAINER_NAME} sudo  /home/dockeruser/venv/bin/pip3 install scipy matplotlib
docker cp /home/ubuntu/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz ${CONTAINER_NAME}:"${GTEX_IRODS_DIR}/RNA-seq/gtex/"

docker exec -ti ${CONTAINER_NAME} git clone https://github.com/unmtransinfo/expression-profiles.git &&
       sh -c "cd '${HOME}/expression-profiles'" &&
       sh -c "echo 'CURRENT DIRECTORY IS: $(pwd)'" &&
       sh -c "git checkout containerize" &&
       sh -c "./Go_rnaseq_prep.sh"



# Copy into container missing GTEx content and notebook
# echo 'Copying required analysis files and UNM Jupyter notebook'
# docker cp /home/data/GTEx/data/biomart_ENSG2NCBI.tsv ${CONTAINER_NAME}:"$INDOCKER_DIR/gtex/"
# docker cp /home/data/GTEx/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz ${CONTAINER_NAME}:"$INDOCKER_DIR/gtex/"
# docker cp gtex_rnaseq_prep.ipynb ${CONTAINER_NAME}:"$INDOCKER_DIR/notebooks/"
