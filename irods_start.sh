#!/usr/bin/env bash

ARGS=2
E_BADARGS=85

GTEX_IRODS_DIR="/renci/irods/gtex/"
DOCKER_DIR="/home/dockeruser/venv"
CONTAINER_NAME="dc_jupyter"

if [$(docker rm -f dc_jupyter)]
then
    echo "Deleting existing iRODS container."
fi
      
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
# docker exec ${CONTAINER_NAME} sudo /home/dockeruser/venv/bin/pip3 install --upgrade pip
docker exec -ti ${CONTAINER_NAME} sudo pip3 install pandas numpy scipy matplotlib
docker cp /home/ubuntu/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz ${CONTAINER_NAME}:"${DOCKER_DIR}"

docker exec ${CONTAINER_NAME} git clone https://github.com/unmtransinfo/expression-profiles.git 
docker exec ${CONTAINER_NAME} sh -c "cd /home/dockeruser/expression-profiles && git checkout containerize"
docker exec -ti ${CONTAINER_NAME} sh -c "cd /home/dockeruser/expression-profiles && ./Go_rnaseq_prep.sh"

# Put generated files outside of container for notebook usage
docker cp ${CONTAINER_NAME}:"/home/dockeruser/expression-profiles/data/output/" ./data/
