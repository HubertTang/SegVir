#Setup scripts to set LD_LIBRARY_PATH
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d

cat <<EOF >$CONDA_PREFIX/etc/conda/activate.d/LD_PATH.sh
export LD_LIBRARY_PATH_CONDA_BACKUP=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
EOF

cat <<EOF >$CONDA_PREFIX/etc/conda/deactivate.d/LD_PATH.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_CONDA_BACKUP
EOF

#Modify the permissions of ORFfinder
chmod a+x ORFfinder
