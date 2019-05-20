#!/bin/bash
## substitute pathnames into config file
#!/bin/bash
SAMTOOLS_BIN=`which samtools`
BOWTIE_BIN=`which bowtie`
BOWTIE_BUILD_BIN=`which bowtie-build`
BLAT_BIN=`which blat`
FATOTWOBIT_BIN=`which faToTwoBit`
GMAP_BIN=`which gmap`
GMAP_BUILD_BIN=`which gmap_build`
GMAP_SETUP_BIN=$GMAP_BUILD_BIN
R_BIN=`which R`
RSCRIPT_BIN=`which Rscript`
DEFUSE_BIN=`which defuse_run.pl`
DEFUSE_PATH=`python -c "import os.path;  print(os.path.dirname(os.path.dirname(os.path.realpath(\"$DEFUSE_BIN\"))))"`
cat $1 | sed "s#__DEFUSE_PATH__#${DEFUSE_PATH}#" | \
sed "s#__SAMTOOLS_BIN__#${SAMTOOLS_BIN}#" | \
sed "s#__BOWTIE_BIN__#${BOWTIE_BIN}#" | \
sed "s#__BOWTIE_BUILD_BIN__#${BOWTIE_BUILD_BIN}#" | \
sed "s#__BLAT_BIN__#${BLAT_BIN}#"| \
sed "s#__FATOTWOBIT_BIN__#${FATOTWOBIT_BIN}#" | \
sed "s#__GMAP_BIN__#${GMAP_BIN}#" | \
sed "s#^gmap_setup_bin#gmap_build_bin#" | \
sed "s#__GMAP_SETUP_BIN__#${GMAP_SETUP_BIN}#" | \
sed "s#__GMAP_BUILD_BIN__#${GMAP_BUILD_BIN}#" | \
sed "s#__R_BIN__#${R_BIN}#" | \
sed "s#__RSCRIPT_BIN__#${RSCRIPT_BIN}#" > $2
export DATASET_DIRECTORY=`grep '^dataset_directory' $1 | awk '{print \$NF}'`
echo "$DATASET_DIRECTORY"
